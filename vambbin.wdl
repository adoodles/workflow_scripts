version 1.0

workflow {
    File inputRead1Files
    String inputRead1Identifier
    String inputRead2Identifier
    String inputExtension
    File humanDB

    Int? MaxMemGB

# 1) Preprocess the reads and check their quality  * QC steps   FASTQ->FASTQ
# 2) Assemble each sample individually and get the contigs out  * megahit / metaspade contigs  FASTQ->FASTA
# 3) Concatenate the FASTA files together while making sure all contig headers stay unique, and filter away small contigs  *vamb script
# 4) Map the reads to the FASTA file to obtain BAM files  * minimap2? (RAM bug) 
# 5) Run Vamb
# 6) Postprocess the results

    # read in a file of the read1 paths
    Array[Array[String]] inputRead = read_tsv(inputReadFiles)
    
    scatter (read1 in inputRead) {
        Array[String] fileSet = [read1[0], sub(read1[0], inputRead1Identifier, inputRead2Identifier), sub(basename(read1[0]), inputRead1Identifier + inputExtension, "")]
    }

    Array[Array[String]] FilePaths = fileSet

    scatter (readPair in FilePaths) {
        call QcAdapters {
            input:
            sampleName = readPair[2],
            read1 = readPair[0],
            read2 = readPair[1],
            dockerImage = 
        }

        call QualityControl {
            input:
            read1 = QcAdapters.outputR1,
            read2 = QcAdapters.outputR2,
            extension = QcAdapters.extension,
            humanDB = humanDB,
            kneaddataDockerImage = kneaddataDockerImage
        }

        call Assemble {
            input:
            read1 = QualityControl.outputR1,
            read2 = QualityControl.outputR2,
            unmatched1 = QualityControl.outputUnmatchedR1,
            unmatched2 = QualityControl.outputUnmatchedR2,
            sampleName = readPair[2]
        }
    }

    call Concatenate {
        input:
        assemblyOutputs = Assemble.
    }

    scatter (readPair in FilePaths){
        call Map {
            input:
            concatCatalogue = Concatenate.concatOutput
            read1 = readPair[0],
            read2 = readPair[1],
            sample = readPair[2],
            MaxMemGB = MaxMemGB
        }
    }

    call Vamb {
        input:
        allBam = Map.mapOutput,
        concatCatalogue = Concatenate.concatOutput
    }

}

task QcAdapters {
    String sampleName
	File read1
	File read2

    String outputFile1 = sampleName + "_trimmed_1.fq.gz"
    String outputFile2 = sampleName + "_trimmed_2.fq.gz"

    # renaming required to catch result files (trim_galore has output dir option but wdl 1.0 cannot use dir as type)
	command {
        mv ${read1} ${sampleName}_1.fq.gz
        mv ${read2} ${sampleName}_2.fq.gz

		trim_galore --paired --phred33 --quality 0 --stringency 5 --length 10 \
		${sampleName}_1.fq.gz ${sampleName}_2.fq.gz

        mv ${sampleName}_1.fq.gz ${outputFile1}
        mv ${sampleName}_2.fq.gz ${outputFile2}
	}
	
	output {
		File outputR1 = "${outputFile1}"
		File outputR2 = "${outputFile2}"
        String extension = ".fq.gz"
	}
	
	runtime {
        # only needs trim_galore
	}
}

task QualityControl {
    File read1
    File read2
    String extension
    File humanDB
    String kneaddataDockerImage

    String read1Basename = basename(read1, extension)

    command <<<
		kneaddata --input ${file1} --input ${file2} -o . \
		-db {humanDB} --trimmomatic-options "HEADCROP:15 SLIDINGWINDOW:1:20 MINLEN:50" -t 4
		rm *trimmed*
		rm *bowtie2*
		
		gzip ${read1Basename}_kneaddata_paired_1.fastq
		gzip ${read1Basename}_kneaddata_paired_2.fastq
		gzip ${read1Basename}_kneaddata_unmatched_1.fastq
		gzip ${read1Basename}_kneaddata_unmatched_2.fastq
    >>>

    output {
        File outputR1 = "${read1Basename}_kneaddata_paired_1.fastq.gz"
        File outputR2 = "${read1Basename}_kneaddata_paired_2.fastq.gz"
        File outputUnmatchedR1 = "${read1Basename}_kneaddata_unmatched_1.fastq.gz"
        File outputUnmatchedR2 = "${read1Basename}_kneaddata_unmatched_2.fastq.gz"
    }
}

task Assemble {
    String sampleName
    File read1
    File read2
    File unmatched1
    File unmatched2

    command <<<
    	rm -f assemble
        spades.py --pe1-1 ${read1} --pe1-2 ${read2} --pe1-s ${unmatched1} --pe1-s ${unmatched2} -t 4 -m 16 -o assemble 
		cat assemble/contigs.fasta | \
		awk -v var="${sampleName}" '
			{if($0 ~ /^>/) {contigName=substr($0, 2,length($0))} 
			else {seq=$0; if(length($0) >= 100) {print ">"var"_"contigName"\n"seq}} }' > assemble/${sampleName}.min100.contigs.fa
    >>>

    output{
        File outputFile = "assemble/${sampleName}.min100.contigs.fa"
    }

    runtime{
        # spades
    }
}

task Concatenate {
    Array[File] assemblyOutputs


    command <<<
        contatenate.py ./concat/catalogue.fna.gz ~{sep=" " assemblyOutputs}
    >>>

    output{
        File concatOutput = "concat/catalogue.fna.gz"
    }
    runtime{
        #spades
    }
}

task Map {
    File concatCatalogue
    File read1
    File read2
    String sample
    Int? MaxMemGB
    
    Int mem = select_first([MaxMemGB, 32])
    String resultCatalogue = "catalogue.mmi"

    command <<<
        minimap2 -I ~{mem} -d ~{resultCatalogue} ~{concatCatalogue} 
        minimap2 -I ~{mem} -t 28 -N 5 -ax sr ~{resultCatalogue} ~{read1} ~{read2} | samtools view -F 3584 -b --threads 8 > ~{sample}.bam
    >>>

    output{
        File mapOutput = "${sample}.bam"
    }

    runtime {
        #minimap, samtools
    }
}

task Vamb {
    Array[File] allBam
    File concatCatalogue

    String resultDir = "vambResults"

    command <<<
        mkdir ~{resultDir}
        vamb -o C --outdir ~{resultDir} --fasta ~{concatCatalogue} --bamfiles ~{sep=" " allBam} --minfasta 200000
    >>>

    output {
        File vambTSV = "${resultDir}/clusters.tsv"
        File vambLogFile = "${resultDir}/log.txt"
    }

    runtime {
        #vamb
    }
}