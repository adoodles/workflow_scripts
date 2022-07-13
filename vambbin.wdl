version 1.0

workflow workflowVamb {
    input {
        File inputRead1Files
        String inputRead1Identifier
        String inputRead2Identifier
        String inputExtension
        File humanDB

        Int? concatMinLength
        Int? vambMinFasta
        Int? assembleMemGB
        Int? mapMemGB
    }

    String vambDockerImage = "diddlydoodles/vamb-workflow:latest"
    String kneaddataDockerImage = "biobakery/kneaddata:0.10.0"

    # read in a file of the read1 paths
    Array[Array[String]] inputRead = read_tsv(inputRead1Files)
    
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
            dockerImage = vambDockerImage
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
            sampleName = readPair[2],
            assembleMemGB = assembleMemGB,
            vambDockerImage = vambDockerImage
        }
    }

    call Concatenate {
        input:
        assemblyOutputs = Assemble.outputFile,
        vambDockerImage = vambDockerImage,
        minLength = concatMinLength
    }

    scatter (readPair in FilePaths){
        call MaptoBam {
            input:
            concatCatalogue = Concatenate.concatOutput,
            read1 = readPair[0],
            read2 = readPair[1],
            sample = readPair[2],
            mapMemGB = mapMemGB,
            vambDockerImage = vambDockerImage
        }
    }

    call Vamb {
        input:
        allBam = MaptoBam.mapOutput,
        concatCatalogue = Concatenate.concatOutput,
        vambDockerImage = vambDockerImage,
        vambMinFasta = vambMinFasta
    }

}

task QcAdapters {
    input {
        String sampleName
        File read1
        File read2
        String dockerImage
    }

    String outputFile1 = sampleName + "_adapterTrim_1.fq.gz"
    String outputFile2 = sampleName + "_adapterTrim_2.fq.gz"

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
        docker: dockerImage
        cpu: 2
        memory: 4 + " GB"
        preemptible: 2
        disks: "local-disk 20 SSD"
    }
}

task QualityControl {
    input {
        File read1
        File read2
        String extension
        File humanDB
        String kneaddataDockerImage
    }
    String humanDatabase = "databases/kneaddata_human/"
    String read1Basename = basename(read1, extension)

    command <<<
    mkdir -p ~{humanDatabase}
    kneaddata_database --download human_genome bowtie2 ~{humanDatabase} --database-location ~{humanDB}

    kneaddata --input ~{read1} --input ~{read2} -o . \
    -db ~{humanDatabase} --trimmomatic-options "HEADCROP:15 SLIDINGWINDOW:1:20 MINLEN:50" -t 4
    rm *trimmed*
    rm *bowtie2*   

    gzip ~{read1Basename}_kneaddata_paired_1.fastq
    gzip ~{read1Basename}_kneaddata_paired_2.fastq
    gzip ~{read1Basename}_kneaddata_unmatched_1.fastq
    gzip ~{read1Basename}_kneaddata_unmatched_2.fastq
    >>>

    output {
        File outputR1 = "${read1Basename}_kneaddata_paired_1.fastq.gz"
        File outputR2 = "${read1Basename}_kneaddata_paired_2.fastq.gz"
        File outputUnmatchedR1 = "${read1Basename}_kneaddata_unmatched_1.fastq.gz"
        File outputUnmatchedR2 = "${read1Basename}_kneaddata_unmatched_2.fastq.gz"
    }

    runtime {
        docker: kneaddataDockerImage
        cpu: 8
        memory: 24 + " GB"
        preemptible: 2
        disks: "local-disk 501 SSD"
    }
}

task Assemble {
    input {
        String sampleName
        File read1
        File read2
        File unmatched1
        File unmatched2
        String vambDockerImage
        Int? assembleMemGB
    }

    Int mem = select_first([assembleMemGB, 64])

    String resultDir = "assemble/"
    String resultFile = resultDir + sampleName + ".contigs.fasta"
    command <<<
    rm -f ~{resultDir}
    spades.py --pe1-1 ~{read1} --pe1-2 ~{read2} --pe1-s ~{unmatched1} --pe1-s ~{unmatched2} -t 4 -m ~{mem} -o assemble 
    mv ~{resultDir}contigs.fasta ~{resultFile}
    >>>

    output{
        File outputFile = "${resultFile}"
    }

	runtime {
		docker: vambDockerImage
		cpu: 4
  		memory: mem + "GB"
  		preemptible: 2
  		bootDiskSizeGb: 50
  		disks: "local-disk 100 SSD"
	}
}

task Concatenate {
    input {
        Array[File] assemblyOutputs
        String vambDockerImage
        Int? minLength
    }

    String resultDirectory = "concat/"
    String catalogueFile = resultDirectory + "catalogue.fna.gz"
    Int minimumLength = select_first([minLength, 2000])

    command <<<
    mkdir -p ~{resultDirectory}
    concatenate.py ~{catalogueFile} ~{sep=" " assemblyOutputs} -m ~{minimumLength}
    >>>

    output{
        File concatOutput = "${catalogueFile}"
    }

	runtime {
		docker: vambDockerImage
		cpu: 8
  		memory: "32GB"
  		preemptible: 2
  		disks: "local-disk 500 SSD"
	}
}

task MaptoBam {
    input {
        File concatCatalogue
        File read1
        File read2
        String sample
        String vambDockerImage
        Int? mapMemGB
    }
    
    Int mem = select_first([mapMemGB, 64])
    String resultCatalogue = "catalogue.mmi"

    command <<<
        minimap2 -I ~{mem}g -d ~{resultCatalogue} ~{concatCatalogue} 
        minimap2 -I ~{mem}g -t 28 -N 5 -a ~{resultCatalogue} ~{read1} ~{read2} | samtools view -F 3584 -b --threads 8 -o ~{sample}.bam
    >>>

    output{
        File mapOutput = "${sample}.bam"
    }

	runtime {
		docker: vambDockerImage
		cpu: 8
  		memory: mem + "GB"
  		preemptible: 2
  		disks: "local-disk 120 SSD"
	}
}

task Vamb {
    input {
        Array[File] allBam
        File concatCatalogue
        String vambDockerImage
        Int? vambMinFasta
    }

    String resultDir = "vambResults/"
    Int minFasta = select_first([vambMinFasta, 200000])

    command <<<
        rm -f ~{resultDir}
        vamb -o C --outdir ~{resultDir} --fasta ~{concatCatalogue} --bamfiles ~{sep=" " allBam} --minfasta ~{minFasta}
    >>>

    output {
        File vambTSV = "${resultDir}clusters.tsv"
        File vambLogFile = "${resultDir}log.txt"
    }

	runtime {
		docker: vambDockerImage
		cpu: 8
  		memory: "24GB"
  		preemptible: 2
  		disks: "local-disk 500 SSD"
	}
}
