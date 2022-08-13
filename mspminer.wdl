version 1.0

workflow workflowMspMiner {
    input {
        File inputRead1Files
        String inputRead1Identifier
        String inputRead2Identifier
        String inputExtension
        File humanDB
		File usearchFile

        Int? assembleMemGB
        Int? mapMemGB
		Int? cdhitCPU
		Int? cdhitMemGB
    }

    String mspminerDockerImage = "diddlydoodles/mspminer-workflow:version1.1"
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
            dockerImage = mspminerDockerImage
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
            dockerImage = mspminerDockerImage
        }

		call PredictGenes {
			input: 
            sampleName = readPair[2],
			fileContigs=Assemble.fileContigs,
			dockerImage=mspminerDockerImage
		}
    }

	call MergeGenePredictions { 
		input: 
		genePredictions=PredictGenes.fileFNA,
		geneAnnotation=PredictGenes.fileGFF,
		dockerImage=mspminerDockerImage
	}

	call CDhit { 
		input: 
		combinedgenepredictions=MergeGenePredictions.out_all,
		usearchFile=usearchFile,
		cdhitCPU=cdhitCPU,
		cdhitMemGB=cdhitMemGB,
		dockerImage=mspminerDockerImage
	}

	Array[Pair[File, File]] QCR1R2 = zip(QualityControl.outputR1, QualityControl.outputR2) 

	scatter (pair in QCR1R2){

		String currentSample = basename(pair.left, "_kneaddata_paired_1.fastq.gz")

		call MaptoCount { 
			input: 
			fileR1=pair.left,
			fileR2=pair.right,
			sample=currentSample,
			nrFa=CDhit.nrFa,
			nrFai=CDhit.nrFai,
			ref1=CDhit.nrRef1,
			ref2=CDhit.nrRef2,
			ref3=CDhit.nrRef3,
			ref4=CDhit.nrRef4,
			ref5=CDhit.nrRef5,
			genomeAnnotation=MergeGenePredictions.out_gff,
			dockerImage=mspminerDockerImage
		}
	}


	call MSP { 
		input: 
		geneCountArray=MaptoCount.fileCount,
		geneNamesArray=MaptoCount.fileGenes[0],
		dockerImage=mspminerDockerImage
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
        maxRetries: 2
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
        maxRetries: 2
    }
}

task Assemble {
	input {
        String sampleName
		File read1
		File read2
		File unmatched1
		File unmatched2
		String dockerImage
		Int? assembleMemGB
	}

    Int mem = select_first([assembleMemGB, 64])

    String resultDir = "assemble/"
    String resultFile = resultDir + sampleName + ".contigs.fasta"

	command <<<
	export PATH="/yes/bin:/yes/condabin:$PATH"
	rm -f ~{resultDir}
	megahit -1 ~{read1} -2 ~{read2} -t 4 -o ~{resultDir}
	cat assemble/final.contigs.fa | \
	awk -v var="~{sampleName}" '
		{if($0 ~ /^>/) {contigName=substr($0, 2,length($0))} 
		else {seq=$0; if(length($0) >= 100) {print ">"var"_"contigName"\n"seq}} }' > assemble/~{sampleName}.min100.contigs.fa
	mv ~{resultDir}final.contigs.fa ~{resultFile}
	>>>
	output {
		File fileContigs = "assemble/${sampleName}.min100.contigs.fa"
		File unfilteredContigs = "${resultFile}"
	}
	runtime {
		docker: dockerImage
		cpu: 4
  		memory: mem + "GB"
  		preemptible: 2
  		bootDiskSizeGb: 50
  		disks: "local-disk 100 SSD"
	}
}

task PredictGenes {
	input {
		File fileContigs
		String sampleName
		String dockerImage
	}

	String resultDir = "prodigal/"
	command <<<
	rm -f ~{resultDir}
	mkdir ~{resultDir}
	if [[ `wc -l ~{fileContigs} | awk '{print $1}'` == "0" ]]; then
		touch ~{resultDir}~{sampleName}.gff
		touch ~{resultDir}~{sampleName}.fna
		touch ~{resultDir}~{sampleName}.faa
	else
		prodigal -p meta -i ~{fileContigs} -f gff \
		-o ~{resultDir}~{sampleName}.gff \
		-d ~{resultDir}~{sampleName}.fna \
		-a ~{resultDir}~{sampleName}.faa \
		2> ~{resultDir}prodigal.stderr > ~{resultDir}prodigal.out
	fi
	>>>
	
	output {
		File fileFNA = "${resultDir}${sampleName}.fna"
		File fileFAA = "${resultDir}${sampleName}.faa"
		File fileGFF = "${resultDir}${sampleName}.gff"
	}

	runtime {
        docker: dockerImage
        cpu: 2
        memory: 8 + " GB"
        preemptible: 2
        disks: "local-disk 100 SSD"
        maxRetries: 2
	}
}

task MergeGenePredictions {
	input {
		Array[File] genePredictions
		Array[File] geneAnnotation
		String dockerImage
	}
	
	command <<<
		cat ~{sep=' ' genePredictions} > combined_genepredictions.fna
		cat ~{sep=' ' geneAnnotation} > combined_geneannotation.gff
	>>>

  	output { 
		File out_all = "combined_genepredictions.fna" 
		File out_gff = "combined_geneannotation.gff"
  	}

  	runtime {
        docker: dockerImage
        cpu: 2
        memory: 8 + " GB"
        preemptible: 2
        disks: "local-disk 100 SSD"
        maxRetries: 2
	}
}

task CDhit {
	input {
		File combinedgenepredictions
		File usearchFile
		Int? cdhitCPU
		Int? cdhitMemGB
		String dockerImage
	}

	Int cpu = select_first([cdhitCPU, 32])
    Int mem = select_first([cdhitMemGB, 64])
	Int memMB = mem * 1000

	command {
	chmod a+rx ${usearchFile}
	${usearchFile} -derep_fulllength ${combinedgenepredictions} -fastaout combinedgenepredictions.derep.fna -minseqlength 1
	${usearchFile} -sortbylength combinedgenepredictions.derep.fna -fastaout combinedgenepredictions.sorted.fna -minseqlength 1
	cd-hit-est -i combinedgenepredictions.sorted.fna -T ${cpu} -aS 0.9 -c 0.95 -M ${memMB} -r 0 -B 0 -d 0 -o nr.fa
	bwa index nr.fa
	samtools faidx nr.fa
	}

  	output { 
  		File nrFa = "nr.fa"
  		File nrFai = "nr.fa.fai"
  		File nrRef1 = "nr.fa.bwt"
  		File nrRef2 = "nr.fa.pac"
  		File nrRef3 = "nr.fa.ann"
  		File nrRef4 = "nr.fa.amb"
  		File nrRef5 = "nr.fa.sa"
        File clusterFile = "nr.fa.clstr"
  	}

  	runtime {
        docker: dockerImage
        cpu: cpu
        memory: mem + " GB"
        preemptible: 2
        disks: "local-disk 500 SSD"
        maxRetries: 2
	}
}

task MaptoCount {
	input {
		File fileR1
		File fileR2
		String sample
		File nrFa
		File nrFai
		File ref1
		File ref2
		File ref3
		File ref4
		File ref5
		File genomeAnnotation
		String dockerImage
	}

	command {
	wget https://raw.githubusercontent.com/lpryszcz/bin/master/bam2counts.py -P /msp
	wget https://raw.githubusercontent.com/lpryszcz/bin/master/python_modules/genome_annotation.py -P /msp
	sed -i "s/commands/subprocess/" /msp/bam2counts.py
	bwa mem -t 8 -M ${nrFa} ${fileR1} ${fileR2} | \
	samtools view -h -Su -F 2308 -q 0 | \
	samtools sort -n -@ 8 -m 2G -O bam -o ${sample}.sort.bam 

	pip install numpy
	python3 /msp/bam2counts.py -v -i ${sample}.sort.bam -g ${genomeAnnotation} > ${sample}.count.txt

	cut -f1 ${sample}.count.txt > gene.names.txt
	cut -f2 ${sample}.count.txt > tmp.count.txt
	
	mv tmp.count.txt ${sample}.count.txt
	}
	
	output {
		File fileBAM = "${sample}.sort.bam"
		File fileCount = "${sample}.count.txt"
		File fileGenes = "gene.names.txt"
	}
	
	runtime {
		docker: dockerImage
		cpu: 8
  		memory: "32GB"
  		preemptible: 2
  		bootDiskSizeGb: 100
  		disks: "local-disk 500 HDD"
        maxRetries: 2
	}
}

task MSP {
	input {
		Array[File] geneCountArray
		File geneNamesArray
		String dockerImage
	}

	String countFile='combined_geneCount.txt'
	String settingsLocation='/msp/mspminer_bin/settings.ini'

	command <<<
	
	mkdir msps
	
	paste ~{geneNamesArray} ~{sep=' ' geneCountArray} > ~{countFile}

	sed -i "s/^count_matrix_file=.*/count_matrix_file=~{countFile}/" ~{settingsLocation}
	sed -i "s/^output_dir=.*/output_dir=msps/" ~{settingsLocation}

	mspminer ~{settingsLocation}
	tar -czf msps.tar.gz msps
	
	>>>
	
	output {
		File fileCountMatrix = "${countFile}"
		File fileMSPs = "msps.tar.gz"	
	}
	
	runtime {
		docker: dockerImage
		cpu: 4
  		memory: "16GB"
  		preemptible: 2
  		disks: "local-disk 100 HDD"
        maxRetries: 2
	}
}


