version 1.0

workflow workflowUnicycler {
	input {
		File inputRead1Files
		String inputRead1Identifier
		String inputRead2Identifier
    String longReadIdentifier
    String projectName
	}
  # Set the docker tags
  String unicyclerDockerImage = "biobakery/kneaddata:0.10.0"

  Array[Array[String]] inputRead1 = read_tsv(inputRead1Files)
  
  # get the sample name and read2 file path
  scatter (read1 in inputRead1) {
     Array[String] pairSet = [read1[0], sub(read1[0], inputRead1Identifier, inputRead2Identifier), sub(read1[0], inputRead1Identifier, longReadIdentifier)]
  }

  Array[Array[String]] PairPaths = pairSet

  # mem settings
  Int UnicyclerMemBase = 8

	scatter (readPair in PairPaths) {
		call Unicycler {
			input:
        read1=readPair[0],
        read2=readPair[1],
        longread=readPair[2],
        MemBase=UnicyclerMemBase,
        dockerImage=unicyclerDockerImage
		}
	}
}

task Unicycler {
  input {
    File read1
    File read2
    File longread
    Int MemBase
    String dockerImage
  }
  String outputDir = "results/"
  Int mem = MemBase
  command {
    mkdir ${outputDir}
    unicycler -1 ${read1} -2 ${read2} -l ${longread} -o ./${outputDir}
  }
  output{
    Array[File] outputs = glob(outputDir + "*")
  }
  runtime{
    docker: dockerImage
    cpu: 8
    memory: mem + " GB"
    preemptible: 2
    disks: "local-disk 501 SSD"
  }
}