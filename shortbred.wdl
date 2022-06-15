version 1.0

workflow workflowShortbred{
    input {
        Boolean skipIdentify
        File? interestProteins
        File? referenceProteins
        File? markerFile
        File inputRead1Files
        String? inputRead1Identifier
        String? inputRead2Identifier
        String inputExtension
        String? tmpResultFolder
        # choice of wgs or genome
        String quantifyInputMethod
    }

    # Set the docker tags
    String shortBredDockerImage = "biobakery/shortbred:0.9.5"

    # Identify task variables
    String identifyMarkerFile = "identifyMarkers.faa"
    
    if (!skipIdentify) {
        call Identify {
            input:
            interestProteins = interestProteins,
            referenceProteins = referenceProteins,
            outputMarkerFile = identifyMarkerFile,
            tmp = tmpResultFolder,
            shortBredDockerImage = shortBredDockerImage
        }
    }

    # read in a file of the read1 paths
    Array[Array[String]] inputRead1 = read_tsv(inputRead1Files)
    
    # get the sample name and read2 file path
    # wdl 1.0 does not allow "" to be File value. Placeholder value read1[0] is used here as "invalid" value temporarily until 1.x allows for similar syntax usage
    scatter (read1 in inputRead1) {
        Array[String] pairSet = [read1[0], if defined(inputRead2Identifier) then sub(read1[0], inputRead1Identifier + "", inputRead2Identifier + "") else read1[0]]
    }

    Array[Array[String]] PairPaths = pairSet

    scatter (Paths in PairPaths){
        call Quantify {
            input: 
            markerFile = select_first([markerFile, Identify.outputMarkerFile]),
            quantifyInputMethod = quantifyInputMethod,
            inputFile = Paths[0],
            inputPairedFile = Paths[1],
            inputExtension = inputExtension,
            shortBredDockerImage = shortBredDockerImage
        }
    }
}

task Identify {
    input {
        File? interestProteins
        File? referenceProteins
        String outputMarkerFile
        String shortBredDockerImage
        String? tmp
    }

    String goiInput = if defined(interestProteins) then 'yes' else 'no'
    String refInput = if defined(interestProteins) then 'yes' else 'no'
    String tmpFolder = if (defined(tmp)) then "--tmp ${tmp}" else ""

    command {
        if [ ${goiInput} == 'no' || ${refInput} == 'no']; then
            echo "interestProteins or referenceProteins not specified"
            exit 1
        fi
        ./shortbred_identify.py --goi  ${interestProteins} --ref ${referenceProteins} --markers ${outputMarkerFile} ${tmpFolder}
    }

    output {
        File outputMarkerFile = outputMarkerFile
    }

    runtime {
        docker: shortBredDockerImage
        cpu: 8
        memory: "4" + " GB"
        disks: "local-disk 501 SSD"
    }
}

task Quantify {
    input {
        File markerFile
        String quantifyInputMethod
        File inputFile
        File inputPairedFile
        String inputExtension
        String shortBredDockerImage
    }

    String compressionExtension = ".gz"
    # call response to "invalid" value read1[0] from pairSet creation
    String pairedFile = if inputFile == inputPairedFile then "" else inputPairedFile
    String fileName = basename(inputFile)
    String resultName = sub(fileName, inputExtension, "")
    String inputFileArg = sub(fileName, compressionExtension, "")

    command <<<
        mkdir -p quantify
        FILENAME=~{inputFile}
        gzip -d ~{inputFile}
        shortbred_quantify.py --markers ~{markerFile} --~{quantifyInputMethod} ${FILENAME%.*} --results quantify/~{resultName}.txt
    >>>

    output {
        File outputFile = "quantify/${resultName}.txt"
    }

    runtime {
        docker: shortBredDockerImage
        cpu: 8
        memory: "4" + " GB"
        disks: "local-disk 501 SSD"
    }
}