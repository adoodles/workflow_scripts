version 1.0

workflow workflowShortbred{
    input {
        Boolean skipIdentify
        File? interestProteins
        File? referenceProteins
        File? markerFile
        File inputReadFiles
        String inputExtension
        # only .gz or .tar.gz or none
        String compressionExtension
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
    Array[Array[String]] inputRead = read_tsv(inputReadFiles)
    
    scatter (read in inputRead) {
        Array[String] fileSet = [read[0]]
    }

    Array[Array[String]] FilePaths = fileSet

    scatter (Paths in FilePaths){
        call Quantify {
            input: 
            markerFile = select_first([markerFile, Identify.outputMarkerFile]),
            quantifyInputMethod = quantifyInputMethod,
            inputFile = Paths[0],
            inputExtension = inputExtension,
            compressionExtension = compressionExtension,
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
        String inputExtension
        String compressionExtension
        String shortBredDockerImage
    }

    String fileName = basename(inputFile)
    String resultName = sub(fileName, inputExtension, "")

    command <<<
        mkdir -p quantify
        FILENAME=~{inputFile}
        if [ ~{compressionExtension} == ".gz" ]; then
            gzip -d ~{inputFile}
        elif [ ~{compressionExtension} == ".tar.gz"]; then
            tar -xzf ~{inputFile}
        fi
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