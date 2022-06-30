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
    String metaphlanDockerImage = "biobakery/metaphlan:3.0.1"

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

    call Collect {
        input:
        quantifyResults = Quantify.outputFile,
        dockerImage = metaphlanDockerImage
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
        EXTRACTEDFILE=${FILENAME}
        if [ ~{compressionExtension} == ".gz" ]; then
            gzip -d ~{inputFile}
            EXTRACTEDFILE=${FILENAME%.gz}
        elif [ ~{compressionExtension} == ".tar.gz" ]; then
            tar -xzf ~{inputFile}
            EXTRACTEDFILE=${FILENAME%.tar.gz}
        fi
        shortbred_quantify.py --markers ~{markerFile} --~{quantifyInputMethod} ${EXTRACTEDFILE} --results quantify/~{resultName}
    >>>

    output {
        File outputFile = "quantify/${resultName}"
    }

    runtime {
        docker: shortBredDockerImage
        cpu: 8
        memory: "4" + " GB"
        preemptible: 2
        disks: "local-disk 100 SSD"
    }
}

task Collect {
    input {
        Array[File] quantifyResults
        String dockerImage
    }

    String resultFilePath = "mergedResults.txt"

    # **check if combineAndWrite's sep delimiter is same as input for listOfFiles.split()**
    command <<<
        python3<<CODE
        import pandas as pd
        import os.path

        def combineAndWrite(listOfFiles):
            mergedTables = pd.DataFrame()
            fileList = listOfFiles.split(" ")
            for f in fileList:
                colName = os.path.basename(f)
                currentTable = pd.read_csv(f, sep = '\t', usecols = ['Family', 'Count'])
                currentTable.rename(columns={'Count' : colName}, inplace=True)
                if mergedTables.empty :
                    mergedTables = currentTable
                else:
                    mergedTables = pd.merge(mergedTables, currentTable, how = 'outer', on = ['Family'])
            mergedTables.to_csv("~{resultFilePath}", sep = '\t', index=False)

        combineAndWrite("~{sep=" " quantifyResults}")
        CODE
    >>>

    output {
        File mergedOutput = resultFilePath
    }

    runtime {
        docker: dockerImage
        cpu: 4
        memory: "8" + " GB"
        preemptible: 2
        disks: "local-disk 50 SSD"
    }
}