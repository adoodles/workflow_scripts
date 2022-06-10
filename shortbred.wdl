version 1.0

workflow workflowShortbred{
    input {
        Boolean skipIdentify
        File? interestProteins
        File? referenceProteins
        File nucleotideReads
        File? markerFile
        String? tmpResultFolder
        String quantifyOutput
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

    call Quantify {
        input: 
        markerFile = select_first([markerFile, Identify.outputMarkerFile]),
        quantifyOutput = quantifyOutput,
        shortBredDockerImage = shortBredDockerImage
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
        String quantifyOutput
        String shortBredDockerImage
    }

    String resultExtension = if ( quantifyOutput == "wgs") then ".fna" else ".faa"
    command {
        ./shortbred_quantify.py --markers ${markerFile} --${quantifyOutput} example/wgs${resultExtension}
    }

    runtime {
        docker: shortBredDockerImage
        cpu: 8
        memory: "4" + " GB"
        disks: "local-disk 501 SSD"
    }
}