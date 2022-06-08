version 1.0

workflow workflowShortbred{
    input {
        File interestProteins
        File referenceProteins
        File nucleotideReads
        String? markerFileName
        String? tmpResultFolder
        String quantifyOutput
    }

    # Set the docker tags
    String shortBredDockerImage = "biobakery/shortbred:0.9.5"

    # Identify task variables
    String outputMarkerFile = select_first([markerFileName, "identifyMarkers.faa"])

    call Identify {
        input:
        interestProteins = interestProteins,
        referenceProteins = referenceProteins,
        outputMarkerFile = outputMarkerFile,
        tmp = tmpResultFolder,
        shortBredDockerImage = shortBredDockerImage
    }

    call Quantify {
        input: 
        File markerFile = Identify.outputMarkerFile,
        String quantifyOutput = quantifyOutput,
        shortBredDockerImage = shortBredDockerImage
    }
}

task Identify {
    input {
        File interestProteins
        File referenceProteins
        String outputMarkerFile
        String shortBredDockerImage
        String? tmp
    }

    String tmpFolder = if (defined(tmp)) then "--tmp ${tmp}" else ""

    command {
        ./shortbred_identify.py --goi  ${interestProteins} --ref ${referenceProteins} --markers ${outputMarkerFile} ${tmpFolder}
    }

    output {
        File outputMarkerFile = outputMarkerFile
    }

    runtime {
        docker: metaWibeleDockerImage
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
        docker: metaWibeleDockerImage
        cpu: 8
        memory: "4" + " GB"
        disks: "local-disk 501 SSD"
    }
}