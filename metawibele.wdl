workflow workflowMetaWibele {
    input {
        Boolean skipPreprocess
        Boolean skipPrioritize
        Directory? fastqDir
        String? extensionPaired
        String? extension
        File? inputSequence
        File? inputCount
        File inputMetadata
    }

    # Set the docker tags
    String metaWibeleDockerImage = "biobakery/metawibele:0.4.4"

    call Preprocess {
        input:
        fastqFolder = fastqDir,
        extensionPaired = extensionPaired,
        extension = extension,
        outputDir = outputDir
    }

    call Categorize {
        input:
        inputSequence = if (!skipPreprocess) then Preprocess.outputSequence else inputSequence,
        inputCount = if (!skipPreprocess) then Preprocess.outputCount else inputCount,
        metadata = metadata,
        outputDir = outputDir
    }

    call Prioritize {
        input:
        inputAnnotation = Categorize.outputAnnotation,
        inputAttribute = Categorize.outputAttribute,
        outputDir = outputDir
    }

    task Preprocess {
        input {
            Directory? fastqFiles
            String? extensionPaired
            String? extension
        }

        
        command {

        }

        output {

        }

        runtime {

        }
    }
}