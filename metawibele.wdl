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
        File uniref90DiamondDatabase
    }

    # Set the docker tags
    String metaWibeleDockerImage = "biobakery/metawibele:0.4.4"

    # Set output file or folder names
    String preprocessOutputDir = "preprocess/"
    String characterizeOutputDir = "characterize/"
    String prioritizeOutputDir = "prioritize/"

    # Basename of project set in .cfg file default value: metawibele
    String basename = "metawibele"

    if (!skipPreprocess) {
        call Preprocess {
            input:
            fastqFolder = fastqDir,
            extensionPaired = extensionPaired,
            extension = extension,
            outputDir = outputDir
        }
    }

    call Characterize {
        input:
        inputSequence = if (!skipPreprocess) then Preprocess.outputSequence else inputSequence,
        inputCount = if (!skipPreprocess) then Preprocess.outputCount else inputCount,
        metadata = metadata,
        outputDir = outputDir
    }

    if (!skipPrioritize) {
        call Prioritize {
            input:
            inputAnnotation = Characterize.outputAnnotation,
            inputAttribute = Characterize.outputAttribute,
            outputDir = outputDir
        }
    }

    task Preprocess {
        input {
            Directory fastqFiles
            String extensionPaired
            String extension
        }

        if (!defined(fastqFiles) || !defined(extensionPaired) || !defined(extension)) {
            command {
                echo 'Input file error preprocess arguments not provided'
                exit 1
            }
        }

        command {
            metawibele preprocess --input ${fastqFiles} --extension-paired ${extensionPaired} --extension ${extension} --output ${preprocessOutputDir}
        }

        output {
            File outputSequence = preprocessOutputDir + "finalized/" + basename + "_genecatalogs.centroid.faa"
            File outputCount = preprocessOutputDir + "finalized/" + basename + "_genecatalogs_counts.all.tsv"
        }

        runtime {
            docker: metaWibeleDockerImage
            cpu: 8
            memory: mem + " GB"
            preemptible: preemptible_attempts
            disks: "local-disk 501 SSD"
        }
    }

    task Characterize {
        input {
            File inputSequence
            File inputCount
            File metadata
            File uniref90DiamondDatabase
        }

        String databaseLocation = "/databases/"
        #install diamond database and run characterize
        command {
            mkdir -p ${databaseLocation}
            metawibele_download_config --config-type global
            sed -i.bak 's/uniref_db =/uniref_db = $(pwd)${databaseLocation}/' metawibele.cfg
            tar -xzf ${uniref90DiamondDatabase} $(pwd)${databaseLocation}

            metawibele characterize --input-sequence inputSequence --input-count inputCount --input-metadata metadata --output ${characterizeOutputDir}
        }

        output{
            File outputAnnotation = characterizeOutputDir + "finalized/" + basename + "_proteinfamilies_annotation.tsv"
            File outputAttribute = characterizeOutputDir + "finalized/" + basename + "_proteinfamilies_annotation.attribute.tsv"
        }
    }
}