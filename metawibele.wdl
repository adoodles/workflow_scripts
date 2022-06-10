version 1.0

workflow workflowMetaWibele {   
    input {
        Boolean skipPreprocess
        Boolean skipPrioritize
        Boolean skipAnnotation
        Boolean? bypass
        # compressed directory of fastq files
        File? fastqFiles
        String? extensionPaired
        String? extension
        File? inputSequence
        File? inputCount
        File inputMetadata
        File? uniref90DiamondDatabase
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
            fastqFiles = fastqFiles,
            extensionPaired = extensionPaired,
            extension = extension,
            basename = basename,
            preprocessOutputDir = preprocessOutputDir,
            metaWibeleDockerImage = metaWibeleDockerImage
        }
    }
    if (skipAnnotation) {
        call Characterize {
            input:
            inputSequence = if (!skipPreprocess) then Preprocess.outputSequence else inputSequence,
            inputCount = if (!skipPreprocess) then Preprocess.outputCount else inputCount,
            metadata = inputMetadata,
            basename = basename,
            characterizeOutputDir = characterizeOutputDir,
            metaWibeleDockerImage = metaWibeleDockerImage
        }
    }
    if (!skipAnnotation) {
        call CharacterizeWithAnnotation {
            input:
            inputSequence = if (!skipPreprocess) then Preprocess.outputSequence else inputSequence,
            inputCount = if (!skipPreprocess) then Preprocess.outputCount else inputCount,
            metadata = inputMetadata,
            basename = basename,
            characterizeOutputDir = characterizeOutputDir,
            metaWibeleDockerImage = metaWibeleDockerImage
        }
        if (!skipPrioritize) {
            call Prioritize {
                input:
                inputAnnotation = CharacterizeWithAnnotation.outputAnnotation,
                inputAttribute = CharacterizeWithAnnotation.outputAttribute,
                basename = basename,
                prioritizeOutputDir = prioritizeOutputDir,
                metaWibeleDockerImage = metaWibeleDockerImage
            }
        }
    }
}

task Preprocess {
    input {
        # compressed directory .tar.gz
        File? fastqFiles
        String? extensionPaired
        String? extension
        String basename
        String preprocessOutputDir
        String metaWibeleDockerImage
    }

    String extractDirectory = "preprocessInput/"

    command {
        mkdir -p ${preprocessOutputDir}
        tar -xzf ${fastqFiles} $(pwd)${extractDirectory}
        metawibele preprocess --input ${extractDirectory} --extension-paired ${extensionPaired} --extension ${extension} --output ${preprocessOutputDir}
    }

    output {
        File outputSequence = preprocessOutputDir + "finalized/" + basename + "_genecatalogs.centroid.faa"
        File outputCount = preprocessOutputDir + "finalized/" + basename + "_genecatalogs_counts.all.tsv"
    }

    runtime {
        docker: metaWibeleDockerImage
        cpu: 1
        memory: "4" + " GB"
        disks: "local-disk 501 SSD"
    }
}

task Characterize {
    input {
        File? inputSequence
        File? inputCount
        File metadata
        String basename
        String characterizeOutputDir
        String metaWibeleDockerImage
    }

    String skipAnnotation = "--bypass-global-homology --bypass-domain-motif --bypass-abundance"
    # skip annotation shouldn't need database
    command {
        metawibele_download_config --config-type global
        mkdir -p ${characterizeOutputDir}
        metawibele characterize --input-sequence inputSequence --input-count inputCount --input-metadata metadata --output ${characterizeOutputDir} ${skipAnnotation}
    }

    output{
        File outputAttribute = characterizeOutputDir + "finalized/" + basename + "_proteinfamilies_annotation.attribute.tsv"
    }

    runtime {
        docker: metaWibeleDockerImage
        cpu: 8
        memory: "4" + " GB"
        disks: "local-disk 501 SSD"
    }
}

# todo: half finished
task CharacterizeWithAnnotation {
    input {
        File? inputSequence
        File? inputCount
        File metadata
        File? uniref90DiamondDatabase
        String basename
        String characterizeOutputDir
        String metaWibeleDockerImage
    }

    String databaseLocation = "/databases/"
    #install diamond database and run characterize
    # download .cfg file and edit uniref_db value
    command {
        mkdir -p ${databaseLocation}
        metawibele_download_config --config-type global
        sed -i.bak 's/uniref_db =/uniref_db = $(pwd)${databaseLocation}/' metawibele.cfg
        tar -xzf ${uniref90DiamondDatabase} $(pwd)${databaseLocation}

        mkdir -p ${characterizeOutputDir}
        metawibele characterize --input-sequence inputSequence --input-count inputCount --input-metadata metadata --output ${characterizeOutputDir}
    }

    output{
        File outputAnnotation = characterizeOutputDir + "finalized/" + basename + "_proteinfamilies_annotation.tsv"
        File outputAttribute = characterizeOutputDir + "finalized/" + basename + "_proteinfamilies_annotation.attribute.tsv"
    }

    runtime {
        docker: metaWibeleDockerImage
        cpu: 8
        memory: "4" + " GB"
        disks: "local-disk 501 SSD"
    }
}

task Prioritize {
    input {
        File inputAnnotation
        File inputAttribute
        String basename
        String prioritizeOutputDir
        String metaWibeleDockerImage
    }
    
    command{
        mkdir -p ${prioritizeOutputDir}
        metawibele prioritize --input-annotation inputAnnotation --input-attribute inputAttribute --output ${prioritizeOutputDir}
    }

    runtime {
        docker: metaWibeleDockerImage
        cpu: 8
        memory: "4" + " GB"
        disks: "local-disk 501 SSD"
    }
}
