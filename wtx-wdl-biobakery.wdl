version 1.0

workflow workflowMTX {

  # Required input variables
  input {
    File inputRead1Files
    String inputRead1Identifier
    String inputRead2Identifier
    String inputExtension
    String ProjectName
    String AdapterType
    
    # Database locations
    File versionSpecifichumanDB
    File versionSpecifictrancriptDB
    File versionSpecificrrnaDB
    File versionSpecificChocophlan
    File versionSpecificUniRef90
    File versionSpecificUtilityMapping
    File versionSpecificMetaphlanDB
    File StrainphlanReferences
    
    # Optional input variables
    Boolean? bypassFunctionalProfiling
    Boolean? bypassStrainProfiling
    Int? MaxStrains
    File? CustomStrainList
    String? dataType
    File? inputMetadataFile
    File? customUtilityMappingEC
    String? customMappingECPath
    File? customUtilityMappingKO
    String? customMappingKOPath
    Int? preemptibleAttemptsOverride
    Int? MaxMemGB_QualityControlTasks
    Int? MaxMemGB_TaxonomicProfileTasks
    Int? MaxDiskGB_TaxonomicProfileTasks
    Int? MaxMemGB_FunctionalProfileTasks
    Int? MaxDiskGB_FunctionalProfileTasks
    
    File? customQCDB1
    File? customQCDB2
    File? customQCDB3
  }
  
  # Set the docker tags
  String kneaddataDockerImage = "biobakery/kneaddata:0.10.0"
  String metaphlanDockerImage = "diddlydoodles/metaphlan:4.0.6"
  String humannDockerImage = "diddlydoodles/humann:3.8"
  
  # for the workflows script tasks use the workflows image without the databases (to reduce data download amount)
  String workflowsDockerImage = "biobakery/workflows:3.0.0.a.6_anadama0.7.9_no_metaphlan_db"
  
  # for the strainphlan tasks use the workflows image with the databases (to include the mapping of references)
  String strainphlanDockerImage = "biobakery/workflows:3.0.0.a.6"
  
  # Set the default data type
  String dataTypeSetting = select_first([dataType, "mtx"])
  
  # Set if metadata is provided
  String metadataSet = if defined(inputMetadataFile) then "yes" else "no"
  
  # Set bypass mode
  Boolean setbypassFunctionalProfiling = select_first([bypassFunctionalProfiling,false])
  Boolean setbypassStrainProfiling = select_first([bypassStrainProfiling,false])  
  
  Int setMaxStrains = select_first([MaxStrains,10])

  # Output file names to match AnADAMA2 workflow output names
  String QCReadCountFileName = "kneaddata_read_count_table.tsv"  

  String JoinedTaxonomicProfilesFileName="metaphlan_taxonomic_profiles.tsv"
  String TaxonomicProfilesCountsFileName="metaphlan_species_counts_table.tsv"
  String JoinGeneFamilesOutFileName="genefamilies.tsv"
  String JoinECsOutFileName="ecs.tsv"
  String JoinKOsOutFileName="kos.tsv"
  String JoinRXNsOutFileName="rxns.tsv"
  String JoinPathwaysOutFileName="pathabundance.tsv"

  String JoinedFeatureCountsFileName="humann_feature_counts.tsv"
  String FunctionalCountFileName = "humann_read_and_species_count_table.tsv"
  
  String StrainPhlAnCladeList = "strainphlan_clade_list.txt"
    
  # mem settings
  # Quality Control Step: Assumes linear space complexity O(n)
  Int QCMemBase = 24
  Int QCMemBaseFileSize = 4
  Int QCMemIncreaseInterval = 8

  # Taxnomic Profile Step: metaphlan4 requires minimum disk size of 15GB
  Int TaxProfileMemBase = 32
  Int TaxProfileDiskBase = 50

  # Functional Profile Step:
  Int FuncProfileMemBase = 32
  Int FuncProfileBaseFileSize = 3
  Int FuncProfileIncreaseInterval = 8
  Int FuncProfileStorage = 500

  Int JoinNormMemDefault = 10
  Int JoinNormMemDefaultGenes = 50

  # read in a file of the read1 paths
  Array[Array[String]] inputRead1 = read_tsv(inputRead1Files)
  
  # get the sample name and read2 file path
  scatter (read1 in inputRead1) {
     Array[String] pairSet = [read1[0], sub(read1[0], inputRead1Identifier, inputRead2Identifier), sub(basename(read1[0]), inputRead1Identifier + inputExtension, "")]
  }

  Array[Array[String]] PairPaths = pairSet

  # add filepaths for all metaphlan db inputs (2 md5 & 2 tar files)
  String baseMetaphlanDB = versionSpecificMetaphlanDB
  String baseMetaphlanMD5 = sub(baseMetaphlanDB, ".tar", ".md5")
  String metaphlanBowtie = sub(baseMetaphlanDB, ".tar", "_bt2.tar")
  String metaphlanBowtieMD5 = sub(baseMetaphlanDB, ".tar", "_bt2.md5")

  # run tasks for each set of read pairs    
  scatter (ReadPair in PairPaths)
  {
    # Part 1: For each sample, run quality control with KneadData
    call QualityControl {
      input:
      rawfile1=ReadPair[0],
      rawfile2=ReadPair[1],
      sample=ReadPair[2],
      adapterType=AdapterType,
      humanDB=versionSpecifichumanDB,
      transcriptDB=versionSpecifictrancriptDB,
      rrnaDB=versionSpecificrrnaDB,
      customDB1=customQCDB1,
      customDB2=customQCDB2,
      customDB3=customQCDB3,
      dataType=dataTypeSetting,
      kneaddataDockerImage=kneaddataDockerImage,
      preemptibleAttemptsOverride=preemptibleAttemptsOverride,
      MaxMemGB=MaxMemGB_QualityControlTasks,
      QCMemBase = QCMemBase,
      QCMemBaseFileSize = QCMemBaseFileSize,
      QCMemIncreaseInterval = QCMemIncreaseInterval
    }
    # Part 2: For each sample, run taxonomic profiling with MetaPhlAn v2
    call TaxonomicProfile {
      input:
      sample=ReadPair[2], 
      QCFastqFile=QualityControl.QCFastqFile,
      chocophlanDB=baseMetaphlanDB,
      chocophlanDBMD5=baseMetaphlanMD5,
      chocophlanIndex=metaphlanBowtie,
      chocophlanIndexMD5=metaphlanBowtieMD5,
      metaphlanDockerImage=metaphlanDockerImage,
      preemptibleAttemptsOverride=preemptibleAttemptsOverride,
      MaxMemGB=MaxMemGB_TaxonomicProfileTasks,
      MaxDiskGB=MaxDiskGB_TaxonomicProfileTasks,
      TaxProfileMemBase = TaxProfileMemBase,
      TaxProfileDiskBase = TaxProfileDiskBase
    }
   }
   
   if (! setbypassFunctionalProfiling ) {
    # Part 3: For each sample, run functional profiling with HUMAnN v3
    scatter (sample_index in range(length(PairPaths))) {
      call FunctionalProfile {
        input: 
        sample=PairPaths[sample_index][2], 
        QCFastqFile=QualityControl.QCFastqFile[sample_index], 
        TaxonomicProfileFile=TaxonomicProfile.TaxonomicProfileFile[sample_index],
        versionSpecificChocophlan=versionSpecificChocophlan,
        versionSpecificUniRef90=versionSpecificUniRef90,
        humannDockerImage=humannDockerImage,
        preemptibleAttemptsOverride=preemptibleAttemptsOverride,
        MaxMemGB=MaxMemGB_FunctionalProfileTasks,
        MaxDiskGB=MaxDiskGB_FunctionalProfileTasks,
        FuncProfileMemBase = FuncProfileMemBase,
        FuncProfileBaseFileSize = FuncProfileBaseFileSize,
        FuncProfileIncreaseInterval = FuncProfileIncreaseInterval,
        FuncProfileStorage = FuncProfileStorage
      }

      # regroup gene families to ECs
      call Regroup as RegroupECs {
        input:
        GeneFamiliesFile=FunctionalProfile.GeneFamiliesFile,
        versionSpecificUtilityMapping=versionSpecificUtilityMapping,
        OutFileName=PairPaths[sample_index][2]+"_ecs.tsv",
        humannDockerImage=humannDockerImage,
        groupName="uniref90_level4ec",
        customUtilityMapping=customUtilityMappingEC,
        customMappingPath=customMappingECPath
      }
      
      # regroup gene families to KOs
      call Regroup as RegroupKOs {
        input:
        GeneFamiliesFile=FunctionalProfile.GeneFamiliesFile,
        versionSpecificUtilityMapping=versionSpecificUtilityMapping,
        OutFileName=PairPaths[sample_index][2]+"_kos.tsv",
        humannDockerImage=humannDockerImage,
        groupName="uniref90_ko",
        customUtilityMapping=customUtilityMappingKO,
        customMappingPath=customMappingKOPath
      }      

      call Regroup as RegroupRXNs {
        input:
        GeneFamiliesFile=FunctionalProfile.GeneFamiliesFile,
        versionSpecificUtilityMapping=versionSpecificUtilityMapping,
        OutFileName=PairPaths[sample_index][2]+"_rxns.tsv",
        humannDockerImage=humannDockerImage,
        groupName="uniref90_rxn",
        specifyGroup="uniref90_rxn"
      }
    }
  
    # count the features during each alignment step
    call FunctionalCount {
      input:
      FunctionalLogFiles=FunctionalProfile.LogFile,
      OutFileName=FunctionalCountFileName,
      workflowsDockerImage=workflowsDockerImage
    }
  
    call JoinTables as JoinGeneFamilies {
      input:
      InFiles=FunctionalProfile.GeneFamiliesFile,
      OutFileName=JoinGeneFamilesOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefaultGenes
    }
    call JoinTables as JoinECs {
      input:
      InFiles=RegroupECs.OutFile,
      OutFileName=JoinECsOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }
    call JoinTables as JoinKOs {
      input:
      InFiles=RegroupKOs.OutFile,
      OutFileName=JoinKOsOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }
    call JoinTables as JoinRXNs {
      input:
      InFiles=RegroupRXNs.OutFile,
      OutFileName=JoinRXNsOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }
    call JoinTables as JoinPathways {
      input:
      InFiles=FunctionalProfile.PathwayAbundanceFile,
      OutFileName=JoinPathwaysOutFileName,
      humannDockerImage=humannDockerImage,
      MaxMemGB=JoinNormMemDefault
    }
  }

  if (! setbypassStrainProfiling ) {
   scatter (sample_index in range(length(PairPaths))) {
    call StrainMarkers {        
      input:
        TaxonomicProfileSam=TaxonomicProfile.TaxonomicProfileSam[sample_index],
        sample=PairPaths[sample_index][2],
        strainphlanDockerImage=strainphlanDockerImage,
        preemptibleAttemptsOverride=preemptibleAttemptsOverride,
        MaxMemGB=MaxMemGB_TaxonomicProfileTasks
      }
    }
    
    if (! defined(CustomStrainList) ) {
      call StrainList {
        input:
          InFiles=StrainMarkers.StrainMarkersOutput,
          OutFileName=StrainPhlAnCladeList,
          strainphlanDockerImage=strainphlanDockerImage,
          preemptibleAttemptsOverride=preemptibleAttemptsOverride,
          MaxMemGB=MaxMemGB_TaxonomicProfileTasks
        }
      }
      
    call PickStrains {
      input:
        StrainList=StrainList.OutFile,
        TaxonomicProfileFile=JoinTaxonomicProfiles.OutFile,
        CustomStrainList=CustomStrainList,
        workflowsDockerImage=workflowsDockerImage
      }
    
	scatter (StrainNumber in range(setMaxStrains)){     
      if (StrainNumber < length(PickStrains.SelectedStrains)) {
      
        call PickReferences {
           input:
           StrainList=PickStrains.SelectedStrains,
           StrainNumber=StrainNumber,
           References=StrainphlanReferences,
           workflowsDockerImage=workflowsDockerImage
        }
        
        call StrainProfile {
           input:
           InFiles=StrainMarkers.StrainMarkersOutput,
           StrainList=PickStrains.SelectedStrains,
           StrainNumber=StrainNumber,
           References=StrainphlanReferences,
           ReferencesOptions=PickReferences.SelectedReferences[0],
           strainphlanDockerImage=strainphlanDockerImage,
           preemptibleAttemptsOverride=preemptibleAttemptsOverride,
           MaxMemGB=MaxMemGB_TaxonomicProfileTasks
        }
      }
    }
  }

  # count the reads during each step of QC
  call QCReadCount {
    input:
    LogFiles=QualityControl.LogFile,
    OutFileName=QCReadCountFileName,
    kneaddataDockerImage=kneaddataDockerImage
  }

  # join all taxonomic profiles, gene families, ecs, and pathways (including relative abundance files) from all samples
  call JoinTaxonomicProfiles {
    input:
    InFiles=TaxonomicProfile.TaxonomicProfileFile,
    OutFileName=JoinedTaxonomicProfilesFileName,
    workflowsDockerImage=workflowsDockerImage,
    MaxMemGB=JoinNormMemDefault
  }

  call Collect {
    input:
    joinTableEC=JoinECs.OutFile,
    joinTableKO=JoinKOs.OutFile,
    joinTableRXN=JoinRXNs.OutFile,
    joinTablePathways=JoinPathways.OutFile,
    joinTableGeneFamilies=JoinGeneFamilies.OutFile,
    joinTaxonomicProfiles=JoinTaxonomicProfiles.OutFile,
    qcReadCount = QCReadCount.OutFile,
    dockerImage = metaphlanDockerImage
  }

  output {
    Array[File] ecOutput = Collect.ecOutput
    Array[File] koOutput = Collect.koOutput
    Array[File] rxnOutput = Collect.rxnOutput
    Array[File] pathwayOutput = Collect.pathwayOutput
    Array[File] geneFamOutput = Collect.geneFamOutput
    File qcOutput = Collect.qcOutput
    File taxoProfileOutput = Collect.taxoProfileOutput
  }
}

task QualityControl {
  input {
    File rawfile1
    File rawfile2
    String sample
    String adapterType
    File humanDB
    File transcriptDB
    File rrnaDB
    File? customDB1
    File? customDB2
    File? customDB3
    String dataType
    String kneaddataDockerImage
    Int? MaxMemGB
    Int? preemptibleAttemptsOverride
    Int QCMemBase
    Int QCMemBaseFileSize
    Int QCMemIncreaseInterval
  }
  Int fileSize = ceil(size(rawfile1, 'GB'))
  Int fileLargerBy = fileSize - QCMemBaseFileSize
  Int scaledMem = if fileLargerBy <= 0 then QCMemBase else QCMemBase + fileLargerBy * QCMemIncreaseInterval
  Int mem = select_first([MaxMemGB, scaledMem])
  Int preemptible_attempts = select_first([preemptibleAttemptsOverride, 2])
  
  String useCustomDB1 = if defined(customDB1) then "yes" else "no"
  String useCustomDB2 = if defined(customDB2) then "yes" else "no"
  String useCustomDB3 = if defined(customDB3) then "yes" else "no"
  
  String humanDatabase = "databases/kneaddata_human/"
  String transcriptDatabase = "databases/kneaddata_transcript/"
  String rrnaDatabase = "databases/kneaddata_rrna/"
  
  # Add additional database to run options depending on data type
  String options = if dataType == "mtx" then "--reference-db ${transcriptDatabase} --reference-db ${rrnaDatabase}" else ""
  
  String customDatabase1 = "databases/db1/"
  String customDatabase2 = "databases/db2/"
  String customDatabase3 = "databases/db3/"
  String custom_options = if defined(customDB2) then "--reference-db ${customDatabase1} --reference-db ${customDatabase2}" else "--reference-db ${customDatabase1}"
  String custom_options_add = if defined(customDB3) then "--reference-db ${customDatabase3} " else " "
  

  # download the two reference databases and then run kneaddata.
  command <<< 
  
    # download second custom database if set
    if [ ~{useCustomDB3} == 'yes' ]; then
        mkdir -p ~{customDatabase3}
        tar xzvf ~{customDB3} -C ~{customDatabase3}
    fi

    # download second custom database if set
    if [ ~{useCustomDB2} == 'yes' ]; then
        mkdir -p ~{customDatabase2}
        tar xzvf ~{customDB2} -C ~{customDatabase2}
    fi
  
    # use custom databases if provided instead of reference
    if [ ~{useCustomDB1} == 'yes' ]; then
        mkdir -p ~{customDatabase1}
        tar xzvf ~{customDB1} -C ~{customDatabase1}
        
        #run kneaddata with custom databases
        kneaddata --input ~{rawfile1} --input ~{rawfile2} --output ./ --serial \
        --threads 8 --output-prefix ~{sample} --cat-final-output --run-fastqc-start ~{custom_options} ~{custom_options_add} --sequencer-source ~{adapterType}
    fi
    
    if [ ~{useCustomDB1} == 'no' ]; then
        # download the human reference
        mkdir -p ~{humanDatabase}
        kneaddata_database --download human_genome bowtie2 ~{humanDatabase} --database-location ~{humanDB}
    
        # if data is of type mtx, then download additional database
        if [ ~{dataType} == 'mtx' ]; then
            #create databases
            mkdir -p ~{transcriptDatabase}
            kneaddata_database --download human_transcriptome bowtie2 ~{transcriptDatabase} --database-location ~{transcriptDB}

            mkdir -p ~{rrnaDatabase}
            kneaddata_database --download ribosomal_RNA bowtie2 ~{rrnaDatabase} --database-location ~{rrnaDB}
        fi
    
        #run kneaddata with two reference databases
        kneaddata --input ~{rawfile1} --input ~{rawfile2} --output ./ --serial --reference-db ~{humanDatabase} \
        --threads 8 --output-prefix ~{sample} --cat-final-output --run-fastqc-start ~{options} --sequencer-source ~{adapterType}
    fi
    
    # gzip outputs to save space
    gzip *.fastq
  >>>
  
  output {
    File QCFastqFile = "${sample}.fastq.gz"
    File LogFile = "${sample}.log"
    Array[File] ContaminateReads = glob("*contam*.fastq.gz") # Keep the intermediate contaminate sequences
    Array[File] FastQCOutputsZip = glob("fastqc/*.zip") # Keep the fastqc output files (zip)
    Array[File] FastQCOutputsHtml = glob("fastqc/*.html") # Keep the fastqc output files (html)
  }
  
  runtime {
    docker: kneaddataDockerImage
    cpu: 8
      memory: mem + " GB"
      preemptible: preemptible_attempts
      disks: "local-disk 501 SSD"
  }
}

task TaxonomicProfile {
  input {
    File QCFastqFile
    String sample
    File chocophlanDB
    File chocophlanDBMD5
    File chocophlanIndex
    File chocophlanIndexMD5
    String metaphlanDockerImage
    Int? MaxMemGB
    Int? MaxDiskGB
    Int? preemptibleAttemptsOverride
    Int TaxProfileMemBase
    Int TaxProfileDiskBase
  }
  
  Int mem = select_first([MaxMemGB, TaxProfileMemBase])
  Int disk = select_first([MaxDiskGB, TaxProfileDiskBase])
  Int preemptible_attempts = select_first([preemptibleAttemptsOverride, 2])
  
  String tmpdir = "tmp/"
  String local_db = tmpdir + "metaphlan_databases/"
  String index = "mpa_vOct22_CHOCOPhlAnSGB_202212"

  # create a temp directory and then run metaphlan
  command {
    mkdir -p ${tmpdir}
    mkdir ${local_db}
    mv ${chocophlanDB} ${chocophlanDBMD5} ${chocophlanIndex} ${chocophlanIndexMD5} -t ${local_db}

    metaphlan ${QCFastqFile} --input_type fastq --nproc 8 --no_map --tmp_dir ${tmpdir} --index ${index} \
    --bowtie2db ${local_db} --output_file ${sample}.tsv --samout ${sample}.sam
  }
    
  output {
    File TaxonomicProfileFile = "${sample}.tsv"
    File TaxonomicProfileSam = "${sample}.sam"
  }
  
  runtime {
    docker: metaphlanDockerImage
    cpu: 8
      memory: mem + " GB"
      preemptible: preemptible_attempts
      disks: "local-disk " + disk + " SSD"
  }
}

task StrainMarkers {
  input {
    File TaxonomicProfileSam
    String sample
    String strainphlanDockerImage
    Int? MaxMemGB
    Int? preemptibleAttemptsOverride
  }
  
  Int mem = select_first([MaxMemGB, 8])
  Int preemptible_attempts = select_first([preemptibleAttemptsOverride, 2])

  # create a temp directory and then run strainphlan
  command {
    sample2markers.py --input ${TaxonomicProfileSam} --input_format sam --output_dir . --nprocs 1
  }
    
  output {
    File StrainMarkersOutput = "${sample}.pkl"
  }
  
  runtime {
    docker: strainphlanDockerImage
    cpu: 1
      memory: mem + " GB"
      preemptible: preemptible_attempts
      disks: "local-disk 10 SSD"
  }
}

task StrainList {
  input {
    Array[File] InFiles
    String OutFileName
    String strainphlanDockerImage
    Int? MaxMemGB
    Int? preemptibleAttemptsOverride
  }
  
  Int mem = select_first([MaxMemGB, 10])

  command {
    for infile in ${sep=' ' InFiles}; do ln -s $infile; done
    
    strainphlan --samples ./*.pkl --output_dir . --print_clades_only > ${OutFileName}
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: strainphlanDockerImage
    cpu: 1
      memory: mem+" GB"
      disks: "local-disk 50 SSD"
  }
}

task PickStrains {
  input {
    File? StrainList
    File TaxonomicProfileFile
    File? CustomStrainList
    String workflowsDockerImage
  }
  
  String UserSelectedStrainList = if defined(CustomStrainList) then "yes" else "no"
  
  command {
    python3 <<CODE
    
    from biobakery_workflows import utilities
    
    if "${UserSelectedStrainList}" == "yes":
        for line in open("${CustomStrainList}"):
            if "s__" in line:
                print(line.rstrip())    
    else:
        species_ranked = utilities.rank_species_average_abundance("${TaxonomicProfileFile}")
    
        clades=set()
        with open("${StrainList}") as file_handle:
            for line in file_handle:
                if "s__" in line:
                    clades.add(line.strip().split("\t")[1].split(": in ")[0])
        
        for taxon in species_ranked:
            if taxon in clades:
                print(taxon)
    CODE
  }
  
  output {
    Array[String] SelectedStrains = read_lines(stdout())
  }
  
  runtime {
    docker: workflowsDockerImage
    cpu: 1
      memory: "5 GB"
      disks: "local-disk 5 SSD"
  }
}

task PickReferences {
  input {
    Array[String] StrainList
    Int StrainNumber
    File References
    String workflowsDockerImage
  }
  
  String CurrentClade = StrainList[StrainNumber]

  command {
    python3 <<CODE
   
    import tarfile

    from biobakery_workflows import data

    refnames = tarfile.open("${References}").getnames()
 
    # get the list of reference genomes
    genomes=set()
    with open(data.get_file("strainphlan_species_gcf.tsv")) as file_handle:
        for line in file_handle:
            if line.startswith("${CurrentClade}"):
                genomes.add(line.rstrip().split("\t")[-1]+".fna.bz2")
        
    genomes=list(filter(lambda x: x in refnames, genomes))
    
    # if more than six references are found, limit total used
    # this resolves e. coli which has so many references it can exceed command line length
    if len(genomes) > 6:
        genomes = genomes[:6]    
        
    # add the reference genome files to the command, if any are found
    command_option =  " "
    if len(genomes):
        command_option = " --references ./references/" + " --references ./references/".join(genomes)

    print(command_option)

    CODE
  }
  
  output {
    Array[String] SelectedReferences = read_lines(stdout())
  }
  
  runtime {
    docker: workflowsDockerImage
    cpu: 1
      memory: "5 GB"
      disks: "local-disk 15 SSD"
  }
}


task StrainProfile {
  input {
    Array[File] InFiles
    Array[String] StrainList
    Int StrainNumber
    File References
    String ReferencesOptions
    String strainphlanDockerImage
    Int? MaxMemGB
    Int? preemptibleAttemptsOverride
  }
  
  Int mem = select_first([MaxMemGB, 10])
  
  String CurrentClade = StrainList[StrainNumber]

  command {
    
    echo ${CurrentClade}
    
    mkdir ${CurrentClade}
    
    mkdir references && tar xzvf ${References} -C references

    for infile in ${sep=' ' InFiles}; do ln -s $infile; done

    extract_markers.py -c ${CurrentClade} -o ./

    strainphlan -s ./*.pkl -o ${CurrentClade} -c ${CurrentClade} -n 8 ${ReferencesOptions}

    tar -czf ${CurrentClade}.tar.gz ${CurrentClade}
  }

  output {
    File OutFile = "${CurrentClade}.tar.gz"
  }

  runtime {
    docker: strainphlanDockerImage
    cpu: 8
      memory: mem+" GB"
      disks: "local-disk 50 SSD"
  }
}

task FunctionalProfile {

  input {
    File QCFastqFile
    File TaxonomicProfileFile
    String sample
    File versionSpecificChocophlan
    File versionSpecificUniRef90
    String humannDockerImage
    Int? MaxMemGB
    Int? MaxDiskGB
    Int? preemptibleAttemptsOverride
    Int FuncProfileMemBase
    Int FuncProfileBaseFileSize
    Int FuncProfileIncreaseInterval
    Int FuncProfileStorage
  }
  Int fileSize = ceil(size(QCFastqFile, 'GB'))
  Int fileLargerBy = fileSize - FuncProfileBaseFileSize
  Int scaledMem = if fileLargerBy <= 0 then FuncProfileMemBase else FuncProfileMemBase + fileLargerBy * FuncProfileIncreaseInterval

  Int mem = select_first([MaxMemGB, scaledMem])
  Int disk = select_first([MaxDiskGB, FuncProfileStorage])
  Int preemptible_attempts = select_first([preemptibleAttemptsOverride, 2])
  
  String databases = "databases/"

  # download the two reference databases and run humann
  command {
    mkdir -p ${databases}
    humann_databases --download chocophlan full ${databases} --database-location ${versionSpecificChocophlan}
    humann_databases --download uniref uniref90_diamond ${databases} --database-location ${versionSpecificUniRef90}

    humann --input ${QCFastqFile} --output ./ --taxonomic-profile ${TaxonomicProfileFile} --threads 8 --o-log ${sample}.log
    }
    
    output {
      File GeneFamiliesFile = "${sample}_genefamilies.tsv"
      File PathwayAbundanceFile = "${sample}_pathabundance.tsv"
      File PathwayCoverageFile = "${sample}_pathcoverage.tsv"
      File LogFile = "${sample}.log"
      Array[File] UnalignedReads = glob("${sample}_humann_temp/*.fa") # Keep the unaligned reads after each mapping step
    }

  runtime {
    docker: humannDockerImage
    cpu: 8
      memory: mem + " GB"
      disks: "local-disk "+ disk +" SSD"
      preemptible: preemptible_attempts
  }
}
task Regroup {
  input {
    File GeneFamiliesFile
    File versionSpecificUtilityMapping
    String OutFileName
    String humannDockerImage
    String groupName
    File? customUtilityMapping
    String? customMappingPath
    String? specifyGroup
  }

  String databases = "databases/"

  String customMapping = if defined(customUtilityMapping) then "-c" else ""
  String groupFlag = if defined(specifyGroup) then "-g" else ""
  # download the utility databases and regroup to ECs
  command {
    mkdir -p ${databases}
    humann_databases --download utility_mapping full ${databases} --database-location ${versionSpecificUtilityMapping}
    humann_regroup_table --input ${GeneFamiliesFile} --output ${OutFileName} ${customMapping} ${customUtilityMapping} ${groupFlag} ${specifyGroup}
  }
    
  output {
      File OutFile = "${OutFileName}"
  }

  runtime {
    docker: humannDockerImage
    cpu: 1
      memory: "10 GB"
      disks: "local-disk 20 SSD"
  }
}

task RenormTable {
  input {
    File InFile
    String OutFileName
    String humannDockerImage
    Int? MaxMemGB
  }
  
  Int mem = select_first([MaxMemGB, 10])

  String databases = "databases/"

  # download the utility databases and renorm tables to relative abundance
  command {
    humann_renorm_table --input ${InFile} --output ${OutFileName} --units relab --special n
  }

  output {
      File OutFile = "${OutFileName}"
  }

  runtime {
    docker: humannDockerImage
    cpu: 1
      memory: mem+" GB"
      disks: "local-disk 20 SSD"
  }
}

task QCReadCount {
  input {
    Array[File] LogFiles
    String OutFileName
    String kneaddataDockerImage
  }
  
  # symlink input files to working directory
  # create a table of read counts during each qc step
  command {
    for logfile in ${sep=' ' LogFiles}; do ln -s $logfile; done
    
    kneaddata_read_count_table --input ./ --output ${OutFileName}
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: kneaddataDockerImage
    cpu: 1
      memory: "5 GB"
      disks: "local-disk 10 SSD"
  }
}

task JoinTaxonomicProfiles {
  input {
    Array[File] InFiles
    String OutFileName
    String workflowsDockerImage
    Int? MaxMemGB
  }
  
  Int mem = select_first([MaxMemGB, 10])

  # symlink input files to working directory
  # join all files into a single file
  command {
    for infile in ${sep=' ' InFiles}; do ln -s $infile; done
    
    join_taxonomic_profiles.py --input ./ --output ${OutFileName} --file_name .tsv
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: workflowsDockerImage
    cpu: 1
      memory: mem+" GB"
      disks: "local-disk 10 SSD"
  }
}

task JoinTables {
  input {
    Array[File] InFiles
    String OutFileName
    String humannDockerImage
    Int? MaxMemGB
  }
  
  Int mem = select_first([MaxMemGB, 10])

  # symlink input files to working directory
  # join all files into a single file
  command {
    for infile in ${sep=' ' InFiles}; do ln -s $infile; done
    
    humann_join_tables --input ./ --output ${OutFileName} --file_name .tsv
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: humannDockerImage
    cpu: 1
      memory: mem+" GB"
      disks: "local-disk 10 SSD"
  }
}

task FunctionalCount {
  input {
    Array[File] FunctionalLogFiles
    String OutFileName
    String workflowsDockerImage
  }

  # symlink logs to working directory
  # compile alignment counts from humann logs
  command {
    for infile in ${sep=' ' FunctionalLogFiles}; do ln -s $infile; done
  
    get_counts_from_humann_logs.py --input ./ --output ${OutFileName}
  }

  output {
    File OutFile = "${OutFileName}"
  }

  runtime {
    docker: workflowsDockerImage
    cpu: 1
      memory: "5 GB"
      disks: "local-disk 10 SSD"
  }
}

task Collect {
  input {
    File? joinTableEC
    File? joinTableKO
    File? joinTableRXN
    File? joinTablePathways
    File? joinTableGeneFamilies
    File qcReadCount
    File joinTaxonomicProfiles
    String dockerImage
  }

  String ecFolder = "ec/"
  String koFolder = "ko/"
  String rxnFolder = "rxn/"
  String pathwayFolder = "pathway/"
  String geneFamilyFolder = "geneFamilies/"

  Array[String] folders = [ecFolder, koFolder, rxnFolder, pathwayFolder, geneFamilyFolder]
  Int ecFilesize = ceil(size(joinTableEC, 'GB')) 
  Int koFilesize = ceil(size(joinTableKO, 'GB'))
  Int largestSize = if (koFilesize > ecFilesize) then koFilesize else ecFilesize
  Int rxnFilesize = if (ceil(size(joinTableRXN, 'GB')) > largestSize) then ceil(size(joinTableRXN, 'GB')) else largestSize
  Int pathwaysFilesize = if (ceil(size(joinTablePathways, 'GB')) > rxnFilesize) then ceil(size(joinTablePathways, 'GB')) else rxnFilesize
  Int geneFilesize = if (ceil(size(joinTableGeneFamilies, 'GB')) > pathwaysFilesize) then ceil(size(joinTableGeneFamilies, 'GB')) else pathwaysFilesize
  Int memoryNeeded = geneFilesize * 10

  command <<<
    for folder in ~{sep=" " folders}; do mkdir -p $folder; done
    python3<<CODE
    import pandas as pd
    import os.path
    import shutil

    def postprocessJoinTable(file, destFolder):
        folderName = os.path.dirname(file)
        fileName = os.path.basename(file)
        resultName = "summarized_" + fileName
        resultFile = destFolder + resultName
        currentTable = pd.read_csv(file, sep = '\t')
        currentTable.columns = currentTable.columns.str.split("_", 1).str.get(0)
        currentTable.to_csv(file, sep = '\t', index=False)
        leftMostColumn = list(currentTable.columns)[0]
        currentTable.drop(currentTable[currentTable[leftMostColumn].str.contains('\|')].index, inplace=True)
        currentTable.to_csv(resultFile, sep = '\t', index=False)
        shutil.move(file, destFolder)
    postprocessJoinTable("~{joinTableEC}", "~{ecFolder}")
    postprocessJoinTable("~{joinTableKO}", "~{koFolder}")
    postprocessJoinTable("~{joinTableRXN}", "~{rxnFolder}")
    postprocessJoinTable("~{joinTablePathways}", "~{pathwayFolder}")
    postprocessJoinTable("~{joinTableGeneFamilies}", "~{geneFamilyFolder}")
    CODE
  >>>

  output{
    Array[File] ecOutput = glob(ecFolder + "*.tsv")
    Array[File] koOutput = glob(koFolder + "*.tsv")
    Array[File] rxnOutput = glob(rxnFolder + "*.tsv")
    Array[File] pathwayOutput = glob(pathwayFolder + "*.tsv")
    Array[File] geneFamOutput = glob(geneFamilyFolder + "*.tsv")
    File qcOutput = qcReadCount
    File taxoProfileOutput = joinTaxonomicProfiles
  }

  runtime {
    docker: dockerImage
    cpu: 4
    memory: memoryNeeded + " GB"
    preemptible: 2
    disks: "local-disk 50 SSD"
  }
}