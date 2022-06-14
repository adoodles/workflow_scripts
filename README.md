Required Files/Directory Structure for each workflow

1. wtx.wdl
    input:
    - fastq.gz 
    humannv3
    - chocophlan full (current version full_chocophlan.v201901_v31.tar.gz seems to be a version mismatch with humannv3)
    - uniref90_annotated_1_1.tar.gz
    - full_utility_mapping_1_1.tar.gz
    - custom map files:
        - ko_uniref.txt
        - level4ec_uniref.txt
    kneaddata
    - Homo_sapiens_hg37_human_contamination_Bowtie2_v0.1.tar.gz
    - Homo_sapiens_hg38_transcriptome_Bowtie2_v0.1.tar.gz
    - SILVA_128_LSUParc_SSUParc_ribosomal_RNA_v0.2.tar.gz


2. metawibele.wdl (*SKIP FIRST*)
    input:
    - fastq.gz DIRECTORY** for Preprocess step
    - fastaa files for Characterize step

    - uniref90 database (only if global homonology is required. atm tentative skip)
    - 

    mspminer (try to make it work)

3. shortbred.wdl
    identify step (creation of markers)
    input:
    - proteins of interest (goi?) in fasta format
    - protein database (uniref90)

    quantify step (using markers from identify/self provided)
    input:
    - markers (fasta)
        try antibiotic markers first
    - sequence reads



task:
create a mspminer .wdl file to get output as sample table
    -how to use can see from https://github.com/biobakery/metawibele/blob/master/metawibele/tasks/characterization.py (ctrl + f "mspminer")