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

2. vambbin.wdl
    input:
    - input reads fastq.gz
    - human DB (hg37? example uses hg19)


3. shortbred.wdl
    input:
    - fastq.gz

    requirements:
    - marker file 

    special note:
    - identify step is usually skipped. in this case. only inputs required are the marker file and the input reads






605ffc3d-2b43-4230-bf79-dbe5b91c6ab7 --> ~1 hour (70 minutes more accurately) blood sample


Questions:

    VAMB results:
        S1... S5 corresponds to each sample
        SxCx is one cluster. (is this too much clusters?) --> postprocess to become cluster 1... cluster n?
        Need abundance report?

    wtx results:
        what to extract?
        latest successful run:
        https://console.cloud.google.com/storage/browser/fc-e92a41e8-a4b7-4424-88e8-b5592bdebc1c/5534632c-6ed2-4a8a-940b-593c04c7dcee/workflowMTX/817b9843-c6ad-4799-b61e-62485c3b308a;tab=objects?authuser=0&prefix=&forceOnObjectsSortingFiltering=false

        Collect stages:
        -QCReadCount
        -JoinGeneFamilies
        -JoinTaxonomicProfile

        -JoinKO, JoinEC, JoinRXN w/post procesing for summation rows

    Remove relab steps.




# 1) Preprocess the reads and check their quality  * QC steps   FASTQ->FASTQ
# 2) Assemble each sample individually and get the contigs out  * megahit / metaspade contigs  FASTQ->FASTA
# 3) Concatenate the FASTA files together while making sure all contig headers stay unique, and filter away small contigs  *vamb script
# 4) Map the reads to the FASTA file to obtain BAM files  * minimap2? (RAM bug) 
# 5) Run Vamb
# 6) Postprocess the results
