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


** Key points while doing sample test runs:
    - sample should consist of the largest sample + smallest sample + 5-10 random ones (limit testing)



605ffc3d-2b43-4230-bf79-dbe5b91c6ab7 --> ~1 hour (70 minutes more accurately) blood sample


Questions:

    - Collect and output results at the end of each scatter task (QC, profiling) (save output to somewhere else?)
    - Output will be checked to remove reliance on call caching. (output will become entry in input table)
    - bypass options provided at every stage/ "stop at"?


mspminer workflow

fastq -> QC (kneaddata) (trim galore? see if kd has trim function) -> raw files -> assembly (megaHIT) w prodigal -> predict (prodigal) -> check evaluation (checkM) (fallback: MetaQUAST) -> 

--> cdhit -> mspminer -> output
--> metaBAT2 -> output

** toolbox image with all required program (including licensed program like usearch/modified settings.ini or can just pull from github)

programs:
    -usearch
    -samtools
    -bwa

**usearch & mspminer settings.ini use wget to download from cloud storage

**** checkM only works on binned stuff: after metabat

# 1) Preprocess the reads and check their quality  * QC steps   FASTQ->FASTQ
# 2) Assemble each sample individually and get the contigs out  * megahit / metaspade contigs  FASTQ->FASTA
# 3) Concatenate the FASTA files together while making sure all contig headers stay unique, and filter away small contigs  *vamb script
# 4) Map the reads to the FASTA file to obtain BAM files  * minimap2? (RAM bug) 
# 5) Run Vamb
# 6) Postprocess the results
