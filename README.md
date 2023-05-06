# FEICP
FEICP (Finding Exon-Intron circular RNAs from paired-end RNA-seq) is for EIciRNAs (Exon-Intron circular RNAs) detection from paired-end RNA-seq.

The whole FEICP pipeline has three main parts: 
  - **1. Detecting circular RNAs from RNA-seq using BWA and CIRI2**;  
  - **2. Getting those BSJs (back splicing reads) which can be mapped within introns**;  
  - **3. Validate the whole intron can be retained within a circular RNA**.  

If you use FEICP in your study, please cite:  
Yinchong Zhong, Yan Yang et al. Genome-wide characterization of EIciRNAs reveals their distinct sequence features and potential function in regulation of gene expression. Manuscript in preparation.  

# Prerequisites
I have tested the pipeline on Red Hat 4.8.5-44 on a 64Bit machine with python 3.7. 

## Softwares and Packages
  - **Software** <br>
    - BWA<br>
    - CIRI2<br>
    - STAR<br>
    - Bowtie<br>
    - STAR<br>
    - Hisat2<br> 
    - samtools<br>
    - bedtools<br>
    - seqtk<br>

  - **Packages**  
    - pandas  
    - scipy  
 
## RNA-seq
The RNase R-treated paired-end RNA-seq is recommended. Total RNA-seq, poly(A)-/ribo- RNA-seq is also acceptable. According to our experience, the RNase R treatment can significantly enrich EIciRNAs, like other circRNAs.

# Installation
All the 3 steps of FEICP is integrated into a PBS script for qsub on on computing cluster and is out of the box so there's no need to install.

# Setup
The FEICP pipeline need the reference genome file and the gene annotation file used to find and annotate the circular RNAs.  
The recommended gene annotation file are gencode, ENSEMBL or NCBI RefSeq GTF files. And the reference genome file should contain all the genome sequences.


# Usage

## A quick start
To get EIciRNAs from multiple paired-end RNA-seq data in multiple processes, clone this repository in local and run:
qsub -v      src_dir="$src_dir",ref_fasta="$ref_fasta",ref_gtf="$ref_gtf",bwa_index="$bwa_index",star_index="$star_index",nth="$nth",ncpus="$ncpus",path="$path",sample_list="$sample_list" FEICP_parallel.pbs

In the above command,
src_dir is the path to the src of this pipeline
ref_fasta is the path to the reference genome fasta file
ref_gtf is the path to the gtf file with gene annotation
bwa_index is the path to the bwa index of reference genome
star_index is the path to the STAR index generated using ref_fasta and ref_gtf
nth is the number of samples from which you want to find EIciRNAs at a time
ncpus is the number of cores used in FEICP for each sample, this parameter will be used manily in mapping stage
path is the path to the fastq file

Build a text file with one sample name each row
For example
```
cat sample_list.txt
```
### test_1
### test_2
and the path directory must contain test_1_1.fastq.gz, test_1_2.fastq.gz; test_2_1.fastq.gz, test_2_2.fastq.gz
Then you can 
After getting the reference genome and gene annotation gtf file, we can run the FEICP pipeline directly without needing to prepeare other annotation file, such as intron, exon-intron, exon-exon file, etc bacause the program itself will get them all.
```
FEICP.py --help  

usage: Getting potential EIciRNAs form RNA-seq [-h] --fastqs FASTQS  
                                               [FASTQS ...] [--paired-end]  
                                               [--sam SAM] [--circ CIRC]  
                                               --genome GENOME --gtf GTF  
                                               [--bwa BWA]  
                                               [--bwa-index BWA_INDEX]  
                                               [--bowtie BOWTIE]  
                                               [--bowtie-build BOWTIE_BUILD]  
                                               [--hisat2 HISAT2]  
                                               [--hisat2-build HISAT2_BUILD]  
                                               [--seqtk SEQTK]  
                                               [--samtools SAMTOOLS]  
                                               [--bedtools BEDTOOLS]  
                                               [--mismatch MISMATCH]  
                                               [--min-overlap MIN_OVERLAP]  
                                               [--src-dir SRC_DIR]  
                                               [--output OUTPUT] --sample  
                                               SAMPLE [--threads THREADS]  
                                               [--log-level {NOTSET,DEBUG,INFO,WARNING,CRITICAL,ERROR,CRITICAL}]  

optional arguments:  
  -h, --help            show this help message and exit  
  --fastqs FASTQS [FASTQS ...], -i FASTQS [FASTQS ...]  
                        List of FASTQs (R1 and R2). FASTQs must be compressed  
                        with gzip (with .gz).  
  --paired-end          Paired-end FASTQs  
  --sam SAM             Sam-formatted file from BWA mapping used to get  
                        circular RNAs  
  --circ CIRC           Circular RNAs file from CIRI2 used to find EIciRNAs  
  --genome GENOME, -f GENOME  
                        The genome fasta file  
  --gtf GTF, -a GTF     The annotation gtf file  
  --bwa BWA             The path to exectuable bwa  
  --bwa-index BWA_INDEX  
                        The path to the bwa index of reference genome fasta  
  --bowtie BOWTIE       The path to exectuable bowtie  
  --bowtie-build BOWTIE_BUILD  
                        The path to exectuable bowtie-build  
  --hisat2 HISAT2       The path to exectuable hisat2  
  --hisat2-build HISAT2_BUILD  
                        The path to exectuable hisat2-build  
  --seqtk SEQTK         The path to exectuable seqtk  
  --samtools SAMTOOLS   The path to exectuable samtools  
  --bedtools BEDTOOLS   The path to exectuable bedtools  
  --mismatch MISMATCH, -m MISMATCH  
                        The maximum number of mismatch permitted during bowtie  
                        mapping  
  --min-overlap MIN_OVERLAP  
                        The minimum overhang between read and reference intron  
  --src-dir SRC_DIR, -c SRC_DIR  
                        The directory storing source code files of CIRI2  
  --output OUTPUT, -o OUTPUT  
                        The output directory  
  --sample SAMPLE, -n SAMPLE  
                        The name of sample used to generate the temp file  
  --threads THREADS, -@ THREADS  
                        The number of CPUs to use  
  --log-level {NOTSET,DEBUG,INFO,WARNING,CRITICAL,ERROR,CRITICAL}  
                        Log level  
```

### Output format ###
| column | name | description |
|--------|------|-------------|
| 1 | chr | chromosome name |
| 2 | start | the start coordinate (0-based) |
| 3 | end | the end coordinate (1-based) |
| 4 | gene_id | the annotated gene name of EIciRNA |
| 5 | EIci_junc_count | the BSJ (back-splice junction) count of EIciRNA |
| 6 | strand | the strand of EIciRNA |
| 7 | intron_start | the start coordinate of retained intron (0-based) |
| 8 | intron_end | the end coordinate of retained intron (1-based) |
| 9 | intron_cov | the reads coverage of retained intron |
| 10 | EI | the junction count of left_exon-intron |
| 11 | IE | the junction count of intron_right_exon |
| 12 | EE | the junction count of left_exon-right_exon |
| 13 | intron_count | the count of reads located within the middle 200 bp of retained intron |


### Hot to filter the outputs ###
Considering the low expression nature of EIciRNAs, according to my experience, I recommend you to filter the outputs according to a relaxed but reasonable condition:  
  - **EIci_junc_count** >= 1  
  - **intron_cov** >= 0.50  
  - **EI** >= 1  
  - **IE** >= 1  
  - **EE** >= 1  
  - **intron_count** >= 1  

You can filter the outputs according to the above condition using your favorite tools. As for me, I am used to processing this kind of task using **pandas** or **tidyverse** or just a single UNIX command:  
```
cat <out_file> | awk 'BEGIN{OFS=FS="\t"}{if(NR==1 || ($9 >= 0.50 && $10*$11*$12*$12*$13)){print}}' > <out_filtered>  
```
