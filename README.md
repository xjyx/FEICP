# FEICP
FEICP (Finding Exon-Intron circular RNAs from paired-end RNA-seq) is for EIciRNAs (Exon-Intron circular RNAs) detection from paired-end RNA-seq.

The whole FEICP pipeline has three main parts: 1. Detecting circular RNAs from RNA-seq using BWA and CIRI2; and 2. Getting those BSJs (back splicing reads) which can be mapped within introns; and 3. Validate the whole intron can be retained within a circular RNA.

If you use FEICP in your study, please cite:
Yan Yang. Characterization of EIciRNAs reveals their distinct sequence features and potential regulatory function in regulation of gene expression. Manuscript in preparation.

# Prerequisites

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
The FEICP pipeline need the reference genome file and the gene annotation file used to find and annotate the circular.  
The recommended gene annotation file are gencode, ENSEMBL or NCBI RefSeq GTF files. And the reference genome file should contain all the genome sequences.


# Usage
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
