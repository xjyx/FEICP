# FEICP
FEICP (Finding Exon-Intron circular RNAs from paired-end RNA-seq) is for EIciRNAs (Exon-Intron circular RNAs) detection from paired-end RNA-seq.

The whole FEICP pipeline has three main parts: 1. Detecting circular RNAs from RNA-seq using BWA and CIRI2; and 2. Getting those BSJs (back splicing reads) which can be mapped within introns; and 3. Validate the whole intron can be retained within a circular RNA.

If you use FEICP in your study, please cite:
Yan Yang. Characterization of EIciRNAs reveals their distinct sequence features and potential regulatory function in regulation of gene expression. Manuscript in preparation.

# Prerequisites

## Softwares and Packages
Software
  BWA
  CIR2
  STAR
  Bowtie
  Hisat2
  Seqtk
Packages
  pandas
  scipy
 
## RNA-seq
The RNase R-trated paired-end RNA-seq is recommended. Total RNA-seq, poly(A)-/ribo- RNA-seq is also acceptable. According to our experience, the RNase R treatment can significantly enrich EIciRNAs, like other circRNAs.

# Installation
All the 3 steps of FEICP is integrated into a PBS script for qsub on on computing cluster and is out of the box so there's no need to install.

