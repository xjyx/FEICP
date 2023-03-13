# FEICP
FEICP (Finding Exon-Intron circular RNAs from paired-end RNA-seq) is for EIciRNAs (Exon-Intron circular RNAs) detection from paired-end RNA-seq.

The whole FEICP pipeline has three main parts: 1. Detecting circular RNAs from RNA-seq using BWA and CIRI2; and 2. Getting those BSJs (back splicing reads) which can be mapped within introns; and 3. Validate the whole intron can be retained within a circular RNA.

If you use FEICP in your study, please cite:
Yan Yang. Characterization of EIciRNAs reveals their distinct sequence features and potential regulatory function in regulation of gene expression. Manuscript in preparation.

# Requirements

Pyhton3
BWA
CIRI2
seqtk
STAR
bowtie
hisat2
bedtools
pandas
