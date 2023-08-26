# FEICP
FEICP (Finding Exon-Intron circular RNAs from paired-end RNA-seq) is for EIciRNAs (Exon-Intron circular RNAs) detection from paired-end RNA-seq.

The whole FEICP pipeline has three main parts: 
  - **1. Detecting circular RNAs from RNA-seq using BWA and CIRI2**;  
  - **2. Getting those BSJs (back splicing reads) which can be mapped within introns**;  
  - **3. Validate the whole intron can be retained within a circular RNA**.  

If you use FEICP in your study, please cite:  
Yinchong Zhong, Yan Yang et al. Systematic identification and characterization of Exon-intron circRNAs.  

# Prerequisites
I have tested the pipeline on Red Hat 4.8.5-44 on a 64Bit machine with python 3.7. 

## Softwares & packages
| Software | Link |
| -------- | ------------------------ |
| Python | https://www.python.org/ |
| BWA | https://github.com/lh3/bwa |
| Bowtie | https://github.com/BenLangmead/bowtie |
| HISAT2 | https://github.com/DaehwanKimLab/hisat2 |
| STAR | https://github.com/alexdobin/STAR |
| CIRI2 | https://sourceforge.net/projects/ciri/files/CIRI2/ |
| Pandas | https://pandas.pydata.org/ |
| Scipy | https://scipy.org/ |
| samtools | https://github.com/samtools/samtools |
| bedtools | https://github.com/arq5x/bedtools2 |
| seqtk | https://github.com/lh3/seqtk |
 
## RNA-seq
The RNase R-treated paired-end RNA-seq is recommended. Total RNA-seq, poly(A)-/ribo- RNA-seq is also acceptable. According to our experience, the RNase R treatment can significantly enrich EIciRNAs, like other circRNAs.

# Installation
All the 3 steps of FEICP is integrated into a PBS script for qsub on computing cluster and is out of the box so there's no need to install.

# Setup
The FEICP pipeline need the reference genome file and the gene annotation file used to find and annotate the circular RNAs.  
The recommended gene annotation file are gencode, ENSEMBL or NCBI RefSeq GTF files. And the reference genome file should contain all the genome sequences.


# Usage

## A quick start
To get EIciRNAs from multiple paired-end RNA-seq data in multiple processes, clone this repository in local and run:
```
qsub -v src_dir="$src_dir",ref_fasta="$ref_fasta",ref_gtf="$ref_gtf",bwa_index="$bwa_index",star_index="$star_index",nth="$nth",ncpus="$ncpus",path="$path",sample_list="$sample_list" FEICP_parallel.pbs
```
In the above command,  
  - **src_dir** is the path to the src of this pipeline  
  - **ref_fasta** is the path to the reference genome fasta file  
  - **ref_gtf** is the path to the gtf file with gene annotation  
  - **bwa_index** is the path to the bwa index of reference genome  
  - **star_index** is the path to the STAR index generated using ref_fasta and ref_gtf  
  - **nth** is the number of samples from which you want to find EIciRNAs at a time  
  - **ncpus** is the number of cores used in FEICP for each sample, this parameter will be used manily in mapping stage  
  - **path** is the path to the fastq file  
  - **sample_list** is the text file with one prefix of gzipped fastq files each row and the names will be used in the output of downstream analysis. For example, is there are fastq files in the path: test_1_1.fastq.gz, test_1_2.fastq.gz; test_2_1.fastq.gz, test_2_2.fastq.gz, then you should build a text file test.txt
  
```
test_1
test_2
```

## Step by step pipeline

### 1. Prepare the necessary reference files

Down load the reference genome and gene annotation gtf file from GENCODE  
```
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
wget -c https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.annotation.gtf.gz
gunzip GRCh38.primary_assembly.genome.fa.gz gencode.v34.annotation.gtf.gz
```
Make BWA index
```
bwa index -p hg38_bwa GRCh38.primary_assembly.genome.fa
```
Make STAR index according to the STAR manual (https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
```
STAR \
        --runThreadN 6 \
        --runMode genomeGenerate \
        --genomeDir /path/to/star_index \
        --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
        --sjdbGTFfile cat gencode.v34.annotation.gtf \
        --sjdbOverhang $(mate_length - 1)
```
Make annotation file used in FEICP
```
cat gencode.v34.annotation.gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="exon"){split($9,A,"\"");print $1,$4-1,$5,A[2],A[4],$7}}' | sort -k1,1 -k2,2n >exon.bed

cat gencode.v34.annotation.gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="transcript"){split($9,A,"\"");print $1,$4-1,$5,A[2],A[4],$7}}' | sort -k1,1 -k2,2n >transcript.bed

cat exon.bed | sort -k1,1 -k2,2n | uniq | python make_EI.py - | sort -k1,1 -k2,2n >exon_intron.bed

python make_intron.py --gtf cat gencode.v34.annotation.gtf --output intron.bed

cat /path/to/star_index/chrNameLength.txt | sort -k1,1 >genome.txt
```

### 2. Finding circular RNAs using CIRI2 pipeline and get circRNAs which are likely EIciRNAs using FEICP.py

```
python FEICP.py --help  

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
  --CIRI2 CIRI2, -c CIRI2  
                        The path to CIRI2.pl
  --output OUTPUT, -o OUTPUT  
                        The output directory  
  --sample SAMPLE, -n SAMPLE  
                        The name of sample used to generate the temp file  
  --threads THREADS, -@ THREADS  
                        The number of CPUs to use  
  --log-level {NOTSET,DEBUG,INFO,WARNING,CRITICAL,ERROR,CRITICAL}  
                        Log level  
```


With the necessary software in the environment, you can run this simply as 
```
FEICP.py \
    --threads 8 \
    --fastqs test_1.fastq.gz test_2.fastq.gz \
    --paired-end \
    --genome GRCh38.primary_assembly.genome.fa \
    --gtf gencode.v34.annotation.gtf \
    --CIRI2 src/CIRI2.pl \
    --output test_FEICP \
    --sample test
```

### 3. Get the circRNAs with a whole intron retained
Map RNA-seq reads against to reference genome using splice-aware aligner STAR
```
STAR 	\
    --genomeLoad NoSharedMemory \
    --runThreadN 8 \
    --genomeDir /path/to/star_index \
    --outFilterMultimapNmax 1 \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix test_STAR_ \
    --outSAMunmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn test_1.fastq.gz test_2.fastq.gz \
    --readFilesCommand gzip -dc
```
Compute genome coverage using deepTools
```
bamCoverage \
    --bam test_STAR_Aligned.sortedByCoord.bam \
    --outFileName test_STAR.bigWig \
    --binSize 10 \
    --numberOfProcessors 8 \
    --normalizeUsing CPM
```
Compute the metrics of introns: Exon_left-intron, intron-Exon_right, Exon_left-Exon_right and reads coverage of intron
```
cat test_EIciRNA.bed | awk 'BEGIN{OFS=FS="\t"}{if(NR>1){split($9,A,",");split($10,B,",");for(i in A){print $1,A[i],B[i],$7,$4,$6}}}' >test_EIci_intron.bed

samtools view -@ 8 -h -q 255 test_STAR_Aligned.sortedByCoord.bam | \
    awk '$0 ~ /^@/ || $6 ~ /N/' | \
    samtools view -@ 8 -b | bedtools bamtobed -i - -split | \
    LC_COLLATE=C sort -k1,1 -k2,2n --parallel=8 -T ./ >test_STAR_Aligned_spliced.bed

samtools view -@ 8 -h -q 255 test_STAR_Aligned.sortedByCoord.bam | \
    bedtools bamtobed -i - -split | \
    LC_COLLATE=C sort -k1,1 -k2,2n --parallel=8 -T ./ >test_STAR_Aligned.bed

bedtools map -a exon_intron.bed -b test_STAR_Aligned.bed -f 1 -nonamecheck -g genome.txt -c 4,4 -o count_distinct,distinct | \
    awk 'BEGIN{OFS=FS="\t"}{if($7!=0){tmp=$4;$4=$5;$5=tmp;print}}' >test_EI_count.bed

bedtools intersect -a test_STAR_Aligned_spliced.bed -b transcript.bed -f 1 -sorted -wa -wb -nonamecheck | \
    bedtools groupby -i - -g 1,2,3,4,5,6 -c 11  -o distinct >test_spliced2tx.bed

python get_EE.py --read-tx test_spliced2tx.bed --output test_EE_count.bed

python get_intron_midpoint.py intron.bed | \
    sort -k1,1 -k2,2n | uniq | \
    bedtools map -a - -b test_STAR_Aligned.bed -f 0.1 -F 0.1 -e -nonamecheck -g genome.txt -c 4,4 -o count_distinct,distinct | \
    awk '$7!=0' >test_intron_count.bed

cat intron.bed | LC_COLLATE=C sort -k1,1 -k2,2n --parallel=8 -T ./ | \
    bedtools coverage -a - -b test_STAR_Aligned.bed -nonamecheck -g genome.txt -sorted | uniq >test_intron_cov.bed

python calculate_PIR.py \
    --ee-count test_EE_count.bed \
    --ei-count test_EI_count.bed \
    --intron-count test_intron_count.bed \
    --query-intron test_EIci_intron.bed \
    --output test_EIci_PIR.bed
    
python merge_PIR_cov.py test_EIci_PIR.bed test_intron_cov.bed test_intron_EIci_intron_stats.bed

python get_tidy_EIciRNA.py test_EIciRNA.bed test_intron_EIci_intron_stats.bed test_intron_EIciRNA_tidy.bed
```



### Output format ###
| column | name            | description                                                            |
| ------ | --------------- | ---------------------------------------------------------------------- |
| 1      | chr             | chromosome name                                                        |
| 2      | circ_start           | the start coordinate (0-based)                                         |
| 3      | circ_end             | the end coordinate (1-based)                                           |
| 4      | gene_id         | the annotated gene name of EIciRNA                                     |
| 5      | EIci_junc_count | the BSJ (back-splice junction) count of EIciRNA                        |
| 6      | strand          | the strand of EIciRNA                                                  |
| 7      | intron_start    | the start coordinate of retained intron (0-based)                      |
| 8      | intron_end      | the end coordinate of retained intron (1-based)                        |
| 9      | coverage      | the reads coverage of retained intron                                  |
| 10     | intron              | the count of reads located within the middle 200 bp of retained intron                                 |
| 11     | E1I              | the junction count of left_exon-intron                                |
| 12     | IE2              | the junction count of intron-right_exon                             |
| 13     | E1E2    | the junction count of left_exon-right_exon |
