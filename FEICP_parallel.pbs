#!/bin/sh
#PBS -N FEICP
#PBS -o FEICP.log
#PBS -e FEICP.err
#PBS -q fat
#PBS -l nodes=1:ppn=30
#PBS -l walltime=96:00:00

# This script implements the FEICP (Find Exon-Intron circRNAs form paired-end RNA sequencing) pipeline in parallel


function find_EIci(){
    sample=$1
    fastq1="$path/${sample}_1.fastq.gz"
    fastq2="$path/${sample}_2.fastq.gz"
    EIci_out_dir="$path/${sample}_FEICP/${sample}_EIci"
    sam_file="$EIci_out_dir/${sample}_bwa.sam"
    mkdir -p $EIci_out_dir && cd $_
    python3 ${find_EIci_py} \
    --fastqs $fastq1 $fastq2 \
    --paired-end \
    --genome $ref_fasta \
    --gtf $ref_gtf \
    --bwa-index $bwa_index \
    --CIRI2 $CIRI2 \
    --output $EIci_out_dir \
    --sample $sample \
    --threads $ncpus
}


function STAR_IR(){
    sample=$1
    fastq1="$path/${sample}_1.fastq.gz"
    fastq2="$path/${sample}_2.fastq.gz"
    STAR_out_dir="$path/${sample}_FEICP/${sample}_STAR_IR"
    mkdir -p $STAR_out_dir && cd $_
    star_prefix="$STAR_out_dir/${sample}_"
    
    STAR 	\
    --genomeLoad NoSharedMemory \
    --runThreadN $ncpus \
    --genomeDir $star_index \
    --outFilterMultimapNmax 1 \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix ${star_prefix} \
    --outSAMunmapped None \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn $fastq1 $fastq2 \
    --readFilesCommand gzip -dc
}


function sam_metrics(){
    sample=$1
    EIci_out_dir="$path/${sample}_FEICP/${sample}_EIci"
    metrics_out_dir="$path/${sample}_FEICP/${sample}_sam_metrics"
    mkdir -p $metrics_out_dir && cd $_
    sam_bwa="$EIci_out_dir/${sample}_bwa.sam"
    bam_bwa="$metrics_out_dir/${sample}_bwa.bam"
    bam_metrics="$metrics_out_dir/${sample}_bwa_flagstat.txt"
    sambamba-0.8.1-linux-amd64-static view -S -f bam -t $ncpus $sam_bwa | \
    sambamba-0.8.1-linux-amd64-static sort -t $ncpus -m 10G -o $bam_bwa -t $ncpus /dev/stdin
    sambamba-0.8.1-linux-amd64-static flagstat -t $ncpus $bam_bwa >$bam_metrics
}

function compute_cov(){
    sample=$1
    STAR_out_dir="$path/${sample}_FEICP/${sample}_STAR_IR"
    metrics_out_dir="$path/${sample}_FEICP/${sample}_sam_metrics"
    cov_out_dir="$path/${sample}_FEICP/${sample}_coverage"
    mkdir -p $cov_out_dir && cd $_
    bwa_bam_file="$metrics_out_dir/${sample}_bwa.bam"
    bwa_cov_file="$cov_out_dir/${sample}_bwa.bw"
    star_bam_file="$STAR_out_dir/${sample}_Aligned.sortedByCoord.out.bam"
    star_cov_file="$cov_out_dir/${sample}_STAR.bw"
    sambamba-0.8.1-linux-amd64-static index -t $ncpus $bwa_bam_file
    sambamba-0.8.1-linux-amd64-static index -t $ncpus $star_bam_file
    bamCoverage \
    --bam $bwa_bam_file \
    --outFileName $bwa_cov_file \
    --binSize 1 \
    --numberOfProcessors $ncpus \
    --normalizeUsing CPM
    
    bamCoverage \
    --bam $star_bam_file \
    --outFileName $star_cov_file \
    --binSize 1 \
    --numberOfProcessors $ncpus \
    --normalizeUsing CPM
}


function count_splice(){
    # input
    sample=$1
    EIci_out_dir="$path/${sample}_FEICP/${sample}_EIci"
    STAR_out_dir="$path/${sample}_FEICP/${sample}_STAR_IR"
    EIci_bed="$EIci_out_dir/${sample}_EIciRNA.bed"
    bam_file="$STAR_out_dir/${sample}_Aligned.sortedByCoord.out.bam"
    
    # output
    splice_out_dir="$path/${sample}_FEICP/${sample}_splice_metrics"
    mkdir -p $splice_out_dir && cd $_
    intron_bed="$splice_out_dir/${sample}_EIci_intron.bed"
    bed_split="$splice_out_dir/${sample}_split.bed" # The bed file composed of reads which are splited mapped (The CIGAR contains 'N')
    bed_file="$splice_out_dir/${sample}.bed" # The bed file composed of all reads
    split_read2tx_bed="$splice_out_dir/${sample}_split_read2tx.bed" # The mapping between splited reads and transcript
    EE_count="$splice_out_dir/${sample}_EE_count.bed" # The count of exon-exon junction for introns
    EI_count="$splice_out_dir/${sample}_EI_count.bed" # The count of exon-intron junction for introns
    intron_count="$splice_out_dir/${sample}_intron_count.bed" # The count of reads within introns
    PIR_bed="$splice_out_dir/${sample}_intron_bed_PIR.bed" # The PIR of introns
    cov_bed="$splice_out_dir/${sample}_intron_coverage.bed" # The coverage of introns
    merge_bed="$splice_out_dir/${sample}_intron_stats.bed" # All the statistics of junction for introns
    EIciRNA_tidy_out="$splice_out_dir/${sample}_EIciRNA_tidy.bed" # The tidy EIciRNAs 
    
    cat $EIci_bed | awk 'BEGIN{OFS=FS="\t"}{if(NR>1){split($9,A,",");split($10,B,",");for(i in A){print $1,A[i],B[i],$7,$4,$6}}}' >$intron_bed
    
    samtools view -@ $ncpus -h -q 255 $bam_file | \
    awk '$0 ~ /^@/ || $6 ~ /N/' | \
    samtools view -@ $ncpus -b | bedtools bamtobed -i - -split | \
    LC_COLLATE=C sort -k1,1 -k2,2n --parallel=$ncpus -T ./ >$bed_split
    bedtools intersect -a $bed_split -b $tx_ref -f 1 -sorted -wa -wb -nonamecheck | \
    bedtools groupby -i - -g 1,2,3,4,5,6 -c 11  -o distinct >$split_read2tx_bed
    python $get_EE --read-tx $split_read2tx_bed --output $EE_count
    
    samtools view -@ $ncpus -h -q 255 $bam_file | \
    bedtools bamtobed -i - -split | \
    LC_COLLATE=C sort -k1,1 -k2,2n --parallel=$ncpus -T ./  >$bed_file
    bedtools map -a $exon_intron_ref -b $bed_file -f 1 -nonamecheck -g $genome -c 4,4 -o count_distinct,distinct | \
    awk 'BEGIN{OFS=FS="\t"}{if($7!=0){tmp=$4;$4=$5;$5=tmp;print}}' >$EI_count
    
    python $get_intron_midpoint $intron_ref | \
    sort -k1,1 -k2,2n | uniq | \
    bedtools map -a - -b $bed_file -f 0.1 -F 0.1 -e -nonamecheck -g $genome -c 4,4 -o count_distinct,distinct | awk '$7!=0' >$intron_count
    
    cat $intron_bed | LC_COLLATE=C sort -k1,1 -k2,2n --parallel=$ncpus -T ./ | \
    bedtools coverage -a - -b $bed_file -nonamecheck -g $genome -sorted | uniq >$cov_bed
    
    python $calculate_PIR \
    --ee-count $EE_count \
    --ei-count $EI_count \
    --intron-count $intron_count \
    --query-intron $intron_bed \
    --output $PIR_bed
    
    python $merge_PIR_cov $PIR_bed $cov_bed $merge_bed

    python $get_tidy $EIci_bed $merge_bed $EIciRNA_tidy_out
}

# Generate all necessary annotation files: transcript; exon; exon-intron junction; intron
function prepare_ref(){
    cat $gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="exon"){split($9,A,"\"");print $1,$4-1,$5,A[2],A[4],$7}}' >$exon_ref
    cat $gtf | awk 'BEGIN{OFS=FS="\t"}{if($3=="transcript"){split($9,A,"\"");print $1,$4-1,$5,A[2],A[4],$7}}' >$tx_ref
    cat $exon_ref | sort -k1,1 -k2,2n | uniq | python $make_EI >$exon_intron_ref
    python $make_intron --gtf $gtf --output $intron_ref
    cat $star_index/chrNameLength.txt | sort -k1,1 >$genome
}


# The directory storing all source code
src_dir="/home/shange/yangyan2015/Galadriel/scripts"
find_EIci_py="$src_dir/find_EIci_PE_v1.4.py"
make_intron="$src_dir/make_intron.py"
make_EI="$src_dir/make_EI.py"
get_EE="$src_dir/get_EE.py"
get_intron_midpoint="$src_dir/get_intron_midpoint.py"
calculate_PIR="$src_dir/calculate_PIR.py"
merge_PIR_cov="$scripts/merge_PIR_cov.py"
get_tidy="$scripts/get_tidy_EIciRNA.py"
CIRI2="$src_dir/CIRI2.pl"

cd $path

# Make necessary reference
exon_ref="$path/exon.bed"
tx_ref="$path/transcript.bed"
exon_intron_ref="$path/exon_intron_junc.bed"
intron_ref="$path/intron.bed"
genome="$path/genome.txt"

prepare_ref


export -f find_EIci
export -f STAR_IR
export -f sam_metrics
export -f compute_cov
export -f count_splice


cat $sample_list | xargs -I {} mkdir -p $path/{}_FEICP
cat $sample_list | xargs -P $nth -I {} bash -c 'find_EIci {}'
cat $sample_list | xargs -P $nth -I {} bash -c 'STAR_IR {}'
cat $sample_list | xargs -P $nth -I {} bash -c 'sam_metrics {}'
cat $sample_list | xargs -P $nth -I {} bash -c 'compute_cov {}'
cat $sample_list | xargs -P $nth -I {} bash -c 'count_splice {}'
