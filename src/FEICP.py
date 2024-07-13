#!/usr/bin/env python

import argparse
import gzip
import logging
import os
import re
import signal
import subprocess
import sys
import time
from collections import defaultdict

logging.basicConfig(
    format='[%(asctime)s %(levelname)s] %(message)s',    stream=sys.stdout
)
log = logging.getLogger(__name__)


# string/system functions
def get_ticks():
    """Returns ticks.
        - Python3: Use time.perf_counter().
        - Python2: Use time.time().
    """
    return getattr(time, 'perf_counter', getattr(time, 'time'))()


def mkdir_p(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

# This function is from ENCODE project to call Linux in Python


def run_shell_cmd(cmd):
    p = subprocess.Popen(
        ['/bin/bash', '-o', 'pipefail'],  # to catch error in pipe
        stdin=subprocess.PIPE,        stdout=subprocess.PIPE,        stderr=subprocess.PIPE,        universal_newlines=True,        preexec_fn=os.setsid)  # to make a new process with a new PGID
    pid = p.pid
    pgid = os.getpgid(pid)
    log.info('run_shell_cmd: PID={}, PGID={}, CMD={}'.format(pid, pgid, cmd))
    t0 = get_ticks()
    stdout, stderr = p.communicate(cmd)
    rc = p.returncode
    t1 = get_ticks()
    err_str = (
        'PID={pid}, PGID={pgid}, RC={rc}, DURATION_SEC={dur:.1f}\n'
        'STDERR={stde}\nSTDOUT={stdo}'
    ).format(
        pid=pid, pgid=pgid, rc=rc, dur=t1 - t0, stde=stderr.strip(), stdo=stdout.strip()
    )
    if rc:
        # kill all child processes
        try:
            os.killpg(pgid, signal.SIGKILL)
        except:
            pass
        finally:
            raise Exception(err_str)
    else:
        log.info(err_str)
    return stdout.strip('\n')


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Getting potential EIciRNAs form RNA-seq')
    parser.add_argument('--fastqs', '-i', nargs='+', required=True,
                        type=str, help='List of FASTQs (R1 and R2). FASTQs must be compressed with gzip (with .gz).')
    parser.add_argument('--paired-end', action="store_true",
                        help="Paired-end FASTQs")
    parser.add_argument(
        '--sam', type=str, help='Sam-formatted file from BWA mapping used to get circular RNAs')
    parser.add_argument(
        '--circ', type=str, help='Circular RNAs file from CIRI2 used to find EIciRNAs')
    parser.add_argument('--genome', '-f', required=True,
                        type=str, help='The genome fasta file')
    parser.add_argument('--gtf', '-a', required=True,
                        type=str, help='The annotation gtf file')
    parser.add_argument('--bwa', type=str, default='bwa',
                        help='The path to exectuable bwa')
    parser.add_argument('--bwa-index', type=str, default='/home/shange/yangyan2015/Galadriel/database/bwaindex/GRCh38/GRCh38_bwa',
                        help='The path to the bwa index of reference genome fasta')
    parser.add_argument('--bowtie', type=str, default='bowtie',
                        help='The path to exectuable bowtie')
    parser.add_argument('--bowtie-build', type=str, default='bowtie-build',
                        help='The path to exectuable bowtie-build')
    parser.add_argument('--hisat2', type=str, default='hisat2',
                        help='The path to exectuable hisat2')
    parser.add_argument('--hisat2-build', type=str, default='hisat2-build',
                        help='The path to exectuable hisat2-build')
    parser.add_argument('--seqtk', type=str, default='seqtk',
                        help='The path to exectuable seqtk')
    parser.add_argument('--samtools', type=str, default='samtools',
                        help='The path to exectuable samtools')
    parser.add_argument('--bedtools', type=str, default='bedtools',
                        help='The path to exectuable bedtools')
    parser.add_argument('--mismatch', '-m', default=2, type=int,
                        help='The maximum number of mismatch permitted during bowtie mapping')
    parser.add_argument('--min-overlap',  default=5, type=int,
                        help='The minimum overhang between read and reference intron')
    parser.add_argument('--src-dir', '-c', type=str, default='/home/shange/yangyan2015/Galadriel/software/CIRI_v2.0.6',
                        help='The directory storing source code files of CIRI2')
    parser.add_argument('--output', '-o', type=str,
                        default='./', help='The output directory')
    parser.add_argument('--sample', '-n', required=True, type=str,
                        help='The name of sample used to generate the temp file')
    parser.add_argument('--threads', '-@', type=int,
                        default=16, help='The number of CPUs to use')
    parser.add_argument('--log-level', default='INFO',
                        choices=['NOTSET', 'DEBUG', 'INFO',
                                 'WARNING', 'CRITICAL', 'ERROR',
                                 'CRITICAL'],
                        help='Log level')
    args = parser.parse_args()

    # check if fastqs have corrrect dimension
    if args.paired_end and len(args.fastqs) != 2:
        raise argparse.ArgumentTypeError('Need 2 fastqs for paired end.')
    if not args.paired_end and len(args.fastqs) != 1:
        raise argparse.ArgumentTypeError('Need 1 fastq for single end.')

    log.setLevel(args.log_level)
    log.info(sys.argv)
    return args


def open_file(infile='', mode='r'):
    if re.search('.gz$', infile):
        mode += 't'
        return gzip.open(infile, mode=mode)
    else:
        return open(infile, mode=mode)


def parse_fasta(infile=''):
    name2seq = {}
    with open_file(infile) as fh:
        for line in fh:
            line = line.strip()
            if re.search('^>', line):
                name = line.split(' ')[0][1:]
                name2seq[name] = []
            else:
                # avoid the poor performance of repeatedly concatenating (immutable) Python strings.
                name2seq[name].append(line)

    for name in name2seq:
        name2seq[name] = ''.join(name2seq[name])

    return name2seq


class Circ(object):
    def __init__(self, chr='chr1', start=0, end=1, gene_name=[], junc_count=0, strand='+', reads_list=[], intron_list=[]):
        self.chr = str(chr)
        self.start = int(start)
        self.end = int(end)
        self.gene_name = gene_name
        self.junc_count = int(junc_count)
        self.strand = strand
        self.reads_list = reads_list
        self.intron_list = intron_list

    def get_seq(self, chr2seq={}):
        # It's uncessary to distinguish the strand because mapping will handle it
        return chr2seq[self.chr][self.start:self.end]


class EIci(Circ):
    """Represent aspects of a circRNA, specific to EIciRNA"""

    def __init__(self, chr='chr1', start=0, end=1, gene_name=[], junc_count=0, strand='+', reads_list=[], intron_list=[], tx_score=0, tx_list=[]):
        super().__init__(chr, start, end, gene_name,
                         junc_count, strand, reads_list, intron_list)
        self.tx_id = tx_list
        self.tx_score = tx_score

    def write_to_file(self, fh):
        fh.write('\t'.join([self.chr, str(self.start), str(self.end), ','.join(self.gene_name), str(
            self.junc_count), self.strand, ','.join(self.tx_id), str(self.tx_score), ','.join([str(i[0]) for i in self.intron_list]), ','.join([str(i[1]) for i in self.intron_list]), ','.join(self.reads_list)]) + '\n')


def get_circ_seq(sample='', name2circ={}, chr2seq={}, out_dir=''):
    prefix = os.path.join(out_dir, sample)
    out_fasta = '{}_circ.fasta'.format(prefix)
    with open(out_fasta, 'w') as output:
        for name in name2circ:
            output.write(
                '\n'.join(['>'+name, name2circ[name].get_seq(chr2seq)]) + '\n')

    return out_fasta


class Block(object):
    def __init__(self, chr='chr1', start=0, end=1, tx_id='', gene_id='', strand='+'):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.tx_id = tx_id
        self.gene_id = gene_id
        self.strand = strand
        self.len = self.get_len()

    def get_len(self):
        return self.end - self.start

    def format_to_string(self):
        return '\t'.join([self.chr, str(self.start), str(self.end), self.tx_id, self.gene_id, self.strand])

    def write_to_file(self, fh):
        fh.write(self.format_to_string() + '\n')


def parse_str(string="", sep='\"'):
    key2value = defaultdict(str)
    str_list = string.split(sep)
    str_len = len(str_list)
    for i in range(0, int(str_len/2)*2-1, 2):
        key2value[str_list[i].strip(';').strip()] = str_list[i+1]

    return key2value


def get_exons(gtf_file):
    tx2exons = defaultdict(list)
    with open_file(gtf_file) as fh:
        for line in fh:
            line = line.strip()
            if re.match('^#', line):
                continue
            lineL = line.split('\t')
            if lineL[2] == "exon":
                bio = parse_str(lineL[8])
                chr, start, end, tx_id, gene_id, strand = lineL[0], int(
                    lineL[3]) - 1, lineL[4], bio["transcript_id"], bio["gene_id"], lineL[6]
                tx2exons[tx_id].append(
                    Block(chr, start, end, tx_id, gene_id, strand))

    for tx in tx2exons:
        tx2exons[tx].sort(key=lambda x: (x.chr, x.start))

    return tx2exons


def get_exons_introns(gtf_file):
    tx2exons = get_exons(gtf_file)
    tx2introns = defaultdict(list)
    for tx in tx2exons:
        for i in range(0, len(tx2exons[tx]) - 1, 1):
            tx2introns[tx].append(Block(tx2exons[tx][i].chr, tx2exons[tx][i].end, tx2exons[tx]
                                        [i+1].start, tx2exons[tx][i].tx_id, tx2exons[tx][i].gene_id, tx2exons[tx][i].strand))

    intron_list = []
    for tx in tx2introns:
        intron_list += tx2introns[tx]
    intron_list.sort(key=lambda x: (x.chr, x.start, x.end))

    return tx2exons, intron_list


def get_all_annotations(ref_gtf='', out_dir=''):
    basename = os.path.basename(ref_gtf)
    prefix = strip_ext_gtf(basename)
    prefix = os.path.join(out_dir, prefix)
    exon_tx_file = '{}_tx_exon.bed'.format(prefix)
    intron_tx_file = '{}_tx_intron.bed'.format(prefix)
    intron_gene_file = '{}_gene_intron.bed'.format(prefix)
    gene2intron = defaultdict(list)
    tx2intron = defaultdict(list)
    gene2tx = defaultdict(list)
    tx2gene = {}

    tx2exons, intron_list = get_exons_introns(ref_gtf)
    with open(intron_tx_file, 'w') as output:
        for intron in intron_list:
            intron.write_to_file(output)

    with open(exon_tx_file, 'w') as output:
        for tx in tx2exons:
            for exon in tx2exons[tx]:
                exon.write_to_file(output)

    cmd0 = 'cat {} | awk \'BEGIN{{OFS=FS="\\t"}}{{$4=".";print}}\' | uniq >{}'
    cmd0 = cmd0.format(intron_tx_file, intron_gene_file)
    run_shell_cmd(cmd0)

    with open_file(intron_tx_file) as fh:
        for line in fh:
            line = line.strip()
            if re.search('^#', line):
                continue
            lineL = line.split('\t')
            tx2intron[lineL[3]].append(
                Block(lineL[0], lineL[1], lineL[2], lineL[3], lineL[4], lineL[5]))
            gene2tx[lineL[4]].append(lineL[3])
            tx2gene[lineL[3]] = lineL[4]

    with open_file(intron_gene_file) as fh:
        for line in fh:
            line = line.strip()
            if re.search('^#', line):
                continue
            lineL = line.split('\t')
            gene2intron[lineL[4]].append(
                Block(lineL[0], lineL[1], lineL[2], lineL[3], lineL[4], lineL[5]))

    return tx2exons, exon_tx_file, intron_tx_file, tx2intron, gene2intron, gene2tx, tx2gene


# The dict indicating complement
COMPLEMENT = {
    'a': 't',
    't': 'a',
    'c': 'g',
    'g': 'c',
    'k': 'm',
    'm': 'k',
    'r': 'y',
    'y': 'r',
    's': 's',
    'w': 'w',
    'b': 'v',
    'v': 'b',
    'h': 'd',
    'd': 'h',
    'n': 'n',
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'K': 'M',
    'M': 'K',
    'R': 'Y',
    'Y': 'R',
    'S': 'S',
    'W': 'W',
    'B': 'V',
    'V': 'B',
    'H': 'D',
    'D': 'H',
    'N': 'N',
}


def complement(s):
    return "".join([COMPLEMENT[x] for x in s])


# return the reversed complement sequences
def rev_comp(seq):
    return complement(seq)[::-1]


def remove_repeat(a_list=[]):
    b_list = sorted(a_list)
    c_list = []
    prev_elem = ''
    for curr_elem in b_list:
        if curr_elem != prev_elem:
            c_list.append(curr_elem)
        prev_elem = curr_elem

    return c_list


# Getting the matched number of bases indicated
# by a CIGAR in SAM-formatted file
def get_match(string):
    match_num = 0
    pattern_char = re.compile('M|I|D|N|S|H|P|=|X')
    pattern_num = re.compile(r'\d')
    list_by_char = [i for i in pattern_char.split(string) if i]
    list_by_num = [i for i in pattern_num.split(string) if i]
    for i in range(0, len(list_by_num)):
        if list_by_num[i] == 'M':
            match_num += int(list_by_char[i])
    return match_num


def get_read_length(fastq):
    line_num = 0
    sample_reads_num = 100000
    analyzed_reads_num = 0
    readLen2count = defaultdict(int)
    with open_file(fastq) as fh:
        for line in fh:
            if line_num % 4 == 1:
                readLen2count[len(line.strip())] += 1
                analyzed_reads_num += 1
            if analyzed_reads_num >= sample_reads_num:
                break
            line_num += 1
    max_occurrence_num = 0
    read_len = 0
    for length in readLen2count:
        if readLen2count[length] > max_occurrence_num:
            max_occurrence_num = readLen2count[length]
            read_len = length
    return read_len


def get_circ(sample, src_dir='', ref_fa='', sam_file='', ref_gtf='', out_dir='./', nThread=8):
    CIRI2 = os.path.join(src_dir, 'CIRI2.pl')
    out_prefix = os.path.join(out_dir, sample)
    out_circ = '{}_CIRI2_circ.txt'.format(out_prefix)
    out_log = '{}_CIRI2.log'.format(out_prefix)
    cmd = 'perl {} -0 -I {} -O {} -F {} -A {} -G {} -T {}'
    cmd = cmd.format(CIRI2, sam_file, out_circ,
                     ref_fa, ref_gtf, out_log, nThread)

    run_shell_cmd(cmd)
    return out_circ


def parse_circ(circ_txt):
    name2circ = {}
    with open_file(circ_txt) as fh_circ:
        circ_num = 0
        for line in fh_circ:
            line = line.strip()
            if line.startswith('circRNA_ID'):
                continue
            lineL = line.split('\t')
            chr, start, end, junc_count, circ_type, gene_id, strand, junc_reads_ID = lineL[
                1], int(lineL[2]) - 1, lineL[3], lineL[4], lineL[8], lineL[9][:-1].split(','), lineL[10], lineL[11][:-1]
            if circ_type == "exon":
                reads_list = junc_reads_ID.split(',')
                circ = Circ(chr, start, end, gene_id,
                            junc_count, strand, reads_list)
                circ_name = 'circ_' + str(circ_num)
                name2circ[circ_name] = circ
                circ_num += 1
    return name2circ


def find_fastq_type(fastq1='', fastq2=''):
    with open_file(fastq1) as fh1, open_file(fastq2) as fh2:
        for line in fh1:
            fq1_name = re.sub(r'^@', '', line.strip().split(' ')[0])
            break
        for line in fh2:
            fq2_name = re.sub(r'^@', '', line.strip().split(' ')[0])
            break

        if fq1_name == fq2_name:
            return 'new'
        elif re.sub(r'/1$', '', fq1_name) == re.sub(r'/2$', '', fq2_name):
            return 'old'


def get_junc_reads(sample='', name2circ={}, fastq1='', fastq2='', out_dir='./', seqtk='seqtk'):
    circ_junc_reads_list = []
    for name in name2circ:
        circ_junc_reads_list += name2circ[name].reads_list
    circ_junc_reads_list = remove_repeat(circ_junc_reads_list)
    prefix = os.path.join(out_dir, sample)

    out_fastq1 = '{}_BSJ_1.fastq.gz'.format(prefix)
    out_fastq2 = '{}_BSJ_2.fastq.gz'.format(prefix)
    cmd0 = '{} subseq {} {} | gzip -nc >{}'
    cmd1 = '{} subseq {} {} | gzip -nc >{}'
    if find_fastq_type(fastq1, fastq2) == "new":
        bsj_reads_file = '{}_BSJ_reads.txt'.format(prefix)
        with open(bsj_reads_file, 'w') as fh:
            fh.write('\n'.join(circ_junc_reads_list))
        cmd0 = cmd0.format(seqtk, fastq1, bsj_reads_file, out_fastq1)
        cmd1 = cmd1.format(seqtk, fastq2, bsj_reads_file, out_fastq2)
    elif find_fastq_type(fastq1, fastq2) == "old":
        bsj_reads_file_r1 = '{}_BSJ_reads_R1.txt'.format(prefix)
        bsj_reads_file_r2 = '{}_BSJ_reads_R2.txt'.format(prefix)
        with open(bsj_reads_file_r1, 'w') as fh1, open(bsj_reads_file_r2, 'w') as fh2:
            fh1.write('\n'.join([i + '/1' for i in circ_junc_reads_list]))
            fh2.write('\n'.join([i + '/2' for i in circ_junc_reads_list]))
        cmd0 = cmd0.format(seqtk, fastq1, bsj_reads_file_r1, out_fastq1)
        cmd1 = cmd1.format(seqtk, fastq2, bsj_reads_file_r2, out_fastq2)

    run_shell_cmd(cmd0)
    run_shell_cmd(cmd1)

    return out_fastq1, out_fastq2


def build_bwa(bwa='bwa', ref_fa='', out_dir='./'):
    prefix = strip_ext_fasta(os.path.basename(ref_fa))
    bwa_index = os.path.join(out_dir, '{}_bwa'.format(prefix))
    cmd = '{} index -p {} {}'
    cmd = cmd.format(bwa, bwa_index, ref_fa)
    run_shell_cmd(cmd)

    return bwa_index


def exist_bwa_index(bwa_index=''):
    suffix_list = ['amb', 'ann', 'bwt', 'pac', 'sa']
    for suffix in suffix_list:
        if not os.path.exists('{}.{}'.format(bwa_index, suffix)):
            return False
    return True


def bwa_pe(sample='', fastq1='', fastq2='', bwa='bwa', bwa_index='hg38_bwa', out_dir='./', nThread=8):
    sam_file = os.path.join(out_dir, '{}_bwa.sam'.format(sample))
    log_file = os.path.join(out_dir, '{}_bwa.log'.format(sample))
    cmd = '{} mem -T 19 -t {} {} {} {} 1> {} 2> {}'
    cmd = cmd.format(bwa, nThread, bwa_index, fastq1,
                     fastq2, sam_file, log_file)

    run_shell_cmd(cmd)
    return sam_file


def make_bt_idx(ref_fasta, bowtie_build='bowtie-build', nthread=8, out_dir='./'):
    basename = os.path.basename(ref_fasta)
    prefix = strip_ext_fasta(basename)
    bt_index = os.path.join(out_dir, '{}_bt'.format(prefix))
    cmd = '{} --threads {} {} {}'.format(bowtie_build,
                                         nthread, ref_fasta, bt_index)
    run_shell_cmd(cmd)

    return bt_index


def make_ht2_idx(ref_fasta, hisat2_build='hisat2-build', nthread=8, out_dir='./'):
    basename = os.path.basename(ref_fasta)
    prefix = strip_ext_fasta(basename)
    ht2_idx = os.path.join(out_dir, '{}_ht2'.format(prefix))
    cmd = '{} -p {} {} {}'
    cmd = cmd.format(hisat2_build, nthread, ref_fasta, ht2_idx)
    run_shell_cmd(cmd)

    return ht2_idx


def strip_ext_fasta(fasta_file):
    return re.sub(r'\.(fa|fasta|Fa|Fasta)$', '', str(fasta_file))


def strip_ext_sam(sam_file):
    return re.sub(r'\.(sam|Sam)$', '', str(sam_file))


def strip_ext_bed(bed_file):
    return re.sub(r'\.(bed|Bed)$', '', str(bed_file))


def strip_ext_fastq(fastq_file):
    return re.sub(r'\.(fq|fastq|Fq|Fastq)\.gz$', '', str(fastq_file))


def strip_ext_gtf(gtf_file):
    return re.sub(r'\.(gtf|GTF)$', '', str(gtf_file))


# Map junction reads file against circular RNA fasta file
# I use bowtie 1.3.0 to map. Make sure that do not use bowtie out of date because it does bad things.
def bowtie_se(fastq='', threads=8, mismatch=2, bt_index='hg38_bt', bowtie='', out_dir='./'):
    prefix = strip_ext_fastq(os.path.basename(fastq))
    out_prefix = os.path.join(out_dir, prefix)
    out_sam_file = '{}_bowtie.sam'.format(out_prefix)
    out_log_file = '{}_bowtie.log'.format(out_prefix)

    cmd = '{} --no-unal -p {} -v {} -x {} {} -S {} 2>{}'
    cmd = cmd.format(bowtie, threads, mismatch, bt_index,
                     fastq, out_sam_file, out_log_file)
    run_shell_cmd(cmd)
    return out_sam_file


# Map junction reads file against circular RNA fasta file
# I use hisat2 to map to be aware of splice junction spanning reads
def hisat2_se(fastq='', threads=8, mismatch=2, ht2_index='hg38_ht2', hisat2='hisat2', out_dir='./'):
    prefix = strip_ext_fastq(os.path.basename(fastq))
    out_prefix = os.path.join(out_dir, prefix)
    out_sam_file = '{}_hisat2.sam'.format(out_prefix)
    out_log_file = '{}_hisat2.log'.format(out_prefix)

    cmd = '{} -p {} -x {} -U {} -S {} 2>{}'
    cmd = cmd.format(hisat2, threads, ht2_index,
                     fastq, out_sam_file, out_log_file)
    run_shell_cmd(cmd)

    return out_sam_file


def sam2bed_smart(sam_file='', samtools='samtools', bedtools='bedtools', threads=10, out_dir='./'):
    basename = os.path.basename(sam_file)
    prefix = strip_ext_sam(basename)
    out_bed = os.path.join(out_dir, '{}.bed'.format(prefix))

    cmd = '{} view -@ {} -Sb {} | '
    cmd += '{} bamtobed -split -i - >{}'
    cmd = cmd.format(samtools, threads, sam_file, bedtools, out_bed)
    run_shell_cmd(cmd)

    return out_bed


class Interval(object):
    def __init__(self, start=0, end=1):
        self.start = int(start)
        self.end = int(end)


def Have_overlap(a=Interval(0, 1), b=Interval(0, 1), min_over=10):
    if min(a.end, b.end) - max(a.start, b.start) >= min_over:
        return True
    return False


def reshape_dict(a_dict={}):
    b_dict = defaultdict(list)
    for key, value in a_dict.items():
        b_dict[value].append(key)

    return b_dict


def map_circ_tx(name2circ={}, tx2gene={}, gene2tx={}, tx2exon={}):
    name2EIci = {}
    for name in name2circ:
        circ = name2circ[name]
        score = 0
        pot_tx_list = []
        pot_gene_list = []
        tx2score = {}
        for gene in circ.gene_name:
            for tx in gene2tx[gene]:
                circ_start_map = False
                circ_end_map = False
                for exon in tx2exon[tx]:
                    if exon.start == circ.start:
                        circ_start_map = True
                    elif exon.end == circ.end:
                        circ_end_map = True
                if circ_start_map and circ_end_map:
                    tx2score[tx] = 2
                elif circ_start_map or circ_end_map:
                    tx2score[tx] = 1
                else:
                    tx2score[tx] = 0
        score2tx = reshape_dict(tx2score)
        if score2tx[2]:
            score = 2
            pot_tx_list = score2tx[2]
            pot_gene_list = remove_repeat([tx2gene[tx] for tx in pot_tx_list])
        elif score2tx[1]:
            score = 1
            pot_tx_list = score2tx[1]
            pot_gene_list = remove_repeat([tx2gene[tx] for tx in pot_tx_list])
        else:
            score = 0
            pot_tx_list = ['NA']
            pot_gene_list = circ.gene_name
        name2EIci[name] = EIci(circ.chr, circ.start, circ.end, pot_gene_list,
                               circ.junc_count, circ.strand, circ.reads_list, [], score, pot_tx_list)

    return name2EIci


def calculate_cov(intron_bed='', bam_file='', out_dir='./', bedtools='bedtools'):
    basename = os.path.basename(intron_bed)
    prefix = strip_ext_bed(basename)
    out_cov = os.path.join(out_dir, '{}_cov.bed'.format(prefix))

    cmd = '{} coverage -a {} -b {} >{}'
    cmd = cmd.format(bedtools, intron_bed, bam_file, out_cov)
    run_shell_cmd(cmd)

    return out_cov


# attempt to find introns of EIciRNAs after identifying potential transcripts
def map_circ_intron(name2EIci={}, gene2intron={}, tx2intron={}):
    for name in name2EIci:
        EIci = name2EIci[name]
        if EIci.tx_score == 2 or EIci.tx_score == 1:
            for tx in EIci.tx_id:
                for intron in tx2intron[tx]:
                    if Have_overlap(Interval(intron.start, intron.end), Interval(EIci.start, EIci.end), intron.len):
                        EIci.intron_list.append(intron)
        else:
            for gene in EIci.gene_name:
                for intron in gene2intron[gene]:
                    if Have_overlap(Interval(intron.start, intron.end), Interval(EIci.start, EIci.end), intron.len):
                        EIci.intron_list.append(intron)

    return name2EIci


def merge_sam(sam_file1='', sam_file2='', out_dir=''):
    basename = os.path.basename(sam_file1)
    prefix = re.sub(r'_1', '', strip_ext_sam(basename))
    out_sam = os.path.join(out_dir, '{}.sam'.format(prefix))
    with open_file(sam_file1) as fh1, open_file(sam_file2) as fh2, open(out_sam, 'w') as output:
        for line in fh1:
            output.write(line)
        for line in fh2:
            if re.search('^@', line):
                continue
            output.write(line)

    return out_sam


def fix_sam(sample='', sam_bt='', sam_ht2='', out_dir='./', mismatch=2, read_len=100):
    prefix = os.path.join(out_dir, sample)
    sam_fix = '{}_fix.sam'.format(prefix)
    bt_map_reads = []
    ht2_map_reads = []
    with open_file(sam_bt) as fh_bt, open_file(sam_ht2) as fh_ht2, open(sam_fix, 'w') as output:
        for line in fh_bt:
            output.write(line)
            if re.search('^@', line):
                continue
            line = line.strip()
            lineL = line.split('\t')
            bt_map_reads.append(lineL[0])
        for line in fh_ht2:
            if re.search('^@', line):
                continue
            line = line.strip()
            lineL = line.split('\t')
            if lineL[1] == 4 or re.search(r'N', lineL[5]) or get_match(lineL[5]) < read_len - mismatch:
                continue
            ht2_map_reads.append(lineL[0])

        ht2_map_uniq_reads = list(set(ht2_map_reads) - set(bt_map_reads))
        fh_ht2.seek(0)
        for line in fh_ht2:
            if re.search('^@', line):
                continue
            line = line.strip()
            lineL = line.split('\t')
            if lineL[1] == 4 or re.search(r'N', lineL[5]) or get_match(lineL[5]) < read_len - mismatch:
                continue
            if lineL[0] in ht2_map_uniq_reads:
                output.write(line + '\n')

    return sam_fix


BIT2DES = {
    1: 'PAIRED',
    2: 'PROPER_PAIR',
    4: 'UNMAP',
    8: 'MUNMAP',
    16: 'REVERSE',
    32: 'MREVERSE',
    64: 'READ1',
    128: 'READ2',
    256: 'SECONDARY',
    512: 'QCFAIL',
    1024: 'DUP',
    2048: 'SUPPLEMENTARY'
}


def decode_flag(flag=99):
    flag = int(flag)
    if flag == 0:
        return ['FORWARD']
    bin_str = re.sub(r'0b', '', bin(flag))
    return [BIT2DES[2**i] for i in range(len(bin_str)) if (1 << i) & flag]


def get_strand(flag=99):
    flag_DES = decode_flag(flag)
    strand = '+'
    if 'REVERSE' in flag_DES:
        strand = '-'
    return strand


def sam2bed(sam_file='', out_dir='./'):
    basename = os.path.basename(sam_file)
    prefix = strip_ext_sam(basename)
    output_bed = os.path.join(out_dir, '{}.bed'.format(prefix))
    with open_file(sam_file) as fh, open(output_bed, 'w') as output:
        for line in fh:
            if re.search('^@', line):
                continue
            line = line.strip()
            lineL = line.split('\t')
            query_name, flag, ref_name, map_pos, mapq, cigar = lineL[
                0], int(lineL[1]), lineL[2], int(lineL[3])-1, lineL[4], lineL[5]
            if 'UNMAP' in decode_flag(flag):
                continue
            strand = get_strand(flag)
            map_len = get_match(cigar)
            output.write('\t'.join([ref_name, str(map_pos), str(map_pos +
                                                                map_len), query_name, mapq, strand]) + '\n')

    return output_bed


def get_genomic_pos(name2circ={}, bed_file='', out_dir=''):
    basename = os.path.basename(bed_file)
    prefix = strip_ext_bed(basename)
    out_bed = os.path.join(out_dir, '{}_genome.bed'.format(prefix))
    with open_file(bed_file) as fh, open(out_bed, 'w') as output:
        bed_list = []
        for line in fh:
            line = line.strip()
            ref_name, map_pos, map_end, query_name, mapq, strand = line.split('\t')[
                :6]
            # because the record of read name in sam file of bowtie and HISAT2
            query_name = query_name.split('/')[0]
            map_pos = int(map_pos)
            map_end = int(map_end)
            circ = name2circ[ref_name]
            bed_list.append('\t'.join([circ.chr, str(
                circ.start+map_pos), str(circ.start+map_end), query_name, mapq, strand]))
        bed_list = remove_repeat(bed_list)
        output.write('\n'.join(bed_list))

    return out_bed


def map_reads_intron(map_bed='', tx_intron_bed='', out_dir='', bedtools='bedtools', min_overlap_base=10):
    basename = os.path.basename(map_bed)
    prefix = strip_ext_bed(basename)
    out_bed = os.path.join(out_dir, '{}_intron_map.bed'.format(prefix))

    cmd1 = '{} intersect -a {} -b {} -wo >{}'
    cmd1 = cmd1.format(bedtools, map_bed, tx_intron_bed, out_bed)
    run_shell_cmd(cmd1)

    read2intron = defaultdict(list)
    with open_file(out_bed) as fh:
        for line in fh:
            line = line.strip()
            lineL = line.split('\t')
            if int(lineL[-1]) >= min_overlap_base:
                read2intron[lineL[3]].append(
                    Block(lineL[6], lineL[7], lineL[8], lineL[9], lineL[10], lineL[11]))

    return read2intron


def map_reads_exon(map_bed='', tx_exon_bed='', out_dir='', bedtools='bedtools'):
    basename = os.path.basename(map_bed)
    prefix = strip_ext_bed(basename)
    out_bed = os.path.join(out_dir, '{}_exon_map.bed'.format(prefix))

    cmd = '{} intersect -a {} -b {} -wo >{}'
    cmd = cmd.format(bedtools, map_bed, tx_exon_bed, out_bed)
    run_shell_cmd(cmd)

    read2exon = defaultdict(list)
    with open_file(out_bed) as fh:
        for line in fh:
            line = line.strip()
            lineL = line.split('\t')
            if abs(int(lineL[-1]) - (int(lineL[2]) - int(lineL[1]))) <= 2:
                read2exon[lineL[3]].append(
                    Block(lineL[6], lineL[7], lineL[8], lineL[9], lineL[10], lineL[11]))

    return read2exon


def have_common_elem(a_list=[], b_list=[]):
    c_list = a_list + b_list
    if len(remove_repeat(c_list)) < len(c_list):
        return True
    return False


def map_EIci_intron_smart(name2EIci_tmp={}, read2intron={}, read2exon={}):
    name2EIci_final = {}
    for name in name2EIci_tmp:
        EIci_tmp = name2EIci_tmp[name]
        EIci_reads_list = []
        EIci_intron_list = []
        EIci_gene_list = []
        EIci_tx_list = []
        if EIci_tmp.tx_score:
            for read in EIci_tmp.reads_list:
                # It is incomprehensible that I made this stupid mistake here
                if have_common_elem(remove_repeat([exon.tx_id for exon in read2exon[read]]), EIci_tmp.tx_id):
                    continue
                for intron in read2intron[read]:
                    for eici_tmp_intron in EIci_tmp.intron_list:
                        if intron.start == eici_tmp_intron.start and intron.end == eici_tmp_intron.end:
                            EIci_reads_list.append(read)
                            EIci_intron_list.append((intron.start, intron.end))
                            EIci_gene_list.append(intron.gene_id)
                            EIci_tx_list.append(intron.tx_id)
        else:
            for read in EIci_tmp.reads_list:
                if read2intron[read]:
                    EIci_reads_list.append(read)
                    EIci_intron_list += [(intron.start, intron.end)
                                         for intron in read2intron[read]]
                    EIci_gene_list += [intron.gene_id for intron in read2intron[read]]
            EIci_tx_list = ['NA']

        if EIci_reads_list:
            EIci_reads_list = remove_repeat(EIci_reads_list)
            EIci_intron_list = remove_repeat(EIci_intron_list)
            EIci_gene_list = remove_repeat(EIci_gene_list)
            EIci_tx_list = remove_repeat(EIci_tx_list)
            EIci_name = 'EIci_{}'.format(name.split('_')[1])
            name2EIci_final[EIci_name] = EIci(EIci_tmp.chr, EIci_tmp.start, EIci_tmp.end, EIci_gene_list, len(
                EIci_reads_list), EIci_tmp.strand, EIci_reads_list, EIci_intron_list, EIci_tmp.tx_score, EIci_tx_list)

    return name2EIci_final


def find_EIci_from_fastq(sample='', fastq1='', fastq2='', bwa_index='', ref_fa='', ref_gtf='', threads=10, mismatch=2, out_dir='./', src_dir='', bwa='', bowtie='bowtie', bowtie_build='bowtie-build', hisat2='', hisat2_build='', seqtk='', samtools='samtools', bedtools='bedtools', min_overlap=10):
    log.info('Mapping {}, {} against {} using bwa mem'.format(
        fastq1, fastq2, bwa_index))
    sam_file = bwa_pe(sample, fastq1, fastq2, bwa, bwa_index, out_dir, threads)
    name2EIci = find_EIci_from_sam(sample, fastq1, fastq2, ref_fa, ref_gtf, sam_file, threads,
                                   mismatch, out_dir, src_dir, bowtie, bowtie_build, hisat2, hisat2_build, seqtk, samtools, bedtools, min_overlap)
    return name2EIci


def find_EIci_from_sam(sample='', fastq1='', fastq2='', ref_fa='', ref_gtf='', sam_file='', threads=10, mismatch=2, out_dir='./', src_dir='', bowtie='bowtie', bowtie_build='bowtie-build', hisat2='', hisat2_build='', seqtk='', samtools='samtools',  bedtools='bedtools', min_overlap=10):
    log.info('Getting circular RNAs from {}, {} using CIRI2'.format(
        fastq1, fastq2))
    circ_file = get_circ(sample, src_dir, ref_fa,
                         sam_file, ref_gtf, out_dir, threads)
    name2EIci = find_EIci_from_CIRI2(sample, fastq1, fastq2, ref_fa, ref_gtf, circ_file,
                                     threads, mismatch, out_dir, bowtie, bowtie_build, hisat2, hisat2_build, seqtk, samtools, bedtools, min_overlap)
    return name2EIci


def find_EIci_from_CIRI2(sample='', fastq1='', fastq2='', ref_fa='', ref_gtf='', circ_file='', threads=10, mismatch=2, out_dir='./', bowtie='bowtie', bowtie_build='bowtie-build', hisat2='hisat2', hisat2_build='', seqtk='', samtools='samtools', bedtools='bedtools', min_overlap=10):
    print('Getting the length of reads')
    read_len = get_read_length(fastq1)

    print('Formatting circRNA file into dictionary')
    name2circ = parse_circ(circ_file)
    
    if not name2circ:
        sys.exit("No circRNA found!")

    print('Getting all the necessary annotations of genome for downstream analysis')
    tx2exons, exon_tx_file, intron_tx_file, tx2intron, gene2intron, gene2tx, tx2gene = get_all_annotations(
        ref_gtf, out_dir)

    print('Identifying all the possible transcripts from which the circRNA is derived')
    name2EIci_tmp = map_circ_tx(name2circ, tx2gene, gene2tx, tx2exons)

    print('Collecting all the potential introns of circRNAs according to the transcript')
    name2EIci_tmp = map_circ_intron(name2EIci_tmp, gene2intron, tx2intron)

    print('Formatting fasta file into dictionary')
    chr2seq = parse_fasta(ref_fa)

    print('Getting sequences of all circular RNAs')
    circ_fasta = get_circ_seq(sample, name2circ, chr2seq, out_dir)

    print('Building bowtie index of circular RNAs')
    bt_index = make_bt_idx(circ_fasta, bowtie_build, threads, out_dir)

    print('Building hisat2 index of circular RNAs')
    ht2_index = make_ht2_idx(circ_fasta, hisat2_build, threads, out_dir)

    print('Getting BSJ reads of all circular RNAs')
    junc_fastq1, junc_fastq2 = get_junc_reads(
        sample, name2circ, fastq1, fastq2, out_dir, seqtk)

    print('Mapping BSJ reads against circular RNAs using bowtie')
    junc_fastq1_bt_sam_file = bowtie_se(
        junc_fastq1, threads, mismatch, bt_index, bowtie, out_dir)
    junc_fastq2_bt_sam_file = bowtie_se(
        junc_fastq2, threads, mismatch, bt_index, bowtie, out_dir)

    print('Mapping BSJ reads against circular RNAs using hisat2')
    junc_fastq1_ht2_sam_file = hisat2_se(
        junc_fastq1, threads, mismatch, ht2_index, hisat2, out_dir)
    junc_fastq2_ht2_sam_file = hisat2_se(
        junc_fastq2, threads, mismatch, ht2_index, hisat2, out_dir)

    print('Merging sam files...')
    junc_sam_file_bt = merge_sam(
        junc_fastq1_bt_sam_file, junc_fastq2_bt_sam_file, out_dir)
    junc_sam_file_ht2 = merge_sam(
        junc_fastq1_ht2_sam_file, junc_fastq2_ht2_sam_file, out_dir)

    print('Fixing the results of bowtie using hisat2 mapping')
    sam_file_final = fix_sam(sample, junc_sam_file_bt,
                             junc_sam_file_ht2, out_dir, mismatch, read_len)

    print('Conerting sam file into bed file')
    junc_bt_bed_file = sam2bed_smart(
        sam_file_final, samtools, bedtools, threads, out_dir)

    print('Getting genomic coordinate for all BSJ reads')
    junc_genome_bed_file = get_genomic_pos(
        name2circ, junc_bt_bed_file, out_dir)

    print('Mapping reads against intron...')
    read2intron = map_reads_intron(
        junc_genome_bed_file, intron_tx_file, out_dir, bedtools, min_overlap)

    print('Mapping reads against exon...')
    read2exon = map_reads_exon(
        junc_genome_bed_file, exon_tx_file, out_dir, bedtools)

    print('Finding EIciRNAs')
    name2EIci = map_EIci_intron_smart(name2EIci_tmp, read2intron, read2exon)

    return name2EIci


def main():
    args = parse_args()
    sample = args.sample
    fastq1, fastq2 = args.fastqs
    paired_end = args.paired_end
    sam_file = args.sam
    circ_file = args.circ
    ref_fasta = args.genome
    ref_gtf = args.gtf
    bwa_index = args.bwa_index
    bwa = args.bwa
    bowtie_build = args.bowtie_build
    bowtie = args.bowtie
    hisat2_build = args.hisat2_build
    hisat2 = args.hisat2
    seqtk = args.seqtk
    samtools = args.samtools
    bedtools = args.bedtools
    threads = args.threads
    mismatch = args.mismatch
    min_overlap = args.min_overlap
    out_dir = args.output
    src_dir = args.src_dir

    mkdir_p(out_dir)

    if not exist_bwa_index(bwa_index):
        log.info('Building BWA index for {}'.format(ref_fasta))
        bwa_index = build_bwa(bwa, ref_fasta, out_dir)
    name2EIci = {}
    if paired_end:
        if circ_file:
            name2EIci = find_EIci_from_CIRI2(sample, fastq1, fastq2, ref_fasta, ref_gtf, circ_file,
                                             threads, mismatch, out_dir, bowtie, bowtie_build, hisat2, hisat2_build, seqtk, samtools, bedtools, min_overlap)
        elif sam_file:
            name2EIci = find_EIci_from_sam(sample, fastq1, fastq2, ref_fasta, ref_gtf, sam_file,
                                           threads, mismatch, out_dir, src_dir, bowtie, bowtie_build, hisat2, hisat2_build, seqtk, samtools, bedtools, min_overlap)
        else:
            name2EIci = find_EIci_from_fastq(sample, fastq1, fastq2, bwa_index, ref_fasta, ref_gtf,
                                             threads, mismatch, out_dir, src_dir, bwa, bowtie, bowtie_build, hisat2, hisat2_build, seqtk, samtools, bedtools, min_overlap)
    output_file = os.path.join(out_dir, '{}_EIciRNA.bed'.format(sample))
    with open(output_file, 'w') as output:
        output.write('\t'.join(['#chr', 'circ_start', 'circ_end', 'gene_name', 'EIci_junc_count',
                                'strand', 'tx_id', 'tx_score', 'intron_start', 'intron_end', 'reads_list']) + '\n')
        for name in name2EIci:
            name2EIci[name].write_to_file(output)


if __name__ == '__main__':
    main()
