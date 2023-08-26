#!/usr/bin/env python

"""This script is used to count the unique reads spanning exon-exon junction from the alignment
======------------------------==========------------------------=========================
                        ======|==========|=========================
The input file can be generated in Linux pipe using the following command:
samtools view -h -q 255 bam_file | \
    awk '$0 ~ /^@/ || $6 ~ /N/' | \
    samtools view -b - | \
    bedtools bamtobed -i - -split | \
    sort -k1,1 -k2,2n | \
    bedtools intersect -a - -b tx_bed -f 1 -wa -wb | \
    bedtools groupby -i - -g 1,2,3,4,5,6 -c 11 -o distinct >read2tx_bed7
"""

import re
import argparse
import gzip
from collections import defaultdict


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Getting exon-exon spanning junction reads number')
    parser.add_argument('--read-tx', '-r', required=True, type=str,
                        help='The read bed file annotated using transcript')
    parser.add_argument('--output', '-o', required=True,
                        type=str, help='The junction count file')
    args = parser.parse_args()
    return args


def open_file(infile='', mode='r'):
    if re.search('.gz$', infile):
        mode += 't'
        return gzip.open(infile, mode=mode)
    return open(infile, mode=mode)


class Block(object):
    def __init__(self, chr='', start='0', end=1, gene_id='', tx_id='', strand='+'):
        self.chr = str(chr)
        self.start = int(start)
        self.end = int(end)
        self.gene_id = gene_id
        self.tx_id = tx_id
        self.strand = str(strand)
        self.len = self.end - self.start


def get_map(read2tx_file=''):
    read2map = defaultdict(list)
    with open_file(read2tx_file) as fh:
        for line in fh:
            line = line.strip()
            lineL = line.split('\t')
            chr, start, end, name, mapq, strand, tx_id = lineL
            read2map[name].append(
                Block(chr, start, end, '', tx_id.split(','), strand))
    # avoid the error: RuntimeError: dictionary changed size during iteration
    for read in list(read2map):
        if len(read2map[read]) == 1:
            read2map.pop(read)
            continue
        # flattening of list
        tx_list_1 = sum([block.tx_id for block in read2map[read]], [])
        tx_list_2 = set(tx_list_1)
        tx_list_2 = [tx for tx in tx_list_2 if tx_list_1.count(
            tx) == len(read2map[read])]  # Make suring that all the block can be mapped to the same transcript
        if not tx_list_2:
            read2map.pop(read)
        else:
            for block in read2map[read]:
                block.tx_id = tx_list_2
    # print('There is {} effective reads'.format(str(len(read2map))))

    return read2map


def get_EE(read2tx_file='', min_overhang=5):
    read2map = get_map(read2tx_file)
    junction2count = defaultdict(int)
    for read in read2map:
        for i in range(len(read2map[read])-1):
            if read2map[read][i].len >= min_overhang and read2map[read][i+1].len >= min_overhang:
                junction2count[(read2map[read][i].chr, read2map[read]
                                [i].end, read2map[read][i+1].start)] += 1

    return junction2count


def main():
    args = parse_args()
    junction2count = get_EE(args.read_tx)
    with open(args.output, 'w') as output:
        for junction in junction2count:
            output.write('\t'.join(
                [str(i) for i in junction] + [str(junction2count[junction])]) + '\n')


if __name__ == '__main__':
    main()
