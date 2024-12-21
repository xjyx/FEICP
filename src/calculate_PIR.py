#!/usr/bin/env python

"""Bedtools does not lie or cheat"""


import argparse
from collections import defaultdict
import re
from FEICP import open_file
from scipy.stats import binomtest
import numpy as np
import sys


EXON_INTRON_JUNC_LEN = 5


def parse_args():
    parser = argparse.ArgumentParser(
        prog='Calculating Percent intron retention (PIR) of specific Introns')
    parser.add_argument('--ee-count', '-ee', required=True,
                        type=str, help='The count of exon-exon junction')
    parser.add_argument('--ei-count', '-ei', required=True,
                        type=str, help='The count of exon-intron junction')
    parser.add_argument('--intron-count', '-i', type=str,
                        required=True, help='The count intron')
    parser.add_argument('--query-intron', '-a', required=True,
                        help='The intron region used to calculate the PIR')
    parser.add_argument('--output', '-o', required=True, type=str,
                        help='The output PIR value of all input introns')

    args = parser.parse_args()
    return args


class Intron(object):
    def __init__(self, chr='chr1', start=0, end=1, tx_id='ENST00', gene_id='ENSG00', strand='',
                 EI_count=-1, IE_count=-1, EE_count=-1, Intron_count=-1):
        self.chr = str(chr)
        self.start = int(start)
        self.end = int(end)
        self.tx_id = str(tx_id.split(',')[0])
        self.gene_id = str(gene_id.split(',')[0])
        self.strand = str(strand)
        self.EE_junction_id = (self.chr, self.start, self.end)
        self.EI_junction_id = (self.chr, self.start-EXON_INTRON_JUNC_LEN, self.start +
                               EXON_INTRON_JUNC_LEN, self.strand)
        self.IE_junction_id = (self.chr, self.end-EXON_INTRON_JUNC_LEN, self.end +
                               EXON_INTRON_JUNC_LEN, self.strand)
        self.Intron_id = (self.chr, self.start, self.end, self.strand) if self.end - self.start <= 200 \
            else (self.chr, int((self.start+self.end)/2-100), int((self.start+self.end)/2+100), self.strand)
        self.EI_count = int(EI_count)
        self.IE_count = int(IE_count)
        self.EE_count = int(EE_count)
        self.Intron_count = int(Intron_count)

    def median_junction_count(self):
        # > 10
        return np.median([self.EI_count, self.IE_count, self.Intron_count]) + self.EE_count

    def PIR(self):
        # PIRmin <= 95
        return 100 * np.average([self.EI_count, self.IE_count]) / (self.EE_count + np.average([self.EI_count, self.IE_count]))

    def pval_binom(self):
        return float(binomtest(k=min(self.EI_count, self.IE_count, self.Intron_count),
                               n=min(self.EI_count, self.IE_count, self.Intron_count)+max(
            self.EI_count, self.IE_count, self.Intron_count),
            p=1/3.5,
            alternative='less').pvalue) \
            if min(self.EI_count, self.IE_count, self.Intron_count)+max(
            self.EI_count, self.IE_count, self.Intron_count) >= 1 else 0  # P-value >= 0.05

    def write_to_file(self, fh):
        fh.write('\t'.join([str(i)
                 for i in [self.chr, self.start, self.end, self.tx_id, self.gene_id, self.strand, self.EI_count, self.IE_count, self.Intron_count, self.EE_count, self.median_junction_count(), self.PIR(), self.pval_binom()]]) + '\n')


def get_count(count_file='', col_num=4):
    if col_num != 4 and col_num != 7:
        sys.exit(
            'The number of column is not right, please check the {}'.format(count_file))
    region2count = defaultdict(int)
    with open_file(count_file) as fh:
        for line in fh:
            if re.search('^#', line):
                continue
            lineL = line.split('\t')
            if col_num == 4:
                region2count[(lineL[0], int(lineL[1]), int(lineL[2]))] = int(
                    lineL[3])
            elif col_num == 7:
                region2count[(lineL[0], int(lineL[1]), int(
                    lineL[2]), lineL[5])] = int(lineL[6])
    return region2count


def get_intron_list(intron_bed=''):
    intron_list = []
    with open_file(intron_bed) as fh:
        for line in fh:
            if re.search('^#', line):
                continue
            line = line.strip()
            lineL = line.split('\t')
            if len(lineL) < 6:
                continue
            intron_list.append(
                Intron(lineL[0], lineL[1], lineL[2], lineL[3], lineL[4], lineL[5]))
    return intron_list


def calculate_PIR(intron_bed='', EE_count='', EI_count='', intron_count=''):
    EE2count = get_count(EE_count, 4)
    EI2count = get_count(EI_count, 7)
    intron2count = get_count(intron_count, 7)
    intron_list = get_intron_list(intron_bed)
    for i in intron_list:
        i.EI_count = EI2count[i.EI_junction_id]
        i.IE_count = EI2count[i.IE_junction_id]
        i.EE_count = EE2count[i.EE_junction_id]
        i.Intron_count = intron2count[i.Intron_id]
    return intron_list


def main():
    args = parse_args()
    intron_list = calculate_PIR(
        args.query_intron, args.ee_count, args.ei_count, args.intron_count)
    with open(args.output, 'w') as output:
        # output.write('\t'.join(['#chr', 'start', 'end', 'tx_id', 'gene_id', 'strand', '#E1I',
        #              '#IE2', '#Intron', '#E1E2', '#median_junction', 'PIR', 'Pval_binomial']) + '\n')
        for intron in intron_list:
            intron.write_to_file(output)


if __name__ == '__main__':
    main()
