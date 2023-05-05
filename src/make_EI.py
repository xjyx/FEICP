#!/usr/bin/env python

""""Getting exon-intron junction from exon bed file
Usage: sort -k1,1 -k2,2n exon.bed | uniq | %prog - >exon_intron.bed"""

# ======------------------------------------====-------------------------------====
# |   |                               |  | |  |                          |   |


import gzip
import sys
from collections import defaultdict
import re


class Exon(object):
    def __init__(self, chr='', start=0, end=1, gene_name='', tx_id='', strand='+'):
        self.chr = str(chr)
        self.start = int(start)
        self.end = int(end)
        self.gene_name = gene_name
        self.tx_id = tx_id
        self.strand = strand

    def output(self):
        sys.stdout.write('\t'.join([self.chr, str(self.start), str(
            self.end), self.gene_name, self.tx_id, self.strand]) + '\n')


def open_file(infile='', mode=''):
    if re.search('.gz$', infile):
        mode += 't'
        return gzip.open(infile, mode=mode)
    else:
        return open(infile, mode=mode)


def parse_exon(exon_bed_file=''):
    tx2exon = defaultdict(list)
    with open_file(exon_bed_file) as fh:
        for line in fh:
            line = line.strip()
            lineL = line.split('\t')
            chr, start, end, gene_name, tx_id, strand = lineL
            tx2exon[tx_id].append(
                Exon(chr, int(start), end, gene_name, tx_id, strand))

    return tx2exon


def parse_exon_stdin():
    tx2exon = defaultdict(list)
    for line in sys.stdin:
        line = line.strip()
        lineL = line.split('\t')
        chr, start, end, gene_name, tx_id, strand = lineL
        tx2exon[tx_id].append(
            Exon(chr, int(start), end, gene_name, tx_id, strand))

    return tx2exon


def get_EI(exon_bed_file='', flank=5):
    """The exon bed file must be sorted according to the chr and coordinate
    """
    if exon_bed_file == '-':
        tx2exon = parse_exon_stdin()
    else:
        tx2exon = parse_exon(exon_bed_file)
    tx2EI = defaultdict(list)
    for tx in tx2exon:
        if len(tx2exon[tx]) == 1:
            continue
        for i in range(len(tx2exon[tx])):
            if i == 0:
                tx2EI[tx].append(
                    Exon(tx2exon[tx][i].chr, tx2exon[tx][i].end-flank, tx2exon[tx][i].end+flank,
                         tx2exon[tx][i].gene_name, tx2exon[tx][i].tx_id, tx2exon[tx][i].strand))
            elif i == len(tx2exon) - 1:
                tx2EI[tx].append(
                    Exon(tx2exon[tx][i].chr, tx2exon[tx][i].start-flank, tx2exon[tx][i].start+flank,
                         tx2exon[tx][i].gene_name, tx2exon[tx][i].tx_id, tx2exon[tx][i].strand))
            else:
                tx2EI[tx].append(
                    Exon(tx2exon[tx][i].chr, tx2exon[tx][i].start-flank, tx2exon[tx][i].start+flank,
                         tx2exon[tx][i].gene_name, tx2exon[tx][i].tx_id, tx2exon[tx][i].strand))
                tx2EI[tx].append(
                    Exon(tx2exon[tx][i].chr, tx2exon[tx][i].end-flank, tx2exon[tx][i].end+flank,
                         tx2exon[tx][i].gene_name, tx2exon[tx][i].tx_id, tx2exon[tx][i].strand))

    return tx2EI


def main():
    try:
        exon_file = sys.argv[1]
    except:
        sys.exit(__doc__)

    tx2EI = get_EI(exon_file)

    for tx in tx2EI:
        for EI in tx2EI[tx]:
            EI.output()


if __name__ == '__main__':
    main()