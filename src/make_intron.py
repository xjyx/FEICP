# /usr/bin/env python

"""This script is used to get intron regions of a genome from gtf file
"""


from collections import defaultdict
import re
import gzip
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser(
        prog='Getting regions of intron from gtf file')
    parser.add_argument('--gtf', '-g', dest='gtf',
                        type=str, required=True, help='The input gtf file')
    parser.add_argument('--output', '-o', dest='output',
                        type=str, required=True, help='The output bed file storing introns')

    args = parser.parse_args()
    return args


def open_file(infile, mode='r'):
    if re.search('.gz$', infile):
        mode = 'rt'
        return gzip.open(infile, mode=mode)
    else:
        return open(infile, mode=mode)


class Block(object):
    def __init__(self, chr, start, end, tx_id, gene_id, strand):
        self.chr = chr
        self.start = int(start)
        self.end = int(end)
        self.tx_id = tx_id
        self.gene_id = gene_id
        self.strand = strand

    def format_to_string(self):
        return '\t'.join([self.chr, str(self.start), str(self.end), self.tx_id, self.gene_id, self.strand])

    def write_to_file(self, fh):
        fh.write(self.format_to_string() + '\n')


def get_exons(gtf_file):
    tx2exons = defaultdict(list)
    with open_file(gtf_file) as fh:
        for line in fh:
            line = line.strip()
            if re.match('^#', line):
                continue
            lineL = line.split('\t')
            if lineL[2] == "exon":
                chr, start, end, tx_id, gene_id, strand = lineL[0], int(
                    lineL[3]) - 1, lineL[4], lineL[8].split('\"')[3], lineL[8].split('\"')[1], lineL[6]
                tx2exons[tx_id].append(
                    Block(chr, start, end, tx_id, gene_id, strand))

    for tx in tx2exons:
        tx2exons[tx].sort(key=lambda x: (x.chr, x.start))

    return tx2exons


def get_introns(gtf_file):
    tx2exons = get_exons(gtf_file)
    tx2introns = defaultdict(list)
    for tx in tx2exons:
        for i in range(0, len(tx2exons[tx]) - 1, 1):
            tx2introns[tx].append(Block(tx2exons[tx][i].chr, tx2exons[tx][i].end, tx2exons[tx]
                                  [i+1].start, tx2exons[tx][i].gene_id, tx, tx2exons[tx][i].strand))

    intron_list = []
    for tx in tx2introns:
        intron_list += tx2introns[tx]
    intron_list.sort(key=lambda x: (x.chr, x.start))

    return intron_list


def main():
    args = parse_args()
    gtf_file = args.gtf
    out_file = args.output
    intron_list = get_introns(gtf_file)
    with open(out_file, 'w') as output:
        for intron in intron_list:
            intron.write_to_file(output)


if __name__ == '__main__':
    main()
