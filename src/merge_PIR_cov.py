#!/usr/bin/env python

"""This script is used to merge the coverage and PIR if introns
Usage: %prog <PIR_bed> <cov_bed> <merged_bed>
"""


import pandas as pd
from sys import argv, exit


try:
    PIR_bed = argv[1]
    cov_bed = argv[2]
    merged_bed = argv[3]
except:
    exit(__doc__)


PIR = pd.read_table(PIR_bed, index_col=None, names=[
    'chr', 'start', 'end', 'tx_id', 'gene_id', 'strand', 'E1I', 'IE2', 'intron', 'E1E2', 'median_junction', 'PIR', 'pval_binominal'])
cov = pd.read_table(cov_bed, index_col=None, names=[
                    'chr', 'start', 'end', 'tx_id', 'gene_id', 'strand', 'feature_num', 'base_num', 'intron_len', 'coverage'])
merge = PIR.merge(cov, how='inner', on=[
                  'chr', 'start', 'end'], suffixes=('', '_cov'))
merge = merge[['chr', 'start', 'end', 'tx_id', 'gene_id', 'strand', 'intron_len',
               'coverage', 'E1I', 'IE2', 'intron', 'E1E2', 'median_junction', 'PIR', 'pval_binominal']]
merge.to_csv(merged_bed, sep='\t', index=False)
