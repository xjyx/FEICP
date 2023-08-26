#!/usr/bin/env python

import pandas as pd
from sys import argv

EIci = pd.read_table(argv[1])
IR = pd.read_table(argv[2])
EIci_out = argv[3]

EIci = EIci.rename(columns={'#chr': 'chr', 'gene_name': 'gene_id'}).query(
    'tx_score >= 1')
EIci['intron_start'] = EIci['intron_start'].str.split(',')
EIci['intron_end'] = EIci['intron_end'].str.split(',')
EIci = EIci.explode(['intron_start', 'intron_end'])
EIci['intron_start'] = EIci['intron_start'].astype('int64')
EIci['intron_end'] = EIci['intron_end'].astype('int64')
IR = IR.query("coverage >= 0.9 and E1I * IE2 * E1E2 * intron > 0").rename(
    {'start': 'intron_start', 'end': 'intron_end'}, axis=1).iloc[:, [0, 1, 2, 7, 8, 9, 10, 11]].drop_duplicates()
EIci_tidy = EIci.merge(IR, on=['chr', 'intron_start', 'intron_end']).iloc[:, [
    i for i in range(6)] + [8, 9, 11, 14, 12, 13, 15]].reset_index(drop=True)
EIci_tidy.to_csv(EIci_out, sep='\t', index=False)
