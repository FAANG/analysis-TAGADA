#!/usr/bin/env python3

from functools import reduce
import pandas as pd
import re
import sys

prefix = sys.argv[1]
files = sys.argv[2:]

dataframes = {
  'cov': [],
  'TPM': []
}

columns = []

for file in files:

  df = pd.read_csv(file, sep = '\t').rename(columns = {
    'Gene ID': 'gene_id',
    't_name': 'transcript_id',
    'Coverage': 'cov',
  })

  if not columns:
    if 'transcript_id' in df.columns:
      columns = ['transcript_id', 'gene_id']
    else:
      columns = ['gene_id']

  file =  re.sub(r'\.(genes|transcripts)\.tsv$', '', file)

  if not 'TPM' in df.columns and 'FPKM' in df.columns:
    df['TPM'] = df['FPKM'].apply(
      lambda x: x / df['FPKM'].sum() * 1e6
    )

  if 'TPM' in df.columns:
    dataframes['TPM'] += [df[columns + ['TPM']].rename(columns = {
      'TPM': 'TPM_' + file
    })]

  if 'cov' in df.columns:
    dataframes['cov'] += [df[columns + ['cov']].rename(columns = {
      'cov': 'cov_' + file
    })]

cov, TPM = [
  reduce(lambda df1, df2:
    pd.merge(
      df1, df2,
      on = columns,
      how = 'outer'
    ),
    dataframes
  )
  for dataframes in [dataframes['cov'], dataframes['TPM']]
]

cov.sort_values(columns).to_csv(prefix+'_coverage.tsv', sep = '\t', index = False)
TPM.sort_values(columns).to_csv(prefix+'_TPM.tsv', sep = '\t', index = False)
