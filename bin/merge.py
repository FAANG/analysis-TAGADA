#!/usr/bin/env python3

from functools import reduce
import pandas as pd
import re
import sys

files = sys.argv[1:]

dataframes = {
  'cov': [],
  'TPM': []
}

keys = []
columns = []
prefix = re.search('\.([^\.]+)\.tsv$', files[0]).group(1)

for file in files:

  df = pd.read_csv(file, sep = '\t').rename(columns = {
    'Gene ID': 'gene_id',
    't_name': 'transcript_id',
    'Coverage': 'cov',
  })

  if not keys:
    if 'transcript_id' in df.columns:
      keys = ['transcript_id', 'gene_id']
    else:
      keys = ['gene_id']

  id = re.sub(r'\.([^\.]+)\.tsv$', '', file)
  columns += [id]

  if not 'TPM' in df.columns and 'FPKM' in df.columns:
    sum_FPKM = df['FPKM'].sum()
    df['TPM'] = df['FPKM'].apply(
      lambda FPKM: 1e6 * FPKM / sum_FPKM
    )

  if 'TPM' in df.columns:
    dataframes['TPM'] += [df[keys + ['TPM']].rename(columns = {
      'TPM': id
    })]

  if 'cov' in df.columns:
    dataframes['cov'] += [df[keys + ['cov']].rename(columns = {
      'cov': id
    })]

cov, TPM = [
  reduce(lambda df1, df2:
    pd.merge(
      df1, df2,
      on = keys,
      how = 'outer'
    ),
    dataframes
  )
  for dataframes in [dataframes['cov'], dataframes['TPM']]
]

cov.reindex(
  columns = keys + sorted(columns, key = str.casefold)
).sort_values(
  keys
).to_csv(prefix+'_coverage.tsv', sep='\t', index=False)

TPM.reindex(
  columns = keys + sorted(columns, key = str.casefold)
).sort_values(
  keys
).to_csv(prefix+'_TPM.tsv', sep='\t', index=False)
