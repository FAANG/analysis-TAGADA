#!/usr/bin/env python3

import pandas as pd
import re
import sys

files = sys.argv[1:]

result = pd.DataFrame()

keys = []
columns = []
prefix = re.search('\.([^\.]+)\.tsv$', files[0]).group(1)

for file in files:

  df = pd.read_csv(file, sep = '\t').rename(columns = {
    'Gene ID': 'gene',
    'gene_id': 'gene',
    'transcript_id': 'transcript',
    't_name': 'transcript'
  })

  if not keys:
    if 'transcript' in df.columns:
      keys = ['transcript', 'gene']
    else:
      keys = ['gene']

  id = re.sub(r'\.([^\.]+)\.tsv$', '', file)
  columns += [id]

  if not 'TPM' in df.columns and 'FPKM' in df.columns:
    sum_FPKM = df['FPKM'].sum()
    df['TPM'] = df['FPKM'].apply(
      lambda FPKM: 1e6 * FPKM / sum_FPKM
    )

  if 'TPM' in df.columns:
    df = df[keys + ['TPM']].rename(columns = {
      'TPM': id
    }).groupby(keys).sum().round(6).astype(str).reset_index()

  elif 'counts' in df.columns:
    df = df[keys + ['counts']].rename(columns = {
      'counts': id
    }).groupby(keys).sum().reset_index()

  result = result.merge(
    df,
    on = keys,
    how = 'outer'
  ) if not result.empty else df

result.reindex(
  columns = keys + sorted(columns, key = str.casefold)
).sort_values(
  keys
).to_csv(prefix + '.tsv', sep = '\t', index = False)
