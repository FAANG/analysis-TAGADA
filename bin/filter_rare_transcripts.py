#!/usr/bin/env python3

import os
import re
import csv
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
  'inputs',
  nargs = '+',
  help = 'Input files in GTF format.'
)
parser.add_argument(
  '-o',
  '--output',
  required = True,
  help = 'Output directory.'
)
parser.add_argument(
  '--min-overlap',
  default = 50,
  type = int,
  help = 'Monoexonic transcripts with at least `--min-overlap` nucleotides' +
         'in common are assumed to be the same transcript.'
)
parser.add_argument(
  '--min-count',
  default = 2,
  type = int,
  help = 'Delete transcripts that appear in less than `--min-count` inputs.'
)
parser.add_argument(
  '--min-tpm',
  default = 0.1,
  type = float,
  help = 'Delete transcripts with all TPM values less than `--min-tpm`.'
)

args = parser.parse_args()

# Parse GTF inputs
columns = [
  'chromosome',
  'source',
  'feature',
  'start',
  'end',
  'score',
  'strand',
  'frame',
  'attribute'
]
inputs = pd.DataFrame()
for file in args.inputs:
  input = pd.read_csv(
    file,
    dtype = {
      i: str if i not in [3, 4] else int for i in range(8)
    },
    comment = '#',
    sep = '\t',
    names = columns
  )
  input['file'] = re.search(r'(.+?)\.[^\.]+$', os.path.basename(file)).group(1)
  input['transcript'] = input['attribute'].str.extract(r'transcript_id "(.+?)"')
  input['tpm'] = input['attribute'].str.extract(r'TPM "(.+?)"').astype(float)
  inputs = inputs.append(input)

# Multiexonic transcripts with identical splice sites are the same
transcripts = inputs[inputs['feature'] == 'exon'].sort_values(
  ['start', 'end']
).groupby(
  ['file', 'transcript', 'chromosome', 'strand']
).apply(
  lambda table: zip(table['start'] - 1, table['end'] + 1)
).apply(
  lambda exons: [str(x) for exon in exons for x in exon][1:-1]
).apply(
  lambda introns: 'multi_' + '_'.join(introns) if len(introns) else pd.NA
).reset_index(
  name = 'splicing'
).merge(
  inputs[inputs['feature'] == 'transcript'],
  on = ['file', 'transcript', 'chromosome', 'strand']
)

# Monoexonic transcripts with overlap of at least `--min-overlap` are the same
transcripts = transcripts.sort_values(['chromosome', 'strand', 'start', 'end'])
monoexonic = transcripts[transcripts['splicing'].isna()]
transcripts.loc[transcripts['splicing'].isna(), 'splicing'] = 'mono_' + (
  (monoexonic['chromosome'] != monoexonic['chromosome'].shift())
  | (monoexonic['strand'] != monoexonic['strand'].shift())
  | (monoexonic['start'] > monoexonic['end'].shift() - args.min_overlap)
).cumsum().astype(str)

# Keep transcripts that appear in at least `--min-count` files
# Keep transcripts with at least one TPM greater or equal to `--min-tpm`
transcripts = transcripts.groupby(['chromosome', 'strand', 'splicing']).agg(**{
  'file': ('file', list),
  'transcript': ('transcript', list),
  'start': ('start', 'first'),
  'end': ('end', 'last'),
  'tpm_max': ('tpm', 'max'),
  'file_count': ('file', 'nunique')
}).reset_index().explode(['file', 'transcript'])

transcripts = transcripts[
  (transcripts['file_count'] >= args.min_count)
  & (transcripts['tpm_max'] >= args.min_tpm)
]

# Output filtered GTF
groups = transcripts[['file', 'transcript']].merge(
  inputs,
  on = ['file', 'transcript']
).groupby('file')

os.makedirs(args.output, exist_ok = True)

for file, group in groups:
  group[columns].to_csv(
    f'{args.output}/{file}.filtered.gtf',
    sep = '\t',
    header = False,
    index = False,
    quoting = csv.QUOTE_NONE
  )
