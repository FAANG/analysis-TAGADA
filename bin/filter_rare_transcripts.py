#!/usr/bin/env python3

import os
import re
import csv
import argparse
import numpy as np
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument(
  'inputs',
  metavar = 'INPUTS',
  nargs = '+',
  help = 'Input files in GTF format.'
)
parser.add_argument(
  '-o',
  '--output',
  metavar = 'OUTPUT',
  required = True,
  help = 'Output directory.'
)
parser.add_argument(
  '--min-overlap',
  metavar = 'MIN OVERLAP',
  default = 0.5,
  type = float,
  help = 'Monoexonic transcripts are assumed to be the same transcript if ' +
         'they both have an overlap percentage of at least `MIN OVERLAP`.'
)
parser.add_argument(
  '--min-count',
  metavar = 'MIN COUNT',
  default = 2,
  type = int,
  help = 'Keep transcripts that appear in at least `MIN COUNT` inputs.'
)
parser.add_argument(
  '--min-tpm',
  metavar = 'MIN TPM',
  default = 0.1,
  type = float,
  help = 'Keep transcripts with at least one TPM greater or equal to `MIN TPM`.'
)

args = parser.parse_args()

# Parse GTF inputs
inputs = pd.DataFrame()

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
  lambda exons: zip(exons['start'] - 1, exons['end'] + 1)
).apply(
  lambda intervals: [str(x) for interval in intervals for x in interval][1:-1]
).apply(
  lambda positions: 'multi_' + '_'.join(positions) if positions else pd.NA
).reset_index(
  name = 'splicing'
).merge(
  inputs[inputs['feature'] == 'transcript'],
  on = ['file', 'transcript', 'chromosome', 'strand']
)

# Monoexonic transcripts with at least `MIN OVERLAP` overlap are the same
mono = transcripts[transcripts['splicing'].isna()]
mono = mono.merge(mono, on = ['chromosome', 'strand'])

mono = mono[
  (mono['start_x'] <= mono['end_y']) &
  (mono['start_y'] <= mono['end_x'])
]

overlap = (
  np.minimum(mono['end_x'], mono['end_y']) -
  np.maximum(mono['start_x'], mono['start_y'])
)

mono = mono[
  (overlap / (mono['end_x'] - mono['start_x']) >= args.min_overlap) &
  (overlap / (mono['end_y'] - mono['start_y']) >= args.min_overlap)
]

mono = mono.rename(columns = {
  'file_x': 'file',
  'transcript_x': 'transcript',
  'start_y': 'start',
  'end_y': 'end',
}).groupby(
  ['file', 'transcript', 'chromosome', 'strand']
).agg(**{
  'start_min': ('start', 'min'),
  'end_max': ('end', 'max'),
}).reset_index().astype({
  'start_min': str,
  'end_max': str
})

transcripts = transcripts.merge(
  mono,
  on = ['file', 'transcript', 'chromosome', 'strand'],
  how = 'left'
)

transcripts['splicing'] = transcripts['splicing'].fillna(
  'mono_' + transcripts['start_min'] + '_' + transcripts['end_max']
)

# Keep transcripts that appear in at least `MIN COUNT` files
# Keep transcripts with at least one TPM greater or equal to `MIN TPM`
transcripts = transcripts.groupby(
  ['chromosome', 'strand', 'splicing']
).agg(**{
  'file': ('file', list),
  'transcript': ('transcript', list),
  'tpm_max': ('tpm', 'max'),
  'file_count': ('file', 'nunique')
}).reset_index().explode(
  ['file', 'transcript']
)

transcripts = transcripts[
  (transcripts['file_count'] >= args.min_count) &
  (transcripts['tpm_max'] >= args.min_tpm)
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
