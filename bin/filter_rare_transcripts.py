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
      0: 'category',
      1: 'category',
      2: 'category',
      3: 'uint32',
      4: 'uint32',
      5: 'float32',
      6: 'category',
      7: 'category',
      8: 'str'
    },
    comment = '#',
    sep = '\t',
    names = columns
  )
  input['file'] = re.search(r'(.+?)\.[^\.]+$', os.path.basename(file)).group(1)
  input['transcript'] = input['attribute'].str.extract(r'transcript_id "(.+?)"')
  input['tpm'] = input['attribute'].str.extract(r'TPM "(.+?)"').astype('float32')
  inputs = inputs.append(input)

inputs = inputs.astype({
  'file': 'category',
  'transcript': 'category'
}).set_index(
  ['file', 'transcript', 'chromosome', 'strand']
)

# Multiexonic transcripts with identical splice sites are the same
transcripts = inputs[inputs['feature'] == 'exon'].sort_values(
  ['start', 'end']
).groupby(
  ['file', 'transcript', 'chromosome', 'strand'],
  observed = True
).apply(
  lambda exons: zip(exons['start'] - 1, exons['end'] + 1)
).apply(
  lambda intervals: [str(x) for interval in intervals for x in interval][1:-1]
).apply(
  lambda positions: 'multi_' + '_'.join(positions) if positions else pd.NA
).reset_index(
  name = 'splicing'
).set_index(
  ['file', 'transcript', 'chromosome', 'strand']
).merge(
  inputs[inputs['feature'] == 'transcript'][['start', 'end', 'tpm']],
  on = ['file', 'transcript', 'chromosome', 'strand']
)

# Monoexonic transcripts with at least `MIN OVERLAP` overlap are the same
mono = transcripts[transcripts['splicing'].isna()][
  ['start', 'end']
].reset_index().sort_values(
  ['chromosome', 'strand', 'start', 'end']
)

mono['group'] = (
  (mono['chromosome'] != mono['chromosome'].shift()) |
  (mono['strand'] != mono['strand'].shift()) |
  (mono['start'] > mono['end'].shift())
).cumsum()

mono = mono.set_index(['chromosome', 'strand', 'group'])

mono = mono.merge(mono, on = ['chromosome', 'strand', 'group'])

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
  ['file', 'transcript', 'chromosome', 'strand'],
  observed = True
).agg(**{
  'start_min': ('start', 'min'),
  'end_max': ('end', 'max'),
}).astype({
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
transcripts = transcripts.reset_index().groupby(
  ['chromosome', 'strand', 'splicing'],
  observed = True
).agg(**{
  'file': ('file', list),
  'transcript': ('transcript', list),
  'tpm_max': ('tpm', 'max'),
  'file_count': ('file', 'nunique')
}).explode(
  ['file', 'transcript']
).reset_index().set_index(
  ['file', 'transcript', 'chromosome', 'strand']
)

transcripts = transcripts[
  (transcripts['file_count'] >= args.min_count) &
  (transcripts['tpm_max'] >= args.min_tpm)
]

# Output filtered GTF
groups = inputs.merge(
  transcripts,
  on = ['file', 'transcript', 'chromosome', 'strand']
).reset_index().groupby(
  'file',
  observed = True
)

os.makedirs(args.output, exist_ok = True)

for file, group in groups:
  group.to_csv(
    f'{args.output}/{file}.filtered.gtf',
    sep = '\t',
    columns = columns,
    header = False,
    index = False,
    quoting = csv.QUOTE_NONE,
    float_format = '%.10g'
  )
