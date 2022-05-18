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
  '--min-monoexonic-overlap',
  metavar = 'MIN_MONOEXONIC_OVERLAP',
  default = 0.5,
  type = float,
  help = 'Monoexonic transcripts are assumed to be the same transcript if ' +
         'they both have an overlap ratio of at least `MIN_MONOEXONIC_' +
         'OVERLAP`. Defaults to `0.5`.'
)
parser.add_argument(
  '--min-transcript-occurrence',
  metavar = 'MIN_TRANSCRIPT_OCCURRENCE',
  default = 2,
  type = int,
  help = 'Filter transcripts that appear in less than `MIN_TRANSCRIPT_' +
         'OCCURRENCE` inputs. Defaults to `2`. This value is lowered to the ' +
         'number of inputs if it exceeds it.'
)
parser.add_argument(
  '--min-monoexonic-occurrence',
  metavar = 'MIN_MONOEXONIC_OCCURRENCE',
  default = None,
  type = int,
  help = 'Filter monoexonic transcripts that appear in less than `MIN_' +
         'MONOEXONIC_OCCURRENCE` inputs. Defaults to `MIN_TRANSCRIPT_' +
         'OCCURRENCE`. This value is lowered to the number of inputs if it ' +
         'exceeds it.'
)
parser.add_argument(
  '--min-transcript-tpm',
  metavar = 'MIN_TRANSCRIPT_TPM',
  default = 0.1,
  type = float,
  help = 'Filter transcripts with no TPM equal to or greater than `MIN_' +
         'TRANSCRIPT_TPM`. Defaults to `0.1`.'
)
parser.add_argument(
  '--min-monoexonic-tpm',
  metavar = 'MIN_MONOEXONIC_TPM',
  default = None,
  type = float,
  help = 'Filter monoexonic transcripts with no TPM equal to or greater than ' +
         '`MIN_MONOEXONIC_TPM`. Defaults to `10 * MIN_TRANSCRIPT_TPM`.'
)

args = parser.parse_args()

args.min_transcript_occurrence = min(
  args.min_transcript_occurrence,
  len(args.inputs)
)

if args.min_monoexonic_occurrence == None:
  args.min_monoexonic_occurrence = args.min_transcript_occurrence
else:
  args.min_monoexonic_occurrence = min(
    args.min_monoexonic_occurrence,
    len(args.inputs)
  )

if args.min_monoexonic_tpm == None:
  args.min_monoexonic_tpm = 10 * args.min_transcript_tpm

# Parse GTF inputs
inputs = pd.DataFrame()

columns = {
  'chromosome': 'category',
  'source': 'category',
  'feature': 'category',
  'start': 'uint32',
  'end': 'uint32',
  'score': 'str',
  'strand': 'category',
  'frame': 'category',
  'attribute': 'str'
}

for gtf in args.inputs:
  input = pd.read_csv(
    gtf,
    sep = '\t',
    comment = '#',
    names = columns.keys(),
    dtype = columns,
    usecols = [
      'chromosome',
      'feature',
      'start',
      'end',
      'strand',
      'attribute'
    ]
  )
  input['file'] = re.search(r'(.+?)\.[^\.]+$', os.path.basename(gtf)).group(1)
  input['transcript'] = input['attribute'].str.extract(r'transcript_id "(.+?)"')
  input['tpm'] = input['attribute'].str.extract(r'TPM "(.+?)"').astype('float32')
  input = input.drop(columns = 'attribute')
  inputs = pd.concat([inputs, input])

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

# Monoexonic transcripts with at least `MIN_MONOEXONIC_OVERLAP` are the same
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
  (overlap / (mono['end_x'] - mono['start_x']) >= args.min_monoexonic_overlap) &
  (overlap / (mono['end_y'] - mono['start_y']) >= args.min_monoexonic_overlap)
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

transcripts['is_monoexonic'] = transcripts['splicing'].isna()

transcripts['splicing'] = transcripts['splicing'].fillna(
  'mono_' + transcripts['start_min'] + '_' + transcripts['end_max']
)

# Filter transcripts that appear in less than `MIN_xxx_OCCURRENCE` files
# Filter transcripts with no TPM equal to or greater than `MIN_xxx_TPM`
transcripts = transcripts.reset_index().groupby(
  ['chromosome', 'strand', 'splicing'],
  observed = True
).agg(**{
  'file': ('file', list),
  'transcript': ('transcript', list),
  'tpm_max': ('tpm', 'max'),
  'file_count': ('file', 'nunique'),
  'is_monoexonic': ('is_monoexonic', 'first')
}).explode(
  ['file', 'transcript']
).reset_index().set_index(
  ['file', 'transcript', 'chromosome', 'strand']
)

monoexonic_count = len(transcripts[transcripts['is_monoexonic'] == True].index)
multiexonic_count = len(transcripts[transcripts['is_monoexonic'] == False].index)
print('All monoexonic transcripts: ' + str(monoexonic_count))
print('All multiexonic transcripts: ' + str(multiexonic_count))

transcripts = transcripts[
  (
    (transcripts['is_monoexonic'] == True) &
    (transcripts['file_count'] >= args.min_monoexonic_occurrence) &
    (transcripts['tpm_max'] >= args.min_monoexonic_tpm)
  ) | (
    (transcripts['is_monoexonic'] == False) &
    (transcripts['file_count'] >= args.min_transcript_occurrence) &
    (transcripts['tpm_max'] >= args.min_transcript_tpm)
  )
]

monoexonic_count = len(transcripts[transcripts['is_monoexonic'] == True].index)
multiexonic_count = len(transcripts[transcripts['is_monoexonic'] == False].index)
print('Kept monoexonic transcripts: ' + str(monoexonic_count))
print('Kept multiexonic transcripts: ' + str(multiexonic_count))

# Output filtered GTFs
os.makedirs(args.output, exist_ok = True)

for gtf in args.inputs:
  input = pd.read_csv(
    gtf,
    sep = '\t',
    comment = '#',
    names = columns.keys(),
    dtype = columns
  )
  file = re.search(r'(.+?)\.[^\.]+$', os.path.basename(gtf)).group(1)
  input['file'] = file
  input['transcript'] = input['attribute'].str.extract(r'transcript_id "(.+?)"')
  input = input.astype({
    'file': 'category',
    'transcript': 'category'
  }).set_index(
    ['file', 'transcript', 'chromosome', 'strand']
  )

  output = input.merge(
    transcripts,
    on = ['file', 'transcript', 'chromosome', 'strand']
  ).reset_index()

  output.to_csv(
    f'{args.output}/{file}.filtered.gtf',
    sep = '\t',
    columns = columns.keys(),
    header = False,
    index = False,
    quoting = csv.QUOTE_NONE,
    float_format = '%.10g'
  )
