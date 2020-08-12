// Parse and check parameters
// ##########################
output = params.containsKey('output') ? params.output : ''
reads = params.containsKey('reads') ? params.reads : ''
genome = params.containsKey('genome') ? params.genome : ''
index = params.containsKey('index') ? params.index : ''
annotation = params.containsKey('annotation') ? params.annotation : ''
metadata = params.containsKey('metadata') ? params.metadata : ''
merge = params.containsKey('merge') ? params.merge.tokenize(',') : ''

error = ''

if (!output) error += 'No --output provided\n'

if (!reads) {
  error += 'No --reads provided\n'
} else {
  number_of_raw_reads = Channel.fromPath(reads).filter { path ->
    filename = path.getName()
    return filename =~ /\.(fastq|fq)(\.gz)?$/
  }.count().get()

  if (number_of_raw_reads > 0 && !genome && !index) {
    error += 'No --genome or --index provided\n'
  }
}

if (!annotation) error += 'No --annotation provided\n'

if (merge && !metadata) {
  error += 'No --metadata provided to execute --merge\n'
}

if (error) exit 1, error

// Check index, genome, annotation, metadata
// #########################################
if (index) {
  Channel.fromPath(index, type: 'dir', checkIfExists: true).set {
    index_to_map
  }
}

if (genome) {
  Channel.fromPath(genome, checkIfExists: true).set {
    genome_to_index
  }
}

Channel.fromPath(annotation, checkIfExists: true).into {
  reference_annotation_to_index
  reference_annotation_to_direction
  reference_annotation_to_assemble
  reference_annotation_to_combine
  reference_annotation_to_quantify
}

if (metadata) {
  Channel.fromPath(
    metadata,
    checkIfExists: true
  ).splitCsv(header: true, sep: '\t').into {
    metadata_to_check
    metadata_to_merge
  }
}

// Split reads into R1/R2/single/mapped
// ####################################
Channel.fromPath(reads, checkIfExists: true).map { path ->

  filename = path.getName()

  return (
    filename =~ /^(.+?)(?:[\._ ][Rr]([12]))?(?:\.(fastq|fq|gz|bam))+$/
  ).with {
    matches() ? [
      'prefix': it[0][1],
      'R': it[0][2] ? it[0][2].toInteger() : null,
      'mapped': it[0][3] == 'bam',
      'path': path,
      'filename': filename
    ] : filename
  }

}.branch {
  invalid: it instanceof String
  mapped: it['mapped']
  R1: it['R'] == 1
  R2: it['R'] == 2
  single: true
}.set {
  reads
}

reads.R1.into {
  r1_reads_to_check
  r1_reads_to_pair
  r1_reads_to_quality
}

reads.R2.into {
  r2_reads_to_check
  r2_reads_to_pair
  r2_reads_to_quality
}

reads.single.into {
  single_reads_to_check
  single_reads_to_append
  single_reads_to_quality
}

reads.mapped.tap {
  mapped_reads_to_check
}.map {
  [it['prefix'], it['path']]
}.set {
  mapped_reads_to_direction
}

error = ''
warning = ''

// Check file names
// ################
invalid = reads.invalid.toList().get()

if (invalid.size() > 0) {
  s = invalid.size() > 1 ? 's' : ''
  error += "Wrong format for file name$s:\n  "
  error += invalid.join('\n  ') + '\n'
}

// Check pairing
// #############
r1 = r1_reads_to_check.map {
  ['prefix': it['prefix'], 'filename': it['filename']]
}.toList().get()

r2 = r2_reads_to_check.map {
  ['prefix': it['prefix'], 'filename': it['filename']]
}.toList().get()

unpaired = r1.findAll {
  !(it['prefix'] in r2.collect{it['prefix']})
} + r2.findAll {
  !(it['prefix'] in r1.collect{it['prefix']})
}

if (unpaired.size() > 0) {
  s = unpaired.size() > 1 ? 's' : ''
  error += "No pair$s found for file$s:\n  "
  error += unpaired.collect{it['filename']}.join('\n  ') + '\n'
}

// Check duplicates
// ################
paired = r1 + r2

single = single_reads_to_check.map {
  ['prefix': it['prefix'], 'filename': it['filename']]
}.toList().get()

mapped = mapped_reads_to_check.map {
  ['prefix': it['prefix'], 'filename': it['filename']]
}.toList().get()

duplicated = (
  paired.intersect(single, { a, b -> a['prefix'] <=> b['prefix'] })
  + paired.intersect(mapped, { a, b -> a['prefix'] <=> b['prefix'] })
  + single.intersect(paired, { a, b -> a['prefix'] <=> b['prefix'] })
  + single.intersect(mapped, { a, b -> a['prefix'] <=> b['prefix'] })
  + mapped.intersect(paired, { a, b -> a['prefix'] <=> b['prefix'] })
  + mapped.intersect(single, { a, b -> a['prefix'] <=> b['prefix'] })
).groupBy {
  it['prefix']
}.collect {
  it.value.collect{it['filename']}
}

if (duplicated.size() > 0) {
  error += 'Duplicates detected:\n  '
  error += duplicated.collect{it.join('  ')}.join('\n  ') + '\n'
}

// Check metadata
// ##############
if (merge && metadata) {

  metadata_to_check = metadata_to_check.toList().get()

  // Check missing metadata rows
  // ###########################
  metadata_rows = metadata_to_check.collect {
    it.values()[0]
  }

  inputs = paired + single + mapped

  missing_rows = inputs.findAll {
    !(it['prefix'] in metadata_rows)
  }

  if (missing_rows.size() > 0) {
    s = missing_rows.size() > 1 ? 's' : ''
    error += "No metadata row$s found for file$s:\n  "
    error += missing_rows.collect{it['filename']}.join('\n  ') + '\n'
  }

  // Check missing metadata columns
  // ##############################
  metadata_columns = metadata_to_check[0].keySet()

  missing_columns = merge.findAll {
    !(it in metadata_columns)
  }

  if (missing_columns.size() > 0) {
    s = missing_columns.size() > 1 ? 's' : ''
    error += "No metadata column$s found for merge factor$s:\n  "
    error += missing_columns.join('\n  ') + '\n'
  }

  // Check missing metadata values
  // #############################
  missing_values = metadata_to_check.findAll {
    it.values()[0] in inputs.collect{it['prefix']}
  }.inject([], { result, row ->
    factors = []
    merge.each { factor ->
      if (!row[factor]) factors += factor
    }
    if (factors) {
      result += [[row.values()[0]] + factors]
    }
    return result
  })

  if (missing_values.size() > 0) {
    s = missing_values.size() > 1 ? 's' : ''
    ss = missing_values.inject([], {result, next ->
      next[1..-1].each {
        if (!(it in result)) result += it
      }
      return result
    }).size() > 1 ? 's' : ''
    warning += "Metadata row$s with missing factor$ss will not be merged:\n  "
    warning += missing_values.collect{it.join('  ')}.join('\n  ')
  }
}

if (error) exit 1, error

if (warning) log.warn warning

// Determine merge groups
// ######################
if (merge && metadata) {

  metadata_to_merge.filter {
    it.values()[0] in inputs.collect{it['prefix']}
  }.map {
    factors = it.subMap(merge).values()
    id = '' in factors ? 'NA_' + it.values()[0] : factors.join('_')
    [id, it.values()[0]]
  }.groupTuple().into {
    groups_to_merge
    groups_to_log
  }

  groups_to_log = groups_to_log.toList().get().sort { a, b -> a[0] <=> b[0] }

  s = groups_to_log.size() > 1 ? 's' : ''
  info = "The following merge group$s will be created:\n  "
  info += groups_to_log.collect{ it[0] + ': ' + it[1].join('  ') }.join('\n  ')

  log.info info
}

// Pair R1/R2 and format reads
// ###########################
r1_reads_to_pair.map {
  [it['prefix'], it['path']]
}.set {
  r1_reads_to_pair
}

r2_reads_to_pair.map {
  [it['prefix'], it['path']]
}.set {
  r2_reads_to_pair
}

single_reads_to_append.map {
  [it['prefix'], '', it['path']]
}.set {
  single_reads_to_append
}

r1_reads_to_pair.join(r2_reads_to_pair).map { it ->
  [it[0], '--paired', [it[1], it[2]]]
}.concat(single_reads_to_append).set {
  reads_to_trim
}

single_reads_to_quality.concat(
  r1_reads_to_quality,
  r2_reads_to_quality
).map {
  [it['prefix'], it['R'] ? '_R' + it['R'] : '', it['path']]
}.set {
  reads_to_quality
}

// Process raw reads
// #################
if (number_of_raw_reads > 0) {

  // Control reads quality
  // #####################
  process quality {

    publishDir "$output/control/quality/raw", mode: 'copy'

    input:
      tuple val(prefix), val(R), path(read) from reads_to_quality

    output:
      path '*_fastqc.html'

    script:
      """
      fastqc $read
      mv *_fastqc.html "$prefix""$R"_raw_fastqc.html
      """
  }

  // Trim adaptators
  // ###############
  process trim {

    label 'high_cpu'
    label 'medium_memory'

    publishDir "$output", mode: 'copy', saveAs: { filename ->
      if (filename.endsWith('trimming_report.txt')) "logs/trim_galore/$filename"
      else if (filename.endsWith('fastqc.html')) "control/quality/trimmed/$filename"
    }

    input:
      tuple val(prefix), val(paired), path(reads) from reads_to_trim

    output:
      path '*_trimming_report.txt'
      path '*_fastqc.html'
      tuple val(prefix), path('*_trimmed.fq.gz') into reads_to_map
      path '*_trimmed.fq.gz' into reads_to_overhang

    script:
      """
      for f in $reads; do
        name=\$(basename "\$f")
        symlink="\${name// /_}"
        target=\$(readlink -m "\$f")
        [[ "\$name" != "\$symlink" ]] && ln -sf "\$target" "\$symlink"
        reads+=("\$symlink")
      done

      trim_galore "\${reads[@]}" \\
                  $paired \\
                  --cores ${task.cpus} \\
                  --fastqc \\
                  --gzip \\
                  --basename "$prefix"

      for f in *_{val_1,val_2,trimmed}.fq.gz; do
        [ -f "\$f" ] || break
        suffix="\${f%.fq.gz}"
        R=""
        [[ "\$suffix" =~ _val_[12]\$ ]] && R=_R"\${suffix: -1}"
        name=\$(basename "\$f")
        rename="$prefix""\$R"_trimmed.fq.gz
        [[ "\$name" != "\$rename" ]] && mv "\$f" "\$rename"
      done

      for f in *_fastqc.html; do
        [ -f "\$f" ] || break
        suffix="\${f%_fastqc.html}"
        R=""
        [[ "\$suffix" =~ _val_[12]\$ ]] && R=_R"\${suffix: -1}"
        name=\$(basename "\$f")
        rename="$prefix""\$R"_trimmed_fastqc.html
        [[ "\$name" != "\$rename" ]] && mv "\$f" "\$rename"
      done

      for f in *_trimming_report.txt; do
        [ -f "\$f" ] || break
        suffix="\${f%_trimming_report.txt}"
        while [[ "\$suffix" =~ \\.(fastq|fq|gz|bam)\$ ]]; do
          suffix="\${suffix%.gz}"
          suffix="\${suffix%.fq}"
          suffix="\${suffix%.fastq}"
          suffix="\${suffix%.bam}"
        done
        R=""
        [[ "\$suffix" =~ ([\\._]|[[:blank:]])[Rr][12]\$ ]] && R=_R"\${suffix: -1}"
        name=\$(basename "\$f")
        rename="$prefix""\$R"_trimming_report.txt
        [[ "\$name" != "\$rename" ]] && mv "\$f" "\$rename"
      done
      """
  }

  // Index genome
  // ############
  if (!index) {

    process overhang {

      input:
        path read from reads_to_overhang.flatten()

      output:
        env overhang into overhangs_to_index

      script:
        """
        \$(gzip -t $read &> /dev/null) && reader=zcat || reader=cat;
        overhang=\$(\$reader $read | head -n 40000 | awk 'NR%4 == 2 {total += length(\$0)} END {print int(total/(NR/4))-1}')
        """
    }

    process index {

      label 'high_cpu'
      label 'high_memory'

      publishDir "$output", mode: 'copy'

      input:
        path genome from genome_to_index
        path annotation from reference_annotation_to_index
        val overhang from overhangs_to_index.max()

      output:
        path 'index' into index_to_map

      script:
        """
        mkdir index
        STAR --runThreadN ${task.cpus} \\
             --runMode genomeGenerate \\
             --genomeDir index \\
             --sjdbGTFfile $annotation \\
             --genomeFastaFiles $genome \\
             --sjdbOverhang $overhang
        """
    }
  }

  // Map reads to genome
  // ###################
  index_to_map.combine(reads_to_map).set {
    reads_to_map
  }

  process map {

    label 'high_cpu'
    label 'high_memory'

    publishDir "$output", mode: 'copy', saveAs: { filename ->
      if (filename.endsWith('.out')) "logs/star/$filename"
      else if (filename.endsWith('.out.tab')) "logs/star/$filename"
      else if (filename.endsWith('.bam')) "maps/$filename"
    }

    input:
      tuple path(index), val(prefix), path(reads) from reads_to_map

    output:
      path '*.out'
      path '*.out.tab'
      tuple val(prefix), path('*.bam') into maps_to_direction

    script:
      """
      STAR --runThreadN ${task.cpus} \\
           --readFilesCommand zcat \\
           --outSAMtype BAM SortedByCoordinate \\
           --genomeDir $index \\
           --readFilesIn $reads \\
           --outFileNamePrefix "$prefix".

      mv *.bam "$prefix".bam
      """
  }

  reference_annotation_to_direction.combine(
    mapped_reads_to_direction.concat(
      maps_to_direction
    )
  ).set {
    maps_to_direction
  }

} else {
  reference_annotation_to_direction.combine(
    mapped_reads_to_direction
  ).set {
    maps_to_direction
  }
}

// Get read directions from maps
// #############################
process direction {

  input:
    tuple path(annotation), val(prefix), path(map) from maps_to_direction

  output:
    tuple val(prefix), env(direction), path(map) into maps_to_length

  script:
    """
    proportions=(\$(infer_library_type.sh $map $annotation))
    difference=\$(awk -v a=\${proportions[0]} -v b=\${proportions[1]} 'BEGIN {print sqrt((a - b)^2)}')
    ratio=\$(awk -v a=\${proportions[0]} -v b=\${proportions[1]} 'BEGIN {print a / b}')
    if [[ \$difference > 50 && \$ratio > 1 ]]; then direction="FR";
    elif [[ \$difference > 50 ]]; then direction="RF";
    else direction="No direction"; fi
    """
}

// Get read lengths from maps
// ##########################
process length {
  input:
    tuple val(prefix), val(direction), path(map) from maps_to_length

  output:
    tuple val(prefix), env(length), val(direction), path(map) into maps_to_merge

  script:
    """
    length=\$(samtools view $map | head -n 10000 | awk '{total += length(\$10)} END {print int(total/NR)}')
    """
}

maps_to_merge.tap {
  maps_to_log
}.map {
  if (it[2] == 'FR') direction = '--fr'
  else if (it[2] == 'RF') direction = '--rf'
  else direction = ''
  [it[0], it[1], direction, it[3]]
}.set {
  maps_to_merge
}

maps_to_log = maps_to_log.toList().get().sort { a, b -> a[0] <=> b[0] }

s = maps_to_log.size() > 1 ? 's' : ''
info = "Proceeding with the following read length$s and direction$s:\n  "
info += maps_to_log.collect{
  it[0] + ' (' + it[1] + ') (' + it[2] + ')'
}.join('\n  ')

log.info info

// Process merge groups
// ####################
if (merge) {

  // Add maps to merge groups
  // ########################
  groups_to_merge = groups_to_merge.toList().get()

  maps_to_merge.map { map ->
    [groups_to_merge.find { group ->
      map[0] in group[1]
    }[0]] + map
  }.groupTuple().tap {
    maps_to_check
  }.map {
    [it[0], it[2][0], it[3][0], it[4]]
  }.set {
    maps_to_merge
  }

  // Check differing read lengths and directions
  // ###########################################
  error = ''

  maps_to_check = maps_to_check.toList().get().sort  { a, b -> a[0] <=> b[0] }

  differing_lengths = maps_to_check.findAll {
    it[2].toUnique().size() > 1
  }.collectEntries { group ->
    [(group[0]): group[1].withIndex().collect { it, i -> it + ' (' + group[2][i] + ')' }]
  }

  if (differing_lengths.size() > 0) {
    error += 'Cannot merge differing read lengths:\n  '
    error += differing_lengths.collect{ it.key + ': ' + it.value.join('  ') }.join('\n  ') + '\n'
  }

  differing_directions = maps_to_check.findAll {
    it[3].toUnique().size() > 1
  }.collectEntries { group ->
    [(group[0]): group[1].withIndex().collect { it, i -> it + ' (' + group[3][i] + ')' }]
  }

  if (differing_directions.size() > 0) {
    error += 'Cannot merge differing read directions:\n  '
    error += differing_directions.collect{ it.key + ': ' + it.value.join('  ') }.join('\n  ') + '\n'
  }

  if (error) exit 1, error

  // Merge maps
  // ##########
  process merge {

    label 'high_cpu'

    input:
      tuple val(prefix), val(length), val(direction), path(maps) from maps_to_merge

    output:
      tuple val(prefix), val(direction), path('*.bam') into maps_to_assemble
      tuple val(prefix), val(length), val(direction), path('*.bam') into maps_to_quantify

    script:
      """
      samtools merge "$prefix".bam $maps --threads ${task.cpus}
      """
  }

} else {
  maps_to_merge.tap {
    maps_to_quantify
  }.map {
    [it[0], it[2], it[3]]
  }.set {
    maps_to_assemble
  }
}

// Assemble transcripts
// ####################
reference_annotation_to_assemble.combine(maps_to_assemble).set {
  maps_to_assemble
}

process assemble {

  input:
    tuple path(annotation), val(prefix), val(direction), path(map) from maps_to_assemble

  output:
    path '*.gff' into assemblies_to_combine

  script:
    """
    stringtie $map $direction -G $annotation -o "$prefix".gff
    """
}

// Combine assemblies
// ##################
process combine {

  publishDir "$output/assembly", mode: 'copy'

  input:
    path annotation from reference_annotation_to_combine
    path assemblies from assemblies_to_combine.collect()

  output:
    path '*.gff' into assembly_annotation_to_quantify

  script:
    """
    stringtie --merge $assemblies -G $annotation -o assembly.gff
    """
}

// Quantify genes and transcripts
// ##############################
reference_annotation_to_quantify.combine(Channel.of('reference')).concat(
  assembly_annotation_to_quantify.combine(Channel.of('assembly'))
).combine(maps_to_quantify).set {
  maps_to_quantify
}

process quantify {

  input:
    tuple path(annotation), val(type), val(prefix), val(length), val(direction), path(map) from maps_to_quantify

  output:
    path('*.reference_genes_TPM.tsv') optional true into reference_genes_TPM_to_format
    path('*.reference_genes_counts.tsv') optional true into reference_genes_counts_to_format
    path('*.reference_transcripts_TPM.tsv') optional true into reference_transcripts_TPM_to_format
    path('*.reference_transcripts_counts.tsv') optional true into reference_transcripts_counts_to_format
    path('*.assembly_genes_TPM.tsv') optional true into assembly_genes_TPM_to_format
    path('*.assembly_genes_counts.tsv') optional true into assembly_genes_counts_to_format
    path('*.assembly_transcripts_TPM.tsv') optional true into assembly_transcripts_TPM_to_format
    path('*.assembly_transcripts_counts.tsv') optional true into assembly_transcripts_counts_to_format

  script:
    """
    stringtie $map \\
              $direction \\
              -e \\
              -B \\
              -G $annotation \\
              -A "$prefix"."$type"_genes_TPM.tsv \\
              -o "$prefix"."$type".gtf

    mv t_data.ctab "$prefix"."$type"_transcripts_TPM.tsv

    mkdir counts
    mv "$prefix"."$type".gtf counts

    prepDE.py -i . \\
              -p counts \\
              -g "$prefix"."$type"_genes_counts.csv \\
              -t "$prefix"."$type"_transcripts_counts.csv \\
              -l $length

    tr ',' '\t' < "$prefix"."$type"_genes_counts.csv > "$prefix"."$type"_genes_counts.tsv

    join <(awk 'BEGIN {FS=","; OFS="\\t"} NR > 1 {print \$1,\$2 | "sort"}' "$prefix"."$type"_transcripts_counts.csv) \\
         <(awk 'BEGIN {OFS="\\t"} NR > 1 {print \$6,\$9 | "sort"}' "$prefix"."$type"_transcripts_TPM.tsv) \\
         -t \$'\\t' \\
         -a 1 \\
         -o "1.1 2.2 1.2" \\
         | cat <(echo "transcript"\$'\\t'"gene"\$'\\t'"counts") - > "$prefix"."$type"_transcripts_counts.tsv
    """
}

// Format results
// ##############
reference_genes_TPM_to_format.tap {
  buffer
}.concat(
  reference_genes_counts_to_format,
  reference_transcripts_TPM_to_format,
  reference_transcripts_counts_to_format,
  assembly_genes_TPM_to_format,
  assembly_genes_counts_to_format,
  assembly_transcripts_TPM_to_format,
  assembly_transcripts_counts_to_format,
).set {
  quantifications_to_format
}

buffer = buffer.count().get()

process format {

  publishDir "$output/quantification", mode: 'copy'

  input:
    path quantifications from quantifications_to_format.buffer(size: buffer)

  output:
    path '*.tsv'

  script:
    """
    format.py $quantifications
    """
}
