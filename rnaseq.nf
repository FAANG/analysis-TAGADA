// Parse and check parameters
// ##########################
output = params.containsKey('output') ? params.output : ''
reads = params.containsKey('reads') ? params.reads : ''
genome = params.containsKey('genome') ? params.genome : ''
index = params.containsKey('index') ? params.index : ''
annotation = params.containsKey('annotation') ? params.annotation : ''
metadata = params.containsKey('metadata') ? params.metadata : ''
merge = params.containsKey('merge') ? params.merge.tokenize(',') : ''
direction = params.containsKey('direction') ? params.'direction' : ''

error = ''

if (!output) error += 'No --output provided.\n'

if (!reads) {
  error += 'No --reads provided.\n'
} else {
  number_of_raw_reads = Channel.fromPath(reads).filter { path ->
    filename = path.getName()
    return filename =~ /\.(fastq|fq)(\.gz)?$/
  }.count().get()

  if (number_of_raw_reads > 0 && !genome && !index) {
    error += 'No --genome or --index provided.\n'
  }
}

if (!annotation) error += 'No --annotation provided.\n'

if (merge && !metadata) {
  error += 'No --metadata provided to execute --merge.\n'
}

if (direction == 'rf') {
  direction = '--rf'
} else if (direction == 'fr') {
  direction = '--fr'
} else if (direction != '') {
  error += 'Invalid --direction. Expected "fr" or "rf".\n'
}

if (error) exit 1, error

// Check index, genome, annotation, metadata
// #########################################
if (index) {
  Channel.fromPath(index, checkIfExists: true).set {
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
  reference_annotation_to_assemble
  reference_annotation_to_combine
  reference_annotation_to_count
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
  r1_reads_to_control
}

reads.R2.into {
  r2_reads_to_check
  r2_reads_to_pair
  r2_reads_to_control
}

reads.single.into {
  single_reads_to_check
  single_reads_to_append
  single_reads_to_control
}

reads.mapped.into {
  mapped_reads_to_check
  mapped_reads_to_merge
}

error = ''
warning = ''

// Check file names
// ################
invalid = reads.invalid.toList().get()

if (invalid.size() > 0) {
  s = invalid.size() > 1 ? 's' : ''
  error += "Wrong format for file name$s:\n  "
  error += invalid.join('\n  ') + '\n\n'
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
  error += unpaired.collect{it['filename']}.join('\n  ') + '\n\n'
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
  error += "Duplicates detected:\n"
  duplicated.each {
    error += '  ' + it.join('  ') + '\n'
  }
  error += '\n'
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
    error += missing_rows.collect{it['filename']}.join('\n  ') + '\n\n'
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
    error += missing_columns.join('\n  ') + '\n\n'
  }

  // Check missing metadata values
  // #############################
  missing_values = metadata_to_check.findAll {
    it.values()[0] in inputs.collect { it['prefix'] }
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
    warning += "Metadata row$s with missing factor$ss will not be merged:\n"
    missing_values.each {
      warning += '  ' + it.join('  ') + '\n'
    }
    warning += '\n'
  }
}

if (error) exit 1, error

if (warning) log.warn '\n' + warning

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

single_reads_to_control.concat(r1_reads_to_control, r2_reads_to_control).map {
  [it['prefix'], it['R'] ? '_R' + it['R'] : '', it['path']]
}.set {
  reads_to_control
}

// Process raw reads
// #################
if (number_of_raw_reads > 0) {

  // Control reads quality
  // #####################
  process control {

    publishDir "$output/quality/raw", mode: 'copy'

    input:
      tuple val(prefix), val(R), path(read) from reads_to_control

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

    publishDir "$output", mode: 'copy', saveAs: { filename ->
      if (filename.endsWith('trimming_report.txt')) "logs/trim_galore/$filename"
      else if (filename.endsWith('fastqc.html')) "quality/trimmed/$filename"
    }

    input:
      tuple val(prefix), val(paired), path(reads) from reads_to_trim

    output:
      path '*_trimming_report.txt'
      path '*_fastqc.html'
      tuple val(prefix), path('*_trimmed.fq.gz') into reads_to_map

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

    process index {

      label 'high_cpu'
      label 'high_memory'

      publishDir "$output", mode: 'copy'

      input:
        path genome from genome_to_index
        path annotation from reference_annotation_to_index

      output:
        path 'index' into index_to_map

      script:
        """
        mkdir index
        STAR --runThreadN ${task.cpus} \\
             --runMode genomeGenerate \\
             --genomeDir index \\
             --sjdbGTFfile $annotation \\
             --genomeFastaFiles $genome
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
      tuple val(prefix), path('*.bam') into maps_to_merge

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

  mapped_reads_to_merge.map {
    [it['prefix'], it['path']]
  }.concat(maps_to_merge).set {
    maps_to_merge
  }

} else {
  mapped_reads_to_merge.map {
    [it['prefix'], it['path']]
  }.set {
    maps_to_merge
  }
}

// Merge maps using metadata
// #########################
if (merge && metadata) {

  metadata_to_merge.cross(maps_to_merge).map {
    missing_values = false
    id = it[0].collectMany { k, v ->
      if (!merge.contains(k)) return []
      else if (!v) missing_values = true
      return [v]
    }.join('_')
    if (missing_values) id = "NA_" + it[1][0]
    [id, it[1][1]]
  }.groupTuple().set {
    maps_to_merge
  }

  process merge {

    input:
      tuple val(prefix), path(maps) from maps_to_merge

    output:
      tuple val(prefix), path('*.bam') into maps_to_assemble
      tuple val(prefix), path('*.bam') into maps_to_count

    script:
      """
      samtools merge "$prefix".bam $maps
      """
  }

} else {
  maps_to_merge.into {
    maps_to_assemble
    maps_to_count
  }
}

// Assemble transcripts
// ####################
reference_annotation_to_assemble.combine(maps_to_assemble).set {
  maps_to_assemble
}

process assemble {

  input:
    tuple path(annotation), val(prefix), path(map) from maps_to_assemble

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

  publishDir "$output/annotation", mode: 'copy'

  input:
    path annotation from reference_annotation_to_combine
    path assemblies from assemblies_to_combine.collect()

  output:
    path '*.gff' into assembly_annotation_to_count

  script:
    """
    stringtie --merge $assemblies -G $annotation -o assembly.gff
    """
}

// Count genes and transcripts
// ###########################
reference_annotation_to_count.combine(Channel.of('reference')).concat(
  assembly_annotation_to_count.combine(Channel.of('assembly'))
).combine(
  maps_to_count
).set {
  maps_to_count
}

process count {

  input:
    tuple path(annotation), val(type), val(prefix), path(map) from maps_to_count

  output:
    tuple path('*.reference_genes.tsv'), path('*.reference_transcripts.tsv') optional true into reference_counts_to_format
    tuple path('*.assembly_genes.tsv'), path('*.assembly_transcripts.tsv') optional true into assembly_counts_to_format

  script:
    """
    stringtie $map \\
              $direction \\
              -e \\
              -B \\
              -G $annotation \\
              -A "$prefix"."$type"_genes.tsv

    mv t_data.ctab "$prefix"."$type"_transcripts.tsv
    """
}

// Format count results
// ####################
reference_counts_to_format.tap {
  buffer
}.concat(assembly_counts_to_format).multiMap {
  genes: it[0]
  transcripts: it[1]
}.set {
  counts_to_format
}

counts_to_format.genes.concat(counts_to_format.transcripts).set {
  counts_to_format
}

buffer = buffer.count().get()

process format {

  publishDir "$output/counts", mode: 'copy'

  input:
    path counts from counts_to_format.buffer(size: buffer)

  output:
    path '*.tsv'

  script:
    """
    merge.py $counts
    """
}
