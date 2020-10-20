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
  reference_annotation_to_get_direction
  reference_annotation_to_assemble
  reference_annotation_to_combine
  reference_annotation_to_quantify
  reference_annotation_to_control_elements
  reference_annotation_to_control_exons
}

if (metadata) {
  Channel.fromPath(
    metadata,
    checkIfExists: true
  ).tap {
    metadata_to_report
  }.splitCsv(header: true, sep: '\t').into {
    metadata_to_check
    metadata_to_merge
  }
} else {
  Channel.empty().into {
    metadata_to_report
    metadata_to_check
    metadata_to_merge
  }
}

// Split reads into R1/R2/single/mapped
// ####################################
Channel.fromPath(reads, checkIfExists: true).map { path ->

  filename = path.getName()

  return (
    filename =~ /^(.+?)(?:[\._ ][Rr]?([12]))?(?:\.(fastq|fq|gz|bam))+$/
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
  r1_reads_to_control_quality
}

reads.R2.into {
  r2_reads_to_check
  r2_reads_to_pair
  r2_reads_to_control_quality
}

reads.single.into {
  single_reads_to_check
  single_reads_to_append
  single_reads_to_control_quality
}

reads.mapped.tap {
  mapped_reads_to_check
}.map {
  [it['prefix'], it['path']]
}.set {
  mapped_reads_to_sort
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
  raw_reads_to_trim
}

single_reads_to_control_quality.concat(
  r1_reads_to_control_quality,
  r2_reads_to_control_quality
).map {
  [it['path'], 'raw']
}.set {
  raw_reads_to_control_quality
}

// Sort input maps
// ###############
process sort {

  label 'cpu_16'
  label 'memory_16'

  publishDir "$output/maps", mode: 'copyNoFollow'

  input:
    tuple val(prefix), path(map, stageAs: 'input') from mapped_reads_to_sort

  output:
    tuple val(prefix), path('*.bam') into maps_to_get_direction
    tuple val(prefix), path('*.bam') into maps_to_control_contigs
    tuple val(prefix), path('*.bam') into maps_to_control_metrics
    tuple val(prefix), path('*.bam') into maps_to_control_flags

  script:
    """
    if [[ \$(grep SO:unsorted <(samtools view -H input | head -n 1)) ]]; then
      samtools sort -@ ${task.cpus} -o "$prefix".bam input
    else
      target=\$(readlink -m input)
      ln -sf "\$target" "$prefix".bam
    fi
    """
}

// Process raw reads
// #################
if (number_of_raw_reads > 0) {

  // Trim adaptators
  // ###############
  process trim {

    label 'cpu_16'
    label 'memory_16'

    input:
      tuple val(prefix), val(paired), path(reads) from raw_reads_to_trim

    output:
      path '*_trimming_report.txt' into trim_to_report
      path '* (trimmed)*' into trimmed_reads_to_get_overhang
      path '* (trimmed)*' into trimmed_reads_to_control_quality
      tuple val(prefix), path('* (trimmed)*') into trimmed_reads_to_map

    script:
      """
      for f in $reads; do
        name=\$(basename "\$f")
        rename="\${name// /_}"
        [[ "\$name" != "\$rename" ]] && mv "\$f" "\$rename"
        reads+=("\$rename")
      done

      trim_galore "\${reads[@]}" \\
                  $paired \\
                  --cores ${task.cpus} \\
                  --gzip \\
                  --basename "$prefix"

      for f in *_trimmed.*; do
        [ -f "\$f" ] || break
        mv "\$f" "$prefix"" (trimmed)""\${f##*_trimmed}"
      done

      for n in 1 2; do
        for f in *_val_\$n.*; do
          [ -f "\$f" ] || break
          i=\$((n - 1))
          prefix="\${reads[i]}"
          while [[ "\$prefix" =~ \\.(fastq|fq|gz)\$ ]]; do
            prefix="\${prefix%.gz}"
            prefix="\${prefix%.fq}"
            prefix="\${prefix%.fastq}"
          done
          mv "\$f" "\$prefix"" (trimmed)""\${f##*_val_\$n}"
        done
      done
      """
  }

  // Control reads quality
  // #####################
  raw_reads_to_control_quality.concat(
    trimmed_reads_to_control_quality.flatten().map {
      [it, null]
    }
  ).set {
    reads_to_control_quality
  }

  process control_quality {

    input:
      tuple path(read), val(suffix) from reads_to_control_quality

    output:
      path '*_fastqc.zip' into control_quality_to_report

    script:
      """
      if [ $suffix != null ]; then
        prefix=$read
        while [[ "\$prefix" =~ \\.(fastq|fq|gz)\$ ]]; do
          prefix="\${prefix%.gz}"
          prefix="\${prefix%.fq}"
          prefix="\${prefix%.fastq}"
        done
        read="\$prefix"" ("$suffix").fq"
        \$(gzip -t $read &> /dev/null) && read="\$read".gz
        [[ $read != "\$read" ]] && mv $read "\$read"
        fastqc "\$read"
      else
        fastqc $read
      fi
      """
  }

  // Index genome
  // ############
  if (!index) {

    process get_overhang {

      input:
        path read from trimmed_reads_to_get_overhang.flatten()

      output:
        env overhang into overhangs_to_index

      script:
        """
        overhang=\$(zcat $read | head -n 40000 | awk 'NR%4 == 2 {total += length(\$0)} END {print int(total/(NR/4))-1}')
        """
    }

    process index {

      label 'cpu_16'
      label 'memory_32'

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
  index_to_map.combine(trimmed_reads_to_map).set {
    trimmed_reads_to_map
  }

  process map {

    label 'cpu_16'
    label 'memory_32'

    publishDir "$output/maps", mode: 'copy', saveAs: { filename ->
      if (filename.endsWith('.bam')) filename
    }

    input:
      tuple path(index), val(prefix), path(reads) from trimmed_reads_to_map

    output:
      path '*.bam'
      path '*.Log.final.out' into map_to_report
      tuple val(prefix), path('*.bam') into mapped_reads_to_get_direction
      tuple val(prefix), path('*.bam') into mapped_reads_to_control_contigs
      tuple val(prefix), path('*.bam') into mapped_reads_to_control_metrics
      tuple val(prefix), path('*.bam') into mapped_reads_to_control_flags

    script:
      """
      STAR --runThreadN ${task.cpus} \\
           --readFilesCommand zcat \\
           --outSAMtype BAM SortedByCoordinate \\
           --outFilterIntronMotifs RemoveNoncanonical \\
           --genomeDir $index \\
           --readFilesIn $reads \\
           --outFileNamePrefix "$prefix".

      mv *.bam "$prefix".bam
      """
  }

  reference_annotation_to_get_direction.combine(
    maps_to_get_direction.concat(
      mapped_reads_to_get_direction
    )
  ).set {
    maps_to_get_direction
  }

  maps_to_control_contigs.concat(
    mapped_reads_to_control_contigs
  ).set {
    maps_to_control_contigs
  }

  maps_to_control_metrics.concat(
    mapped_reads_to_control_metrics
  ).set {
    maps_to_control_metrics
  }

  maps_to_control_flags.concat(
    mapped_reads_to_control_flags
  ).set {
    maps_to_control_flags
  }

} else {
  reference_annotation_to_get_direction.combine(
    maps_to_get_direction
  ).set {
    maps_to_get_direction
  }

  Channel.empty().into {
    control_quality_to_report
    trim_to_report
    map_to_report
  }
}

// Control mapped reads per contig
// ###############################
process control_contigs {

  input:
    tuple val(prefix), path(map) from maps_to_control_contigs

  output:
    path('*.idxstats') into control_contigs_to_report

  script:
    """
    samtools index -@ \$((${task.cpus} - 1)) $map
    samtools idxstats $map > "$prefix".idxstats
    """
}

// Control mapping metrics
// #######################
process control_metrics {

  input:
    tuple val(prefix), path(map) from maps_to_control_metrics

  output:
    path('*.stats') into control_metrics_to_report

  script:
    """
    samtools stats -@ \$((${task.cpus} - 1)) $map > "$prefix".stats
    """
}

// Control mapping flags
// #####################
process control_flags {

  input:
    tuple val(prefix), path(map) from maps_to_control_flags

  output:
    path('*.flagstat') into control_flags_to_report

  script:
    """
    samtools flagstat -@ \$((${task.cpus} - 1)) $map > "$prefix".flagstat
    """
}


// Get read directions from maps
// #############################
process get_direction {

  input:
    tuple path(annotation), val(prefix), path(map) from maps_to_get_direction

  output:
    tuple val(prefix), env(direction), path(map) into maps_to_get_length
    tuple val(prefix), env(direction), path(map) into maps_to_coverage

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
process get_length {

  input:
    tuple val(prefix), val(direction), path(map) from maps_to_get_length

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
  [it[0], it[1].toInteger(), direction, it[3]]
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

// Get read coverage from maps
// ###########################
process coverage {

  label 'memory_4'

  publishDir "$output/coverage", mode: 'copy'

  input:
    tuple val(prefix), val(direction), path(map) from maps_to_coverage

  output:
    path '*.bed'

  script:
    """
    if [[ "$direction" == "RF" || "$direction" == "FR" ]]; then
      bedtools genomecov -ibam $map -bg -strand + > +.tsv
      bedtools genomecov -ibam $map -bg -strand - > -.tsv
      cat <(awk 'BEGIN {OFS="\\t"} {print \$0,"+"}' +.tsv) \\
          <(awk 'BEGIN {OFS="\\t"} {print \$0,"-"}' -.tsv) | sort -T "." -k1,3 -k5 \\
          > "$prefix".bed
    else
      bedtools genomecov -ibam $map -bg > "$prefix".bed
    fi
    """
}

// Process merge groups
// ####################
if (merge) {

  // Add maps to merge groups
  // ########################
  groups_to_merge = groups_to_merge.toList().get()

  maps_to_merge.map { map ->
    [
      groups_to_merge.find { group ->
        map[0] in group[1]
      }[0]
    ] + map
  }.groupTuple().tap {
    maps_to_check
  }.map {
    [it[0], (it[2].sum()/it[2].size()).toInteger(), it[3][0], it[4]]
  }.set {
    maps_to_merge
  }

  // Check differing read lengths and directions
  // ###########################################
  error = ''

  maps_to_check = maps_to_check.toList().get().sort { a, b -> a[0] <=> b[0] }

  differing_lengths = maps_to_check.findAll {
    it[2].subsequences().findAll {
      it.size() == 2
    }.collect {
      (it[0] - it[1]).abs()
    }.findAll {
      it > 10
    }.size() > 0
  }.collectEntries { group ->
    [
      (group[0]): group[1].withIndex().collect { it, i ->
        it + ' (' + group[2][i] + ')'
      }
    ]
  }

  if (differing_lengths.size() > 0) {
    error += 'Cannot merge differing read lengths:\n  '
    error += differing_lengths.collect {
      it.key + ': ' + it.value.join('  ')
    }.join('\n  ') + '\n'
  }

  differing_directions = maps_to_check.findAll {
    it[3].toUnique().size() > 1
  }.collectEntries { group ->
    [
      (group[0]): group[1].withIndex().collect { it, i ->
        it + ' (' + group[3][i] + ')'
      }
    ]
  }

  if (differing_directions.size() > 0) {
    error += 'Cannot merge differing read directions:\n  '
    error += differing_directions.collect{ it.key + ': ' + it.value.join('  ') }.join('\n  ') + '\n'
  }

  if (error) exit 1, error

  // Merge maps
  // ##########
  process merge {

    label 'cpu_16'

    input:
      tuple val(prefix), val(length), val(direction), path(maps) from maps_to_merge

    output:
      tuple val(prefix), val(direction), path('*.bam') into maps_to_assemble
      tuple val(prefix), val(length), val(direction), path('*.bam') into maps_to_quantify
      path '*.bam' into maps_to_control_exons

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
  }.tap {
    maps_to_assemble
  }.map {
    it[2]
  }.set {
    maps_to_control_exons
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
    path '*.gff' into assembly_annotation_to_control_elements
    path '*.gff' into assembly_annotation_to_control_exons

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
    path 'assembly_transcripts_TPM.tsv' optional true into formatted_transcripts_to_control_elements
    path 'assembly_genes_TPM.tsv' optional true into formatted_genes_to_control_elements
    path 'reference_genes_TPM.tsv' optional true into formatted_genes_tpm_to_control_expression
    path 'reference_genes_counts.tsv' optional true into formatted_genes_counts_to_control_expression

  script:
    """
    format.py $quantifications
    """
}

// Control elements detected
// #########################
process control_elements {

  publishDir "$output/control/elements", mode: 'copy'

  input:
    path reference_annotation from reference_annotation_to_control_elements
    path assembly_annotation from assembly_annotation_to_control_elements
    path formatted_transcripts from formatted_transcripts_to_control_elements
    path formatted_genes from formatted_genes_to_control_elements

  output:
    path '*.png' into control_elements_to_report

  script:
    """
    detected_elements_sumstats.sh \\
      $reference_annotation \\
      $assembly_annotation \\
      $formatted_transcripts \\
      $formatted_genes

    for f in \\
      Plots/ExonPerTranscript/*.png \\
      Plots/ExonLength/*.png \\
      Plots/DistinctExonLength/*.png \\
      Plots/5pExonLength_Tr/*.png \\
      Plots/InternalExonLength/*.png \\
      Plots/5pExonLength_Gn/*.png \\
      Plots/DistinctInternalExonLength/*.png \\
      Plots/TrLength/*.png \\
      Plots/cDNALength/*.png; \\
    do
      [ -f "\$f" ] || break
      convert "\$f" -crop 1400x2100+0+0 +repage "\${f%.png}"_cropped.png
    done

    convert Plots/ExonPerTranscript/*_cropped.png \\
            Plots/TranscriptPerGene/prediction_nbtringn_forggplot.png \\
            +append \\
            +repage \\
            expertr_trpergn.png

    convert Plots/ExonLength/*_cropped.png \\
            Plots/DistinctExonLength/*_cropped.png \\
            Plots/MonoExTrExLength/*.png \\
            +append \\
            +repage \\
            exlg_distinctexlg_monoextrlg.png

    convert Plots/5pExonLength_Tr/*_cropped.png \\
            Plots/InternalExonLength/*_cropped.png \\
            Plots/3pExonLength_Tr/*.png \\
            +append \\
            +repage \\
            5pexlgpertr_internexlg_3pexlgpertr.png

    convert Plots/5pExonLength_Gn/*_cropped.png \\
            Plots/DistinctInternalExonLength/*_cropped.png \\
            Plots/3pExonLength_Gn/*.png \\
            +append \\
            +repage \\
            5pexlgpergn_distinctinternexlg_3pexlgpergn.png

    convert Plots/TrLength/*_cropped.png \\
            Plots/cDNALength/*_cropped.png \\
            Plots/Exact_tr_dist_to_Genc_TSS/*.png \\
            +append \\
            +repage \\
            trlg_cdnalg_exacttrdisttoreftss.png
    """
}

// Control reference gene expression detected
// ##########################################
process control_expression {

  publishDir "$output/control/expression", mode: 'copy'

  input:
    path formatted_genes_tpm from formatted_genes_tpm_to_control_expression
    path formatted_genes_counts from formatted_genes_counts_to_control_expression

  output:
    path '*.png' into control_expression_to_report

  script:
    """
    awk 'BEGIN{OFS="\\t"; print "labExpId"} NR==1{for(i=2; i<=NF; i++){print \$i}}' $formatted_genes_tpm > metadata.tsv
    plot_gene_expression.sh \\
      $formatted_genes_tpm \\
      $formatted_genes_counts \\
      metadata.tsv labExpId .
    """
}

// Control exonic reads counts
// ###########################
reference_annotation_to_control_exons.combine(Channel.of('reference')).concat(
  assembly_annotation_to_control_exons.combine(Channel.of('assembly'))
).set {
  annotations_to_control_exons
}

process control_exons {

  label 'cpu_16'

  publishDir "$output/control/exons", mode: 'copy'

  input:
    tuple path(annotation), val(type) from annotations_to_control_exons
    path(maps) from maps_to_control_exons.collect()

  output:
    path '*.tsv'
    path '*.summary' into control_exons_to_report

  script:
    """
    featureCounts -t exon \\
                  -g gene_id \\
                  -s 0 \\
                  --primary \\
                  -T ${task.cpus} \\
                  -a $annotation \\
                  -o "$type"_exons_counts.tsv \\
                  $maps

    mv "$type"_exons_counts.tsv.summary "$type"_exons_counts.summary
    sed -i -e '2s/\\.bam\\b//g' "$type"_exons_counts.tsv
    sed -i -e '1s/\\.bam\\b//g' "$type"_exons_counts.summary
    """
}

// Create multiqc report
// #####################
config_to_report = Channel.fromPath("$baseDir/multiqc.yaml", checkIfExists: true)

process report {

  publishDir "$output/control", mode: 'copy'

  input:
    path config from config_to_report
    path '*' from metadata_to_report.flatten().collect().ifEmpty([])
    path '*' from control_quality_to_report.flatten().collect().ifEmpty([])
    path '*' from trim_to_report.flatten().collect().ifEmpty([])
    path '*' from map_to_report.flatten().collect().ifEmpty([])
    path '*' from control_elements_to_report.flatten().collect().ifEmpty([])
    path '*' from control_expression_to_report.flatten().collect().ifEmpty([])
    path '*' from control_exons_to_report.flatten().collect().ifEmpty([])
    path '*' from control_contigs_to_report.flatten().collect().ifEmpty([])
    path '*' from control_metrics_to_report.flatten().collect().ifEmpty([])
    path '*' from control_flags_to_report.flatten().collect().ifEmpty([])

  output:
    path 'multiqc_report.html'

  script:
    """
    for f in *.tsv *.png; do
      [ -f "\$f" ] || break
      mv "\$f" "\${f%.*}"_mqc."\${f##*.}"
    done
    multiqc --config $config .
    """
}