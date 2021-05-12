// Parse and check parameters
// ##########################
output = params.containsKey('output') ? params.output : ''
reads = params.containsKey('reads') ? params.reads : ''
genome = params.containsKey('genome') ? params.genome : ''
index = params.containsKey('index') ? params.index : ''
annotation = params.containsKey('annotation') ? params.annotation : ''
metadata = params.containsKey('metadata') ? params.metadata : ''
merge = params.containsKey('merge') ? params.merge.tokenize(',') : ''
feelnc_args = params.containsKey('feelnc-args') ? params.'feelnc-args' : ''
skip_feelnc = params.containsKey('skip-feelnc') ? true : false

error = ''

if (!output) error += 'No --output provided\n'

if (!reads) error += 'No --reads provided\n'

if (!genome) error += 'No --genome provided\n'

if (!annotation) error += 'No --annotation provided\n'

if (merge && !metadata) {
  error += 'No --metadata provided to execute --merge\n'
}

if (error) exit 1, error

// Check index, genome, annotation, metadata
// #########################################
if (index) {
  Channel.fromPath(
    index,
    type: 'dir',
    checkIfExists: true
  ).map { path ->
    filename = path.getName()
    if (filename.endsWith('.gz'))
      filename = filename.substring(0, filename.length() - 3)
    if (filename.endsWith('.tar'))
      filename = filename.substring(0, filename.length() - 4)
    return [path, filename]
  }.set {
    index_to_decompress
  }
} else {
  Channel.empty().set {
    index_to_decompress
  }
}

Channel.fromPath(
  genome,
  checkIfExists: true
).map { path ->
  filename = path.getName()
  if (filename.endsWith('.gz'))
    filename = filename.substring(0, filename.length() - 3)
  if (filename.endsWith('.tar'))
    filename = filename.substring(0, filename.length() - 4)
  return [path, filename]
}.set {
  genome_to_decompress
}

Channel.fromPath(
  annotation,
  checkIfExists: true
).map { path ->
  filename = path.getName()
  if (filename.endsWith('.gz'))
    filename = filename.substring(0, filename.length() - 3)
  if (filename.endsWith('.tar'))
    filename = filename.substring(0, filename.length() - 4)
  return [path, filename]
}.set {
  reference_annotation_to_decompress
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
if (reads.endsWith('.txt')) {
  Channel.fromPath(reads, checkIfExists: true).splitText().map {
    file(it.strip(), checkIfExists: true)
  }.into {
    reads_to_count
    reads_to_split
  }
} else {
  Channel.fromPath(reads, checkIfExists: true).into {
    reads_to_count
    reads_to_split
  }
}

number_of_raw_reads = reads_to_count.filter { path ->
  filename = path.getName()
  return filename =~ /\.(fastq|fq)(\.gz)?$/
}.count().get()

reads_to_split.map { path ->

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

// Decompress index, genome, annotation
// ####################################
process decompress {

  input:
    tuple path(index), val(index_name) from index_to_decompress.ifEmpty([file('   '),''])
    tuple path(genome), val(genome_name) from genome_to_decompress.ifEmpty([file('  '),''])
    tuple path(annotation), val(annotation_name) from reference_annotation_to_decompress.ifEmpty([file(' '),''])

  output:
    path "$index_name" optional true into index_to_map
    path "$genome_name" optional true into (
      genome_to_index,
      genome_to_detect_lncRNA
    )
    path "$annotation_name" optional true into (
      reference_annotation_to_index,
      reference_annotation_to_direction,
      reference_annotation_to_assemble,
      reference_annotation_to_combine,
      reference_annotation_to_quantify,
      reference_annotation_to_detect_lncRNA,
      reference_annotation_to_control_elements,
      reference_annotation_to_control_exons
    )

  script:
    """
    for f in $index $genome $annotation; do
      if [[ \${f: -7} == .tar.gz ]]; then
        tar -xf "\$(readlink -m "\$f")"
      elif [[ \${f: -3} == .gz ]]; then
        gzip -c -d "\$(readlink -m "\$f")" > "\${f%.gz}"
      fi
    done
    """
}

// Sort input maps
// ###############
process sort {

  label 'cpu_16'
  label 'memory_16'

  publishDir "$output/maps", mode: 'copyNoFollow', overwrite: false

  input:
    tuple val(prefix), path(map, stageAs: 'input') from mapped_reads_to_sort

  output:
    tuple val(prefix), path('*.bam') into maps_to_direction
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
      path '* (trimmed)*' into trimmed_reads_to_overhang
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
        [ -f "\$f" ] || continue
        mv "\$f" "$prefix"" (trimmed)""\${f##*_trimmed}"
      done

      for n in 1 2; do
        for f in *_val_\$n.*; do
          [ -f "\$f" ] || continue
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

    process overhang {

      input:
        path read from trimmed_reads_to_overhang.flatten()

      output:
        env overhang into overhangs_to_index

      script:
        """
        overhang=\$(zcat $read | head -n 40000 | awk 'NR%4 == 2 {total += length(\$0)} END {print int(total/(NR/4))-1}')
        """
    }

    process index {

      label 'cpu_16'
      label 'memory_64'

      publishDir "$output", mode: 'copy'

      input:
        path genome from genome_to_index
        path annotation from reference_annotation_to_index
        val overhang from overhangs_to_index.max()

      output:
        path 'index' into indexed_genome_to_map

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

    indexed_genome_to_map.set {
      index_to_map
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
      tuple val(prefix), path('*.bam') into mapped_reads_to_direction
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

  reference_annotation_to_direction.combine(
    maps_to_direction.concat(
      mapped_reads_to_direction
    )
  ).set {
    maps_to_direction
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
  reference_annotation_to_direction.combine(
    maps_to_direction
  ).set {
    maps_to_direction
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
process direction {

  input:
    tuple path(annotation), val(prefix), path(map) from maps_to_direction

  output:
    tuple val(prefix), env(direction), path(map) into maps_to_length
    tuple val(prefix), env(direction), path(map) into maps_to_coverage

  script:
    """
    proportions=(\$(infer_library_type.sh $map $annotation))
    difference=\$(awk -v a=\${proportions[0]} -v b=\${proportions[1]} 'BEGIN {print sqrt((a - b)^2)}')
    if [[ \$difference > 50 && \$proportions[0] > \$proportions[1] ]]; then direction="FR";
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

  publishDir "$output/annotation", mode: 'copy'

  input:
    path annotation from reference_annotation_to_combine
    path assemblies from assemblies_to_combine.collect()

  output:
    path 'novel.gff' into (
      novel_annotation_to_detect_lncRNA,
      novel_annotation_to_quantify,
      novel_annotation_to_control_elements,
      novel_annotation_to_control_exons
    )

  script:
    """
    stringtie --merge $assemblies -G $annotation -o novel_no_biotype.gff

    # Add lines for genes
    awk -f /usr/local/src/Scripts/compute_boundaries.awk -v toadd=gene -v fldno=10 -v keys=gene_name,ref_gene_id novel_no_biotype.gff > genes.gff

    cat genes.gff novel_no_biotype.gff | sort -k1,1 -k4,4n -k5,5rn > novel.genes.gff

    # Write transcript biotype in merged assembly
    awk 'BEGIN { FS = "\t" }
        NR==FNR {
          match(\$9,/transcript_id "([^;]*)";*/,tId)
          match(\$9,/transcript_biotype "([^;]*)";*/,biotype)
          biotypes[tId[1]]=biotype[1]
          next
        }
        {
          if (substr(\$1,1,1)!="#") {
            match(\$9,/transcript_id "([^;]*)";*/,tId)
            if (tId[1] in biotypes) {
              print \$0 "transcript_biotype \\""biotypes[tId[1]]"\\";"
            } else {
              print \$0
            }
          } else {
            print \$0
          }
        }  ' $annotation novel.genes.gff > novel.gff


    """
}

// Detect long non-coding RNAs
// #############################
process detect_lncRNA {

  label 'memory_16'

  publishDir "$output", mode: 'copy', saveAs: { filename ->
    if (filename == 'lncRNA_classes.txt') "assembly/lnc_classification/$filename"
    else if (filename.startsWith('exons.')) "assembly/lnc_classification/$filename"
    else if (filename.endsWith('.gtf')) "assembly/$filename"
    else if (filename.endsWith('.gff')) "assembly/$filename"
    else "control/lnc/$filename"
  }

  when:
    !skip_feelnc

  input:
    path reference_annotation from reference_annotation_to_detect_lncRNA
    path novel_annotation from novel_annotation_to_detect_lncRNA
    path genome from genome_to_detect_lncRNA

  output:
    path '*classes.txt' into feelnc_classes_to_report
    path 'exons.*.gtf'
    path 'exons_RF_summary.txt' into feelnc_rf_summary_to_report
    path '*.feelncclassifier.log' into feelnc_classifier_log_to_report
    path 'assembly.feelnc_biotype.gff' into assembly_annotation_with_coding_potential

  script:
    """
    script=\$(which FEELnc_codpot.pl)
    export FEELNCPATH=\${script%/*}/..

    FEELnc_filter.pl --mRNAfile $reference_annotation \\
                     --infile $novel_annotation \\
                     --biotype transcript_biotype=protein_coding \\
                     > candidate_transcripts.gtf

    FEELnc_codpot.pl --genome $genome \\
                     --mRNAfile $reference_annotation \\
                     --infile candidate_transcripts.gtf \\
                     --biotype transcript_biotype=protein_coding \\
                     --numtx 5000,5000 \\
                     --kmer 1,2,3,6,9,12 \\
                     --outdir . \\
                     --outname exons \\
                     --mode shuffle \\
                     --spethres=0.98,0.98 \\
                     $feelnc_args

    FEELnc_classifier.pl --mrna $reference_annotation \\
                         --lncrna exons.lncRNA.gtf \\
                         > lncRNA_classes.txt

    # Enrich assembled annotation with new biotypes
    cp $novel_annotation assembly.feelnc_biotype.gff
    for biotype in lncRNA mRNA noORF TUCp
      do if [ -f exons.\$biotype.gtf ] ; then
        awk -v biotype=\$biotype '
        BEGIN { FS = "\t" }
        NR==FNR {
        match(\$9,/transcript_id "([^;]*)";*/,tId)
              transcripts[tId[1]]=0
              next
        }
        #(substr(\$1,1,1)=="#"){print \$0}
        {match(\$9,/transcript_id "([^;]*)";*/,tId)
        if (tId[1] in transcripts) {
              # Check if there is already a biotype in the annotation
              match(\$9,/biotype=([^;]*)*/,oldBiotype)
              if (oldBiotype[1]){
                      print \$0
              } else {
                      print \$0 " feelnc_biotype \\"" biotype "\\";"
              }
        } else { print \$0 }
      }' exons."\$biotype".gtf assembly.feelnc_biotype.gff > tmp.gff
      mv tmp.gff assembly.feelnc_biotype.gff
      fi
    done

    """
}

// Quantify genes and transcripts
// ##############################
reference_annotation_to_quantify.combine(Channel.of('reference')).concat(
  novel_annotation_to_quantify.combine(Channel.of('novel'))
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
    path('*.novel_genes_TPM.tsv') optional true into novel_genes_TPM_to_format
    path('*.novel_genes_counts.tsv') optional true into novel_genes_counts_to_format
    path('*.novel_transcripts_TPM.tsv') optional true into novel_transcripts_TPM_to_format
    path('*.novel_transcripts_counts.tsv') optional true into novel_transcripts_counts_to_format

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
  novel_genes_TPM_to_format,
  novel_genes_counts_to_format,
  novel_transcripts_TPM_to_format,
  novel_transcripts_counts_to_format,
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
    path 'novel_transcripts_TPM.tsv' optional true into formatted_transcripts_to_control_elements
    path 'novel_genes_TPM.tsv' optional true into formatted_genes_to_control_elements
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

  publishDir "$output/control/elements", mode: 'copy', saveAs: { filename ->
    if (
      filename.endsWith('.png') ||
      filename.endsWith('.txt') ||
      filename == 'Plots' ||
      filename == 'Tables'
    ) filename
  }

  input:
    path reference_annotation from reference_annotation_to_control_elements
    path novel_annotation from novel_annotation_to_control_elements
    path formatted_transcripts from formatted_transcripts_to_control_elements
    path formatted_genes from formatted_genes_to_control_elements

  output:
    path '*.png'
    path 'Plots'
    path 'Tables'
    path '*_expressed_*.txt'
    path '*_annotation.tsv' into control_elements_to_report

  script:
    """
    detected_elements_sumstats.sh \\
      $reference_annotation \\
      $novel_annotation \\
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
      Plots/Exact_tr_dist_to_Genc_TSS/*.png \\
      Plots/TrLength/*.png \\
    do
      [ -f "\$f" ] || continue
      convert "\$f" -crop 1400x2100+0+0 +repage "\${f%.png}"_cropped.png
    done

    convert Plots/ExonPerTranscript/*_cropped.png \\
            Plots/TranscriptPerGene/prediction_nbtringn_forggplot.png \\
            +append \\
            +repage \\
            exon_per_transcript_and_transcript_per_gene.png

    convert Plots/ExonLength/*_cropped.png \\
            Plots/DistinctExonLength/*_cropped.png \\
            Plots/MonoExTrExLength/*.png \\
            +append \\
            +repage \\
            all_distinct_exon_length_and_monoexonic_transcript_length.png

    convert Plots/5pExonLength_Tr/*_cropped.png \\
            Plots/InternalExonLength/*_cropped.png \\
            Plots/3pExonLength_Tr/*.png \\
            +append \\
            +repage \\
            5p_internal_3p_exon_length_per_transcript.png

    convert Plots/5pExonLength_Gn/*_cropped.png \\
            Plots/DistinctInternalExonLength/*_cropped.png \\
            Plots/3pExonLength_Gn/*.png \\
            +append \\
            +repage \\
            5p_internal_3p_exon_length_per_gene.png

    convert Plots/Exact_tr_dist_to_Genc_TSS/*_cropped.png \\
            Plots/TrLength/*_cropped.png \\
            Plots/cDNALength/*.png \\
            +append \\
            +repage \\
            transcript_cdna_length_and_TSStorefgeneTSS_distance_for_exact_transcripts.png

    reference_genes=\$(awk 'NR == 3 {print \$2}' detected_transcripts_genes_numbers.tsv)
    novel_genes=\$(awk 'NR == 3 {print \$2 + \$7}' detected_transcripts_genes_numbers.tsv)
    reference_transcripts=\$(awk 'NR == 2 {print \$2}' detected_transcripts_genes_numbers.tsv)
    novel_transcripts=\$(awk 'NR == 2 {print \$2 + \$7}' detected_transcripts_genes_numbers.tsv)

    perl -pe 's/^"([^"]+)".+\$/\$1/g' \\
      ref_expr/ref.annot.tpm0.1.2samples.exons_complete_gnid_nbtr.txt \\
      > reference_expressed_genes.txt
    reference_expressed_genes=\$(wc -l reference_expressed_genes.txt | awk '{print \$1}')

    perl -pe 's/^"([^"]+)".+\$/\$1/g' \\
      string_expr/stringtie.annot.tpm0.1.2samples.exons_complete_gnid_nbtr.txt \\
      > novel_expressed_genes.txt
    novel_expressed_genes=\$(wc -l novel_expressed_genes.txt | awk '{print \$1}')

    perl -pe 's/^"([^"]+)".+\$/\$1/g' \\
      ref_expr/ref.annot.tpm0.1.2samples.exons_complete_trid_nbex.txt \\
      > reference_expressed_transcripts.txt
    reference_expressed_transcripts=\$(wc -l reference_expressed_transcripts.txt | awk '{print \$1}')

    perl -pe 's/^"([^"]+)".+\$/\$1/g' \\
      string_expr/stringtie.annot.tpm0.1.2samples.exons_complete_trid_nbex.txt \\
      > novel_expressed_transcripts.txt
    novel_expressed_transcripts=\$(wc -l novel_expressed_transcripts.txt | awk '{print \$1}')

    percent_novel_expressed_genes=\$(
      echo | \\
      awk -v expressed=\$novel_expressed_genes \\
          -v all=\$novel_genes \\
          '{print 100 * expressed / all}'
    )
    percent_reference_expressed_genes=\$(
      echo | \\
      awk -v expressed=\$reference_expressed_genes \\
          -v all=\$reference_genes \\
          '{print 100 * expressed / all}'
    )
    percent_novel_expressed_transcripts=\$(
      echo | \\
      awk -v expressed=\$novel_expressed_transcripts \\
          -v all=\$novel_transcripts \\
          '{print 100 * expressed / all}'
    )
    percent_reference_expressed_transcripts=\$(
      echo | \\
      awk -v expressed=\$reference_expressed_transcripts \\
          -v all=\$reference_transcripts \\
          '{print 100 * expressed / all}'
    )

    echo -e "Category\tTotal\tPercentage" > reference_annotation.tsv
    echo -e "Genes\t\$reference_genes\t" >> reference_annotation.tsv
    echo -e "Expressed genes\t\$reference_expressed_genes\t\$percent_reference_expressed_genes" >> reference_annotation.tsv
    echo -e "Transcripts\t\$reference_transcripts\t" >> reference_annotation.tsv
    echo -e "Expressed transcripts\t\$reference_expressed_transcripts\t\$percent_reference_expressed_transcripts" >> reference_annotation.tsv

    echo -e "Category\tTotal\tPercentage" > novel_annotation.tsv
    echo -e "Genes\t\$novel_genes\t" >> novel_annotation.tsv
    echo -e "Expressed genes\t\$novel_expressed_genes\t\$percent_novel_expressed_genes" >> novel_annotation.tsv
    echo -e "Transcripts\t\$novel_transcripts\t" >> novel_annotation.tsv
    echo -e "Expressed transcripts\t\$novel_expressed_transcripts\t\$percent_novel_expressed_transcripts" >> novel_annotation.tsv

    awk '
      BEGIN {OFS = "\\t"}
      NR == 1 {
        print "Annotation subset", "Exact spliced transcripts",
        "Extended spliced transcripts", "Shortened spliced transcripts",
        "Other spliced transcripts", "Monoexonic transcripts";
      }
      NR >= 4 {
        if (\$1 == "string") \$1 = "All transcripts";
        if (\$1 == "string_expr") \$1 = "Expressed transcripts";
        print \$1, \$6, \$11-\$6, \$16-\$11, \$3-\$16, \$2-\$3;
      }
    ' Tables/prediction_sets_eval_wrt_ref_for_table.txt > transcripts_comparison_annotation.tsv
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
    path '*.png'

  script:
    """
    awk 'BEGIN{OFS="\\t"; print "labExpId", "Name"} NR==1{for(i=2; i<=NF; i++){print \$i, \$i}}' $formatted_genes_tpm > metadata.tsv
    plot_gene_expression.sh \\
      $formatted_genes_tpm \\
      $formatted_genes_counts \\
      metadata.tsv Name .
    mv histogram.log_T.psd_0.genes_TPM.png refgenes_log10TPM_distribution_nozero.png
    mv histogram.log_T.psd_1e-04.genes_TPM.png refgenes_log10TPM_distribution_withzero.png
    mv histogram.log_T.psd_0.genes_readcount.png refgenes_log10readcount_distribution_nozero.png
    mv histogram.log_T.psd_1e-04.genes_readcount.png refgenes_log10readcount_distribution_withzero.png
    mv TPM_fraction.genes_totalTPM_captured_by_top_genes.png cumulative_fraction_of_refgeneTPMsum_captured_by_N_most_expr_refgenes.png
    """
}

// Control exonic reads counts
// ###########################
reference_annotation_to_control_exons.combine(Channel.of('reference')).concat(
  novel_annotation_to_control_exons.combine(Channel.of('novel'))
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
    path '*' from control_exons_to_report.flatten().collect().ifEmpty([])
    path '*' from control_contigs_to_report.flatten().collect().ifEmpty([])
    path '*' from control_metrics_to_report.flatten().collect().ifEmpty([])
    path '*' from control_flags_to_report.flatten().collect().ifEmpty([])
    path '*' from feelnc_classifier_log_to_report.flatten().collect().ifEmpty([])
    path '*' from feelnc_classes_to_report.flatten().collect().ifEmpty([])
    path '*' from feelnc_rf_summary_to_report.flatten().collect().ifEmpty([])

  output:
    path 'multiqc_report.html'

  script:
    """
    for f in *.tsv *.png; do
      [ -f "\$f" ] || continue
      mv "\$f" "\${f%.*}"_mqc."\${f##*.}"
    done
    multiqc --config $config .
    """
}
