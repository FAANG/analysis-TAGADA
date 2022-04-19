process STAR_index_genome {

  label 'cpu_16'
  label 'memory_64'

  publishDir = [
    path: params.output,
    mode: 'copy',
    overwrite: true
  ]

  input:
    path(genome)
    path(annotation)
    path(fastqs)

  output:
    path('index')

  shell:
    '''
    for fastq in !{fastqs}; do
      echo $(
        zcat "$fastq" \\
        | head -n 40000 \\
        | awk 'NR%4 == 2 {total += length($0)} END {print int(total/(NR/4))-1}'
      ) >> overhangs
    done

    max_overhang=$(
      awk '$1 > max || NR == 1 {max = $1} END {print max}' overhangs
    )

    mkdir index

    STAR \\
      --runThreadN !{task.cpus} \\
      --runMode genomeGenerate \\
      --genomeDir index \\
      --sjdbGTFfile !{annotation} \\
      --genomeFastaFiles !{genome} \\
      --sjdbOverhang "$max_overhang"
    '''
}

process STAR_align_reads {

  label 'cpu_16'
  label 'memory_32'

  publishDir = [
    path: params.output,
    mode: 'copy',
    overwrite: true,
    saveAs: { filename ->
      if (filename.endsWith('.bam'))
        'alignment/' + filename
      else if (filename.endsWith('.splicing.tsv'))
        'control/splicing/' + filename
    }
  ]

  input:
    tuple val(prefix), path(fastqs), path(index)

  output:
    path(splicing)
    path('*.Log.final.out'), emit: reports
    tuple val(prefix), path(bam), emit: results

  shell:
    bam = prefix + '.bam'
    splicing = prefix + '.splicing.tsv'
    '''
    STAR \\
      --runThreadN !{task.cpus} \\
      --readFilesCommand zcat \\
      \\
      `# output sorted BAM file` \\
      --outSAMtype BAM SortedByCoordinate \\
      \\
      `# keep unmapped reads` \\
      --outSAMunmapped Within \\
      \\
      `# filter spurious junctions` \\
      --outFilterType BySJout \\
      \\
      `# filter non-canonical junctions` \\
      --outFilterIntronMotifs RemoveNoncanonical \\
      \\
      `# tags to specify` \\
      `# NH: number of reported alignments` \\
      `# HI: query hit index` \\
      `# AS: alignment score` \\
      `# MD: mismatching positions` \\
      `# NM: number of mismatches in each mate` \\
      `# nM: number of mismatches per (paired) alignment` \\
      `# XS: strand` \\
      --outSAMattributes NH HI AS MD NM nM XS \\
      \\
      --genomeDir !{index} \\
      --readFilesIn !{fastqs} \\
      --outFileNamePrefix '!{prefix}'.

    mv *.Aligned.sortedByCoord.out.bam '!{bam}'
    mv *.SJ.out.tab '!{splicing}'
    '''
}
