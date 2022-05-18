process SAMTOOLS_sort_reads {

  input:
    tuple val(prefix), path(bam)

  output:
    tuple val(prefix), path(sorted)

  shell:
    sorted = prefix + '.bam'
    '''
    source=$(readlink -m !{bam})
    if [[ $(grep SO:unsorted <(samtools view -H "$source" | head -n 1)) ]]; then
      samtools sort -@ !{task.cpus} -o '!{sorted}' "$source"
    elif [[ !{bam} != '!{sorted}' ]]; then
      ln -s "$source" '!{sorted}'
    fi
    '''
}

process SAMTOOLS_control_flags {

  input:
    tuple val(prefix), path(bam)

  output:
    path('*.flagstat')

  shell:
    '''
    samtools flagstat -@ $((!{task.cpus} - 1)) !{bam} > '!{prefix}'.flagstat
    '''
}

process SAMTOOLS_control_alignments {

  input:
    tuple val(prefix), path(bam)

  output:
    path('*.stats')

  shell:
    '''
    samtools stats -@ $((!{task.cpus} - 1)) !{bam} > '!{prefix}'.stats
    '''
}

process SAMTOOLS_control_contigs {

  input:
    tuple val(prefix), path(bam)

  output:
    path('*.idxstats')

  shell:
    '''
    samtools index -@ $((!{task.cpus} - 1)) !{bam}
    samtools idxstats !{bam} > '!{prefix}'.idxstats
    '''
}

process SAMTOOLS_merge_reads {

  input:
    tuple val(id), path(bams), val(length), val(direction)

  output:
    tuple val(id), path(merged), val(length), val(direction)

  shell:
    merged = id + '.bam'
    '''
    samtools merge '!{merged}' !{bams} --threads !{task.cpus}
    '''
}
