process FEATURECOUNTS_control_exons {

  input:
    path(bams)
    tuple path(annotation), val(type)

  output:
    path('*.tsv')
    path('*.summary'), emit: reports

  shell:
    '''
    featureCounts \\
      -t exon \\
      -g gene_id \\
      -s 0 \\
      --primary \\
      -T !{task.cpus} \\
      -a !{annotation} \\
      -o '!{type}'_exons_counts.tsv \\
      !{bams}

    mv '!{type}'_exons_counts.tsv.summary '!{type}'_exons_counts.summary
    sed -i -e '2s/\\.bam\\b//g' '!{type}'_exons_counts.tsv
    sed -i -e '1s/\\.bam\\b//g' '!{type}'_exons_counts.summary
    '''
}
