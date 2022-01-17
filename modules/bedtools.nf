process BEDTOOLS_compute_coverage {

  label 'memory_4'

  publishDir = [
    path: params.output + '/coverage',
    mode: 'copy'
  ]

  input:
    tuple val(prefix), path(bam), val(direction)

  output:
    path(bed)

  shell:
    bed = prefix + '.bed'
    '''
    if [[ '!{direction}' == RF || '!{direction}' == FR ]]; then
      bedtools genomecov -ibam !{bam} -bg -strand + > +.tsv
      bedtools genomecov -ibam !{bam} -bg -strand - > -.tsv
      cat \\
        <(awk 'BEGIN {OFS="\\t"} {print $0,"+"}' +.tsv) \\
        <(awk 'BEGIN {OFS="\\t"} {print $0,"-"}' -.tsv) \\
        | sort -T "." -k1,3 -k5 \\
        > '!{bed}'
    else
      bedtools genomecov -ibam !{bam} -bg > '!{bed}'
    fi
    '''
}
