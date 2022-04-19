process FASTQC_control_reads {

  publishDir = [
    path: params.output + '/control/reads',
    mode: 'copy',
    overwrite: true
  ]

  input:
    tuple val(prefix), path(fastq), val(suffix)

  output:
    path('*_fastqc.zip'), emit: reports

  shell:
    '''
    rename=!{fastq}

    if [[
      !{fastq} =~ (([\\._ ][Rr]?[12])( \\([^\\)]+\\))?)?((\\.(fastq|fq|gz))+)$
    ]]; then
      rename='!{prefix}'"${BASH_REMATCH[2]}"
      if [[ -n '!{suffix}' ]]; then
        rename="$rename ("'!{suffix}'")"
      else
        rename="$rename${BASH_REMATCH[3]}"
      fi
      rename="$rename${BASH_REMATCH[4]}"
    fi

    [[ !{fastq} != "$rename" ]] && mv !{fastq} "$rename"

    fastqc "$rename"
    '''
}
