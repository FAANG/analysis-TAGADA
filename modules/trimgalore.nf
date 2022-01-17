process TRIMGALORE_trim_adapters {

  label 'cpu_16'
  label 'memory_16'

  input:
    tuple val(prefix), path(fastqs)

  output:
    path('*_trimming_report.txt'), emit: reports
    tuple val(prefix), path('* (trimmed)*'), emit: results

  shell:
    paired = fastqs instanceof List ? '--paired' : ''
    '''
    for fastq in !{fastqs}; do
      nospaces=${fastq// /_}
      [[ "$fastq" != "$nospaces" ]] && mv "$fastq" "$nospaces"
      fastqs+=("$fastq")
      renamed+=("$nospaces")
    done

    trim_galore \\
      "${renamed[@]}" \\
      !{paired} \\
      --cores !{task.cpus} \\
      --gzip \\
      --basename '!{prefix}'

    for output in *_trimmed.*; do
      [[ -f "$output" ]] || continue
      mv "$output" '!{prefix}'" (trimmed)${output##*_trimmed}"
    done

    for n in 1 2; do
      for output in *_val_$n.*; do
        [[ -f "$output" ]] || continue
        R=_R$n
        i=$((n - 1))
        [[ "${fastqs[i]}" =~ ([\\._ ][Rr]?[12])(\\.(fastq|fq|gz))+$ ]] && \\
          R="${BASH_REMATCH[1]}"
        mv "$output" '!{prefix}'"$R (trimmed)${output##*_val_$n}"
      done
    done
    '''
}
