process STRINGTIE_assemble_transcripts {

  input:
    tuple val(id), path(bam), val(direction), path(annotation)

  output:
    path(assembly)

  shell:
    assembly = id + '.gtf'
    direction = direction == 'FR' ? '--fr' : direction == 'RF' ? '--rf' : ''
    '''
    stringtie !{bam} !{direction} -G !{annotation} -o '!{assembly}'
    '''
}

process STRINGTIE_coalesce_transcripts {

  input:
    path(assemblies)
    path(annotation)

  output:
    path('novel.gtf')

  shell:
    '''
    mkdir results

    stringtie \\
      --merge !{assemblies} \\
      -G !{annotation} \\
      -o results/novel.gtf

    # Add genes
    awk \\
      -f $(which compute_boundaries.awk) \\
      -v toadd=gene \\
      -v fldno=10 \\
      -v keys=gene_name,ref_gene_id \\
      results/novel.gtf \\
      > results/novel.genes.gtf

    cat \\
      results/novel.genes.gtf \\
      results/novel.gtf \\
      | sort -k1,1 -k4,4n -k5,5rn \\
      > results/novel.all.gtf

    # Add transcript biotypes
    awk \\
      'BEGIN {
        FS = "\\t"
      }
      NR == FNR {
        match($9, /transcript_id "([^;]*)";*/, tId)
        match($9, /!{params.transcript_biotype_field} "([^;]*)";*/, biotype)
        biotypes[tId[1]] = biotype[1]
        next
      }
      !/^#/ && $3 != "gene" {
        match($9, /transcript_id "([^;]*)";*/, tId)
        if (tId[1] in biotypes) {
          print $0 " !{params.transcript_biotype_field} \\""biotypes[tId[1]]"\\";"
          next
        }
      }
      {
        print $0
      }' \\
      !{annotation} \\
      results/novel.all.gtf \\
      > results/novel.done.gtf

    mv results/novel.done.gtf novel.gtf
    '''
}

process STRINGTIE_quantify_expression {

  input:
    tuple val(id), path(bam), val(length), val(direction), path(annotation), val(type)

  output:
    path('*.tsv')

  shell:
    direction = direction == 'FR' ? '--fr' : direction == 'RF' ? '--rf' : ''
    '''
    stringtie \\
      !{bam} \\
      !{direction} \\
      -e \\
      -B \\
      -G !{annotation} \\
      -A '!{id}'.'!{type}'_genes_TPM.tsv \\
      -o '!{id}'.'!{type}'.gtf

    mv t_data.ctab '!{id}'.'!{type}'_transcripts_TPM.tsv

    mkdir counts
    mv '!{id}'.'!{type}'.gtf counts

    prepDE.py3 \\
      -i . \\
      -p counts \\
      -g '!{id}'.'!{type}'_genes_counts.csv \\
      -t '!{id}'.'!{type}'_transcripts_counts.csv \\
      -l '!{length}'

    tr ',' $'\\t' \\
      < '!{id}'.'!{type}'_genes_counts.csv \\
      > '!{id}'.'!{type}'_genes_counts.tsv

    join \\
      <(awk 'BEGIN {FS=","; OFS="\\t"} NR > 1 {print $1,$2 | "sort"}' \\
      '!{id}'.'!{type}'_transcripts_counts.csv) \\
      <(awk 'BEGIN {OFS="\\t"} NR > 1 {print $6,$9 | "sort"}' \\
      '!{id}'.'!{type}'_transcripts_TPM.tsv) \\
      -t $'\\t' \\
      -a 1 \\
      -o '1.1 2.2 1.2' \\
      | cat <(echo "transcript"$'\\t'"gene"$'\\t'"counts") - \\
      > '!{id}'.'!{type}'_transcripts_counts.tsv
    '''
}
