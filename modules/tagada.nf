process TAGADA_estimate_reads {

  input:
    tuple val(prefix), path(bam), path(annotation)

  output:
    tuple val(prefix), path(bam), env(length), env(direction)

  shell:
    '''
    length=$(
      samtools view !{bam} \\
      | head -n 10000 \\
      | awk '{total += length($10)} END {print int(total/NR)}'
    )

    proportions=($(infer_library_type.sh !{bam} !{annotation}))

    direction=$(
      awk \\
        -v fr=${proportions[0]} \\
        -v rf=${proportions[1]} \\
        'BEGIN {
          difference = int(sqrt((fr - rf)^2))
          fr = int(fr)
          rf = int(rf)
          if (difference > 50 && fr > rf) {
            print "FR"
          } else if (difference > 50) {
            print "RF"
          } else {
            print "Undirected"
          }
        }'
    )
    '''
}

process TAGADA_format_quantifications {

  input:
    path(quantifications)

  output:
    path('*.tsv')

  shell:
    '''
    format_count_matrices.py !{quantifications}
    '''
}

process TAGADA_filter_transcripts {

  input:
    path(assemblies)

  output:
    path('filtered_transcripts.txt')
    path('filtered/*.gtf'), emit: results

  shell:
    min_transcript_occurrence =
      params.min_transcript_occurrence ?
      '--min-transcript-occurrence ' + params.min_transcript_occurrence : ''

    min_monoexonic_occurrence =
      params.min_monoexonic_occurrence ?
      '--min-monoexonic-occurrence ' + params.min_monoexonic_occurrence : ''

    min_transcript_tpm =
      params.min_transcript_tpm ?
      '--min-transcript-tpm ' + params.min_transcript_tpm : ''

    min_monoexonic_tpm =
      params.min_monoexonic_tpm ?
      '--min-monoexonic-tpm ' + params.min_monoexonic_tpm : ''

    '''
    filter_rare_transcripts.py \\
      !{assemblies} \\
      -o filtered \\
      !{min_transcript_occurrence} \\
      !{min_monoexonic_occurrence} \\
      !{min_transcript_tpm} \\
      !{min_monoexonic_tpm} \\
      > filtered_transcripts.txt
    '''
}

process TAGADA_cluster_expression {

  input:
    path(genes_tpm_quantification)
    val(quantification_metadata)

  output:
    path('*_correlation.{pdf,tsv}')

  shell:
    rows = quantification_metadata.unique({ row -> row['id'] })

    columns = rows.first()['columns'].keySet()
    columns = columns.collectEntries({ column -> [(column): [] as Set] })

    metadata = 'labExpId' + '\t' + 'Name'

    if (columns.size() > 0) metadata += '\t' + columns.keySet().join('\t')
    else columns = ['Name': rows.collect({ row -> row['id'] })]

    rows.each({ row ->
      metadata += '\n' + ([row['id']] * 2 + row['columns'].values()).join('\t')
      row['columns'].each({ column, value -> columns[column].add(value) })
    })

    factors = columns.keySet().join(',')
    sizes = columns.values().collect({ values -> values.size() })
    palettes = sizes.collect({ size -> 'rainbow.' + size + '.txt' }).join(',')
    sizes = sizes.join(' ')

    '''
    echo -e '!{metadata}' > metadata.tsv

    for size in !{sizes}; do make.rainbow.palette.sh $size .; done

    cat \\
      <(head -n 1 !{genes_tpm_quantification} | cut -f 2-) \\
      <(tail -n +2 !{genes_tpm_quantification}) \\
      > reference_genes_TPM.reformatted.tsv

    matrix_to_dist.R \\
      -i reference_genes_TPM.reformatted.tsv \\
      --log10 \\
      -c pearson \\
      -p 0.1 \\
      -v \\
      -o reference_genes_TPM_log_pearson_correlation.tsv
      
    args_dendro=""
    if [ "$(wc -l reference_genes_TPM_log_pearson_correlation.tsv)" -gt 2 ]; then
        args_dendro=" --col_dendro --row_dendro"
    fi

    ggheatmap.R \\
      -i reference_genes_TPM_log_pearson_correlation.tsv \\
      -d i \\
      --row_metadata metadata.tsv \\
      --col_metadata metadata.tsv \\
      $args_dendro \\
      --rowSide_by '!{factors}' \\
      --matrix_legend_title \\
        'reference genes TPM pearson correlation\n(log10 pseudocount 0.1)' \\
      -B 10 \\
      --matrix_palette <(echo -e '#00A600FF\n#ECB176FF\n#F2F2F2FF') \\
      --rowSide_palette '!{palettes}' \\
      --col_labels '!{factors}' \\
      --row_labels '!{factors}' \\
      -v \\
      -o reference_genes_TPM_log_pearson_correlation.pdf
    '''
}

process TAGADA_control_expression {

  input:
    path(genes_tpm_quantification)
    path(genes_counts_quantification)

  output:
    path('*.png'), optional: true, emit: reports

  shell:
    '''
    samples_count=$(awk 'NR == 1 {print NF - 1}' !{genes_tpm_quantification})
    (( $samples_count > 40 )) && exit 0

    awk \\
      'BEGIN {
        OFS = "\\t"
        print "labExpId", "Name"
      }
      NR == 1 {
        for (i = 2; i <= NF; i++) {
          print $i, $i
        }
      }' \\
      !{genes_tpm_quantification} \\
      > metadata.tsv

    plot_gene_expression.sh \\
      !{genes_tpm_quantification} \\
      !{genes_counts_quantification} \\
      metadata.tsv Name .

    mv \\
      histogram.log_T.psd_0.genes_TPM.png \\
      refgenes_log10TPM_distribution_nozero.png

    mv \\
      histogram.log_T.psd_1e-04.genes_TPM.png \\
      refgenes_log10TPM_distribution_withzero.png

    mv \\
      histogram.log_T.psd_0.genes_readcount.png \\
      refgenes_log10readcount_distribution_nozero.png

    mv \\
      histogram.log_T.psd_1e-04.genes_readcount.png \\
      refgenes_log10readcount_distribution_withzero.png

    mv \\
      TPM_fraction.genes_totalTPM_captured_by_top_genes.png \\
      cumulative_fraction_of_refgeneTPMsum_captured_by_N_most_expr_refgenes.png
    '''
}

process TAGADA_control_annotation {

  input:
    path(reference_annotation)
    path(novel_annotation)
    path(transcripts_tpm_quantification)
    path(genes_tpm_quantification)

  output:
    path('Plots')
    path('Tables')
    path('*_expressed_*.txt')
    path('{*.png,*_annotation.tsv}'), emit: reports

  shell:
    '''
    detected_elements_sumstats.sh \\
      !{reference_annotation} \\
      !{novel_annotation} \\
      !{transcripts_tpm_quantification} \\
      !{genes_tpm_quantification}

    for f in \\
      Plots/ExonPerTranscript/*.png \\
      Plots/ExonLength/*.png \\
      Plots/DistinctExonLength/*.png \\
      Plots/5pExonLength_Tr/*.png \\
      Plots/InternalExonLength/*.png \\
      Plots/5pExonLength_Gn/*.png \\
      Plots/DistinctInternalExonLength/*.png \\
      Plots/Exact_tr_dist_to_Genc_TSS/*.png \\
      Plots/TrLength/*.png
    do
      [[ -f "$f" ]] || continue
      convert "$f" -crop 1400x2100+0+0 +repage "${f%.png}"_cropped.png
    done

    convert \\
      Plots/ExonPerTranscript/*_cropped.png \\
      Plots/TranscriptPerGene/prediction_nbtringn_forggplot.png \\
      +append \\
      +repage \\
      exon_per_transcript_and_transcript_per_gene.png

    convert \\
      Plots/ExonLength/*_cropped.png \\
      Plots/DistinctExonLength/*_cropped.png \\
      Plots/MonoExTrExLength/*.png \\
      +append \\
      +repage \\
      all_distinct_exon_length_and_monoexonic_transcript_length.png

    convert \\
      Plots/5pExonLength_Tr/*_cropped.png \\
      Plots/InternalExonLength/*_cropped.png \\
      Plots/3pExonLength_Tr/*.png \\
      +append \\
      +repage \\
      5p_internal_3p_exon_length_per_transcript.png

    convert \\
      Plots/5pExonLength_Gn/*_cropped.png \\
      Plots/DistinctInternalExonLength/*_cropped.png \\
      Plots/3pExonLength_Gn/*.png \\
      +append \\
      +repage \\
      5p_internal_3p_exon_length_per_gene.png

    convert \\
      Plots/Exact_tr_dist_to_Genc_TSS/*_cropped.png \\
      Plots/TrLength/*_cropped.png \\
      Plots/cDNALength/*.png \\
      +append \\
      +repage \\
      transcript_cdna_length_and_TSStorefgeneTSS_distance_for_exact_transcripts.png

    # Novel annotation stats

    novel_genes=$(
      wc -l string/novel_gnid_nbtr.txt \\
      | awk '{print $1}'
    )

    novel_transcripts=$(
      wc -l string/novel_trid_nbex.txt \\
      | awk '{print $1}'
    )

    perl -pe 's/^"([^"]+)".+$/$1/g' \\
      string_expr/stringtie.annot.tpm0.1.2samples.exons_complete_gnid_nbtr.txt \\
      > novel_expressed_genes.txt

    novel_expressed_genes=$(
      wc -l novel_expressed_genes.txt \\
      | awk '{print $1}'
    )

    perl -pe 's/^"([^"]+)".+$/$1/g' \\
      string_expr/stringtie.annot.tpm0.1.2samples.exons_complete_trid_nbex.txt \\
      > novel_expressed_transcripts.txt

    novel_expressed_transcripts=$(
      wc -l novel_expressed_transcripts.txt \\
      | awk '{print $1}'
    )

    percent_novel_expressed_genes=$(
      echo | awk \\
        -v expressed=$novel_expressed_genes \\
        -v all=$novel_genes \\
        '{print 100 * expressed / all}'
    )

    percent_novel_expressed_transcripts=$(
      echo | awk \\
        -v expressed=$novel_expressed_transcripts \\
        -v all=$novel_transcripts \\
        '{print 100 * expressed / all}'
    )

    awk \\
      'BEGIN {
        OFS = "\\t"
      }
      NR == 1 {
        print "Annotation subset", "Exact spliced transcripts",
        "Extended spliced transcripts", "Shortened spliced transcripts",
        "Other spliced transcripts", "Monoexonic transcripts";
      }
      NR >= 4 {
        if ($1 == "string") $1 = "All transcripts";
        if ($1 == "string_expr") $1 = "Expressed transcripts";
        print $1, $6, $11-$6, $16-$11, $3-$16, $2-$3;
      }' \\
      Tables/prediction_sets_eval_wrt_ref_for_table.txt \\
      > transcripts_comparison_annotation.tsv

    # Reference annotation stats

    reference_genes=$(
      awk 'NR == 3 {print $2}' detected_transcripts_genes_numbers.tsv
    )

    reference_transcripts=$(
      awk 'NR == 2 {print $2}' detected_transcripts_genes_numbers.tsv
    )

    perl -pe 's/^"([^"]+)".+$/$1/g' \\
      ref_expr/ref.annot.tpm0.1.2samples.exons_complete_gnid_nbtr.txt \\
      > reference_expressed_genes.txt

    reference_expressed_genes=$(
      wc -l reference_expressed_genes.txt \\
      | awk '{print $1}'
    )

    perl -pe 's/^"([^"]+)".+$/$1/g' \\
      ref_expr/ref.annot.tpm0.1.2samples.exons_complete_trid_nbex.txt \\
      > reference_expressed_transcripts.txt

    reference_expressed_transcripts=$(
      wc -l reference_expressed_transcripts.txt \\
      | awk '{print $1}'
    )

    percent_reference_expressed_genes=$(
      echo | awk \\
        -v expressed=$reference_expressed_genes \\
        -v all=$reference_genes \\
        '{print 100 * expressed / all}'
    )

    percent_reference_expressed_transcripts=$(
      echo | awk \\
        -v expressed=$reference_expressed_transcripts \\
        -v all=$reference_transcripts \\
        '{print 100 * expressed / all}'
    )

    echo \\
      -e "Category\\tTotal\\tPercentage" > novel_annotation.tsv
    echo \\
      -e "Genes\\t$novel_genes\\t" >> novel_annotation.tsv
    echo \\
      -e "Expressed genes\\t$novel_expressed_genes\\t$percent_novel_expressed_genes" \\
      >> novel_annotation.tsv
    echo \\
      -e "Transcripts\\t$novel_transcripts\\t" \\
      >> novel_annotation.tsv
    echo \\
      -e "Expressed transcripts\\t$novel_expressed_transcripts\\t$percent_novel_expressed_transcripts" \\
      >> novel_annotation.tsv

    echo \\
      -e "Category\\tTotal\\tPercentage" \\
      > reference_annotation.tsv
    echo \\
      -e "Genes\\t$reference_genes\\t" \\
      >> reference_annotation.tsv
    echo \\
      -e "Expressed genes\\t$reference_expressed_genes\\t$percent_reference_expressed_genes" \\
      >> reference_annotation.tsv
    echo \\
      -e "Transcripts\\t$reference_transcripts\\t" \\
      >> reference_annotation.tsv
    echo \\
      -e "Expressed transcripts\\t$reference_expressed_transcripts\\t$percent_reference_expressed_transcripts" >> \\
      reference_annotation.tsv
    '''
}
