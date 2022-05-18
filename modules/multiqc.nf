process MULTIQC_generate_report {

  input:
    path(reports)
    path(config)
    path(images, stageAs: 'custom_images.tsv')
    path(metadata, stageAs: 'metadata_mqc.tsv')

  output:
    path('multiqc_report.html')

  shell:
    '''
    for report in !{reports}; do
      [[ -f "$report" && "${report##*.}" == tsv ]] || continue
      mv "$report" "${report%.*}"_mqc."${report##*.}"
    done

    multiqc --config !{config} .
    '''
}
