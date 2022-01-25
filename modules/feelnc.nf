process FEELNC_classify_transcripts {

  label 'memory_16'

  publishDir = [
    path: params.output,
    mode: 'copy',
    overwrite: true,
    saveAs: { filename ->
      if (filename == 'novel.gtf')
        'annotation/' + filename
      else if (filename == 'lncRNA_classes.txt')
        'annotation/lnc_classification/' + filename
      else if (filename.endsWith('.gtf'))
        'annotation/lnc_classification/' + filename
      else if (filename == 'feelnc_classification_summary.txt')
        'control/lnc/' + filename
      else if (filename.endsWith('.feelncclassifier.log'))
        'control/lnc/' + filename
    }
  ]

  input:
    path(genome, stageAs: 'genome.input.fa')
    path(reference_annotation, stageAs: 'reference.input.gtf')
    path(novel_annotation, stageAs: 'novel.input.gtf')

  output:
    path('novel.gtf')
    path('exons.*.gtf')
    path('*.{txt,log}'), emit: reports

  shell:
    '''
    script=$(which FEELnc_codpot.pl)
    export FEELNCPATH=${script%/*}/..

    FEELnc_filter.pl \\
      --mRNAfile !{reference_annotation} \\
      --infile !{novel_annotation} \\
      --biotype transcript_biotype=protein_coding \\
      > candidate_transcripts.gtf

    FEELnc_codpot.pl \\
      --genome !{genome} \\
      --mRNAfile !{reference_annotation} \\
      --infile candidate_transcripts.gtf \\
      --biotype transcript_biotype=protein_coding \\
      --numtx 5000,5000 \\
      --kmer 1,2,3,6,9,12 \\
      --outdir . \\
      --outname exons \\
      --mode shuffle \\
      --spethres=0.98,0.98 \\
      '!{params.feelnc_args}'

    # Update annotation with new biotypes
    cp "$(readlink -m !{novel_annotation})" updated.gtf
    for biotype in lncRNA mRNA noORF TUCp; do
      [[ -f exons."$biotype".gtf ]] || continue

      awk \\
        -v biotype=$biotype \\
        'BEGIN {
          FS = "\\t"
        }
        NR == FNR {
          match($9, /transcript_id "([^;]*)";*/, tId)
          transcripts[tId[1]] = 0
          next
        }
        {
          match($9, /transcript_id "([^;]*)";*/, tId)
          if (tId[1] in transcripts) {
            # Check if there is already a biotype in the annotation
            match($9, /biotype=([^;]*)*/, oldBiotype)
            if (oldBiotype[1]) {
              print $0
            } else {
              print $0 " feelnc_biotype \\"" biotype "\\";"
            }
          } else {
            print $0
          }
        }' \\
        exons."$biotype".gtf \\
        updated.gtf \\
        > temp.gtf

      mv temp.gtf updated.gtf
    done

    # Write FEELnc classification summary
    awk \\
      'BEGIN {
        FS = OFS = "\\t"
        feelnc_classes["lncRNA"] = 0
        feelnc_classes["noORF"] = 0
        feelnc_classes["mRNA"] = 0
        feelnc_classes["TUCp"] = 0
        feelnc_classes[""] = 0
      }
      $3 == "transcript" {
        ++nb_transcripts
        match($9, /feelnc_biotype "([^;]*)";*/, feelnc_biotype)
        ++feelnc_classes[feelnc_biotype[1]]
      }
      END {
        print "Lnc transcripts",feelnc_classes["lncRNA"]
        print "Coding transcripts from FEELnc classification",feelnc_classes["mRNA"]
        print "Transcripts with no ORF",feelnc_classes["noORF"]
        print "Transcripts of unknown coding potential",feelnc_classes["TUCp"]
      }' \\
      updated.gtf \\
      > feelnc_classification_summary.txt

    # Filter coding transcripts for lnc-messenger interactions
    grep \\
      -E '#|transcript_biotype "protein_coding"|feelnc_biotype "mRNA"' \\
      !{novel_annotation} \\
      > coding_transcripts.gtf

    FEELnc_classifier.pl \\
      --mrna coding_transcripts.gtf \\
      --lncrna exons.lncRNA.gtf \\
      > lncRNA_classes.txt

    mv updated.gtf novel.gtf
    '''
}
