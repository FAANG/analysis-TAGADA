process TMERGE_coalesce_transcripts {

  label 'memory_16'

  publishDir = [
    path: params.output + '/annotation',
    mode: 'copy'
  ]

  input:
    path(assemblies)
    path(annotation)

   output:
    path('novel.gtf')

  shell:
    '''
    mkdir results

    # Add file names to transcript ids
    # Start with the reference
    awk \\
      '!/^#/ && $1 != "." && $3 == "exon" {
        sub(/transcript_id "/,"transcript_id \\"ref:")
        print $0
      }' \\
      !{annotation} \\
      > results/all_exons_from_assemblies_and_ref.gtf

    for assembly in !{assemblies}; do
      name=$(basename "$assembly" .gtf)
      awk \\
        -v name=$name \\
        '!/^#/ && $1 != "." && $3 == "exon" {
          sub(/transcript_id "/,"transcript_id \\""name":")
          print $0
        }' \\
        "$assembly" \\
        >> results/all_exons_from_assemblies_and_ref.gtf
    done

    # Sort the file
    sort \\
      -k1,1 -k4,4n -k5,5n \\
      results/all_exons_from_assemblies_and_ref.gtf \\
      > results/all_exons_from_assemblies_and_ref_sorted.gtf

    # Run tmerge
    tmerge \\
      results/all_exons_from_assemblies_and_ref_sorted.gtf \\
      --endFuzz 10000 \\
      --exonOverhangTolerance 10 \\
      '!{params.tmerge_args}' \\
      > results/tmerged.gtf

    # Add gene ids
    bedtools intersect \\
      -s \\
      -wao \\
      -a results/tmerged.gtf \\
      -b results/tmerged.gtf \\
      | awk '$NF != 0' \\
      | buildLoci.pl - \\
      | sort -k1,1 -k4,4n -k5,5n \\
      > results/tmerged.gene_id.gtf

    # Add ref_gene_id and gene/transcript rows
    cat \\
      results/tmerged.gene_id.gtf \\
      | awk \\
        -f $(which make_gff_ok.awk) \\
      | awk \\
        -v fileRef=!{annotation} \\
        -f $(which tmerge.buildLoci2okgff.awk) \\
      | awk \\
        -f $(which gff2gff.awk) \\
      | \\
      sort -k1,1 -k4,4n -k5,5n \\
      > results/novel.gtf

    # Add transcript biotypes
    awk \\
      'BEGIN {
        FS = "\\t"
      }
      NR == FNR {
        match($9, /transcript_id "([^;]*)";*/, tId)
        match($9, /transcript_biotype "([^;]*)";*/, biotype)
        biotypes[tId[1]] = biotype[1]
        next
      }
      !/^#/ && $3 != "gene" {
        match($9, /transcript_id "([^;]*)";*/, tId)
        if (tId[1] in biotypes) {
          print $0 " transcript_biotype \\""biotypes[tId[1]]"\\";"
          next
        }
      }
      {
        print $0
      }' \\
      !{annotation} \\
      results/novel.gtf \\
      > results/novel.done.gtf

    mv results/novel.done.gtf novel.gtf
    '''
}
