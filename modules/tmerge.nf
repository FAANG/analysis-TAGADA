process TMERGE_merge_assemblies {

  label 'memory_16'

  input:
    path(assemblies)
    path(annotation)
  
   output:
    path('novel.gtf')

  when:
    !stringtie_merge

  script:
    """
    # Add file names to transcript ids
    # Start with the reference 

    awk '!/^#/ && \$1!="." && \$3=="exon" {
        sub(/transcript_id "/,"transcript_id \\"ref:")
        print \$0
    }' > all_exons_from_assemblies_and_ref.gtf

    for assembly in $assemblies
    do
        name=`basename \$assembly .gtf`
        awk -v name=\$name \
            '!/^#/ && \$1!="." && \$3=="exon" {
                sub(/transcript_id "/,"transcript_id \\""name":")
                print \$0
            }' \$assembly >> all_exons_from_assemblies_and_ref.gtf
    done

    # Sort the file
    sort -k1,1 -k4,4n -k5,5n all_exons_from_assemblies_and_ref.gtf > all_exons_from_assemblies_and_ref_sorted.gtf
    
    # Run tmerge
    tmerge --endFuzz $params.tmerge_endfuzz --exonOverhangTolerance $params.tmerge_overhang all_exons_from_assemblies_and_ref_sorted.gtf > tmerged.gtf

    # Add gene ids
    bedtools intersect -s -wao -a tmerged.gtf -b tmerged.gtf | \
        awk '\$NF!=0' | \
        /usr/local/src/buildLoci/buildLoci.pl - | \
        sort -k1,1 -k4,4n -k5,5n > tmerged.gene_id.gtf  

    # Add ref_gene_id and gene/transcript rows
    cat tmerged.gene_id.gtf | \\
        awk -f /usr/local/bin/make_gff_ok.awk | \\
        awk -v fileRef=$annotation -f /usr/local/bin/tmerge.buildLoci2okgff.awk | \\
        awk -f /usr/local/bin/gff2gff.awk | \\
        sort -k1,1 -k4,4n -k5,5n \\
        > novel.gtf
    """
}