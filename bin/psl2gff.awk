# psl2gff.awk

# usage
# awk -f psl2gff.awk input.psl > output.gff
# makes an exon gff file out of an alignment psl file
# the gene and transcript id will be the alignment query sequence plus a counter

BEGIN{OFS="\t";
}
NF==21{
    split($21,tbeg,",");  # The begining positions of the blocks on the target sequence (0-based)
    split($19,size,",");  # The sizes of the blocks (there are $18 of them)
    split($20,qbeg, ","); # The begining positions of the blocks on the query sequence (0-based)

    nb[$10]++;
    for(k=1; k<=$18; k++)
    {
	print $14, "p2g", "exon", tbeg[k]+1, tbeg[k]+size[k], ".", $9, ".", "gene_id \""$10"."nb[$10]"\"\; transcript_id \""$10"."nb[$10]"\"\;";
    }
}
