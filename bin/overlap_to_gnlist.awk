# ~/Awk/overlap_to_gnlist.awk
# takes as input the output of overlap on exons from predictions and exons from annot and reports
# for each gene of the prediction file the list of annot genes it connects and . if none and also the
# number of such annot genes

# cd /work/project/fragencode/workspace/sdjebali/rnaseq/fragencode/analysis/transcript_models/bos_taurus/all/all
# time awk -v fld1=10 -v fld2=16 -f ~/Awk/overlap_to_gnlist.awk cuffex_over_annotex.gff > cuffgn_annotgnlist.tsv
# real	0m5.095s

# input
# 1	Cufflinks	exon	19774	19899	.	-	.	 gene_id "XLOC_001013"; transcript_id "TCONS_00002595"; nb_ov_ex: 1 list_ex: "ENSBTAG00000046619";,
# 1	Cufflinks	exon	449750	452465	.	-	.	 gene_id "XLOC_001018"; transcript_id "TCONS_00002600"; nb_ov_ex: 0 list_ex: .
# 1092854 (16 fields)

# output
# XLOC_039637	.	0
# XLOC_029124	ENSBTAG00000019696,	1
# 44912 (3 fields)

BEGIN{
    if(fld1=="")
    {
	fld1=10
    }
    if(fld2=="")
    {
	fld2=16;
    }
}

{
    split($fld1,a,"\"");
    ok[a[2]]=1;
    split($fld2,b,",");
    k=1;
    while(b[k]!="")
    {
	split(b[k],c,"\"");
	seen[a[2],c[2]]++;
	if(c[2]!=""&&seen[a[2],c[2]]==1)
	{
	    gnlist[a[2]]=(gnlist[a[2]])(c[2])(",");
	}
	k++;
    }
}

END{
    OFS="\t";
    for(g in ok)
    {
	print g, (gnlist[g]!="" ? gnlist[g] : "."), (gnlist[g]!="" ? (split(gnlist[g],a,",")-1) : 0);
    }
}
