# ~sdjebali/Awk/intersectBed_pebam_with_exon_to_mateid_with_nbmappings_okorientation.awk
# takes as input the result of intersectBed on a pe rnaseq bam file and an exon gff file with gene id in field no 10
# (surrounded by quotes and ending with semicon) as well as a mate_strand parameter for the rnaseq (among the following:
# MATE2_SENSE, MATE1_SENSE, MATE_STRAND_CSHL and NONE), and produces a 2 column file with the mate id (like name/1 
# or name/2) and the number of mappings that (strandedly if data is stranded) overlap it on the genome, according to 
# the mate_strand parameter

# usage
# CUTGFF=~sdjebali/Awk/cutgff.awk
# INTER=intersectBed
# REMOVEREDUND=~sdjebali/Awk/remove_redund_better.awk
# awk '$3=="exon"' $annot | awk -v to=12 -f $CUTGFF | $INTER -abam $bamfile -b stdin -split -bed -wo | awk -v mate_strand=$mate_strand -f ~sdjebali/Awk/intersectBed_pebam_with_exon_to_mateid_with_nbmappings_okorientation.awk | awk -v fldlist=gnlist:2 -f $REMOVEREDUND | awk '{split($4,a,","); k=1; s=""; while(a[k]!=""){split(a[k],b,":"); s=(s)(b[1])(","); k++} print $1, s}' | gzip > $outdir/readid_gnlist_whoseexoverread_noredund.txt.gz

# input
# chr1	12656	12706	HWI-ST227:223:D1D4EACXX:8:2205:13723:170780/2	254	+	12656	12706	0,0,0	1	50,	0,	chr1	HAVANA	exon	12613	12721	.	+	.gene_id "ENSG00000223972.4"; transcript_id "ENST00000456328.2";	50
# for some reason intersectBed puts together . and gene_id!!! this is why I am localising the column of gene id first

{
    split($4,a,"/");        # read id, only used for MATE2_SENSE and MATE1_SENSE
    if(NR==1)               # determine the position of the gene_id value for the first row (assumes it is the same for all rows)
    {
	i=13;
	while($i!~/gene_id/)
	{
	    i++;
	}
	if($i~/gene_id/)
	{
	    j=(i+1);
	}
    }
    split($j,b,"\"");  # gene, is used for all mate_strand
    if((mate_strand=="MATE2_SENSE")||(mate_strand=="MATE1_SENSE")||(mate_strand=="MATE_STRAND_CSHL")||(mate_strand=="NONE")) # only possible values to produce anything
    {
	if($6==$19)       # same strand between the mapping and the exon
	{
	    if(((mate_strand=="MATE2_SENSE")&&(a[2]==2))||((mate_strand=="MATE1_SENSE")&&(a[2]==1))||(mate_strand=="MATE_STRAND_CSHL")||(mate_strand=="NONE"))
	    {
		gnlist[$4]=(gnlist[$4])(b[2])(",");
	    }
	}
	else             # diff strand between the mapping and the exon, in the case only MATE2_SENSE, MATE1_SENSE and NONE are considered, not MATE_STRAND_CSHL
	{
	    if(((mate_strand=="MATE2_SENSE")&&(a[2]==1))||((mate_strand=="MATE1_SENSE")&&(a[2]==2))||(mate_strand=="NONE"))
	    {
		gnlist[$4]=(gnlist[$4])(b[2])(",");
	    }
	}
    }
}



END{
    for(m in gnlist)
    {
	print m, gnlist[m];
    }
}