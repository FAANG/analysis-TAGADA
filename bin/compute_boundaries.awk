#~/Awk/compute_boundaries.awk
#takes as input a gff files and a field number fldno in this file where the features for which we want 
#to compute the boundaries from the gff features are, and compute the most5' and the most3' boundary 
#of it, toadd is the name of the feature in fldnoth field

#usage
#awk -v toadd=gene -v fldno=12 -f ~/Awk/compute_boundaries.awk /projects/encode/scaling_up/chr21_22/annotations/new_annotations/hg17/HAVANA_chr22_hg17_exons_2cat_CDS.gff

#/projects/encode/scaling_up/chr21_22/annotations/new_annotations/hg17/HAVANA_chr22_hg17_exons_2cat_CDS.gff
#chr22   VEGA_Novel_CDS  exon    14636441        14636677        .       -       .       transcript_id "LA16c-3G11.6-001"; gene_id "LA16c-3G11.6"; gene_alias ""; exon_id "LA16c-3G11.6-001-exon11";


$1!~/#/{
  seen[$fldno]++;
  if(seen[$fldno]==1)
    {
      chr[$fldno]=$1;
      strand[$fldno]=$7;
      cat[$fldno]=$2;
      fstbeg[$fldno]=$4;
      lstend[$fldno]=$5
    }
  else
    {
	if($7!=strand[$fldno])
	{
	    strand[$fldno]=".";
	}
	if($2!=cat[$fldno])
	{
	    cat[$fldno]=".";
	}
	if($4<fstbeg[$fldno])
	{
	    fstbeg[$fldno]=$4;
	} 
	if($5>lstend[$fldno])
	{
	    lstend[$fldno]=$5;
	} 
    }
}

END{
  OFS="\t";
  for(k in seen)
    {
      print chr[k], cat[k], toadd, fstbeg[k], lstend[k], ".", strand[k], ".", toadd"_id "(k);
    }
}

