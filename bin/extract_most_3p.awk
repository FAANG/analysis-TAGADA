#~/Awk/extract_most_3p.awk
#For all the features of a gff input file that have the field no fldno in common,
#extract the most 3' features, taking the strand into account
#be careful: the features must be stranded
#when two features are equally upstream then we take the longer one
# on 12/18/2013 edited so that only + and - strand features are considered


#usage
#Chr21=/projects/encode/scaling_up/chr21_22/annotations/hg18/Havana_chr21_exons_nopseudo.gff
#awk -v fldno=10 -f ~/Awk/extract_most_3p.awk $Chr21 

#$Chr21
#chr21   VEGA_Novel_CDS  exon    9884493 9884538 .       +       0       transcript_id "AF254982.1-001"; gene_id "AF254982.1"; mRNA_start_not_found "0"; mRNA_end_not_found "0"; start_codon_not_found "0"; stop_codon_not_found "0";


(($7=="+")||($7=="-")){
  seen[$fldno]++;
  if(seen[$fldno]==1)   #initialization
    {
      chr[$fldno]=$1;
      strand[$fldno]=$7;
      most3p_beg[$fldno]=$4;
      most3p_end[$fldno]=$5;
      most3p_all[$fldno]=$0;
    }

  if(((strand[$fldno]=="+")&&(after($4,$5,most3p_beg[$fldno],most3p_end[$fldno])==1))||((strand[$fldno]=="-")&&(after($4,$5,most3p_beg[$fldno],most3p_end[$fldno])==0)))
    {
      most3p_beg[$fldno]=$4;
      most3p_end[$fldno]=$5;
      most3p_all[$fldno]=$0;
    }
}

END{
  for(k in seen)
    {
      print most3p_all[k];
    }
}

#the after function takes as input two objects and returns 1 if the first one is 
#after the second one. Note that when two features have the same start then it returns 1 
#if the first one is longer than the second one
function after(beg1,end1,beg2,end2)
{
  return ((beg1>beg2)||((beg1==beg2)&&((end1-beg1)>(end2-beg2)))) ? 1 : 0;
}

