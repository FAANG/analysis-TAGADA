# ~sdjebali/Awk/reconstruct_chimtr.awk
# From an annotation gtf file with at least exons that have transcript in $12 and a tsv file where each row is the specification for a chimeric transcript
# - donor gene id
# - acceptor gene id
# - donor transcript id
# - acceptor transcript id
# - donor chromosome
# - donor strand
# - acceptor chromosome
# - acceptor strand
# - donor coord of the chimeric junction
# - acceptor coord of the chimeric junction
# outputs the same tsv file as the input but with two additional columns:
# - comma separated list of exons for the donor part of the chimeric transcript
# - comma separated list of exons for the acceptor part of the chimeric transcript

# !!!!! Important note !!!!!
# !!!!! the annotation file in fileRef HAS to be sorted according to transcript id ($12) and then gbeg ($4) and gend ($5), otherwise will not work !!!!!

# Example
#########
# cd /no_backup/rg/sdjebali/Chimeras/Benchmark/Data/Make_Chimeras_from_Annot
# awk -v fileRef=gencode.v19.annotation.long.pcg.splicedgn.gtf -f ~/Awk/reconstruct_chimtr.awk genepairs_nonoverlap_samechr_samestr_okgxorder_150_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc.tsv > genepairs_nonoverlap_samechr_samestr_okgxorder_150_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist.tsv


# Inputs
# fileRef
#########
# chrX	HAVANA	gene	99883667	99894988	.	-	.	gene_id "ENSG00000000003.10"; transcript_id "ENSG00000000003.10"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "TSPAN6"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "TSPAN6"; level 2; havana_gene "OTTHUMG00000022002.1";
# 242 (26 fields)
# 16504 (28 fields)
# 130007 (30 fields)
# 53809 (32 fields)
# 495999 (34 fields)
# 336699 (36 fields)
# 384954 (38 fields)
# 360265 (40 fields)
# 508411 (42 fields)
# 94980 (44 fields)
# 9233 (46 fields)
# 848 (48 fields)
# 116 (50 fields) 
# normal input
##############
# ENSG00000188368.5	ENSG00000183207.8	ENST00000499536.2	ENST00000598768.1	chr19	+	chr19	+	42814337	49502580
# 150 (10 fields)

# output
########
# ENSG00000188368.5	ENSG00000183207.8	ENST00000499536.2	ENST00000598768.1	chr19	+	chr19	+	42814337	49502580	chr19_42812926_42814337_+,	chr19_49502580_49502621_+,
# ENSG00000159173.14	ENSG00000021852.8	ENST00000555948.1	ENST00000535057.1	chr1	-	chr1	-	201384341	57397551	chr1_201391147_201391149_-,chr1_201390801_201390854_-,chr1_201386911_201386940_-,chr1_201386244_201386247_-,chr1_201384341_201384382_-,	chr1_57397483_57397551_-,chr1_57395026_57395231_-,
# 150 (12 fields) *** real	0m18.214s


BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	if($3=="exon")
	{
	    split($12,a,"\"");
	    nbex[a[2]]++; 
	    ex[a[2],nbex[a[2]]]=$1"_"$4"_"$5"_"$7;
	}
    }
}

{
    # Making the list of exons for the donor side of the chimeric transcript
    sdon="";
    found=0;
    # if the donor transcript is on the + then we want to go from 1st exon to the one with gend=don coord
    if($6=="+")
    {
	i=1;
	while((found==0)&&(i<=nbex[$3]))
	{
	    sdon=(sdon)(ex[$3,i])(",");
	    split(ex[$3,i],a,"_");
	    if(a[3]==$9)
	    {
		found=1;
	    }
	    i++;
	}
    }
    # if the donor transcript is on the - then we want to go from the last exon to the one with gbeg=don coord (in this order)
    else
    {
	i=nbex[$3];
	while((found==0)&&(i>=1))
	{
	    sdon=(sdon)(ex[$3,i])(",");
	    split(ex[$3,i],a,"_");
	    if(a[2]==$9)
	    {
		found=1;
	    }
	    i--;
	}
    }

    # Making the list of exons for the acceptor side of the chimeric transcript
    sacc="";
    found=0;
    # if the acceptor transcript is on the + then we want to go from the exon with gebg=acc coord to the last one 
    if($8=="+")
    {
	i=1;
	while((found==0)&&(i<=nbex[$4]))
	{
	    split(ex[$4,i],a,"_");
	    if(a[2]==$10)
	    {
		found=1;
	    }
	    i++;
	}
	if(found==1)
	{
	    for(j=(i-1); j<=nbex[$4]; j++)
	    {
		sacc=(sacc)(ex[$4,j])(",");
	    }
	}
    }
    # if the acceptor transcript is on the - then we want to go from the one with gend=acc coord to the 1st one (in this order)
    else
    {
	i=nbex[$4];
	while((found==0)&&(i>=1))
	{
	    split(ex[$4,i],a,"_");
	    if(a[3]==$10)
	    {
		found=1;
	    }
	    i--;
	}
	if(found==1)
	{
	    for(j=(i+1); j>=1; j--)
	    {
		sacc=(sacc)(ex[$4,j])(",");
	    }
	}
    }
    print $0, sdon, sacc;
}

