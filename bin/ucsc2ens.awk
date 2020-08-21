# ~/Awk/ucsc2ens.awk
# works for bed, gff and psl files (and possibly others provided there is no space in a field)
# on dec 20th 2016 I add the fld param so that the chr name can be anywhere in the file
# note that for psl and gff files a retabulation might be needed afterwards

# example
# cd ~/fragencode/workspace/sdjebali/comparativegx/projections/genes_same_species/pig/susCr4.susCr3
# awk -v fld=14 -f ~/Awk/ucsc2ens.awk pig10.2_ensembl_annot_from_ucsc.psl > pig10.2_ensembl_annot_from_ucsc_enschrnames.psl

# input
# 1430	0	0	0	0	0	8	364401	-	ENSSSCT00000005917.2	1430	0	1430	chr1	315321322	268374985	268740816	9	104,61,237,201,66,102,171,138,350,	0,104,165,402,603,669,771,942,1080,	268374985,268389243,268432047,268445275,268451820,268547181,268585500,268616531,268740466,
# 30585 (21 fields)
# output
# 1430	0	0	0	0	0	8	364401	-	ENSSSCT00000005917.2	1430	0	1430	1	315321322	268374985	268740816	9	104,61,237,201,66,102,171,138,350,	0,104,165,402,603,669,771,942,1080,	268374985,268389243,268432047,268445275,268451820,268547181,268585500,268616531,268740466,
# 30585 (21 fields)

BEGIN{
    OFS="\t";
    if(fld=="")
    {
	fld=1;
    }
}

{
    if($fld~/^chr/)
    {
	if($fld=="chrM")
	{
	    $fld="MT";
	}
	else
	{
	    split($fld,a,"chr");
	    $fld=a[2];
	}
    }
    print;
}
