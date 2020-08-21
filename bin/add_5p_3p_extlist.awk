# add_5p_3p_extlist.awk
# takes in a bed file of annotated genes (fileRef1) and a bed file of new genes (fileRef2)
# and as input a 4 column tsv file with annot gene id, whether it is spliced, the comma separated
# list of new genes connected to it and the corresponding list of booleans for splicing status
# and gives as output a 6 column tsv file with the same info as in the input plus 2 columns
# which are 2 comma separated lists of the same size as the one in column 4 and with booleans
# indicating whether each new gene extends the annot gene in 5' and in 3' respectively

# example
# resdir=/work/project/fragencode/results/rnaseq
# andir=/work/project/fragencode/data/species
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/add_5p_3p_extlist.awk
# sp=bos_taurus
# cd $resdir/$sp/assembled
# annot=$andir/$sp/`tail -1 $andir/$sp/versions.list`/$sp\_gene_sorted.bed
# awk -v fileRef1=$annot -v fileRef2=fragn.bed -f $pgm annotgnid.spliced.samestrexonicoverfraggnlist.splicestatuslist.tsv > annotgnid.spliced.samestrexonicoverfraggnlist.splicestatuslist.5pextlist.3pextlist.tsv

# $annot
# 1	19773	19899	ENSBTAG00000046619	0	-
# 1	34626	35558	ENSBTAG00000006858	0	+
# 24616 (6 fields)

# fragn.bed
# 19	14349895	14354999	XLOC_016441	0	+
# 15	27932199	27934085	XLOC_011083	0	-
# 34296 (6 fields)

# annotgnid.spliced.samestrexonicoverfraggnlist.splicestatuslist.tsv 
# ENSBTAG00000032047	1	XLOC_013599,	1,
# ENSBTAG00000003922	1	XLOC_006961,	1,
# 24616 (4 fields)

BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	str1[$4]=$6;
	gbeg1[$4]=$2;
	gend1[$4]=$3;
    }
    while (getline < fileRef2 >0)
    {
	str2[$4]=$6;
	gbeg2[$4]=$2;
	gend2[$4]=$3;
    }
}

{
    s1="";
    s2="";
    if($3!=".")
    {
	split($3,a,",");
	k=1;
	while(a[k]!="")
	{
	    if(str1[$1]=="+")
	    {
		if(gbeg2[a[k]]<gbeg1[$1])
		{
		    s1=(s1)("1,");
		}
		else
		{
		    s1=(s1)("0,");  
		}
		if(gend2[a[k]]>gend1[$1])
		{
		    s2=(s2)("1,");
		}
		else
		{
		    s2=(s2)("0,");
		}
	    }
	    else
	    {
		if(gend2[a[k]]>gend1[$1])
		{
		    s1=(s1)("1,");
		}
		else
		{
		    s1=(s1)("0,");
		}
		if(gbeg2[a[k]]<gbeg1[$1])
		{
		    s2=(s2)("1,");
		}
		else
		{
		    s2=(s2)("0,");
		}
	    }
	    k++;
	}
    }
    print $0, (s1!="" ? s1 : "."), (s2!="" ? s2 : ".");
}
