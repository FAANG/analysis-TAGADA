# ens2ucsc.awk

# usage
# awk -f ens2ucsc.awk features_in_ens.bed > features_in_ucsc.bed
# but also works for gff like files
# for gff or gtf files it is better to use gff2gff afterwards otherwise tabs everywhere


BEGIN{
    OFS="\t";
}

{
    if($1~/^[0-9]+$/||$1=="X"||$1=="Y"||$1=="MT"||$1=="W"||$1=="Z")
    {
	if($1!="MT")
	{
	    print "chr"$0;
	}
	else
	{
	    gsub(/MT/,"M",$1)
	    print "chr"$0;
	}
    }
    else
    {
	print;
    }
}
