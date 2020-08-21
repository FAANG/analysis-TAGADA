# ensgenome2ucsc.awk

# usage
# awk -f ensgenome2ucsc.awk ens_genome.fa > ucsc_genome.fa


BEGIN{
    OFS="\t";
}

{
    if($1~/^>[0-9]+$/||$1==">X"||$1==">Y"||$1==">MT"||$1==">W"||$1==">Z")
    {
	new=substr($0,2,length($0)-1)
	if($1!=">MT")
	{
	    print ">chr"new;
	}
	else
	{
	    gsub(/MT/,"M",new)
	    print ">chr"new;
	}
    }
    else
    {
	print;
    }
}
