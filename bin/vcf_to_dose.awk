# vcf_to_dose.awk

BEGIN{
    OFS="\t";
}

$1=="#CHROM"{
    split($1,a,"#");
    str=a[2]"\t"$2"\tSNP\tA1\tA2\t"$6"\t"$7"\t"$8"\t"$9"\t";
    for(i=10; i<=NF-1; i++)
    {
	str=(str)((i-10+1)" "$i)("\t");
    }
    print (str)((i-10+1)" "$i);
}

$1!~/\#/{
    str=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t";
    for(i=9; i<=NF-1; i++)
    {
	split($i,a,":");
	str=(str)(a[3])("\t");
    }
    split($i,a,":");
    print (str)(a[3]);
}
