
# ~/Awk/gff2gff.awk

$1!~/#/{
    for (i=1;i<=7;i++)
    {
	printf $i"\t";
    }
    printf $8;
    if(NF>8)
    {
	printf "\t"$9;
	for (i=10;i<=NF;i++)
	{
	    printf " "$i;
	}
    }
    print "";
}

$1~/#/{print}