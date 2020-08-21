
BEGIN{
    OFS="\t";
}

NR>=2{
    if($2==1)
    {
	annot[1]+=$6;
	annot[2]+=$7;
	annot[3]+=$8;
    }
    else
    {
	if($3==1)
	{
	    extens[1]+=$6;
	    extens[2]+=$7;
	    extens[3]+=$8;
	}
	else
	{
	    if($4==1)
	    {
		alternative[1]+=$6;
		alternative[2]+=$7;
		alternative[3]+=$8;
	    }
	    else
	    {
		if($5==1)
		{
		    novel[1]+=$6;
		    novel[2]+=$7;
		    novel[3]+=$8;
		}
	    }
	}
    }
}

END{
    print sp, "1.annot", "1.mRNA", annot[1];
    print sp, "1.annot", "2.lncRNA", annot[2];
    print sp, "1.annot", "3.otherRNA", annot[3];
    print sp, "2.extension", "1.mRNA", extens[1];
    print sp, "2.extension", "2.lncRNA", extens[2];
    print sp, "2.extension", "3.otherRNA", extens[3];
    print sp, "3.alternative", "1.mRNA", alternative[1];
    print sp, "3.alternative", "2.lncRNA", alternative[2];
    print sp, "3.alternative", "3.otherRNA", alternative[3];
    print sp, "4.novel", "1.mRNA", novel[1];
    print sp, "4.novel", "2.lncRNA", novel[2];
    print sp, "4.novel", "3.otherRNA", novel[3];
}
