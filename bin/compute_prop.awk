#compute_prop.awk

#awk -v comp=1 -v fld=no -f ~/Awk/compute_prop.awk file 


{if(fld=="")
  {
    fld=NF;
  }

 n++;
 if($fld!=0)
   n1++;
}
END{
    if(n1=="")
    {
	n1=0;
    }
    if(comp=="")
    {
	if(n>=1)
	    print "#", n, n1, n1/n*100;
	else
	    print "#", 0, 0, "NA";
    }
    else
    {
	if(n>=1)
	    print "#", n, n-n1, (n-n1)/n*100;
	else
	    print "#", 0, 0, "NA";
    }
}
