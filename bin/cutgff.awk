# ~/Awk/cutgff.awk
# this script cuts a gff file to the toth field (outpuf: from field no 1 to no to
# which is specified as an argument)

# awk -v to=10 -f ~/Awk/cutgff.awk in.gff > out.gff


{
    s="";
    if(to<=9)
    {
	for(i=1; i<=to-1; i++)
	{
	    s=(s)($i)("\t");
	}  
	s=(s)($i);
	print s;
    }
    else
    {
	for(i=1; i<=8; i++)
	{
	    s=(s)($i)("\t");
	}  
	for(i=9; i<=to-1; i++) 
	{
	    s=(s)($i)" ";
	}	
	s=(s)($i);
	print s;
    }
}