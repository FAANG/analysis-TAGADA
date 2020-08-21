
# ~sdjebali/Awk/chim_txt_to_gff_ssext.awk

# Takes as input a file with header that contains chimeric junctions in chimpipe's format with their ids in 1st column 
# (donchr_donpos_donstrand:accchr_accpos_accstrand) and outputs a gff file with the donor and the acceptor parts on
# different rows and a common junction id, extending each part by the ext bp (given as parameter using -v)
# More general than chim_txt_to_bedpe_ss.awk since chim_txt_to_bedpe_ss.awk is the present script with ext=0
# TODO: put a default of 0 here for ext and try this script by default and chim_txt_to_bedpe_ss.awk on the same input
# if we get the same then rename this one chim_txt_to_bedpe_ss.awk  and replace by this name in all shell scripts using it

# example
# awk -v ext=25 -f $ChimToGff $ref | awk -f $GFF2GFF > ref.gff

NR>=2{
    pos11=""; 
    pos12=""; 
    pos21=""; 
    pos22=""; 
    split($1,a,":"); 
    split(a[1],a1,"_"); 
    split(a[2],a2,"_"); 
    if(a1[3]=="+")
    {
	pos11=a1[2]-ext; 
	pos12=a1[2];
    }
    else
    {
	if(a1[3]=="-")
	{
	    pos11=a1[2]; 
	    pos12=a1[2]+ext;
	}
    } 

    if(a2[3]=="+")
    {
	pos21=a2[2]; 
	pos22=a2[2]+ext;
    }
    else
    {
	if(a2[3]=="-")
	{
	    pos21=a2[2]-ext; 
	    pos22=a2[2];
	}
    } 

    if((pos11!="")&&(pos21!=""))
    {
	print a1[1], ".", ".", pos11, pos12, ".", a1[3], ".", "junc_id", "\""($1)"\"\;"; 
	print a2[1], ".", ".", pos21, pos22, ".", a2[3], ".", "junc_id", "\""($1)"\"\;";
    }
}