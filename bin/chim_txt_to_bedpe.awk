# ~sdjebali/Awk/chim_txt_to_bedpe.awk
# improved on nov 3rd 2015 in order to be able to have the junction beg and end anywhere in the file
# it has to be specified with -v jbeg= and -v jend= but by default it will be 4 and 5 to be backward
# compatible

# takes as input a chimeric junction file in txt format (as produced by ChimPipe) 
# and an optional jbeg and jend argument for the junction beg and end (default 4 and 5)
# and returns a bedpe file with the two exonic parts of the junction for the junctions
# where the two blocks are stranded, the donor first and then the acceptor. 
# This is useful to extract the exact sequence of the chimeric junction, either for 
# blatting it against the genome to discard FP or to translate into peptide ...

# input
# chr14_23346027_+:chr14_23354472_+ 4 4 23345958 23354548 1 1 8445 GT AG LRP10, REM2, LRP10, REM2, . . 1-1:44, . .
# chr7_112462026_-:chr7_112424834_- 13 13 112462101 112424750 1 1 37192 GT AG C7orf60, TMEM168, C7orf60, TMEM168, . . 1-1:28, . .
# 1984 (19 fields)  
 
# output
# chr14	23345957	23346027	chr14	23354471	23354548	chr14_23346027_+:chr14_23354472_+	.	+	+
# chr7	112462025	112462101	chr7	112424749	112424834	chr7_112462026_-:chr7_112424834_-	.	-	-
# 1984 (10 fields)


BEGIN{OFS="\t";
    if(jbeg=="")
    {
	jbeg=4;
    }
    if(jend=="")
    {
	jend=5;
    }
} 

$1!~/unc/{
    split($1,a,":"); 
    split(a[1],a1,"_"); 
    split(a[2],a2,"_"); 
    if((a1[3]!=".")&&(a2[3]!="."))
    {
	beg1=((a1[2]<$jbeg) ? a1[2] : $jbeg);
	end1=((a1[2]<$jbeg) ? $jbeg : a1[2]);
	beg2=((a2[2]<$jend) ? a2[2] : $jend);
	end2=((a2[2]<$jend) ? $jend : a2[2]);
	print a1[1], beg1-1, end1, a2[1], beg2-1, end2, $1, ".", a1[3], a2[3];
    }
}