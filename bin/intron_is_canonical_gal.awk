# ~/Awk/intron_is_canonical_gal.awk
# same as intron_is_canonical.awk except that it takes as input any chimeric intron, not nec on same chr same str and ok gx order
# so the input file has to be a bedpe file like the one duplseq_gal takes as input
# with:
# - the 24 mer corresponding to the donor site in field no 11 (this can be changed)
# - the 24 mer corresponding to the acc site in field no 12 (this can be changed)
# note that independtly of the combination of strands (+/+, -/-, +/-, -/+) and since the 24mers are revcomp when on -
# the donor splice site is always in pos 13 and 14 and the acceptor splice site is always on pos 11 and 12
# This program will output the same bedpe as input with added information about whether this intron is canonical wrt Havana criteria = 
# GT-AG or GC-AG or AT-AC (even if the strand is - we have this pattern here, due to the extraction)

# cd ~/bin/OcamlSrc/DuplSeq_General/Test/All_version_0.7.0
# awk -v flddon=11 -v fldacc=12 -f ~/Awk/intron_is_canonical_gal.awk junc_stag_total_ss1_ss2_LID16627_knownstrand_with24mer_donacc.bedpe_withduplseq.bedpe > junc_stag_total_ss1_ss2_LID16627_knownstrand_with24mer_donacc.bedpe_withduplseq_can.bedpe

# input
# chr3	156260387	156260388	chr11	71806062	71806063	chr3_156260388_+:chr11_71806063_+	.	+	+	 CACTATTTTTAGGTTACTACCTTG CGCCTGGGGGAGGTGAATAAGCTG . . . .
# 37703 (16 fields)

# output


BEGIN{
    if(flddon=="")
    {
	flddon=11;
    }
    if(fldacc=="")
    {
	fldacc=12;
    }
}

{
    canonical=0;
    split($flddon,d,"");   # the piece of interest is in d[13]d[14] whatever the strand
    split($fldacc,a,"");   # the piece of interest is in a[11]a[12] whatever the strand   

    don=(toupper(a[11]))(toupper(a[12]));
    acc=(toupper(d[13]))(toupper(d[14]));
    
    if(don=="AG")
    {
	if((acc=="GT")||(acc=="GC"))
	{
	    canonical=1;
	}	
    }
    else
    {
	if(acc=="AC")
	{
	    if(don=="AT")
	    {
		canonical=1;
	    }	
	}
    }
    print $0, canonical;
}
