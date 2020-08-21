

# ~sdjebali/Awk/compare_chim_junc.awk

# This script takes as input a file of chimeric junctions in chimpipe format = donchr_donpos_donstr:accchr_accpos_accstr
# at column 1 with a header and as variable file another file of chimeric junctions (reference) also in the same format
# and with a header and outputs the input file with additional information when compared to the second file, which is,
# in case the junction has at least one reference junction with the same chr and strand as it on both the don and the acc
# side:
# - the comma separated list of reference junctions that have the same chr and strand as the predicted one on both sides, 
# - the distance between the predicted junction and the reference at the donor site for each of them (comma separated list)
# - the distance between the predicted junction and the reference at the acceptor site for each of them (comma separated list)
# - the sum of those distances for each of them (comma separated list)
# - the comma separately list of minimal sums of don and acc distances 
# - the comma separated list of reference junctions corresponding to them


# IMPORTANT NOTE:
# this script proceeds in a very naive way by comparing each junction in one file to all the junctions in the other
# so it assumes that the files are not too big, with a max of 1000 items per file


# example
#########
# cd /no_backup/rg/sdjebali/Chimeras/ChimPipe/benchmark/ChimPipe-0.8.8/Berger/K-562_lib1
# ref=~/Chimeras/Benchmark/Data/Berger/K-562_junctions.txt 
# awk '{print $1}' chimeric_junctions_K-562_lib1.txt | awk -v fileRef=$ref -f ~sdjebali/Awk/compare_chim_junc.awk > chimeric_junctions_K-562_lib1_vs_K-562_junctions.txt

# inputs
########
# juncId nbstag nbtotal maxbeg maxEnd samechr samestr dist ss1 ss2 gnlist1 gnlist2 gnname1 gnname2 bt1 bt2 PEsupport maxSim maxLgal
# chr7_6077055_-:chr7_6075896_- 2 2 6077097 6075872 1 1 1159 GT AG ENSG00000086232.8, ENSG00000157999.5, EIF2AK1, ANKRD61, protein_coding, protein_coding, 1-1:4, . .
# 283 (19 fields)
# juncid	gnpair
# chr6_31620201_-:chr6_31833172_-	BAT3-SLC44A4
# 4 (2 fields)

# output
########
# juncId	refjunc	dondist	accdist	sumdist	bestdist	bestref
# chr6_41250075_-:chr6_30692222_-	chr6_31620201_-:chr6_31833172_-,	9629874,	1140950,	10770824,	10770824,	chr6_31620201_-:chr6_31833172_-,
# 283 (7 fields)


BEGIN{
    while (getline < fileRef >0)
    {
	n++;
	if(n>=2)
	{
	    split($1,a,":");
	    split(a[1],a1,"_");
	    split(a[2],a2,"_");
	    reflist[a1[1]"_"a1[3]":"a2[1]"_"a2[3]]=(reflist[a1[1]"_"a1[3]":"a2[1]"_"a2[3]])($1)(",");
	}
    }
}


{
    if(NR==1)
    {
	refjunc="refjunc";
	dondist="dondist";
	accdist="accdist";
	sumdist="sumdist";
	bestdist="bestdist";
	bestref="bestref";
    }
    else
    {
	split($1,a,":");
	split(a[1],a1,"_");
	split(a[2],a2,"_");
	if(reflist[a1[1]"_"a1[3]":"a2[1]"_"a2[3]]=="")
	{
	    refjunc=".";
	    dondist=".";
	    accdist=".";
	    sumdist=".";
	    bestdist=".";
	    bestref="."; 
	}
	else
	{
	    refjunc=reflist[a1[1]"_"a1[3]":"a2[1]"_"a2[3]];
	    dondist="";
	    accdist="";
	    sumdist="";
	    bestdist="";
	    bestref="";
	    minsum=4000000000;   # may have to be increased for genomes with largest chromosome bigger than 4 billion nt
	    split(refjunc,a,",");
	    k=1;
	    while(a[k]!="")
	    {
		split(a[k],b,":");
		split(b[1],b1,"_");
		split(b[2],b2,"_");
		dd=abs(b1[2]-a1[2]);
		ad=abs(b2[2]-a2[2]);
		dondist=(dondist)(dd)(",");
		accdist=(accdist)(ad)(",");
		sumdist=(sumdist)(dd+ad)(",");
		if((dd+ad)<minsum)
		{
		    minsum=dd+ad;
		    bestdist=(bestdist)(minsum)(",");
		    bestref=(bestref)(a[k])(",");
		}
		k++;
	    }
	}
    }
    print $0"\t"refjunc"\t"dondist"\t"accdist"\t"sumdist"\t"bestdist"\t"bestref;
}


function abs(x)
{
    return x >= 0 ? x : -x;
}
