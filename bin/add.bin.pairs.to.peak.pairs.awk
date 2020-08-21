# add.bin.pairs.to.peak.pairs.awk
# from a bed like file of peak pairs on a given chr (picked from last row)
# add the bin pairs corresponding to them according to a bin size bsz and a chr size sz
# provided as parameters with -v
# this is provided on the main output file but another file named after the chr, its size and the bin size
# is also provided that only had chr, beg, end, bin number for each bin
# this is to be able to use hitc for plotting heatmaps out of correlations

BEGIN{
    OFS="\t";
    m=1000000000000000;
}

NR>=2{
    chr=$1;
    mid1=($2+$3)/2;
    mid2=($6+$7)/2;
    if(mid1<m)
    {
	m=mid1;
    }
    if(mid2<m)
    {
	m=mid2;
    }

    n++;
    line[n]=$0"\t"mid1"\t"mid2;
}

END{
    for(i=1; i<=n; i++)
    {
	split(line[i],a,"\t");
	if((a[10]-m)<=sz&&(a[11]-m)<=sz)
	{
	    fst=int((a[10]-m)/bsz)+1;
	    snd=int((a[11]-m)/bsz)+1;
	    print line[i], fst, snd;
	    if(snd>M)
	    {
		M=snd;
	    }
	}
    }
    for(i=1; i<=M; i++)
    {
	print chr, (i-1)*bsz+1, i*bsz, i > chr"_"sz"_"bsz"_bins.bin.coord.tsv";
    }
}
