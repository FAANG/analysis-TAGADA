# atacsumstat2barplotinput3.awk

# converts a peak summary stat file into a tsv file for ggplot2 where we have access for each species
# each da status and each unicity class (1_unique, 2_multi), to the number of peaks in the da class
# the number of peaks in the unicity class and its % wrt the total number of peaks in the da class

# Example
# cd /work/project/fragencode/results/atacseq/multi/tissue.peaks.merged/project.peaks.to.human
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/atacsumstat2barplotinput3.awk
# time awk -f $pgm sp.peak.coord.size.align.palign.psim.nborth.over.dhs.enh.avgnormaccess.da.class.tsv > species.peakda.nb.unicityclass.peaks.nb.pcent.tsv
# real	0m0.512s

# input
# species	peak.coord	peak.size	peak.align	pcent.align	pcent.sim	unique	nborth	over.dhs	over.enh	peak.avg.normaccess	peak.da	peak.class
# sus_scrofa	1:14949-15262	313	aligned	36.1022	74.3363	1	0	1	1	16.3673	non.da	4_exon
# 449018 (13 fields)

# output
# species	peak.da	nb.in.da.class	unicity.class	nb.in.class	pcent.in.class
# bos_taurus	da	1631	1_unique	1583	97.057
# 17 (6 fields)

BEGIN{
    OFS="\t";
    sp[1]="bos_taurus";
    sp[2]="capra_hircus";
    sp[3]="gallus_gallus";
    sp[4]="sus_scrofa";
    corr[1]="da";
    corr[2]="non.da";
    c[1]="1_unique";
    c[2]="2_multi";
}

NR==1{
    print $1, $12, "nb.in.da.class", "unicity.class", "nb.in.class", "pcent.in.class";
}

NR>=2&&$4=="aligned"{
    class=($7==1 ? "1_unique" : "2_multi");
    nb[$1,$12,class]++;
}

END{
    for(i=1; i<=4; i++)
    {
	for(j=1; j<=2; j++)
	{
	    for(k=1; k<=2; k++)
	    {
		tot[i,j]+=nb[sp[i],corr[j],c[k]];
	    }
	    for(k=1; k<=2; k++)
	    {
		print sp[i], corr[j], tot[i,j], c[k], nb[sp[i],corr[j],c[k]], nb[sp[i],corr[j],c[k]]*100/tot[i,j];
	    }
	}
    }
}
