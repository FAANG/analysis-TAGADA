# atacsumstat2barplotinput2.awk

# converts a peak summary stat file into a tsv file for ggplot2 where we have access for each species
# each da status and each conservation class (not aligned, 0, 1, 2, 3 other livestock species), to the
# number of peaks in the conservation class and its % wrt the total number of peaks in the da class

# Example
# cd /work/project/fragencode/results/atacseq/multi/tissue.peaks.merged/project.peaks.to.human
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/atacsumstat2barplotinput2.awk
# time awk -f $pgm sp.peak.coord.size.align.palign.psim.nborth.over.dhs.enh.avgnormaccess.da.class.tsv > species.peakda.nb.conservclass.peaks.nb.pcent.tsv
# real	0m0.686s

# input
# species	peak.coord	peak.size	peak.align	pcent.align	pcent.sim	unique	nborth	over.dhs	over.enh	peak.avg.normaccess	peak.da	peak.class
# sus_scrofa	1:14949-15262	313	aligned	36.1022	74.3363	1	0	1	1	16.3673	non.da	4_exon
# 449018 (13 fields)

# output
# species	peak.da	nb.in.da.class	conserv.class	nb.in.class	pcent.in.class
# bos_taurus	da	1971	1_not.aligned	340	17.2501
# 41 (6 fields)

BEGIN{
    OFS="\t";
    sp[1]="bos_taurus";
    sp[2]="capra_hircus";
    sp[3]="gallus_gallus";
    sp[4]="sus_scrofa";
    corr[1]="da";
    corr[2]="non.da";
    corr2[0]="zero";
    corr2[1]="one";
    corr2[2]="two";
    corr2[3]="three";
    c[1]="1_not.aligned";
    c[2]="2_zero";
    c[3]="3_one";
    c[4]="4_two";
    c[5]="5_three";
}

NR==1{
    print $1, $12, "nb.in.da.class", "conserv.class", "nb.in.class", "pcent.in.class";
}

NR>=2{
    class=($4=="not.aligned" ? "1_"$4 : ($8+2)"_"corr2[$8]);
    nb[$1,$12,class]++;
}

END{
    for(i=1; i<=4; i++)
    {
	for(j=1; j<=2; j++)
	{
	    for(k=1; k<=5; k++)
	    {
		tot[i,j]+=nb[sp[i],corr[j],c[k]];
		# print sp[i], corr[j], c[k], i,j,tot[i,j];
	    }
	    for(k=1; k<=5; k++)
	    {
		print sp[i], corr[j], tot[i,j], c[k], nb[sp[i],corr[j],c[k]], nb[sp[i],corr[j],c[k]]*100/tot[i,j];
	    }
	}
    }
}
