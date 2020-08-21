# atacsumstat2barplotinput2-2.awk

# converts a peak summary stat file into a tsv file for ggplot2 where we have access for each species
# each conservation class (0, 1, 2, 3 other livestock species) and each da class, to the total
# number of peaks in the conservation class, the nb of peaks in the conservation class and da class
# and the corresponding %

# Example
# cd /work/project/fragencode/results/atacseq/multi/tissue.peaks.merged/project.peaks.to.human
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/atacsumstat2barplotinput2-2.awk
# time awk 'NR==1||($4=="aligned"&&$7==1)' sp.peak.coord.size.align.palign.psim.nborth.over.dhs.enh.avgnormaccess.da.class.tsv | awk -f $pgm | awk 'NF==6' > species.aligned.unique.conservclass.nb.da.peaks.nb.pcent.tsv

# input
# species	peak.coord	peak.size	peak.align	pcent.align	pcent.sim	unique	nborth	over.dhs	over.enh	peak.avg.normaccess	peak.da	peak.class
# sus_scrofa	1:14949-15262	313	aligned	36.1022	74.3363	1	0	1	1	16.3673	non.da	4_exon
# 251007 (13 fields)

# output
# species	conservclass	nb.in.conserv.class	da.class	nb.in.class	pcent.in.class
# bos_taurus	0_zero	49250	da	1120	2.27411
# 33 (6 fields)

BEGIN{
    OFS="\t";
    sp[1]="bos_taurus";
    sp[2]="capra_hircus";
    sp[3]="gallus_gallus";
    sp[4]="sus_scrofa";
    corr[0]="0_zero";
    corr[1]="1_one";
    corr[2]="2_two";
    corr[3]="3_three";
    c[1]="da";
    c[2]="non.da";
}

NR==1{
    print $1, "conservclass", "nb.in.conserv.class", "da.class", "nb.in.class", "pcent.in.class";
}

NR>=2{
    nb[$1,$8,$12]++;
}

END{
    for(i=1; i<=4; i++)
    {
	for(j=0; j<=3; j++)
	{
	    for(k=1; k<=2; k++)
	    {
		tot[i,j]+=nb[sp[i],j,c[k]];
	    }
	    for(k=1; k<=2; k++)
	    {
		print sp[i], corr[j], tot[i,j], c[k], nb[sp[i],j,c[k]], nb[sp[i],j,c[k]]*100/tot[i,j];
	    }
	}
    }
}
