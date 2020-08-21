# atacsumstat2barplotinput1.awk

# converts a peak summary stat file into a tsv file for ggplot2 where we have access for each species
# each peak number of orthologs and each peak class (wrt gene annot), to the number of peaks with such
# number of orthologs, the number of peaks in the gx class and its % wrt the total number of peaks with
# this number of orthologs

# Example
# cd /work/project/fragencode/results/atacseq/multi/tissue.peaks.merged/project.peaks.to.human
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/atacsumstat2barplotinput1.awk
# time awk -f $pgm 4species.besthit.merged.max1eachspecies.sp.nbsp.avgnormaccess.pcentalign.pcentsim.class.peaksize.tsv > species.nbsp2.nbinorthclass.class.peaks.nb.pcent.tsv

# input
# species	nbsp	nbsp2	avg.norm.access	pcent.aligned	pcent.sim	class	peaksize
# gallus_gallus	1	1_one	1492.32	3.63436	76.7677	5_intergenic	193
# 251007 (8 fields)

# output
# species	nbsp2	nb.in.orth.class	class	nb.in.class	pcent.in.class
# bos_taurus	0_zero	49250	1_tss	2849	5.78477
# 81 (6 fields)  *** OK format

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
    c[1]="tss";
    c[2]="tts";
    c[3]="intron";
    c[4]="exon";
    c[5]="intergenic";
}

NR==1{
    print $1, $3, "nb.in.orth.class", $7, "nb.in.class", "pcent.in.class";
}

NR>=2{
    nb[$1,$2,$7]++;
}

END{
    for(i=1; i<=4; i++)
    {
	for(j=0; j<=3; j++)
	{
	    for(k=1; k<=5; k++)
	    {
		tot[i,j]+=nb[sp[i],j,k"_"c[k]];
	    }
	    for(k=1; k<=5; k++)
	    {
		print sp[i], corr[j], tot[i,j], k"_"c[k], nb[sp[i],j,k"_"c[k]], nb[sp[i],j,k"_"c[k]]*100/tot[i,j];
	    }
	}
    }
}
