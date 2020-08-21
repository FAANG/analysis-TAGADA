# gather.atac.and.villar.peaks.awk
# this is to get a tsv file with header of atac-seq peaks with information about
# - their avg accessibility in liver
# - their avg accessibility in T cells
# - their overlap with villar h3k4me3 peaks
# - their overlap with villar h3k27ac peaks
# - whether they are over accessible in liver
# - whether they are over accessible in T cells

# cd /work/project/fragencode/workspace/sdjebali/fragencode/atacseq/compare_to_other_data/villar_et_all_cell_2015/sus_scrofa
# head=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed.readme.idx
# atac=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# de=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.cdvsliver/data/mergedpeaks.peaknb.allexp.readnb.normcounts.diff.cd.liver.bed
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/gather.atac.and.villar.peaks.awk
# time cat $head $atac | awk -v fileRef1=fragencode.atacpeaks.H3K4me3peakscorelist.tsv -v fileRef2=fragencode.atacpeaks.H3K27Acpeakscorelist.tsv -v fileRef3=$de -f $pgm > atacpeaks.normaccessavg.liver.tcell.overlap.h3k4me3.h3k27ac.overexpr.liv.tcell.tsv

# inputs
########
# 9	127400576	127401482	.	2	+	.
# 1	270235515	270236031	.	1	+	15.7466666667,
# 149333 (7 fields)

# 9	127400576	127401482	.	2	+	.
# 149333 (7 fields)

# AEMK02000693.1	8381	9237	.	3	+	57.74	121.51	77.25	95.75	62.05	114.16	64.83	94.62	158.97	385.57	267.37	1.9472750971068e-278	4.47191866578495	4.43018532978147
# 9145 (20 fields)  *** when $(NF-1) is positive there is overexpression in liver
# and the corresponding header was
# work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.cdvsliver/data/mergedpeaks.peaknb.allexp.readnb.normcounts.diff.readme.idx 
# chr	beg	end	id	peaknb	str	pig1.cd4	pig2.cd4	pig3.cd4	pig4.cd4	pig1.cd8	pig2.cd8	pig3.cd8	pig4.cd8	pig1.liver	pig2.liver	pig3.liver	adj.pval	logFC	normLogFC
# 1 (20 fields)


# chr	beg	end	id	peaknb	str	pig1.cd4	pig2.cd4	pig3.cd4	pig4.cd4	pig1.cd8	pig2.cd8	pig3.cd8	pig4.cd8	pig1.liver	pig2.liver	pig3.liver
# 1	14949	15262	.	1	+	9.34	16.9	20.49	10.74	16.77	18.46	19.43	14.08	19.61	13.18	21.04
# 149334 (17 fields)


# output
########
# chr	beg	end	id	peaknb	str	avg.liver	avg.tcells	h3k4me3	h3k27ac	over.in.liv	over.in.tcel
# 1	14949	15262	.	1	+	17.9433	15.7762	0	0	0	0
# 1	274421	274867	.	1	+	51.96	23.095	0	0	1.11536758943911	0
# 1	349361	349711	.	1	+	32.14	25.775	74.39	33.1333	0	0
# 149334 (12 fields)  *** some checks done for different columns but not exhausively


BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	# in case there is overlap with a H3K4me3 peak, we do the average of the intensities of this mark
	if($NF!=".")
	{
	    n=split($NF,a,",");
	    k=1;
	    while(a[k]!="")
	    {
		ok1[$1":"$2":"$3]+=a[k];
		k++;
	    }
	    ok1[$1":"$2":"$3]=ok1[$1":"$2":"$3]/(n-1);
	}
    }
    while (getline < fileRef2 >0)
    {
	# in case there is overlap with a H3K27ac peak, we do the average of the intensities of this mark
	if($NF!=".")
	{
	    n=split($NF,a,",");
	    k=1;
	    while(a[k]!="")
	    {
		ok2[$1":"$2":"$3]+=a[k];
		k++;
	    }
	    ok2[$1":"$2":"$3]=ok2[$1":"$2":"$3]/(n-1);
	}
    }
    while (getline < fileRef3 >0)
    {
	# if there is overexpr in liver we put the corresponding logFC
	if($(NF-1)>0)
	{
	    overliv[$1":"$2":"$3]=$(NF-1);
	}
	# otherwise we do the same but for tcell
	else
	{
	    overtcel[$1":"$2":"$3]=(-1*$(NF-1));
	}
    }
}

NR==1{
    for(i=7; i<=NF; i++)
    {
	split($i,a,".");
	if(a[2]=="liver")
	{
	    okliv[i]=1;
	    nliv++;
	}
	else
	{
	    ntcel++;
	}
    }
    print $1, $2, $3, $4, $5, $6, "avg.liver", "avg.tcells", "h3k4me3", "h3k27ac", "over.in.liv", "over.in.tcel";
}

NR>=2{
    liv=0;
    tcel=0;
    for(i=7; i<=NF; i++)
    {
	if(okliv[i]==1)
	{
	    liv+=$i;
	}
	else
	{
	    tcel+=$i;
	}
    }
    print $1, $2, $3, $4, $5, $6, liv/nliv, tcel/ntcel, (ok1[$1":"$2":"$3]!="" ? ok1[$1":"$2":"$3] : 0), (ok2[$1":"$2":"$3]!="" ? ok2[$1":"$2":"$3] : 0), (overliv[$1":"$2":"$3]!="" ? overliv[$1":"$2":"$3] : 0), (overtcel[$1":"$2":"$3]!="" ? overtcel[$1":"$2":"$3] : 0);
}
