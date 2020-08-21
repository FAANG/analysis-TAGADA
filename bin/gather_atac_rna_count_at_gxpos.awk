# !!! need to make sure that fileRef5 is indeed a file that has the common labels between fileRef1 and fileRef3 !!!
# !!! in particular if we have a dot as delimitor it should be in the 3 files !!!
# fld indicate the kind of gx position we want the peak to be in
# this position could be tss, tss1kb, tss5kb or others
# by default it will be tss1kb

# example
# cd /work/project/fragencode/workspace/sdjebali/fragencode/multi/atacseq.rnaseq
# mkdir -p Aux
# corr=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/mergedpeaks_allinfo.tsv
# metatac=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed.readme.idx
# atac=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# metarna=/work/project/fragencode/results/rnaseq/sus_scrofa/diffcounts.nominsum/data/refgenes.counts.min2tpm0.1.normcounts.bed.readme.idx
# rna=/work/project/fragencode/results/rnaseq/sus_scrofa/diffcounts.nominsum/data/refgenes.counts.min2tpm0.1.normcounts.bed
# time awk -v fld=7 -v dir=Aux -v fileRef1=$metatac -v fileRef2=$atac -v fileRef3=$metarna -v fileRef4=$rna -v fileRef5=an.tiss.common.atac.rna.tsv -f gather_atac_rna_count_at_gxpos.awk $corr
# real	0m11.646s

# fileRef1 = header of the atacseq peak count file, experiments must start from the 7th column
##############################################################################################
# chr	beg	end	id	peaknb	str	pig1.cd4	pig2.cd4	pig3.cd4	pig4.cd4	pig1.cd8	pig2.cd8	pig3.cd8	pig4.cd8	pig1.liver	pig2.liver	pig3.liver
# 1 (17 fields) 
# fileRef2 = atacseq peak count file, 6 first fields must be proper bed
######################################################################
# 1	697	1143	.	7	+	847.05	-222.5	355.78	124.33	102.79	158.41	219.62	90.97	-233.59	-134.81	-34.1
# 120914 (17 fields)
# fileRef3 = header of the gene count file, experiments must start from the 7th column
######################################################################################
# chr	beg	end	id	score	str	pig1.cd8	pig1.liver	pig2.cd4	pig2.cd8	pig2.liver	pig3.cd4	pig3.cd8	pig3.liver	pig4.cd4	pig4.cd8	pig4.liver
# 1 (17 fields)
# fileRef4 = gene count file, must have the gene id in 4th column
# 1	228293	268388	ENSSSCG00000029697	11	-	21.12	25.55	29.1	27.64	23.87	33.64	26.67	26.58	29.17	20.32	30.46
# 15928 (17 fields)
# fileRef5 = labels of the experiments in common between atacseq and rnaseq (in case there are more experiments in the two files)
#################################################################################################################################
# pig1.cd8
# 10 (1 fields)
# input file (corr)
###################
# chr	beg	end	exon	intron	tss	tss1000	tss5000	tts	tts1000	tts5000
# 1	697	1143	NA	NA	NA	NA	NA	NA	NA	NA
# ...
# 1	267883	269193	ENSSSCG00000029697,	NA	ENSSSCG00000029697,	ENSSSCG00000029697,	ENSSSCG00000029697,	NA	NA	NA
# ...
# 120915 (11 fields)  *** 7th column indicate the gene list for which the tss+-1kb overlaps the peaks, usually one gene in the list but could be two or more

BEGIN{
    if(fld=="")
    {
	fld=7;
    }
    OFS="\t";
    # reads the header corresponding to the atacseq peak count file
    while (getline < fileRef1 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    atacexp[i]=$i;
	}
	# check the dir variable
	if(dir=="")
	{
	    dir="Aux";
	}
    }
    
    # reads the atacseq peak count file (norm or not) and remember for each peak its count for each experiment in a double key hashtable
    while (getline < fileRef2 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    atacount[$1":"$2":"$3,atacexp[i]]=$i;
	}
    }

    # reads the header corresponding to the rnaseq gene count file
    while (getline < fileRef3 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    rnaexp[i]=$i;
	}
    }

    # reads the rnaseq gene count file (norm or not)
    while (getline < fileRef4 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    rnacount[$4,rnaexp[i]]=$i;
	}
    }

    # reads the common experiment ids (for example $animal.$tissue)
    while (getline < fileRef5 >0)
    {
	l++;
	expt[l]=$1;
    }
}

# reads the correspondance file between the atac-seq peaks and the gene ids (contains more info than just
# the genes at the gx position in question and the user can choose to look at any gx pos, not just tss or tts
# but also exonic or intronic peaks by changing fld)
# for each atac-seq peak that overlaps the gx position in question and that corresponds to a list of genes, writes a file with for each experiment
# common to atac and rna, the counts (norm or not) of atac for the peak and of rnaseq for the gene
# on which we then want to compute a correlation
NR>=2&&$fld!="NA"{
    split($fld,a,",");
    k=1;
    while(a[k]!="")
    {
	if((atacount[$1":"$2":"$3,expt[1]]!="")&&(rnacount[a[k],expt[1]]!=""))  # writes something only if both the peak and the gene were present in the matrices
	{
	    for(m=1; m<=l; m++)
	    {
		print expt[m], atacount[$1":"$2":"$3,expt[m]], rnacount[a[k],expt[m]] > dir"/"$1":"$2":"$3"."a[k]".tsv";
	    }
	}
	k++;
    }
}
