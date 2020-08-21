# starts from the paste of 8 4 column tsv files with gene lists of the same peak wrt exons, introns
# tss, tss 1kb window, tss 5kb window, tts, tts 1kb window, tts 5kb window and outputs a 11 column
# tsv file with header and the following information:
# - chromosome of the atac-seq peak
# - beg of the atac-seq peak (in bed coord)
# - end of the atac-seq peak (in bed coord)
# - list of genes with at least one exon overlapping the atac-seq peak (by at least 1 bp)
# - list of genes with at least one intron encompassing the atac-seq peak
# - list of genes with their most 5' bp overlapping the atac-seq peak
# - list of genes with the 1kb window around their most 5' bp overlapping the atac-seq peak
# - list of genes with the 5kb window around their most 5' bp overlapping the atac-seq peak
# - list of genes with their most 3' bp overlapping the atac-seq peak
# - list of genes with the 1kb window around their most 3' bp overlapping the atac-seq peak
# - list of genes with the 5kb window around their most 3' bp overlapping the atac-seq peak
# The main difficulty here is to know what is really intronic because a peak can be totally included
# in an intron and yet overlapping an exon, so we want to subtract from the intronic gene list
# the genes that are in the exonic gene list
# we also want to replace the . by NA

# example
#########
# cd /work/project/fragencode/workspace/sdjebali/atacseq/fragencode/peaks/peakdistrib/pergene_andrea/test
# paste mergedpeaks_over_exons.bed mergedpeaks_over_introns.tsv mergedpeaks_over_tss.bed mergedpeaks_over_tss_ext1000.bed mergedpeaks_over_tss_ext5000.bed mergedpeaks_over_tts.bed mergedpeaks_over_tts_ext1000.bed mergedpeaks_over_tts_ext5000.bed | awk -f ../peakoverlap2classif.awk > mergedpeaks_allinfo.tsv

# input
#######
# 1	161924	162248	.	1	161924	162248	ENSSSCG00000030218,	1	161924	162248	.	1	161924	162248	.	1	161924	162248	ENSSSCG00000030218,	1161924	162248	.	1	161924	162248	.	1	161924	162248	ENSSSCG00000030218,
# ...
# X	127252267	127252660	ENSSSCG00000012699,	X	127252267	127252660	ENSSSCG00000012699,	X	127252267	127252660	ENSSSCG00000012699,	X127252267	127252660	ENSSSCG00000012699,	X	127252267	127252660	ENSSSCG00000012699,	X	127252267	127252660	.	X	127252267	127252660	.	X	127252267	127252660	.
# ... 

# output
########
# chr	beg	end	exon	intron	tss	tss1000	tss5000	tts	tts1000	tts5000
# ...
# 1	155874	156202	NA	NA	NA	NA	NA	NA	NA	ENSSSCG00000030218,    *** at tss5000 of a gene only (possible), looked at ucsc and ok
# 1	159261	159521	ENSSSCG00000030218,	NA	NA	NA	NA	NA	ENSSSCG00000030218,	ENSSSCG00000030218, *** exonic but also at the tts1000 and tts5000 but not at the tts, possible if far from the tts itself, looked at ucsc and ok
# 1	161924	162248	NA	ENSSSCG00000030218,	NA	NA	ENSSSCG00000030218,	NA	NA	ENSSSCG00000030218, *** intronic and at tss5000 and tts5000, possible, looked at ucsc and ok
# 1	164463	165057	ENSSSCG00000030218,	NA	NA	NA	ENSSSCG00000030218,	NA	NA	NA  *** exonic and tss5000, possible and looked at ucsc and ok
# 1	168067	168491	NA	NA	NA	NA	ENSSSCG00000030218,	NA	NA	NA
# 1	173226	173592	NA	NA	NA	NA	NA	NA	NA	NA
# 1	235658	235907	NA	ENSSSCG00000029697,	NA	NA	NA	NA	NA	NA
# 1	267912	269030	ENSSSCG00000029697,	NA	ENSSSCG00000029697,	ENSSSCG00000029697,	ENSSSCG00000029697,	NA	NA	NA
# 1	280058	280799	ENSSSCG00000027726,	NA	ENSSSCG00000027726,	ENSSSCG00000027726,	ENSSSCG00000027726,	NA	NA	NA
# ...
# X	127252267	127252660	ENSSSCG00000012699,	NA	ENSSSCG00000012699,	ENSSSCG00000012699,	ENSSSCG00000012699,	NA	NA	NA
# ...
# 56365 (11 fields)


BEGIN{
    OFS="\t";
    print "chr", "beg", "end", "exon", "intron", "tss", "tss1000", "tss5000", "tts", "tts1000", "tts5000";
}

{
    s="";
    split($4,a,",");
    split($8,b,",");
    if($4==".")   # if the peak is not exonic then there is nothing to remove from the intronic list
    {
	s=na($8);
    }
    else
    {
	if($8==".")  # if the peak is not intronic then neither
	{
	    s=na($8);
	}
	else     # here the peak is both exonic and intronic so we need to subtract from the intronic gene list the genes that are in the exonic list
	{
	    k=1;
	    while(b[k]!="")   # go over the intronic gene list
	    {
		found=0;      # consider that the current gene is not in the exonic list
		l=1;
		while(found==0&&a[l]!="")  # go over the exonic gene list
		{
		    if(b[k]!=a[l])   # if we find the current intronic gene in the exonic gene list we put found to 1 
		    {
			found=1;
		    }
		    l++;
		}
		k++;
		if(found==0&&b[k]!="")    # if we have not found the current intronic gene in the whole exonic gene list then we print it
		{
		    s=(s)(b[k])(",");
		}
	    }
	}
	s=na(s);
    }
    print $1, $2, $3, na($4), s, na($12), na($16), na($20), na($24), na($28), na($32);
}

function na(x){
    return ((x=="."||x=="") ? "NA" : x);
}
