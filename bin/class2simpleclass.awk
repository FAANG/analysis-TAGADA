# class2simpleclass.awk
# from a file of peaks output by peak_distrib_matrix.sh with 8 classes makes a file with the same number of peaks
# but with 5 classes only and associated to the corresponding gene list:
# - tss (TSS1kb)
# - tts (TTS1kb and not tss)
# - intron (intronic but not tss and not tts)
# - exon (exonic but not intronic or tss or tts)
# - intergenic

# example
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/class2simpleclass.awk
# cd /work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged
# awk -f $pgm mergedpeaks_allinfo.tsv > mergedpeaks_simpleclass.tsv
# real	0m0.379s

# input
# chr	beg	end	exon	intron	tss	tss1000	tss5000	tts	tts1000	tts5000
# 1	14949	15262	ENSSSCG00000037372,	NA	NA	NA	ENSSSCG00000037372,	NA	NA	NA
# 149334 (11 fields)

# output
# chr	beg	end	class	elt.list
# chr1	10264	10433	intergenic	NA
# 215621 (5 fields) 


BEGIN{
    OFS="\t";
}

NR==1{
    print $1, $2, $3, "class", "elt.list";
}

NR>=2{
    if($7!="NA")
    {
	class="tss";
	list=$7;
    }
    else
    {
	if($10!="NA")
	{
	    class="tts";
	    list=$10;
	}
	else
	{
	    if($5!="NA")
	    {
		class="intron";
		list=$5;
	    }
	    else
	    {
		if($4!="NA")
		{
		    class="exon";
		    list=$4;
		}
		else
		{
		    class="intergenic";
		    list="NA";
		}
	    }
	}
    }
    print $1, $2, $3, class, list;
}
