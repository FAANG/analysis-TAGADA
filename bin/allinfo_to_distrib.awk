# allinfo_to_distrib.awk
# takes as input a segment matrix with indication of 9 genomic domains (not nec a partition of the whole)
# which are the ones given by the script peak_distrib_matrix.sh plus the ig one which are peaks not in any
# of the 8 ctegories and provides as output a tsv file to make a barplot with two pieces of information
# that could be species and kind of peaks (from the 3 we had for atac peaks), or species
# and all, da or not da for atac peaks so that we can make a barplot split by 1st information (~ in ggplot)
# and then with as many gx distributions as we have the 2nd information
# note: 9th category (intergenic) added on May 2nd 2018

# example
# chr	beg	end	exon	intron	tss	tss1000	tss5000	tts	tts1000	tts5000
# 1	5009	6313	NA	NA	NA	NA	NA	NA	NA	NA
# 58385 (11 fields)
# meta=/work/project/fragencode/data/metadata/atacseq/atacseq_all_narrowpeak_metadata.tsv
# meta2=/work/project/fragencode/workspace/sdjebali/fragencode/atacseq/peaks/pooling_all/strategy_fromcol.tsv
# dir=/work/project/fragencode/results/atacseq
# cd $dir/allpeaks/compare_3_peak_calling_methods
# for each species
# awk 'NR>=2{print $2}' $meta | sort | uniq | while read sp
# do
    # for each peak calling method
#    awk '{print NR, $1}' $meta2 | while read no strat
#    do
#	awk -v info1=$sp -v info2=$no"_"$strat -f ~/save/Awk/allinfo_to_distrib.awk $dir/$sp/$strat/mergedpeaks_allinfo.tsv
#    done
# done | awk 'BEGIN{OFS="\t"; print "species", "peakcallmeth", "peaks_in_class_nb", "peaks_in_class_pcent", "class"} {print}' > atacseqpeaks.sp.peakcallmeth.inclass.nb.pcent.class.tsv2

# Make a barplot with the different categories
##############################################
# R
# library(reshape2)
# library(ggplot2)
# set variables
###############
# theme_set(theme_bw(base_size = 16))
# read the input file
#####################
# data = read.delim("/work/project/fragencode/results/atacseq/allpeaks/compare_3_peak_calling_methods/atacseqpeaks.sp.peakcallmeth.inclass.nb.pcent.class.tsv", sep="\t", h=TRUE)
# Make the plot
###############
# plot = ggplot(data=data, aes(x=peakcallmeth, y=peaks_in_class_pcent, fill=class)) + scale_fill_manual(values = c("#FF0000FF", "#FFBF00FF", "#80FF00FF", "#00FF40FF", "#00FFFFFF", "#0040FFFF", "#8000FFFF", "#FF00BFFF")) + geom_bar(position = 'dodge', stat="identity") + facet_wrap(~species)
# plot = plot + scale_x_discrete(labels=gsub("[0-9]_","",levels(as.factor(data$peakcallmeth))))
# plot = plot + theme(axis.text.x=element_text(angle=35, hjust=1, vjust=1))
# plot = plot + labs(title = '', x='Peak calling method', y='% ATAC-seq peaks in class')
# ggsave(filename="atacseqpeaks.sp.peakcallmeth.inclass.nb.pcent.class.pdf", path="/work/project/fragencode/results/atacseq/allpeaks/compare_3_peak_calling_methods")
# ggsave(filename="atacseqpeaks.sp.peakcallmeth.inclass.nb.pcent.class.png", path="/work/project/fragencode/results/atacseq/allpeaks/compare_3_peak_calling_methods")
# dev.off()

BEGIN{
    OFS="\t";
}

NR>=2{
    n++;
    nb0=0;
    if($4!="NA")
    {
	nex++;
    }
    else
    {
	nb0++;
    }

    if($5!="NA")
    {
	nintr++;
    }
    else
    {
	nb0++;
    }
    
    if($6!="NA")
    {
	ntss++;
    }
    else
    {
	nb0++;
    }

    if($7!="NA")
    {
	ntss1++;
    }
    else
    {
	nb0++;
    }

    if($8!="NA")
    {
	ntss5++;
    }
    else
    {
	nb0++;
    }

    if($9!="NA")
    {
	ntts++;
    }
    else
    {
	nb0++;
    }

    if($10!="NA")
    {
	ntts1++;
    }
    else
    {
	nb0++;
    }

    if($11!="NA")
    {
	ntts5++;
    }
    else
    {
	nb0++;
    }
    if(nb0==8)
    {
	nig++;
    }
}

END{
    print info1, info2, nex, nex/n*100, "1_exonic";
    print info1, info2, nintr, nintr/n*100, "2_intronic";
    print info1, info2, ntss, ntss/n*100, "3_tss";
    print info1, info2, ntss1, ntss1/n*100, "4_tss1kb";
    print info1, info2, ntss5, ntss5/n*100, "5_tss5kb";
    print info1, info2, ntts, ntts/n*100, "6_tts";
    print info1, info2, ntts1, ntts1/n*100, "7_tts1kb";
    print info1, info2, ntts5, ntts5/n*100, "8_tts5kb";
    print info1, info2, nig, nig/n*100, "9_intergenic";
}
