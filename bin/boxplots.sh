#!/bin/bash

# boxplots.sh
# make boxplots for values (y) belonging to several categories (x) from an input tsv file with ggplot2 
# - takes as input
#   * absolute path to input tsv file
#   * header key for x values in the tsv file
#   * header key for y values in the tsv file
#   * label for y axis
#   * min y value to be plotted
#   * max y value to be plotted
#   * absolute path to output file (with extension), could be pdf, png, eps.
# Note1: needs reshape2 and ggplot2 libraries to be installed
# Note2: in order to have the boxplots in a given order one can number the different sets for which we want boxplots

# example
# cd ~/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Stringtie/WithoutAnnot/Test
# input=~/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Plots/ExonPerTranscript/annotation_nbexintr_forggplot.tsv
# output=annotation_nbexintr_forggplot.png
# time boxplots.sh $input annotation nb_ex_in_transcript "Number of exons per transcript" 20 $output
# real    0m12.362s

# input is like this
####################
# annotation    nb_ex_in_transcript
# 1_Gencv19_all   2
# 1220427 (2 fields)

if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ] || [ ! -n "$6" ] || [ ! -n "$7" ]
then
   echo "" >&2
    echo Usage: boxplots.sh input.tsv headkey_xval headkey_yval yaxis_label min_yval max_yval output_file >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- absolute path to input tsv file" >&2
    echo "- header key for x values in the tsv file" >&2
    echo "- header key for y values in the tsv file" >&2
    echo "- label for y axis" >&2
    echo "- min y value to be plotted" >&2
    echo "- max y value to be plotted" >&2
    echo "- absolute path to output file (with extension), could be pdf, png, eps." >&2
    echo "" >&2
    echo "produces as output:" >&2
    echo "- the file specified as the last argument with as many boxplots as categories for x and with the values indicated in y" >&2
    echo "" >&2
    echo "Note1: needs reshape2 and ggplot2 libraries to be installed" >&2
    echo "Note2: in order to have the boxplots in a given order one can number the different sets for which we want boxplots" >&2
    exit 1
fi

echo '
library(ggplot2)
theme_set(theme_bw(base_size = 16))
data = read.delim("'$1'", sep="\t", h=TRUE)
gp = ggplot(data) + geom_boxplot(aes(y='$3',x=factor('$2'),fill=factor('$2')), varwidth = TRUE, notch=T) + scale_fill_brewer(palette="Set1") 
gp = gp + scale_x_discrete(labels=gsub("[0-9][0-9]_","",levels(as.factor(data$'$2'))))
gp = gp + theme(axis.text.x=element_text(angle=35, hjust=1, vjust=1))
gp = gp + labs(y='\'$4\'') 
gp = gp + ylim(c('$5','$6'))   
ggsave(filename="'$7'")
' | R --vanilla

