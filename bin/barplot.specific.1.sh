#!/bin/bash

# barplot.specific.1.sh
# This is to make a barplot from two columns of a tsv file, for two types of elements (for example defined by vertex.type)
# and dividing/faceting them according to a factor (for example score.quantile)
# It takes as input:
# - a tsv file with header that contains the two columns necesary for the barplot (x and y) as well as the type of element and the factor
# - the id in the header for the x value to be plotted
# - the id in the header for the y value to be plotted (this column should be numeric)
# - the id in the header for the type of element
# - the id in the header for the factor by which to facet the plot
# - a one word label for the title of the plot
# - an optional parameter saying whether the y axis has to be transformed in log10 scale
# And it produces as output in the same directory as the input file:
# - a png file named the same way as the input file but ending in png instead of tsv, with the barplot in question

# Example
#########
# cd ~/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019/sumstats
# module load system/R-3.6.1
# pgm=~/fragencode/tools/multi/Scripts/Bash/barplot.specific.1.sh
# time $pgm refelt.scorequantile.nbconn.nbtimes.tsv degree number vertex.type score.quantile ylog 2> barplot.specific.1.err

# input is like this
####################
# vertex.type	score.quantile	degree	number
# elt1	1	1	0
# 2791 (4 fields)

# output is a file called refelt.scorequantile.nbconn.nbtimes.png with the barplot in question (here 4 rows and 2 columns since 4 score.quantile and 2 vertex.type)

# Check all the obligatory inputs are indeed provided (should also check the input file exists and is not empty and that the ids are fine, for later)
#####################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ] || [ ! -n "$6" ]
then
   echo "" >&2
    echo "Usage: barplot.specific.1.sh input.tsv xcolid ycolid elttypeid factorid title <ylog10>" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a tsv file with header that contains the two columns necesary for the barplot (x and y) as well as the type of element and the factor" >&2
    echo "- the id in the header for the x value to be plotted (e.g degree)" >&2
    echo "- the id in the header for the y value to be plotted (this column should be numeric) (e.g number)" >&2
    echo "- the id in the header for the type of element" >&2
    echo "- the id in the header for the factor by which to facet the plot" >&2
    echo "- a one word label for the title of the plot" >&2
    echo "- an optional parameter saying whether the y axis has to be transformed in log10 scale" >&2
    echo "" >&2
    echo "produces as output and in the same directory as the input:" >&2
    echo "- a png file named the same way as the input file but ending in png instead of tsv, with the barplot in question" >&2
    echo "" >&2
    echo "Note: needs reshape2 and ggplot2 libraries to be installed" >&2
    exit 1
fi

# Check if the optional argument about the y log scale and in this case set str to "scale_y_log10() +"
#######################################################################################################
if [ -n "$7" ]
then
    str="+ scale_y_log10()"
fi

# Set the basename of the input file and the output directory as the directory of the input tsv file and go there
#################################################################################################################
input=$1
inbase=`basename ${input%.tsv}`
outdir=`dirname $input`
cd $outdir

# Script content
#################
echo '
library(reshape2)
library(ggplot2)
theme_set(theme_bw(base_size = 16))
data = read.delim("'$input'",sep="\t", h=TRUE)
gp = ggplot(data, aes(x='$2',y='$3',fill='$4')) + geom_bar(stat="identity") '$str' + facet_grid('$5' ~ '$4') 
gp = gp + labs(title = "'$6'", x="'$2'", y="'$3'") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=24))
w=7
h=10
ggsave(filename="'$inbase'.png", h=h, w=w)
' | R --vanilla
