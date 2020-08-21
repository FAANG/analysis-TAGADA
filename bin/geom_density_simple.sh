#!/bin/bash

# geom_density_simple.sh (following model of boxplots.sh)
# Make a simple density plot for a numeric column of a tsv file with a shade inside and a vertical dash line for the mean
# - takes as input
#   * the input tsv file with header, from which the column will be taken, column that must be numeric
#   * the id of the column we want to plot in the header of the input file
#   * a one word label for the x axis which is the meaning of the column in question
#   * a string with minimum and maximum value for the column in question (min-max)
#   * a string to put as a title of the plot
# - produces as output and in the same directory as the input:
#   * a png file (named after the one word label of the column), with the density plot of this column
# Note: needs reshape2 and ggplot2 libraries to be installed

# example
#########
# cd ~/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019/sumstats
# module load system/R-3.6.1
# pgm=~/fragencode/tools/multi/Scripts/Bash/geom_density_simple.sh
# time $pgm dist.score.small.tsv score scoresmall 1-6 2> geom_density_simple.err
# real	0m5.510s

# input is like this
####################
# distance	score
# 989202	0.634330229095194
# 1001 (2 fields)

# Check all the obligatory inputs are indeed provided (should also check the input file exists and is not empty and that the ids are fine, for later)
#####################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ]
then
   echo "" >&2
    echo "Usage: geom_density_simple.sh input.tsv colidtoplot xaxislabel min-max title" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- the input tsv file with header, from which the numeric column will be taken" >&2
    echo "- the id of the column we want to plot in the header of the input file" >&2
    echo "- a one word label for the x axis which is the meaning of the column in question" >&2
    echo "- a string with <min>-<max> values for the column to be plotted" >&2
    echo "- a string that will be the title for the plot" >&2
    echo "" >&2
    echo "produces as output and in the same directory as the input:" >&2
    echo "- a png file (named after the one word label of the column), with the density plot of this column" >&2
    echo "" >&2
    echo "Note: needs reshape2 and ggplot2 libraries to be installed" >&2
    exit 1
fi

# Set the output directory as the directory of the input tsv file and go there
##############################################################################
outdir=`dirname $1`
cd $outdir

# Set the min and max for the plot in case proper integers are given
#####################################################################
mM=$4
arg41=`echo $mM | awk '{split($1,a,"-"); print a[1]}'`
arg42=`echo $mM | awk '{split($1,a,"-"); print a[2]}'`
if [[ "$arg41" =~ ^[0-9]+$ ]] && [[ "$arg42" =~ ^[0-9]+$ ]] 
then
    m=$arg41
    M=$arg42
else
    echo "min and max values for the column to be plotted should be proper integers and be provided separated by a dash on the command line" >&2
    exit 1
fi

# Script content
#################
echo '
library(reshape2)
library(ggplot2)
theme_set(theme_bw(base_size = 16))
data = read.delim("'$1'", sep="\t", h=TRUE)
gp = ggplot(data, aes(x='$2')) + geom_density(alpha=.2, fill="#FF6666", size=1) + geom_vline(aes(xintercept=mean('$2')),color="blue", linetype="dashed", size=1) + coord_cartesian(xlim = c('$m', '$M'))
gp = gp + labs(title = "'$5'", x="'$3'", y="Density") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=24))
w=5
h=5
ggsave(filename="'$3'.png", h=h, w=w)
' | R --vanilla

