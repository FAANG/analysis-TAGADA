#!/bin/bash

# geom_density_by_two_factors.sh
# Make a density plot for a numeric column of a tsv file according to two factors (for example type of element and score quantile)
# - takes as input
#   * the input tsv file with header from which the numeric column will be taken
#   * the column id in the header of the numeric column to be plotted
#   * the column id in the header for the 1st factor (for example refelt)
#   * the column id in the header for the 2nd factor (for example score.quantile)
#   * a string with <min>-<max> values for the column to be plotted
#   * a one word label for the title of the plot
# - produces as output in the same directory as the input:
#   * a png file named the same way as the input tsv file but ending in png instead of tsv with the density plot in question but faceted by the two factors
#     (number of rows given by the 2nd factor and number of columns given by the 1st factor)

# Example
#########
# cd ~/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019/sumstats
# module load system/R-3.6.1
# pgm=~/fragencode/tools/multi/Scripts/Bash/geom_density_by_two_factors.sh
# time $pgm refelt.scorequantile.fraglength.tsv frag.length refelt score.quantile 0-20000 2> geom_density_by_two_factors.err
# real	1m7.598s

# Input is like this
####################
# refelt	score.quantile	frag.length
# elt1	1	16656
# 30367151 (3 fields)

# Check all the obligatory inputs are indeed provided (should also check the input file exists and is not empty and that the ids are fine, for later)
#####################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ] || [ ! -n "$6" ]
then
   echo "" >&2
    echo "Usage: geom_density_by_two_factors.sh input.tsv colidtoplot factor1id factor2id min-max title" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- the input tsv file with header, from which the numeric column will be taken" >&2
    echo "- the id of the column we want to plot in the header of the input file" >&2
    echo "- the id of the column for the 1st factor in the header of the input file" >&2
    echo "- the id of the column for the 2nd factor in the header of the input file" >&2
    echo "- a string with <min>-<max> values for the column to be plotted" >&2
    echo "- a string for the title of the plot" >&2
    echo "" >&2
    echo "produces as output and in the same directory as the input:" >&2
    echo "- a png file named the same way as the input tsv file but ending in png, with the density plot of the column to be plotted" >&2
    echo "  facetted by the 1st and the 2nd factor (1st factor defining columns and 2nd factors defining rows)" >&2
    echo "" >&2
    echo "Note: needs reshape2 and ggplot2 libraries to be installed" >&2
    exit 1
fi


# Set the basename of the input file and the output directory as the directory of the input tsv file and go there
#################################################################################################################
input=$1
inbase=`basename ${input%.tsv}`
outdir=`dirname $input`
cd $outdir

# Set the min and max for the plot in case proper integers are given
#####################################################################
mM=$5
arg51=`echo $mM | awk '{split($1,a,"-"); print a[1]}'`
arg52=`echo $mM | awk '{split($1,a,"-"); print a[2]}'`
if [[ "$arg51" =~ ^[0-9]+$ ]] && [[ "$arg52" =~ ^[0-9]+$ ]] 
then
    m=$arg51
    M=$arg52
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
data = read.delim("'$input'",sep="\t", h=TRUE)
gp = ggplot(data, aes('$2',fill='$3')) + geom_density(alpha=.5, size=1) + facet_grid('$4' ~ '$3') + coord_cartesian(xlim = c('$m', '$M'))
gp = gp + labs(title = "'$6'", x="'$2'", y="Density") + theme(plot.title = element_text(hjust = 0.5)) + theme(plot.title = element_text(size=24))
w=7
h=10
ggsave(filename="'$inbase'.png", h=h, w=w)
' | R --vanilla

