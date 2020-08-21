#!/bin/bash

# make_hist_rna_in_ref_better_onlyplotting_gal.sh 
#################################################
# same as make_hist_rna_in_ref_better_onlyplotting.sh except that it will print
# as many bins as it has in the txt input file and even if they are negative
# will take min and max from input file
# takes as imput 2 arguments
# - 1st = name of the distance bin file from which to make the plot
# - 2nd = label for the main title for the plot. Typically cell line and replicate number instead of 
#   exp short name since too cryptic

# Usage:
########
# make_hist_rna_in_ref_better_onlyplotting_gal.sh distbinfile.txt cellline_biorep


# Programs
###########
EPSTOJPG=eps2jpg.sh

# If the user does not provide any reference segment file or any nr/distinct RNA feature file or any strandedness boolean
#################################################################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] 
then
    echo  >&2
    echo Usage:    make_hist_rna_in_ref_better_onlyplotting_gal.sh distbinfile.txt label >&2
    echo "         where:" >&2
    echo "         - distbinfile.txt is the file of distance bins from which to make the plot" >&2
    echo "         - label is the main title for the plot (choose something telling)" >&2
    exit 0
fi

tmp=`basename $1`
output=${tmp%.txt}

echo I am making the plot... >&2 
echo 'postscript(file="'$output'.eps"); data=read.table("'$1'"); m=min(data$V1); M=max(data$V1); hist(data[data$V1>=m&data$V1<=M,],breaks=seq(0,100,1),col="blue",xlab="",ylab="",main="'$2'",cex.main=4, cex.axis=3, cex.names=3); dev.off()' | R --vanilla   # ,breaks=seq(m,M,1)
$EPSTOJPG $output.eps





