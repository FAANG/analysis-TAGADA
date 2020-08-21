#!/bin/bash

# bedpe.sumstats.sh
###################
# This script takes as input:
#############################
# - the absolute path of a bedpe.gz file of interactions (either predicted from 3D data in a cell type or from 1D data across cell lines) where the 1st element has a special meaning (for example promoter)
# - a min-max thresholds for the score plot
# - a min-max threshold for the fragment length plot
# and that produces as output in a sumstats directory that it creates where the bedpe file is:
##############################################################################################
# - the same bedpe.gz file but with score quantiles
# - several (5) png plot files with their associated input tsv files and as a title the basename of the bedpe.gz file
# - several tsv or tsv.gz files that represent summary statistics tables and that are supposed to complete the messages from the plots

# TODO:
# - change the inter fragment distance from middle to middle to 3' end of most 5' to 5' end of most 3'


# Example (on genologin toulouse)
#########
# cd ~/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019
# module load system/R-3.6.1
# pgm=~/fragencode/tools/multi/Scripts/Bash/bedpe.sumstats.sh
# input=~/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019/GM/pall/GM.pall.bedpe.gz
# time $pgm $input 0-5 0-20000 2> bedpe.sumstats.err
# real	12m39.221s

# Note: for this precise example, plots and table have also been put in the following documents
# google ppt here
# https://docs.google.com/presentation/d/1kCzcPIgFrVm-i-3R091ruGRY-R2DXaqZMQok9OMET_8/edit?folder=1iIN-thsZ0yR6tyQqUwXgM6UMMKZ04zC6#slide=id.g6d5bf7fc51_1_28
# google xls here
# https://docs.google.com/spreadsheets/d/1_wRYm78legumZ9h_tuRO9NYXkTwO5rgpLKZp7Cv0EcM/edit#gid=0

# Input file
############
# chr1	915520	932176	chr1	1041678	1048651	chr1:915520-932176,chr1:1041678-1048651	1.65928294258351	.	.	GM.po	9969	5226	12	121316	4.2857	4.2857	1.9875	1.7686
# 6074049 (19 fields)
# with no header but meaning of the different fields would be the following:
# chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	type	frag1_cov	frga2_cov	freq	dist	all_capture_res	all_capture_res	exp_value_dist	dist_res
# 1 (19 fields)

# Produces all these output files
#################################
# - GM.pall.scorequantile.bedpe.gz
# - Score.png
# - Distance.png
# - Distance_by_score.quantile.png
# - refelt.scorequantile.nbconn.nbtimes.png
# - refelt.scorequantile.fraglength.png
# - dist.score.quantiles.tsv.gz
# - GM.pall.nbconn.acc.to.score.tsv
# - refelt.scorequantile.fraglength.tsv.gz
# - refelt.scorequantile.nbconn.nbtimes.tsv


# Check obligatory inputs are indeed provided otherwise exit with error (should also check the input file exists, for later)
#######################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ]
then
    echo "" >&2
    echo "Usage: bedpe.sumstats.sh abs.path.to.file.bedpe.gz scoremin-scoremax lengthmin-lengthmax" >&2
    echo "" >&2
    echo "produces as output:" >&2
    echo "- the same bedpe.gz file but with score quantiles" >&2
    echo "- several (5) png plot files with their associated input tsv files" >&2
    echo "- several tsv or tsv.gz files that represent summary statistics tables and that are supposed to complete the messages from the plots" >&2
    echo "" >&2
    echo "Note1:" >&2 
    echo "- This script writes files with fixed names and therefore should not be run twice in the same directory for different inputs" >&2
    echo "- It will also not work if some gzipped files from a previous run were made so erase the sumstats dir before running it" >&2
    echo "Note2:" >&2 
    echo "- Quality controls to check the input file is indeed in bedpe format and in particular checking that" >&2
    echo "  gend is bigger than gbeg for both all first and all second elements need to be done beforehand" >&2
    echo "Note3:" >&2 
    echo "- Plots will get as title the basename of the input bedpe.gz file" >&2
    echo "Note4:" >&2 
    echo "- For intra-chromosomal connections there is a hard-coded limit of 2 Mb for the plots" >&2
    exit 1
fi

# Variable assignment
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
input=$1
inbase=`basename ${input%.bedpe.gz}`
outbase=`dirname $input`
outdir=$outbase/sumstats

# Programs
##########
ADDQUANT=$rootDir/add_quantile.R
DENSITY1=$rootDir/geom_density_simple.sh
DENSITY2=$rootDir/geom_density_by_factor.sh
BARPLOT=$rootDir/barplot.specific.1.sh
DENSITY3=$rootDir/geom_density_by_two_factors.sh

# 0. Create the output directory and go there
#############################################
mkdir -p $outdir
cd $outdir

# 1. Add score quantile to bedpe file, keep only compulsory columns, add a header and place in sumstat dir
##########################################################################################################
zcat $input | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10}' | Rscript $ADDQUANT -i stdin -c 8 -q 4 -m quantiles -s -o stdout | awk 'BEGIN{OFS="\t"; print "chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score", "strand1", "strand2", "quantile_score", "quant_index_score", "score.quantile"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $12"_"$11}' > $inbase.scorequantile.bedpe
# chrom1	start1	end1	chrom2	start2	end2	name	score	strand1	strand2	quantile_score	quant_index_score	score.quantile
# chr1	915520	932176	chr1	1041678	1048651	chr1:915520-932176,chr1:1041678-1048651	1.65928294258351	.	.	(0.557,79.7]	4	4_(0.557,79.7]
# 6074050 (13 fields)

# 2. Score distribution
########################
$DENSITY1 $inbase.scorequantile.bedpe score Score $2 $inbase
# produces the file Score.png with the density plot for the score

# 3. Total number of connections, of intra and of inter and up and down if the first segment has a meaning (here prom)
######################################################################################################################
# and when asking for a score larger than a quartile (quartile 1 will be the complete network), make tsv file with header for xls
#################################################################################################################################
for i in $(seq 1 4); do awk -v i=$i 'BEGIN{OFS="\t"} NR>=2&&$NF>=i{ntot++; if($1==$4){nintra++; mid1=($2+$3)/2; mid2=($5+$6)/2; if(mid1<mid2){nup++}}} END{ninter=ntot-nintra; ndown=nintra-nup; print i, ntot, nintra, nintra/ntot*100, (ninter=="" ? 0 : ninter), ninter/ntot*100, nup, nup/nintra*100, ndown, ndown/nintra*100}' $inbase.scorequantile.bedpe; done | awk 'BEGIN{OFS="\t"; print "score.quantile", "tot_nb", "intra_nb", "intra_pcent", "inter_nb", "inter_pcent", "up_nb", "up_pcent", "down_nb", "down_pcent"} {print}' > $inbase.nbconn.acc.to.score.tsv
# score.quantile	tot_nb	intra_nb	intra_pcent	inter_nb	inter_pcent	up_nb	up_pcent	down_nb	down_pcent
# 1	6074049	6074049	100	0	0	3029189	49.871	3044860	50.129
# 5 (10 fields)  *** exact same file as when done step by step

# 4. Connection distance distribution (for intra chr) when taking middle to middle, globally and for each score quantile
########################################################################################################################
# a. input tsv file with header
###############################
awk 'BEGIN{OFS="\t"} NR==1{print "distance", $8, $11, $12, $13} $1==$4{mid1=($2+$3)/2; mid2=($5+$6)/2; print abs(mid1-mid2), $8, $11, $12, $13} function abs(x){return (x>=0 ? x : (-1)*x)}' $inbase.scorequantile.bedpe > dist.score.quantiles.tsv
# distance	score	quantile_score	quant_index_score	score.quantile
# 121316	1.65928294258351	(0.557,79.7]	4	4_(0.557,79.7]
# 6074050 (5 fields)  *** exact same file as when done step by step
# b. plot the distance of all interactions with ggplot
######################################################
$DENSITY1 dist.score.quantiles.tsv distance Distance 0-2000000 $inbase
# exact same plot as when done step by step
# c. plot the distance of the interactions belonging to each quantile of score
##############################################################################
$DENSITY2 dist.score.quantiles.tsv distance score.quantile Distance $inbase 0-2000000
# exact same plot as when done step by step

# 5. Degree of first and second elements (vertex degree), according to score
############################################################################
# a. tsv file with header (first 4 files and then gather)
#########################################################
for i in $(seq 1 4); do awk -v i=$i 'NR>=2&&$NF>=i' $inbase.scorequantile.bedpe | awk -v i=$i '{nb1[$1":"$2":"$3]++; nb2[$4":"$5":"$6]++} END{OFS="\t"; for(s in nb1){nb1super[nb1[s]]++; if(nb1[s]>m1){m1=nb1[s]}} for(s in nb2){nb2super[nb2[s]]++; if(nb2[s]>m2){m2=nb2[s]}} for(j=1; j<=m1; j++){print "elt1", i, j, nn(nb1super[j])} for(j=1; j<=m2; j++){print "elt2", i, j, nn(nb2super[j])}} function nn(x){return (x!="" ? x : 0)}' ; done | awk 'BEGIN{OFS="\t"; print "vertex.type", "score.quantile", "degree", "number"} {print}' > refelt.scorequantile.nbconn.nbtimes.tsv
# vertex.type	score.quantile	degree	number
# elt1	1	1	0
# 2791 (4 fields)  *** exact same file as when doing step by step
# b. plot
##########
$BARPLOT refelt.scorequantile.nbconn.nbtimes.tsv degree number vertex.type score.quantile $inbase ylog
# exact same plot as when done step by step

# 6. Fragment length distribution and distinguishing the first and the second elt, for interactions with a score higher than a quantile
#######################################################################################################################################
# a. tsv file with header
#########################
for i in $(seq 1 4); do awk -v i=$i 'BEGIN{OFS="\t"} NR>=2&&$NF>=i{print "elt1", i, $3-$2; print "elt2", i, $6-$5}' $inbase.scorequantile.bedpe; done | awk 'BEGIN{OFS="\t"; print "refelt", "score.quantile", "frag.length"} {print}' > refelt.scorequantile.fraglength.tsv
# refelt	score.quantile	frag.length
# elt1	1	16656
# 30367151 (3 fields)  *** exact same file as when doing step by step
# b. plot with ggplot2
######################
$DENSITY3 refelt.scorequantile.fraglength.tsv frag.length refelt score.quantile $3 $inbase
# real	1m7.598s  *** exact same plot as when done step by step

# c. make a tsv file to see the complete distrib (not truncated like here)
###########################################################################
for i in $(seq 1 4); do for e in elt1 elt2; do awk -v i=$i -v e=$e '$1==e&&$2==i{print $3}' refelt.scorequantile.fraglength.tsv > tmp; stats.sh tmp | awk -v i=$i -v e=$e 'BEGIN{OFS="\t"} NR==1{tot=$3}NR==3{print i, e, tot, $2, $3, $4, $5, $6, $7}' ; done; done > scorequantile.refelt.nb.fraglength.distrib.tsv 
# 1	elt1	6074049	582	7285	10764	12318	15627	122914
# 1	elt2	6074049	321	3968	5150	6340	7125	1664250
# 8 (9 fields)  *** exact same file as when doing step by step

# 7. remove temporary files and gzip big files
##############################################
rm tmp
gzip $inbase.scorequantile.bedpe
gzip dist.score.quantiles.tsv
gzip refelt.scorequantile.fraglength.tsv

