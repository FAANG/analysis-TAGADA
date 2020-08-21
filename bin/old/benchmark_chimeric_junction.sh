
# benchmark_chimeric_junction.sh
# Script that provides several numbers to assess sensitivitu (sn) and precision (pr) of a program that 
# produces chimeric junctions from rnaseq of a sample with respect to a list of reference junctions known 
# to be present in the same sample

# AN improvement would be to detect the reference and non ref cases with NA in one of their coordinates
# and also not to run the program if either the ref or the file to assess does not contain anything
# Also from BR's comments I should add
# - reporting of nb cases in wrong order
# - reporting of nb cases in wrong strand (on both sides? on one?)
# - reporting of TP in different flavours such as exact, exact+25, exact+25+100nt as was done in Berger and Edgren for 0.8.8

# Takes as input:
#################
# - a file with header which 1st column corresponds to the ids of the chimeric junctions that need to be detected (reference) 
#   (in the chimpipe format: donchr_donpos_donstrand:accchr_accpos_accstrand)
# - a file with header which 1st column corresponds to the chimeric junctions that are predicted by a program 
#   (in the chimpipe format: donchr_donpos_donstrand:accchr_accpos_accstrand)
# Provides:
###########
# - on the standard output tabulated information with the number junctions in each file, the number of junctions in common, the number of junctions that are 
#   in one and not in the other and a sensitivity and a precision measure (although for positive sets this precision is an underestimate of the true one 
#   since we do not know whether there are other chimeras to be found in this set)
# - a 1 column file called common.txt with the coordinates of the common chimeric junctions 
# - a 2 column tsv file called refjunc_closestpred.tsv with the reference junctions which have a close predicted junction (1st col = ref; 2nd col = close)
# - a file with the predicted junctions (same as 1st input file) but with information about all the reference junctions sharing the same chr and strand
#   for the two parts of the junction, their donor distance to the predicted junction don, their acceptor distance to the predicted junction acc
#   the sum of those and the subset of sums that are minimum together with their associated ref junctions

# usage
#######
# benchmark_chimeric_junction.sh ref_junctions.txt pred_junctions.txt > report.tsv


# example
#########
# cd ~/Chimeras/ChimPipe/benchmark/ChimPipe-0.8.8/Edgren
# ref=~/Chimeras/Benchmark/Data/Edgren/edgren_juncid_withheader.txt
# time benchmark_chimeric_junction.sh $ref chimeric_junctions_candidates_Edgren.txt > benchmark_report.tsv
# real	0m0.967s

# output1
# chr3_33055548_-:chr3_32483332_+		chr3_33055548_-:chr3_32483332_+
# 18 (2 fields)

# output2
# ref	5
# pred	195
# common	1
# ref_not_in_common	4
# pred_not_in_common	194
# sn	20
# pr	0.512821
# close25nt_not_exact	0
# samechrstr	39

# Check the input files and if OK assign to variables
#####################################################
if [ ! -e "$1" ] || [ ! -e "$2" ]
then
    echo "" >&2
    echo USAGE: >&2
    echo "" >&2
    echo "benchmark_chimeric_junction.sh ref_junctions.txt pred_junctions.txt" >&2
    echo "" >&2
    echo "Takes as input:" >&2
    echo "- a file with header which 1st column corresponds to ids of chimeric junctions that need to be detected (reference) (in the chimpipe format: donchr_donpos_donstrand:accchr_accpos_accstrand)" >&2
    echo "- a file with header which 1st column corresponds to chimeric junctions that are predicted by a program (in the chimpipe format: donchr_donpos_donstrand:accchr_accpos_accstrand)" >&2
    echo "" >&2
    echo "Provides:" >&2
    echo "- on the standard output tabulated information with the number junctions in each file, the number of junctions in common, the number of junctions that are" >&2
    echo "  in one and not in the other and a sensitivity and a precision measure (although for positive sets this precision is an underestimate of the true one" >&2
    echo "  since we do not know whether there are other chimeras to be found in this set)" >&2
    echo "- a 1 column file called common.txt with the coordinates of the common chimeric junctions" >&2
    echo "" >&2
    exit 1
else
    ref=$1
    pred=$2
fi
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

# Programs and scripts
ChimToGff=$rootDir/chim_txt_to_gff_ssext.awk
GFF2GFF=$rootDir/gff2gff.awk
OVER=$rootDir/../bin/overlap
RMRND=$rootDir/remove_redund_better.awk
COMP=$rootDir/compare_chim_junc.awk

# Compute basic numbers and print junctions that are exactly identical
######################################################################
refno=`awk 'NR>=2{print $1}' $ref | sort | uniq | wc -l | awk '{print $1}'`
predno=`awk 'NR>=2{print $1}' $pred | sort | uniq | wc -l | awk '{print $1}'`

cat $ref $pred | awk '{print $1}' | sort | uniq -c | awk '$1==2&&$2~/:/{print $2}' > common.txt
common=`wc -l common.txt | awk '{print $1}'`

sn=`echo $common $refno | awk '{print ($1/$2)*100}'`
pr=`echo $common $predno | awk '{print ($1/$2)*100}'`


# Find the junctions which two sides overlap but that are not exact
###################################################################
awk 'NR>=2' $ref | awk -v ext=25 -f $ChimToGff | awk -f $GFF2GFF > ref.gff
awk 'NR>=2' $pred | awk -v ext=25 -f $ChimToGff | awk -f $GFF2GFF > pred.gff
$OVER ref.gff pred.gff -st 1 -f pred -m 10 -nr -o ref_over_pred.gff
awk '$NF!="."{split($10,a,"\""); gsub(/\"/,"",$NF); gsub(/\;/,"",$NF); list[a[2]]=(list[a[2]])($NF)}END{for(j in list){print j, list[j]}}' ref_over_pred.gff | awk -v fldlist=junc:2 -f $RMRND | awk '{found=""; split($(NF-2),a,","); k=1; while(a[k]!=""){split(a[k],b,":"); if(b[3]==2){found=(found)","(b[1]":"b[2])} k++} if(found!=""){gsub(/\,/,"\t",found); print $1"\t"found}}' | awk '$1!=$2' > refjunc_closepred.tsv
close=`wc -l refjunc_closepred.tsv | awk '{print $1}'`

# For each predicted junction, report the list of reference junctions that have the same chr and strand
# as the predicted one on both sides, together with the distances between the predicted and the reference
# splice sites for each part of the junction, as well as the sum of those distances, and at the end
# put either the reference junction that minimizes both the donor and the acceptor distance, or if it
# does not exist, put the one that minimizes the sum of those distances (if there are many put the list)
awk -v fileRef=$ref -f $COMP $pred > pred_vs_ref.tsv
close2=`awk '$NF!="."' pred_vs_ref.tsv | wc -l`

printf "ref\t$refno\n"
printf "pred\t$predno\n"
printf "common\t$common\n"
printf "ref_not_in_common\t"$((refno-common))"\n"
printf "pred_not_in_common\t"$((predno-common))"\n"
printf "sensitivity\t$sn\n"
printf "precision\t$pr\n"
printf "close25nt_not_exact\t"$close"\n"
printf "samechrstr\t"$close2"\n"

# Clean
# rm ref.gff pred.gff
