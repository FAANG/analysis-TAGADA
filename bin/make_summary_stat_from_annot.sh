#!/bin/bash

# make_summary_stat_from_annot.sh
#################################
# takes as input:
################
# - a gff file of which it only considers the exon rows and assumes the gnid is in field no 10 and trid in field no 12
# it outputs:
#############
# - a two row files with some numbers of annotated elements (nbex nbdistinctex nbtr nbgn nbintrons nbdistinctintrons)
# - some intermediate useful files:
#   * $b2\_complete.gff: the complete annotation file with exon, intron, tr and genes
#   * $b2\_distinct_exons.txt: the file of distinct exon coordinates (1/gff based)
#   * $b2\_distinct_introns.txt: the file of distinct intron coordinates (1/gff based)
#   * $b2\_distinct_introns_bed.txt: the file of distinct intron coordinates (bed based)
#   * $b2\_distinct_introns_minus1nteachside.txt: the file of distinct exon-exon junction (1/gff based); 
#     junction coordinates can be obtained by going extending the intron coord by 1 nt on each side, which explains the current name (should be changed at one point?)
#   * $b2\_trids.txt: the file of transcript ids
#   * $b2\_gnids.txt: the file of gene ids
#   * $b2\_trid_nbex.txt: the file of tr id with nb of exons
#   * $b2\_gnid_nbtr.txt: the file of gene id with nb of transcripts
# from this it is then easy to retrieve additionally
####################################################
# - nbexpertr
# - nbexpergn
# - nbdistinctexpergn
# - nbtrpergn
# - nbintronpertr
# - nbintronpergn
# - nbdistinctintronpergn
# by simply running this command:
#################################
# awk 'NR==2' $output_from_this_script | awk 'BEGIN{OFS="\t"; print "nbex\tnbdistinctex\tnbtr\tnbgn\tnbintrons\tnbdistinctintrons\tnbexpertr\tnbexpergn\tnbdistinctexpergn\tnbtrpergn\tnbintronpertr\tnbintronpergn\tnbdistinctintronpergn"}{print $0, $1/$3, $1/$4, $2/$4, $3/$4, $5/$3, $5/$4, $6/$4}'

# example:
##########
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version10/Long/Element/gen10.long.exon.gtf
# time make_summary_stat_from_annot.sh $annot > gen10.long.sumstat.txt 2> gen10.long.sumstat.err
# real	0m47.351s

# NOTE: on Dec 13th 2016 added generation of files with nb of tr for each gene and nb of ex for each tr
# in order to be able to compute the distribution of those things
# On Dec 15th 2016 replaced the _ delimitor by : in order to be able to accept ncbi annotation

if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: make_summary_stat_from_annot.sh annot.gff >&2
    echo "" >&2
    exit 1
fi

path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

annot=$1
b=`basename ${annot%.gff}`
b2=${b%.gtf} 

# Programs
##########
MAKEOK=$rootDir/make_gff_ok.awk
INTRONS=$rootDir/make_introns.awk
GFF2GFF=$rootDir/gff2gff.awk
BOUNDARIES=$rootDir/compute_boundaries.awk
STATS=$rootDir/stats.sh

# Make necesary gff files for the stats
echo Making the necesary gff files for the stats >&2
awk '$3=="exon"' $annot | awk -f $MAKEOK | awk -f $GFF2GFF | sort -k12,12 -k4,4n -k5,5n > $b2\_exons_sorted_by_tr.gff
awk -v fldgn=10 -v fldtr=12 -f $INTRONS $b2\_exons_sorted_by_tr.gff | awk -f $GFF2GFF > $b2\_introns.gff
awk -v toadd=transcript -v fldno=12 -f $BOUNDARIES $b2\_exons_sorted_by_tr.gff | awk -v fileRef=$b2\_exons_sorted_by_tr.gff 'BEGIN{while (getline < fileRef >0){gnid[$12]=$10}} {print $1, $2, $3, $4, $5, $6, $7, $8, "gene_id", gnid[$10], $9, $10}' | awk -f $GFF2GFF > $b2\_transcripts.gff
awk -v toadd=gene -v fldno=10 -f $BOUNDARIES $b2\_exons_sorted_by_tr.gff | awk '{print $0, "transcript_id", $NF}' | awk -f $GFF2GFF > $b2\_genes.gff
echo done >&2

# Make necesary txt files for the simple stats and distributions
echo Making the necesary txt files for the simple stats >&2
awk '{print $1":"$4":"$5":"$7}' $b2\_exons_sorted_by_tr.gff | sort | uniq > $b2\_distinct_exons.txt
awk '{print $1":"$4":"$5":"$7}' $b2\_introns.gff | sort | uniq > $b2\_distinct_introns.txt
awk '{print $1":"($4-1)":"$5":"$7}' $b2\_introns.gff | sort | uniq > $b2\_distinct_introns_bed.txt
awk '{print $1":"($4-1)":"($5+1)":"$7}' $b2\_introns.gff | sort | uniq > $b2\_distinct_introns_minus1nteachside.txt
awk '{split($12,a,"\""); print a[2]}' $b2\_transcripts.gff | sort | uniq > $b2\_trids.txt
awk '{split($10,a,"\""); print a[2]}' $b2\_genes.gff | sort | uniq > $b2\_gnids.txt
echo done >&2
echo Making the necesary txt files for the distributions >&2
awk '{nbex[$12]++}END{for(t in nbex){print t, nbex[t]}}' $b2\_exons_sorted_by_tr.gff > $b2\_trid_nbex.txt
awk '{seen[$12,$10]++; if(seen[$12,$10]==1){nbtr[$10]++}}END{for(g in nbtr){print g, nbtr[g]}}' $b2\_exons_sorted_by_tr.gff > $b2\_gnid_nbtr.txt
echo done >&2

# Make the simple stats and distributions
echo Making the simple stats >&2
nbex=`wc -l $b2\_exons_sorted_by_tr.gff | awk '{print $1}'`
nbdistinctex=`wc -l $b2\_distinct_exons.txt | awk '{print $1}'`
nbtr=`wc -l $b2\_trids.txt | awk '{print $1}'`
nbgn=`wc -l $b2\_gnids.txt | awk '{print $1}'`
nbintrons=`wc -l $b2\_introns.gff | awk '{print $1}'`
nbdistinctintrons=`wc -l $b2\_distinct_introns.txt | awk '{print $1}'`
printf "nbex\tnbdistinctex\tnbtr\tnbgn\tnbintrons\tnbdistinctintrons\n"
printf $nbex"\t"$nbdistinctex"\t"$nbtr"\t"$nbgn"\t"$nbintrons"\t"$nbdistinctintrons"\n"
echo done >&2
echo Making the distributions >&2
printf "distribution of the number of exons per transcript\n"
$STATS $b2\_trid_nbex.txt 2 
printf "distribution of the number of transcripts per gene\n"
$STATS $b2\_gnid_nbtr.txt 2
echo done >&2

# Clean
echo Cleaning >&2
cat $b2\_exons_sorted_by_tr.gff $b2\_introns.gff $b2\_transcripts.gff $b2\_genes.gff > $b2\_complete.gff
rm $b2\_exons_sorted_by_tr.gff $b2\_introns.gff $b2\_transcripts.gff $b2\_genes.gff
rm $b2\_distinct_exons.txt $b2\_distinct_introns.txt $b2\_distinct_introns_bed.txt $b2\_distinct_introns_minus1nteachside.txt $b2\_trids.txt $b2\_gnids.txt 
echo done >&2
