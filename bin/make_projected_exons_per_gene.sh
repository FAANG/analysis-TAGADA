#!/bin/bash

# make_projected_exons_per_gene.sh
# this script takes as input an annotation gff file with gene_id in column no 10 and transcript_id in column no 12
# and an output directory, and makes a gff file of projected exons per gene

# usage
#######
# make_projected_exons_per_gene.sh annot.gff outdir

if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: make_projected_exons_per_gene.sh annot.gff outdir >&2
    echo This script takes as input an annotation gff file with gene_id \in column no 10 and transcript_id \in column no 12, and an output directory >&2
    echo and makes a gff file of projected exons per gene, \in the directory where it is launched called after the input file >&2
    echo "" >&2
    exit 1
fi

# Variable from input
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
annot=$1
b=`basename $annot`
b2tmp=${b%.gff}
b2=${b2tmp%.gtf}
outdir=$2

# Programs
##########
CUTGFF=$rootDir/cutgff.awk 
GFF2GFF=$rootDir/gff2gff.awk
MAKESP=$rootDir/makeSP

# Call the projected exons of each gene 
#######################################
echo I am making the projected exons of each gene >&2
echo "    first make a proper exon file for makeSP..." >&2
awk '$3=="exon"{gn=$10; tr=$12; $10=tr; $12=gn; $9="transcript_id"; $11="gene_id"; print}' $annot | awk -v to=12 -f $CUTGFF | awk -f $GFF2GFF > $outdir/$b2.formakeSP.gff
echo "    second apply makeSP..." >&2
$MAKESP $outdir/$b2.formakeSP.gff -f gff -v > $outdir/$b2.formakeSP.gff_segproj.gff
echo "    third make the projected exon file..." >&2
awk '{split($10,a,":"); print a[1], "makeSP", "projex", a[2], a[3], ".", (a[4]=="p" ? "+" : "-"), ".", "gene_id", "\""$12"\"\;";}' $outdir/$b2.formakeSP.gff_segproj.gff | sort | uniq | awk -f $GFF2GFF > $outdir/$b2.exonproj.gff
echo done >&2

# Cleaning
##########
echo I am cleaning >&2
rm $outdir/$b2.formakeSP.gff $outdir/$b2.formakeSP.gff_segproj.gff
echo done >&2
