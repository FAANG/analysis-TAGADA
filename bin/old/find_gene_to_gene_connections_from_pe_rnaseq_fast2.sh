#!/bin/bash

# find_gene_to_gene_connections_from_pe_rnaseq_fast2.sh
# Works for pe data, stranded or not.
# Improvement of find_gene_to_gene_connections_from_pe_rnaseq_fast.sh to be able to
# deal with standard PE strabded bam files = where the two mates point towards each other
# (will be MATE1_SENSE or MATE2_SENSE depending on which mates is on the tr strand) 
# but it can also deal with star v1 bam files which were not standard since the two 
# mates where on the same strand as the transcript.
# It uses intersectBed instead of overlap to be faster. 
# Its can use any kind of element from the gff file, eg exon, gene...
# It has the same arguments as for the 2 subscripts of chimsplice, except that instead of 
# a txt file with the gff files we now directly take a bam file. 
# It considers the mappings that just (strandedly if data is strand) overlap the elemtns, and 
# not be totally included in it like in find_gene_to_gene_connections_from_pe_rnaseq_gal.sh

# usage:
########
# find_gene_to_gene_connections_from_pe_rnaseq_fast2.sh mapping.bam annot.gff outputdir mate_strand elt
# where:
# - mapping.bam is the pe mapping file (compulsory)
# - annot.gff is the annotation file, which should at least have exons (default v19 on hg19)
# - outputdir is the directory in which the user wants the results to be stored (default cwd)
# - mate_strand is a string that specifies the directionality convention for the dataset among those: 
#   MATE1_SENSE, MATE2_SENSE, MATE_STRAND_CSHL, NONE (default NONE)
# - elt is the kind of element from the gff file that we want to overlap with the mappings (eg exon, gene ...) (default exon)
# NOTE: the mapping must be paired end and the read names must end by /1 and /2 (should be more general in the future)
# NOTE: we consider all exons of the annotation and the mappings only have to overlap them by 1 bp (not inclusion)
# NOTE: assumes that the geneid in the gff file is in column 10 (with quotes and semilcolon)
# NOTE: assumes that the geneid string does not contain any dash (-)


# In case the user does not provide any input file, an error message is raised
##############################################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: find_gene_to_gene_connections_from_pe_rnaseq_fast2.sh mapping.bam annot.gff outputdir mate_strand elt >&2
    echo "" >&2
    echo Example:  find_gene_to_gene_connections_from_pe_rnaseq_fast2.sh >&2 
    echo /users/rg/dgonzalez/Projects/TripleMama/hg19/TM05B9085PE/SAM/05B9085.merged.bam >&2  
    echo /users/rg/dgonzalez/ReferenceAnnotations/H.sapiens/gencode.v6.annot.plus.tRNAs.female.gtf >&2  
    echo ~/BreastCancer/Selection_for_RTPCR/Confirm_by_Solexa/05B9085_map.PE.40 NONE exon >&2
    echo "" >&2
    echo "Takes a PE RNAseq mapping file in bam format of an experiment, an annotation in gff format," >&2
    echo "an output directory, the directionality convention (mate_strand) of the dataset (among" >&2 
    echo "MATE1_SENSE, MATE2_SENSE,MATE_STRAND_CSHL, NONE (default NONE), and the kind of elements" >&2
    echo "we want to intersect with the mappings (eg exon, gene ...) (default exon)." >&2 
    echo "It produces the file of (directed if data is stranded) gene to gene connections" >&2  
    echo "found by the pe mappings with the number of mappings supporting the connection." >&2 
    echo "For a connection g1 to g2 to exist there must be at least one read where the first" >&2 
    echo "mate is (strandedly if data is stranded) overlapping an exon of g1 and the second" >&2  
    echo "mate is (strandedly if data is stranded) overlapping an exon of g2." >&2
    echo "" >&2
    exit 1
fi

bamfile=$1

# In case the user does not provide any annotation file
########################################################
# or output directory or mate_strand, default values are provided
#################################################################
if [ ! -n "$2" ]
then
    annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/gencode.v19.annotation.gtf
    outdir=.
    mate_strand=NONE
    elt=exon
else
    annot=$2
    if [ ! -n "$3" ]
    then
	outdir=.
	mate_strand=NONE
	elt=exon
    else
	outdir=$3
	mkdir -p $outdir 
	if [ ! -n "$4" ]
	then
	    mate_strand=NONE
	    elt=exon
	else
	    if [ "$4" != "MATE1_SENSE" ] &&  [ "$4" != "MATE2_SENSE" ] && [ "$4" != "MATE_STRAND_CSHL" ] && [ "$4" != "NONE" ]
	    then 
		echo "" >&2
		echo Usage: find_gene_to_gene_connections_from_pe_rnaseq_fast2.sh mapping.bam annot.gff outputdir mate_strand elt >&2
		echo Wrong value for mate_strand: it can only be MATE1_SENSE, MATE2_SENSE, MATE_STRAND_CSHL or NONE, not anything else >&2
		echo "" >&2
		exit 1
	    else
		mate_strand=$4
		if [ ! -n "$5" ]
		then
		    elt=exon
		else
		    elt=$5
		fi
	    fi
	fi
    fi
fi


path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

# Will exit if there is an error in a pipe
###########################################
set -e -o pipefail

# Programs
##########
CUTGFF=$rootDir/cutgff.awk
INTER=$rootDir/../bin/intersectBed
INTER2GNLIST=$rootDir/intersectBed_pebam_with_elt_to_mateid_with_geneidlist_okorientation.awk
REMOVEREDUND=$rootDir/remove_redund_better.awk
GFF2GFF=$rootDir/gff2gff.awk

# Intersect the bam file with the exons of the annotation (longest step)
########################################################################
# treat the split mappings as separate blocks
#############################################
echo I am intersecting the bam file with the exons of the annotation and providing >&2
echo for each mapping of a read the read id with the non redundant list of genes whose >&2
echo exons are overlapped by the read \(longest step, overlap is done unstrandedly but >&2
echo type of stranded bam file is considered to filter the reads overlapping the annotation\) >&2
awk -v elt=$elt '$3==elt' $annot | awk -v to=12 -f $CUTGFF | $INTER -abam $bamfile -b stdin -split -bed -wo | awk -v mate_strand=$mate_strand -f $INTER2GNLIST | awk -v fldlist=gnlist:2 -f $REMOVEREDUND | awk '{split($4,a,","); k=1; s=""; while(a[k]!=""){split(a[k],b,":"); s=(s)(b[1])(","); k++} print $1, s}' | gzip > $outdir/readid_gnlist_whoseexoverread_noredund.txt.gz
echo done >&2

# For each read that has both mates with a gene list, provide all the pairs 
###########################################################################
# of different genes where the first one is in the gene list of /1 and the 
#########################################################################
# second one is in the gene list of /2
#######################################
echo For each read that has both mates associated to a gene list >&2
echo I am making all the pairs of different genes where the first one is >&2
echo \in the list of the first mate and where the second one >&2
echo is \in the list of the second mate >&2
zcat $outdir/readid_gnlist_whoseexoverread_noredund.txt.gz | awk '{split($1,a,"/"); exists[a[1]]++; gnlist[a[1],a[2]]=$2}END{for(r in exists){s=""; if((gnlist[r,1]!="")&&(gnlist[r,2]!="")){split(gnlist[r,1],a,","); split(gnlist[r,2],b,","); k=1; while(a[k]!=""){l=1; while(b[l]!=""){if(a[k]!=b[l]){s=(s)(a[k]"-"b[l])(",");} l++} k++}} if(s!=""){print r, s}}}' | gzip > $outdir/readid_twomateswithgnlist_alldiffgnpairs_where_1stassociatedto1stmate_and2ndto2ndmate.txt.gz
echo done >&2

# Gather this information in order to report for each pair of different genes (in alphabetical order)
####################################################################################################
# the number of reads where the two mates support the pair
##########################################################
echo I am gathering this information \in order to report all the pairs of diffent genes >&2
echo \in alphabetical order\, that are supported by pe reads together with the number of reads >&2
echo supporting them >&2
zcat $outdir/readid_twomateswithgnlist_alldiffgnpairs_where_1stassociatedto1stmate_and2ndto2ndmate.txt.gz | awk '{split($2,a,","); k=1; while(a[k]!=""){split(a[k],b,"-"); if(b[1]<b[2]){print b[1], b[2]}else{print b[2], b[1]} k++}}' | sort | uniq -c | awk '{print $2, $3, $1}' > $outdir/pairs_of_diff_gn_supported_by_pereads_nbpereads.txt
echo done >&2


