#!/bin/bash

# find_retained_introns.sh
# takes as input an annotation file in gff2 or gtf format with at least exon rows and with gene_id and transcript_id information
# and outputs two 1 column file with the ids of the transcripts that have a retained intron not at the 3' utr and the ids of their genes

# the strategy si the following
# if an exon that is not the last exon of a tr strandedly overlaps 2 or more introns
# from the same other tr and if its beg corresponds to the beg of an exon from this 
# other tr and its end correspond to the end of an exon from this other tr
# then I will extract the id of the tr with this big exon and consider it a retained intron tr

# example
# cd ~/faang/projects/piglncRNA/data/lncRNA_analysis/FEELnc_winter20152016/Retained_introns/Test
# annot=/home/projects/faang/projects/piglncRNA/data/RNA_seq/STAR_cuffmerge/merged.gtf
# time ~/Scripts/Bash/find_retained_introns.sh $annot 2> find_retained_introns.err &
# real	6m00.869s 

# Checks an input file is given otherwise exits
###############################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo Usage: find_retained_introns.sh annot.gff2 >&2
    echo "Should not be run twice in the same directory" >&2
    echo "" >&2
    exit 1
else
    annot=$1
    b=`basename ${annot%.gff}`
    b2=${b%.gtf} 
fi

# Paths for programs
####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path


# Programs
##########
MAKEOK=$rootDir/make_gff_ok.awk
INTRONS=$rootDir/make_introns.awk
MOST3P=$rootDir/extract_most_3p.awk
OVER=$rootDir/overlap
REDUND=$rootDir/remove_redund_better.awk
GFF2GFF=$rootDir/gff2gff.awk

# Make an exon gff file with gene id and transcript id as the first two key,value pairs
#######################################################################################
echo Making an exon gff file with gene id and transcript id as the first two key,value pairs >&2
awk '$3=="exon"' $annot | awk -f $MAKEOK | sort -k12,12 -k4,4n -k5,5n > $b2\_exons_sorted_by_tr.gff
echo done >&2

# Make an intron gff file
#########################
echo Making an intron gff file >&2
awk -v fldgn=10 -v fldtr=12 -f $INTRONS $b2\_exons_sorted_by_tr.gff > $b2\_introns.gff
echo done >&2

# Extract most 3' exons
#######################
echo Extracting most 3\' exons >&2
awk -v fldno=12 -f $MOST3P $b2\_exons_sorted_by_tr.gff > $b2\_exons_most3p.gff
echo done >&2
# chr7  Cufflinks       exon    43026302        43026919        .       +       .       gene_id "XLOC_040584"; transcript_id "TCONS_00119233";
# 145880 (12 fields)  

# Make the file of exons that are not most 3'
############################################
echo Making a file of exons that are not the most 3\' >&2
awk -v fileRef=$b2\_exons_most3p.gff 'BEGIN{while (getline < fileRef >0){ko[$0]=1}} ko[$0]!=1' $b2\_exons_sorted_by_tr.gff > $b2\_exons_notmost3p.gff
echo done >&2
# GL892166.1    Cufflinks       exon    446     624     .       -       .       gene_id "XLOC_000014"; transcript_id "TCONS_00000014";
# 1210433 (12 fields) 

# Compute the stranded overlap between them and the introns
###########################################################
echo Computing the stranded overlap between the exons that are not the most 3\' and the introns >&2
$OVER $b2\_exons_notmost3p.gff $b2\_introns.gff -st 1 -m 12 -f intr -o $b2\_exons_notmost3p_over_intr.gff
echo done >&2
# Command line has been read
# Names of sequences retrieved
# Hashtable of feature arrays created
# Hashtable of feature arrays filling - 1st step
# Hashtable of feature arrays filled
# Hashtable of feature arrays sorted
# hashtable of feature arrays made for file2
# I have created the hashtable of indices for file2
# I have treated (sorted and put in temp file) file1 according to what the user wants
# Overlap did its work ! Now writing output file
# I have removed the temporary sorted file if it has been created
# chr1  Cufflinks       exon    14428   14605   .       +       .        gene_id "XLOC_006535"; transcript_id "TCONS_00012427"; nb_ov_intr: 0 list_intr: .
# chr1    Cufflinks       exon    252724  252794  .       -       .        gene_id "XLOC_008649"; transcript_id "TCONS_00018808"; nb_ov_intr: 12 list_intr: "TCONS_00018798";,"TCONS_00018799";,"TCONS_00018800";,"TCONS_00018801";,"TCONS_00018802";,"TCONS_00018803";,"TCONS_00018804";,"TCONS_00018805";,"TCONS_00018806";,"TCONS_00018807";,"TCONS_00018809";,"TCONS_00018810";,
# 1210433 (16 fields) *** real    0m27.099s

# Compute redundancy of these tr to see how many appear at least twice
######################################################################
echo Eliminating the corresponding transcript redundancy >&2
awk -v fldlist=overtr:16 -f $REDUND $b2\_exons_notmost3p_over_intr.gff > $b2\_exons_notmost3p_over_intr_redund.gff
echo done >&2
# chr1  Cufflinks       exon    14428   14605   .       +       .        gene_id "XLOC_006535"; transcript_id "TCONS_00012427"; nb_ov_intr: 0 list_intr: . overtr_nr: .:1,  totno: 1
# chr1    Cufflinks       exon    252724  252794  .       -       .        gene_id "XLOC_008649"; transcript_id "TCONS_00018808"; nb_ov_intr: 12 list_intr: "TCONS_00018798";,"TCONS_00018799";,"TCONS_00018800";,"TCONS_00018801";,"TCONS_00018802";,"TCONS_00018803";,"TCONS_00018804";,"TCONS_00018805";,"TCONS_00018806";,"TCONS_00018807";,"TCONS_00018809";,"TCONS_00018810";, overtr_nr: "TCONS_00018798";:1,"TCONS_00018799";:1,"TCONS_00018800";:1,"TCONS_00018801";:1,"TCONS_00018802";:1,"TCONS_00018803";:1,"TCONS_00018804";:1,"TCONS_00018805";:1,"TCONS_00018806";:1,"TCONS_00018807";:1,"TCONS_00018809";:1,"TCONS_00018810";:1,  totno: 12
# 1210433 (20 fields)

# Add the list of transcripts that have an exon (not most 3') overlapping at least two introns from the same transcript
#######################################################################################################################
echo Adding the list of transcript with an exon \(not most 3\'\) overlapping at least two introns from the same transcript >&2
awk '$16=="."{print $0, "trwith2intr:", "."}$16!="."{split($18,a,","); k=1; s=""; while(a[k]!=""){split(a[k],b,":"); if(b[2]>=2){split(b[1],c,"\""); s=(s)(c[2])(",")} k++} print $0, "trwith2intr:", (s=="" ? "." : s)}' $b2\_exons_notmost3p_over_intr_redund.gff | awk -f $GFF2GFF > $b2\_exons_notmost3p_over_intr_redund_tr2intr.gff
echo done >&2
# chr1  Cufflinks       exon    14428   14605   .       +       .       gene_id "XLOC_006535"; transcript_id "TCONS_00012427"; nb_ov_intr: 0 list_intr: . overtr_nr: .:1, totno: 1 trwith2intr: .
# 1210433 (22 fields)

# Now form each such big exon, and for each of the transcripts for which it overlaps at least 2 introns
#######################################################################################################
# I need to know whether there is one exon from this tr which beg is the same as its beg and another
####################################################################################################
# exon from this tr which ends is the same as its end
#####################################################
echo Retaining the exons \(not most 3\'\) overlapping at least two introns from the same transcript and with one from this transcript starting where it starts and another exon from this same transcript ending where it ends >&2
awk -v fileRef=$b2\_exons_sorted_by_tr.gff 'BEGIN{while (getline < fileRef >0){split($12,a,"\""); nbex[a[2]]++; beg[a[2],nbex[a[2]]]=$4; end[a[2],nbex[a[2]]]=$5}} $NF!="."{s=""; split($NF,a,","); k=1; while(a[k]!=""){b1=0; b2=0; i=1; while(b1==0&&i<=nbex[a[k]]){if(beg[a[k],i]==$4){b1=1} i++} i=1; while(b2==0&&i<=nbex[a[k]]){if(end[a[k],i]==$5){b2=1} i++} k++; if(b1==1&&b2==1){s=(s)(a[k])(",")}} print $0, "trok:", (s=="" ? "." : s)}' $b2\_exons_notmost3p_over_intr_redund_tr2intr.gff > $b2\_exons_notmost3p_over_intr_redund_tr2introk_trexok.gff
echo done >&2
# chr1  Cufflinks       exon    253254  253435  .       -       .       gene_id "XLOC_008649"; transcript_id "TCONS_00018804"; nb_ov_intr: 11 list_intr: "TCONS_00018798";,"TCONS_00018799";,"TCONS_00018800";,"TCONS_00018800";,"TCONS_00018801";,"TCONS_00018802";,"TCONS_00018802";,"TCONS_00018803";,"TCONS_00018806";,"TCONS_00018807";,"TCONS_00018809";, overtr_nr: "TCONS_00018798";:1,"TCONS_00018799";:1,"TCONS_00018800";:2,"TCONS_00018801";:1,"TCONS_00018802";:2,"TCONS_00018803";:1,"TCONS_00018806";:1,"TCONS_00018807";:1,"TCONS_00018809";:1, totno: 11 trwith2intr: TCONS_00018800,TCONS_00018802, trok: .
# 12888 (24 fields) *** real    0m10.642s

# Make the list of transcript ids with retained introns and a list of their gene ids
####################################################################################
echo Making the list of transcript id and the list of gene ids with retained introns as two 1 column files >&2
awk '$NF!="."{split($12,a,"\""); print a[2]}' $b2\_exons_notmost3p_over_intr_redund_tr2introk_trexok.gff | sort | uniq > $b2\_trid_with_retainedintron.txt
awk '$NF!="."{split($10,a,"\""); print a[2]}' $b2\_exons_notmost3p_over_intr_redund_tr2introk_trexok.gff | sort | uniq > $b2\_gnid_with_retainedintron.txt
echo done >&2
# TCONS_00000880
# 627 (1 fields) 

# Clean
#######
echo I am cleaning >&2
rm $b2\_exons_sorted_by_tr.gff
rm $b2\_introns.gff 
rm $b2\_exons_most3p.gff
rm $b2\_exons_notmost3p.gff 
rm $b2\_exons_notmost3p_over_intr.gff
rm $b2\_exons_notmost3p_over_intr_redund.gff
rm $b2\_exons_notmost3p_over_intr_redund_tr2intr.gff
echo done >&2
