#!/bin/bash
# annot_to_introns_with_ds_and_can.sh
######################################
# !!! can be generalized easily to any kind of junction since duplseq_gal exists !!!
# !!! needs fastalength to be installed !!!
# !!! should not be run twice in the same directory since uses files with fixed names !!!
# !!! TODO: do not crash even if a sequence to extract goes beyond the boundaries of a chromosome !!!

# From an annotation makes its introns with information about duplicated sequence (as computed by duplseq) and canonicity

# Difference with annot_to_introns_with_ds_and_can_info.sh is that it uses a single genome fasta file and the bedtoold
# getfasta program to extract sequences from it (rather than the perl script from Tyler which was very slow and was starting
# from genome files split by chr)

# To be more general this script will have to start from junctions in chimpipe's format (donchr_donpos_donstr:accchr_accpos_accstr)
# with header, as is the case for the program ~/bin/junc_txt_to_canonical.sh, and will produce the same tsv file with additional information

# More formally, this script takes as input:
############################################
# - an annotation in gtf or gff2 format with at least exon rows, and with gene id and transcript id as the first two keys
# - a fasta file for the chromosome sequences 
# and it provides as output:
############################
# - the introns from the annotation associated to both the duplicated sequence (DS) and the canonicity information
# (the DS is computed with duplseq)
# - on the standard output some statistics

# Dependencies:
###############
# - fastalength from exonerate
# - fastaFromBed from bedtools
# - duplseq from Sarah Djebali
# - 5 awk scripts (see below)


# example:
##########
# cd ~sdjebali/DuplicatedSequences/InIntrons/v20
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version20/gencode.v20.annotation.gtf
# genome=/users/rg/projects/references/Genome/H.sapiens/hg38/hg38.regChr.fa
# time annot_to_introns_with_ds_and_can.sh $annot $genome > annot_to_introns_with_ds_and_can.out 2> annot_to_introns_with_ds_and_can.err

# outputs
# an intron gff file with all info
##################################
# gencode.v20.annotation_introns_24mer_don_acc_seq_ds_can.gff
# chr10 HAVANA  intron  14300   14496   .       -       .       gene_id "ENSG00000260370.1"; transcript_id "ENST00000562162.1"; 24merdonbed: "chr10_14484_14508_-"; 24meraccbed: "chr10_14287_14311_-"; 24merdonseq ATTCCTCCAAAGGTAAAGTGTCGC 24meraccseq AATGTTTTCCAGGTGGATATTGAG duplseq: . length: . startdon: . startacc: . can 1
# 962475 (30 fields) *** real    0m21.772s
# and a report in standard output
#################################
# introns with DS >= 5nt
# 5 962475 84450 8.77425
# Read 84450 items
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   5.000   5.000   5.000   5.909   6.000  24.000 
# canonical introns with DS
# 5 778759 76014 9.76091
# Read 76014 items
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   5.000   5.000   5.000   5.673   6.000  24.000 
# non canonical introns with DS
# 5 183716 8436 4.59187
# Read 8436 items
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   5.000   5.000   6.000   8.029   8.000  24.000 

if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: annot_to_introns_with_ds_and_can.sh annot.gtf genome.fa >&2
    echo "where:" >&2
    echo "- annot.gtf is an annotation in gtf or gff2 format with at least exon rows, and with gene id and transcript id as the first two keys" >&2
    echo "- genome.fa is a fasta file containing the chromosome sequences for this annotation" >&2
    echo "will provide as output an intron file made from the annotation with both duplicated sequence (DS, computed by duplseq) and canonicity information" >&2
    echo "as well as some DS statistics depending on the canonicity on the standard output" >&2
    echo "!!! should not be run twice in the same directory since uses files with fixed names !!!" >&2
    echo "" >&2
    exit 1
fi

# Assign variables
##################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

annot=$1
genome=$2
basetmp=`basename $annot`
basetmp2=${basetmp%.gtf}
base=${basetmp2%.gff}

# Programs
##########
MAKEINTRONS=$rootDir/make_introns.awk
CUTGFF=$rootDir/cutgff.awk
GFF2GFF=$rootDir/gff2gff.awk
DUPLSEQ=$rootDir/duplseq
CANONICAL=$rootDir/intron_is_canonical_gal.awk
PROP=$rootDir/compute_prop.awk

# Makes the introns from the annotation gtf file
################################################
# !!! the exons have to be sorted according to tr id and then coord !!!
echo I am making the introns from the annotation gtf file >&2
awk '$3=="exon"' $annot | sort -k12,12 -k4,4n -k5,5n | awk -v fldgn=10 -v fldtr=12 -f $MAKEINTRONS | awk -v to=12 -f $CUTGFF > $base\_introns.gff
# chr1  HAVANA  intron  12058   12178   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2";
# 962475 (12 fields) *** real    1m4.374s
echo done >&2

# Makes the chromosome length file from the genome fasta file
#############################################################
echo I am making the chromosome length file from the genome fasta file >&2
fastalength $genome > genome.len
# 196202544 1
# 23475 (2 fields)
echo done >&2

# Make a file of 24 mers for the donor and acceptor sites of each intron
########################################################################
# the 24don seq at 5' end of the intron is composed of 12 nt exon and then 12 nt introns which start with the GT
# the 24acc seq at 3' end of the intron is composed of 12 nt intron which end with the AG and then 12 nt exon
# !!! here I should check that each feature does not go beyond the boundaries of a chromosome otherwise fastaFromBed fails !!!
echo I am making a file of coordinates for the 24mers around donor and acceptor >&2
awk -v fileRef=genome.len 'BEGIN{while (getline < fileRef >0){size[$2]=$1}} {intrbeg=$4; intrend=$5; if($7=="+"){donbeg=intrbeg-12; donend=intrbeg+12-1; accbeg=intrend-12+1; accend=intrend+12;} else{if($7=="-"){donbeg=intrend-12+1; donend=intrend+12; accbeg=intrbeg-12; accend=intrbeg+12-1;}} if((donbeg>=1)&&(donend<=size[$1])&&(accbeg>=1)&&(accend<=size[$1])){print $0, "24merdonbed", "\""$1":"(donbeg-1)":"donend":"$7"\"\;", "24meraccbeg", "\""$1":"(accbeg-1)":"accend":"$7"\"\;"}}' $base\_introns.gff | awk -f $GFF2GFF > $base\_introns_24mer_don_acc.gff
# chr1  HAVANA  intron  12058   12178   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; 24merdonbed "chr1_12045_12069_+"; 24meraccbed "chr1_12166_12190_+";
# 962475 (16 fields)
awk '{split($(NF-2),a,"\""); split(a[2],a1,":"); split($NF,b,"\""); split(b[2],b1,":"); print a1[1], a1[2], a1[3], ".", ".", a1[4]; print b1[1], b1[2], b1[3], ".", ".", b1[4]}' $base\_introns_24mer_don_acc.gff | sort | uniq | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > $base\_introns_24mers_don_acc_nr.bed
# chr10 100009825       100009849       .       .       -
# 652433 (6 fields)  *** real    0m19.786s
echo done >&2

# Extract the sequences using fastaFromBed from bedtools and put this information back to the intron file
##########################################################################################################
# !! the gff produced is not standard because there are no double quotes around the sequences but it can be used for duplseq !!!
echo I am extracting the sequences of 24nt around donor and acceptor >&2
fastaFromBed -tab -s -fi $genome -bed $base\_introns_24mers_don_acc_nr.bed -fo $base\_introns_24mer_don_acc_nr.tsv
# chr10:100009825-100009849(-)  GGCGGAAAGCAGGTCAGAGgccgg
# 652433 (2 fields)  *** real    0m32.324s
awk -v fileRef=$base\_introns_24mer_don_acc_nr.tsv 'BEGIN{while (getline < fileRef >0){split($1,a,":"); split(a[2],b,"\("); split(b[1],c,"-"); split(b[2],d,"\)"); seq["\""a[1]":"c[1]":"c[2]":"d[1]"\"\;"]=$2}} {print $0, "24merdonseq", seq[$(NF-2)], "24meraccseq",seq[$NF]}' $base\_introns_24mer_don_acc.gff | awk -f $GFF2GFF > $base\_introns_24mer_don_acc_seq.gff
# chr1  HAVANA  intron  12058   12178   .       +       .       gene_id "ENSG00000223972.5"; transcript_id "ENST00000450305.2"; 24merdonbed: "chr1_12045_12069_+"; 24meraccbed: "chr1_12166_12190_+"; 24merdonseq GTGCAAGCTGAGCACTGGAGTGGA 24meraccseq ATAGGGGAAAGATTGGAGGAAAGA
# 962475 (20 fields) *** real    0m17.557s
echo done >&2

# Find the DS
#############
echo I am finding the DS around the splice junctions >&2
$DUPLSEQ $base\_introns_24mer_don_acc_seq.gff -d 18 -a 20 -o $base\_introns_24mer_don_acc_seq_ds.gff
# chr10 HAVANA  intron  14300   14496   .       -       .        gene_id "ENSG00000260370.1"; transcript_id "ENST00000562162.1"; 24merdonbed: "chr10_14484_14508_-"; 24meraccbed: "chr10_14287_14311_-"; 24merdonseq ATTCCTCCAAAGGTAAAGTGTCGC 24meraccseq AATGTTTTCCAGGTGGATATTGAG duplseq: . length: . startdon: . startacc: .
# 962475 (28 fields)  *** real    0m25.069s
echo done >&2

# Find the canonicity
#####################
echo I am finding the canonicity of the introns >&2
awk -v flddon=18 -v fldacc=20 -f $CANONICAL $base\_introns_24mer_don_acc_seq_ds.gff | awk '{s=""; can=$NF; for(i=1; i<=(NF-1); i++){s=(s)($i)(" ")} print (s)("can ")(can)}' | awk -f $GFF2GFF > $base\_introns_24mer_don_acc_seq_ds_can.gff
# chr10 HAVANA  intron  14300   14496   .       -       .       gene_id "ENSG00000260370.1"; transcript_id "ENST00000562162.1"; 24merdonbed: "chr10_14484_14508_-"; 24meraccbed: "chr10_14287_14311_-"; 24merdonseq ATTCCTCCAAAGGTAAAGTGTCGC 24meraccseq AATGTTTTCCAGGTGGATATTGAG duplseq: . length: . startdon: . startacc: . can 1
# 962475 (30 fields) *** real    0m21.772s
echo done >&2

# Compute some statistics
##########################
echo I am computing some statistics >&2
echo "# number and % of canonical introns"
awk -f $PROP $base\_introns_24mer_don_acc_seq_ds_can.gff
# in how many introns do we find a DS >= 5nt?
##############################################
echo "# introns with DS >= 5nt"
for i in $(seq 5 15)
do
cat $base\_introns_24mer_don_acc_seq_ds_can.gff | awk -v minlg=$i '{n++; found=0; if($22!="."){split($24,a,","); k=1; while((found==0)&&(a[k]!="")){if(a[k]>=minlg){found=1;}k++;}} if(found==1){n1++;}}END{print "#", minlg, n, (n1!=""?n1:0), n1/n*100}' 
done
# 5 962475 84450 8.77425  *** 9% which is what we had before
# 6 962475 31811 3.30512
# 7 962475 13058 1.35671
# 8 962475 6447 0.669836
# 9 962475 4086 0.424531
# 10 962475 3250 0.337671
# 11 962475 2984 0.310034
# 12 962475 2675 0.277929
# 13 962475 2194 0.227954
# 14 962475 1942 0.201771
# 15 962475 1702 0.176836
cat $base\_introns_24mer_don_acc_seq_ds_can.gff | awk -v minlg=5 '{found=0; maxlg=0; if($22!="."){split($24,a,","); k=1; while(a[k]!=""){if(a[k]>maxlg){maxlg=a[k]} if(a[k]>=minlg){found=1;} k++;}} if(found==1){print maxlg}}' > $base.ds.lg.tmp
stats.sh $base.ds.lg.tmp 
# Read 84450 items
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   5.000   5.000   5.000   5.909   6.000  24.000   *** length is usually small, as we found before
# can with ds?
##############
echo "# canonical introns with DS" 
for i in $(seq 5 15)
do
cat $base\_introns_24mer_don_acc_seq_ds_can.gff | awk -v minlg=$i '$30==1{n++; found=0; if($22!="."){split($24,a,","); k=1; while((found==0)&&(a[k]!="")){if(a[k]>=minlg){found=1;}k++;}} if(found==1){n1++;}}END{print "#", minlg, n, (n1!=""?n1:0), n1/n*100}' 
done
# 5 778759 76014 9.76091
# 6 778759 27523 3.53421
# 7 778759 10234 1.31414
# 8 778759 4139 0.531487
# 9 778759 2017 0.259002
# 10 778759 1275 0.163722
# 11 778759 1079 0.138554
# 12 778759 978 0.125584
# 13 778759 696 0.089373
# 14 778759 613 0.078715
# 15 778759 555 0.0712672
cat $base\_introns_24mer_don_acc_seq_ds_can.gff | awk -v minlg=5 '$30==1{found=0; maxlg=0; if($22!="."){split($24,a,","); k=1; while(a[k]!=""){if(a[k]>maxlg){maxlg=a[k]} if(a[k]>=minlg){found=1;} k++;}} if(found==1){print maxlg}}' > $base.ds.lg.tmp
stats.sh $base.ds.lg.tmp 
# Read 76014 items
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   5.000   5.000   5.000   5.673   6.000  24.000 
# non can with ds?
##################
echo "# non canonical introns with DS" 
for i in $(seq 5 15)
do
cat $base\_introns_24mer_don_acc_seq_ds_can.gff | awk -v minlg=$i '$30==0{n++; found=0; if($22!="."){split($24,a,","); k=1; while((found==0)&&(a[k]!="")){if(a[k]>=minlg){found=1;}k++;}} if(found==1){n1++;}}END{print "#", minlg, n, (n1!=""?n1:0), n1/n*100}' 
done
# 5 183716 8436 4.59187
# 6 183716 4288 2.33404
# 7 183716 2824 1.53716
# 8 183716 2308 1.25629
# 9 183716 2069 1.12619
# 10 183716 1975 1.07503
# 11 183716 1905 1.03693
# 12 183716 1697 0.923708
# 13 183716 1498 0.815389
# 14 183716 1329 0.723399
# 15 183716 1147 0.624333
cat $base\_introns_24mer_don_acc_seq_ds_can.gff | awk -v minlg=5 '$30==0{found=0; maxlg=0; if($22!="."){split($24,a,","); k=1; while(a[k]!=""){if(a[k]>maxlg){maxlg=a[k]} if(a[k]>=minlg){found=1;} k++;}} if(found==1){print maxlg}}' > $base.ds.lg.tmp
stats.sh $base.ds.lg.tmp 
# Read 8436 items
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   5.000   5.000   6.000   8.029   8.000  24.000 
echo done >&2

# Cleaning
##########
rm genome.len
rm $base\_introns.gff
rm $base\_introns_24mer_don_acc.gff
rm $base\_introns_24mers_don_acc_nr.bed 
rm $base\_introns_24mer_don_acc_nr.tsv
rm $base\_introns_24mer_don_acc_seq.gff
rm $base\_introns_24mer_don_acc_seq_ds.gff
rm $base.ds.lg.tmp 
