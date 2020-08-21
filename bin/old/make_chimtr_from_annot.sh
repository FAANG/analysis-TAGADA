#!/bin/bash

# make_chimtr_from_annot.sh
###########################
# Given as input:
#################
# - a gene annotation (1st arg, gt file with at least exon and gene rows and where the gene id is in column 10 and transcript id is in column 12), 
# - a genome index produced by gem (2nd argument)
# - 4 numbers (3th to 6th arguments) for the number of chimeric transcripts to make on
#   - same chr, same str, ok gx order (n1)
#   - same chr, same str, ko gx order (n2)
#   - same chr, diff str (n3)
#   - diff chr (n4)
# - a set of gene biotypes to consider (7th argument) (if not provided then no filter is done)
# Provides as output in the working directory:
##############################################
# - a single fasta file with n1+n2+n3+n4 chimeric transcripts randomly selected from the gene pairs of the 4 classes
# - an Aux directory with some intermediate files
# !!! be careful: the names of the 4 classes is hard-coded here and has to match what is produced by the script !!!
# !!! make_gene_pairs_from_annot_better.awk !!!
# !!! be careful: cannot be run twice in the same directory without loosing previous outputs since uses fixed names for outputs !!!

# example
#########
# cd /no_backup/rg/sdjebali/Chimeras/Benchmark/Data/Make_Chimeras_from_Annot/Test
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/Long/gencode.v19.annotation.long.gtf
# genome=/users/rg/projects/references/Genome/H.sapiens/hg19/gemtools1.7.1-i3/Homo_sapiens.GRCh37.chromosomes.chr.M.gem
# time make_chimtr_from_annot.sh $annot $genome 150 50 50 50 ../wanted_gnbt.txt 2> make_chimtr_from_annot.err

# inputs
# protein_coding
# 1 (1 fields)

# Check if the arguments are correct
####################################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ] || [ ! -n "$5" ] || [ ! -n "$6" ]
then
    echo "" >&2
    echo "Usage: make_chimtr_from_annot.sh annot.gtf genome.gem n1 n2 n3 n4 [list_gn_bt.txt]" >&2
    echo "where" >&2
    echo "- annot.gtf is a gene annotation (with at least exon and gene rows and where gene id and transcript id are in column 10 and 12) (mandatory)" >&2
    echo "- genome.gem is a genome index produced by gem (mandatory)" >&2
    echo "- n1 is the number of chimeric transcripts to make on same chr, same str, ok gx order (mandatory)" >&2
    echo "- n2 is the number of chimeric transcripts to make on same chr, same str, ko gx order (mandatory)" >&2
    echo "- n3 is the number of chimeric transcripts to make on same chr, diff str (mandatory)" >&2
    echo "- n4 is the number of chimeric transcripts to make on diff chr (mandatory)" >&2
    echo "- list_gn_bt.txt is a list of gene biotypes to consider from the annotation (optional)" >&2
    echo "Will produce in the working directory:" >&2
    echo "- a single fasta file with n1+n2+n3+n4 chimeric transcripts randomly selected from the gene pairs of the 4 classes" >&2 
    echo "- an Aux directory with some intermediate files" >&2
    echo "NOTE1: the names of the 4 classes is hard-coded here and has to match what is produced by the script ~sdjebali/Awk/make_gene_pairs_from_annot_better.awk" >&2 
    echo "NOTE2: cannot be run twice in the same directory without loosing previous outputs since uses fixed names for outputs" >&2 
    echo "NOTE3: needs gem-retriever to be installed !!!" >&2 
    echo "" >&2
    exit 1
fi

# Assign variables
##################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
annot=$1
anntmp=`basename $annot`
anntmp2=${anntmp%.gtf}
annbase=${anntmp2%.gff}
genome=$2
n1=$3
n2=$4
n3=$5
n4=$6

# Programs
##########
GFF2GFF=$rootDir/gff2gff.awk
MAKEGNPAIRS=$rootDir/make_gene_pairs_from_annot_better.awk
MAKEINTRONS=$rootDir/make_introns.awk
RECONSTRUCT=$rootDir/reconstruct_chimtr.awk
RETRIEVER=$rootDir/../bin/gem-retriever 


# Make a file with the name of the 4 categories and the numbers of chimeric transcripts wanted (useful for the loops)
#####################################################################################################################
# !!! be careful: the names of the 4 classes have to match what is produced by the script ~sdjebali/Awk/make_gene_pairs_from_annot_better.awk !!!
printf nonoverlap_samechr_samestr_okgxorder"\t"$n1"\n" > gnpaircategory_nbwanted.tsv
printf nonoverlap_samechr_samestr_kogxorder"\t"$n2"\n" >> gnpaircategory_nbwanted.tsv
printf nonoverlap_samechr_diffstr"\t"$n3"\n" >> gnpaircategory_nbwanted.tsv
printf diffchr"\t"$n4"\n" >> gnpaircategory_nbwanted.tsv
# nonoverlap_samechr_samestr_okgxorder 150
# 4 (2 fields)

# Create an aux directory for the aux files that we want to keep
################################################################
mkdir -p Aux


# Start with the 5 actions
##########################
# 1. filter the annotation to only get the rows corresponding to a list of biotypes provided by the user
########################################################################################################
# and corresponding to genes with at least one spliced transcript. Sort the file by transcript and then beg and end
###################################################################################################################
# for many downstream steps
###########################
if [ -n "$7" ] 
then
awk -v fileRef1=$7 -v fileRef2=$annot 'BEGIN{while (getline < fileRef1 >0){ok["\""$1"\"\;"]=1;} while (getline < fileRef2 >0){if((ok2[$10]!=1)&&($3=="exon")){nbex[$10,$12]++; if(nbex[$10,$12]>=2){ok2[$10]=1;}}}} (ok2[$10]==1){k=9; while($k!=""){if($k=="gene_type"){if(ok[$(k+1)]==1){print}} k+=2}}' $annot | sort -k12,12 -k4,4n -k5,5n | awk -f $GFF2GFF > $annbase.filt.splicedgn.gtf
else
awk -v fileRef2=$annot 'BEGIN{while (getline < fileRef2 >0){if((ok2[$10]!=1)&&($3=="exon")){nbex[$10,$12]++; if(nbex[$10,$12]>=2){ok2[$10]=1;}}}} (ok2[$10]==1)' $annot | sort -k12,12 -k4,4n -k5,5n | awk -f $GFF2GFF > $annbase.filt.splicedgn.gtf
fi
# chrX	HAVANA	gene	99883667	99894988	.	-	.	gene_id "ENSG00000000003.10"; transcript_id "ENSG00000000003.10"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "TSPAN6"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "TSPAN6"; level 2; havana_gene "OTTHUMG00000022002.1";
# 242 (26 fields)
# 16504 (28 fields)
# 130007 (30 fields)
# 53809 (32 fields)
# 495999 (34 fields)
# 336699 (36 fields)
# 384954 (38 fields)
# 360265 (40 fields)
# 508411 (42 fields)
# 94980 (44 fields)
# 9233 (46 fields)
# 848 (48 fields)
# 116 (50 fields) *** real	3m0.973s

# 2. generate all the possible gene pairs from the 4 categories
###############################################################
#    1) genepairs_nonoverlap_samechr_samestr_okgxorder.tsv 
#    2) genepairs_nonoverlap_samechr_samestr_kogxorder.tsv 
#    3) genepairs_nonoverlap_samechr_diffstr.tsv 
#    4) genepairs_diffchr.tsv
awk -v fileRef=$annbase.filt.splicedgn.gtf -f $MAKEGNPAIRS $annbase.filt.splicedgn.gtf
# real	9m46.792s  *** 3 same chr files of same size (eg v19 pcg 165M) and the diff chr one much bigger (eg v19 pcg 5.9G)
# gzip all the files at the end to save space
#############################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
gzip genepairs_$cat.tsv
done
# 730M in total when zipped, instead of 6.4G

# 3. randomly sample n1 elements from 1, n2 elements from 2, n3 elements from 3 and n4 elements from 4 as the user wants
########################################################################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
zcat genepairs_$cat.tsv.gz | shuf | head -n $nb > genepairs_$cat\_$nb\_sampled.tsv
done
# ENSG00000188368.5	ENSG00000183207.8
# 150 (2 fields)
# ENSG00000123836.10	ENSG00000203857.5
# 50 (2 fields)
# ENSG00000115946.3	ENSG00000144161.8
# 50 (2 fields)
# ENSG00000116871.11	ENSG00000165507.8
# 50 (2 fields) *** real	1m32.214s

# 4. for each randomly selected gene pair, randomly pick 1 spliced transcript for each gene and then 1 donor from 1st tr and
############################################################################################################################
#    one acceptor from 2nd transcript, also in a random way
###########################################################
# Make gene pair file but with spliced transcript list information
###################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v fileRef=$annbase.filt.splicedgn.gtf 'BEGIN{OFS="\t"; while (getline < fileRef >0){if($3=="exon"){split($10,a,"\""); split($12,b,"\""); nbex[a[2],b[2]]++; if(nbex[a[2],b[2]]==2){trlist[a[2]]=(trlist[a[2]])(b[2])(",")}}}} {print $1, $2, trlist[$1], trlist[$2]}' genepairs_$cat\_$nb\_sampled.tsv > Aux/genepairs_$cat\_$nb\_sampled_splicedtrlist_eachgn.tsv
done
# ENSG00000188368.5	ENSG00000183207.8	ENST00000598490.1,ENST00000341747.3,ENST00000595750.1,ENST00000499536.2,	ENST00000595090.1,ENST00000601968.1,ENST00000596837.1,ENST00000595811.1,ENST00000594017.1,ENST00000598768.1,ENST00000596247.1,ENST00000221413.6,ENST00000413176.2,ENST00000593570.1,ENST00000594338.1,ENST00000595002.1,
# 150 (4 fields)
# ENSG00000123836.10	ENSG00000203857.5	ENST00000411990.2,ENST00000367080.3,ENST00000464777.1,ENST00000468857.1,ENST00000367079.2,ENST00000545806.1,ENST00000541914.1,ENST00000473310.1,ENST00000483688.1,	ENST00000531340.1,ENST00000369413.3,ENST00000235547.6,ENST00000492140.1,ENST00000487520.1,ENST00000528909.1,
# 50 (4 fields)
# ENSG00000115946.3	ENSG00000144161.8	ENST00000263657.2,ENST00000430742.1,ENST00000488728.1,	ENST00000409573.2,ENST00000272570.5,ENST00000476902.1,ENST00000474234.1,ENST00000466259.1,ENST00000464305.1,ENST00000495264.1,
# 50 (4 fields)
# ENSG00000116871.11	ENSG00000165507.8	ENST00000429533.2,ENST00000316156.4,ENST00000527764.1,ENST00000373150.4,ENST00000373151.2,ENST00000530729.1,ENST00000462118.1,ENST00000474796.1,ENST00000530975.1,ENST00000487131.2,ENST00000373148.4,ENST00000532131.1,ENST00000487114.1,	ENST00000496638.1,ENST00000298295.3,ENST00000448778.1,
# 50 (4 fields) *** real	0m19.871s  *** OK
# Then the gn pair files but with 1 transcript randomly selected from each list
###############################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
cat Aux/genepairs_$cat\_$nb\_sampled_splicedtrlist_eachgn.tsv | while read g1 g2 trl1 trl2
do
tr1=`echo $trl1 | awk '{split($1,a,","); k=1; while(a[k]!=""){print a[k]; k++}}' | shuf | head -1`
tr2=`echo $trl2 | awk '{split($1,a,","); k=1; while(a[k]!=""){print a[k]; k++}}' | shuf | head -1`
echo $g1 $g2 $tr1 $tr2
done > genepairs_$cat\_$nb\_sampled_rndtr_eachgn.tsv
done
# ENSG00000188368.5 ENSG00000183207.8 ENST00000499536.2 ENST00000598768.1
# 150 (4 fields)
# ENSG00000123836.10 ENSG00000203857.5 ENST00000411990.2 ENST00000235547.6
# 50 (4 fields)
# ENSG00000115946.3 ENSG00000144161.8 ENST00000488728.1 ENST00000409573.2
# 50 (4 fields)
# ENSG00000116871.11 ENSG00000165507.8 ENST00000373151.2 ENST00000496638.1
# 50 (4 fields)  *** real	0m3.427s
# Then the donor list of the first tr and the acceptor list of the second tr in addition to these 4 columns
###########################################################################################################
# !!! note: they have to exist since here we only considered spliced tr !!!
# first make all introns from the genes which biotypes are provided by the user
###############################################################################
# !!! note: the annot file is already sorted by transcript id and then beg and end so no need to redo here !!!
awk -v fldgn=10 -v fldtr=12 -f $MAKEINTRONS $annbase.filt.splicedgn.gtf > $annbase.filt.splicedgn.introns.gff
# chr7	HAVANA	intron	127228620	127229136	.	+	.	gene_id "ENSG00000004059.6"; transcript_id "ENST00000000233.5"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "ARF5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "ARF5-001"; level 2; tag "basic"; tag "appris_principal"; tag "CCDS"; ccdsid "CCDS34745.1"; havana_gene "OTTHUMG00000023246.5"; havana_transcript "OTTHUMT00000059567.2";
# 155 (28 fields)
# 250303 (30 fields)
# 132014 (32 fields)
# 134495 (34 fields)
# 146352 (36 fields)
# 216592 (38 fields)
# 40781 (40 fields)
# 4018 (42 fields)
# 373 (44 fields)
# 53 (46 fields)  *** real	0m9.160s
# then make the list of donors and acceptors of the 1st and 2nd tr respectively
###############################################################################
# !!! be careful introns boundaries are not exonic while here we want donor and acceptor coord to be exonic !!!
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v fileRef=$annbase.filt.splicedgn.introns.gff 'BEGIN{OFS="\t"; while (getline < fileRef >0){split($12,a,"\""); chr[a[2]]=$1; str[a[2]]=$7; donlist[a[2]]=(donlist[a[2]])($7=="+" ? ($4-1) : ($5+1))(","); acclist[a[2]]=(acclist[a[2]])($7=="+" ? ($5+1) : ($4-1))(",");}} {print $1, $2, $3, $4, chr[$3], str[$3], chr[$4], str[$4], (donlist[$3]!="" ? donlist[$3] : "."), (acclist[$4]!="" ? acclist[$4] : ".")}' genepairs_$cat\_$nb\_sampled_rndtr_eachgn.tsv > Aux/genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_donlist_acclist.tsv
done
# ENSG00000188368.5	ENSG00000183207.8	ENST00000499536.2	ENST00000598768.1	chr19	+	chr19	+	42814337,	49497762,49502580,
# 150 (10 fields)
# ENSG00000123836.10	ENSG00000203857.5	ENST00000411990.2	ENST00000235547.6	chr1	+	chr1	+	207222962,207228147,207235423,207236061,207236554,207236766,207237174,207238505,207241051,207241654,207242873,207243754,207244595,207244918,	120050015,120054126,120056457,
# 50 (10 fields)
# ENSG00000115946.3	ENSG00000144161.8	ENST00000488728.1	ENST00000409573.2	chr2	+	chr2	-	68400549,	112974045,112988527,112989524,112990948,112991813,112994272,112996105,113007849,
# 50 (10 fields)
# ENSG00000116871.11	ENSG00000165507.8	ENST00000373151.2	ENST00000496638.1	chr1	+	chr10	-	36622064,36636916,36637182,36638228,36639079,36640609,36642182,36642443,36643802,36644197,36644424,36644626,36644916,36645158,36645357,36645668,	45467123,
# 50 (10 fields)  **** real	0m15.431s  *** OK
# finally randomly select one donor and one acceptor from those lists
#####################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
cat Aux/genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_donlist_acclist.tsv | while read a b c d e f g h i j
do
don=`echo $i | awk '{split($1,a,","); k=1; while(a[k]!=""){print a[k]; k++}}' | shuf | head -1`
acc=`echo $j | awk '{split($1,a,","); k=1; while(a[k]!=""){print a[k]; k++}}' | shuf | head -1`
printf $a"\t"$b"\t"$c"\t"$d"\t"$e"\t"$f"\t"$g"\t"$h"\t"$don"\t"$acc"\n"
done > genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc.tsv
done
# ENSG00000188368.5	ENSG00000183207.8	ENST00000499536.2	ENST00000598768.1	chr19	+	chr19	+	42814337	49502580
# 150 (10 fields)
# ENSG00000123836.10	ENSG00000203857.5	ENST00000411990.2	ENST00000235547.6	chr1	+	chr1	+	207235423	120056457
# 50 (10 fields)
# ENSG00000115946.3	ENSG00000144161.8	ENST00000488728.1	ENST00000409573.2	chr2	+	chr2	-	68400549	112991813
# 50 (10 fields)
# ENSG00000116871.11	ENSG00000165507.8	ENST00000373151.2	ENST00000496638.1	chr1	+	chr10	-	36644197	45467123
# 50 (10 fields) *** real	0m4.671s


# 5. make the complete fasta sequence of all these transcripts while recording in each header all info about genes, 
###################################################################################################################
#    transcripts and donor and acceptor coordinates
###################################################
# !!! note: the annot file is already sorted by transcript id and then beg and end so no need to redo here !!!
# - For the donor side
#   * if tr is on + take exons from 1st to the one with gend=don coord (in THIS order)
#     otherwise take exons from last to the one with gbeg=don coord (in THIS order)
# - For the acc side
#   * if tr is on + take exons from the one with gbeg=acc coord to the last (in THIS order)
#     otherwise take exons from the one with gend=acc coord to the 1st (in THIS order)
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v fileRef=$annbase.filt.splicedgn.gtf -f $RECONSTRUCT genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc.tsv > genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist.tsv
done
# ENSG00000188368.5	ENSG00000183207.8	ENST00000499536.2	ENST00000598768.1	chr19	+	chr19	+	42814337	49502580	chr19_42812926_42814337_+,	chr19_49502580_49502621_+,
# 150 (12 fields)
# ENSG00000123836.10	ENSG00000203857.5	ENST00000411990.2	ENST00000235547.6	chr1	+	chr1	+	207235423	120056457	chr1_207222801_207222962_+,chr1_207228046_207228147_+,chr1_207235298_207235423_+,	chr1_120056457_120057681_+,
# 50 (12 fields)
# ENSG00000115946.3	ENSG00000144161.8	ENST00000488728.1	ENST00000409573.2	chr2	+	chr2	-	68400549	112991813	chr2_68400239_68400549_+,chr2_112991697_112991813_-,chr2_112990837_112990948_-,chr2_112989415_112989524_-,chr2_112988480_112988527_-,chr2_112969102_112974045_-,
# 50 (12 fields)
# ENSG00000116871.11	ENSG00000165507.8	ENST00000373151.2	ENST00000496638.1	chr1	+	chr10	-	36644197	45467123	chr1_36621803_36622064_+,chr1_36636572_36636916_+,chr1_36637114_36637182_+,chr1_36638065_36638228_+,chr1_36638965_36639079_+,chr1_36640499_36640609_+,chr1_36641800_36642182_+,chr1_36642298_36642443_+,chr1_36643474_36643802_+,chr1_36644020_36644197_+,	chr10_45466429_45467123_-,
# 50 (12 fields) *** real	0m18.986s

# Extract the exons using gem retriever for all categories
##########################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk 'BEGIN{OFS="\t"} {split($11,a,","); k=1; while(a[k]!=""){split(a[k],b,"_"); print b[1], b[4], b[2], b[3]; k++} split($12,a,","); k=1; while(a[k]!=""){split(a[k],b,"_"); print b[1], b[4], b[2], b[3]; k++}}' genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist.tsv 
done | sort | uniq > genepairs_distinct_exon_coord.tsv 
# chr10	+	115312825	115312949
# 1965 (4 fields) *** real	0m0.056s 
cat genepairs_distinct_exon_coord.tsv | $RETRIEVER $genome > genepairs_distinct_exon_coord.seq
# ATAGACAACAAAAGAAATTTTATTGAGAGGAAAACACAAGTCCTTAAACTGCAAAGATGTTTGCCAGGATGTCTGATCTCCATGTTCTGCTGTTAATGGCTCTGGTGGGAAAGACAGCCTGTGGG
# 1965 (1 fields) *** real	0m19.801s
paste genepairs_distinct_exon_coord.tsv genepairs_distinct_exon_coord.seq | awk '{print $1"_"$3"_"$4"_"$2, $5}' > genepairs_distinct_exon_coord_seq.txt
# chr10_115312825_115312949_+ ATAGACAACAAAAGAAATTTTATTGAGAGGAAAACACAAGTCCTTAAACTGCAAAGATGTTTGCCAGGATGTCTGATCTCCATGTTCTGCTGTTAATGGCTCTGGTGGGAAAGACAGCCTGTGGG
# 1965 (2 fields) *** real	0m0.063s
# And then concatenate the exons in the order they appear in the donor and acceptor list
########################################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v fileRef=genepairs_distinct_exon_coord_seq.txt 'BEGIN{OFS="\t"; while (getline < fileRef >0){seqex[$1]=$2}} {donseq=""; split($11,a,","); k=1; while(a[k]!=""){donseq=(donseq)(seqex[a[k]]); k++} accseq=""; split($12,a,","); k=1; while(a[k]!=""){accseq=(accseq)(seqex[a[k]]); k++} print $0, donseq, accseq}' genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist.tsv > Aux/genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist_donseq_accseq.tsv
done
# ENSG00000188368.5	ENSG00000183207.8	ENST00000499536.2	ENST00000598768.1	chr19	+	chr19	+	42814337	49502580	chr19_42812926_42814337_+,	chr19_49502580_49502621_+,	CCTAG...CCAG	CCACAACCAAAGTCCCGGAGATCCGTGATGTAACAAGGATTG
# 150 (14 fields)
# ENSG00000123836.10	ENSG00000203857.5	ENST00000411990.2	ENST00000235547.6	chr1	+	chr1	+	207235423	120056457	chr1_207222801_207222962_+,chr1_207228046_207228147_+,chr1_207235298_207235423_+,	chr1_120056457_120057681_+,	TTT...AAG	GTAC...TGTG
# 50 (14 fields)
# ENSG00000115946.3	ENSG00000144161.8	ENST00000488728.1	ENST00000409573.2	chr2	+	chr2	-	68400549	112991813	chr2_68400239_68400549_+,	chr2_112991697_112991813_-,chr2_112990837_112990948_-,chr2_112989415_112989524_-,chr2_112988480_112988527_-,chr2_112969102_112974045_-,	TAAG...CTTGG	GATG...TAAA
# 50 (14 fields)
# ENSG00000116871.11	ENSG00000165507.8	ENST00000373151.2	ENST00000496638.1	chr1	+	chr10	-	36644197	45467123	chr1_36621803_36622064_+,chr1_36636572_36636916_+,chr1_36637114_36637182_+,chr1_36638065_36638228_+,chr1_36638965_36639079_+,chr1_36640499_36640609_+,chr1_36641800_36642182_+,chr1_36642298_36642443_+,chr1_36643474_36643802_+,chr1_36644020_36644197_+,	chr10_45466429_45467123_-,	GGG...ACAA	GAGT...GAGG
# 50 (14 fields) *** real	0m0.057s

# Make a proper fasta file with all info in header and then a nt sequence correctly formatted (60nt per row)
############################################################################################################
# !!! include the coord of the junction in transcript space !!!
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
awk -v cat=$cat '{print ">"$3"-"$4, cat, $1, $2, $5, $6, $7, $8, $9, $10, length($13); s=($13)($14); n=length(s); n2=int(n/60); for(i=0; i<=(n2-1); i++){print substr(s,i*60+1,60)} if(n>n2*60){print substr(s,n2*60+1,n-n2*60)}}' Aux/genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist_donseq_accseq.tsv > genepairs_$cat\_$nb.fasta
done
# >ENST00000354289.4-ENST00000476862.1 nonoverlap_samechr_samestr_okgxorder ENSG00000105835.7 ENSG00000085563.10 chr7 - chr7 - 105917469 87230394 251
# 4553 (1 fields)
# 150 (11 fields)
# >ENST00000514078.1-ENST00000504807.1 nonoverlap_samechr_samestr_kogxorder ENSG00000169777.5 ENSG00000183258.7 chr5 - chr5 - 9712284 176939877 207
# 1843 (1 fields)
# 50 (11 fields)
# >ENST00000422000.1-ENST00000555554.1 nonoverlap_samechr_diffstr ENSG00000120647.5 ENSG00000118307.14 chr12 + chr12 - 514730 25267804 191
# 2581 (1 fields)
# 50 (11 fields)
# >ENST00000407793.2-ENST00000546806.1 diffchr ENSG00000013573.12 ENSG00000196511.9 chr12 + chr7 - 31247568 144380071 1665
# 1689 (1 fields)
# 50 (11 fields) *** real	0m0.065s  *** seems to be fine  (here output from another test)
# Makes a single fasta file with all these transcripts (the chimeric transcript class is in the header)
#######################################################################################################
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
cat genepairs_$cat\_$nb.fasta
done > genepairs_allcat.fasta
# >ENST00000354289.4-ENST00000476862.1 nonoverlap_samechr_samestr_okgxorder ENSG00000105835.7 ENSG00000085563.10 chr7 - chr7 - 105917469 87230394 251
# 10666 (1 fields)
# 300 (11 fields) *** final output with all needed information (also from other test than the previous steps)


# Clean
#######
rm $annbase.filt.splicedgn.gtf
rm $annbase.filt.splicedgn.introns.gff
cat gnpaircategory_nbwanted.tsv | while read cat nb
do
rm genepairs_$cat\_$nb\_sampled.tsv
rm genepairs_$cat.tsv.gz
rm genepairs_$cat\_$nb\_sampled_rndtr_eachgn.tsv
rm genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc.tsv
rm genepairs_$cat\_$nb\_sampled_rndtr_eachgn_chr_str_eachtr_rnddon_rndacc_donexlist_accexlist.tsv
rm genepairs_$cat\_$nb.fasta
done
rm genepairs_distinct_exon_coord.tsv genepairs_distinct_exon_coord.seq