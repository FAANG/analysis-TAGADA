#!/bin/bash

# find_orf_preservation_for_chimeras.sh
# takes as input 4 mandatory arguments:
# - a chimeric junction file made from stranded data with a header and with junction id in col no 1, and external beg and end in col 2 and 3
# - an annotation gtf file including at least CDS rows
# - two integers indicating the positions of the gene biotype lists for gene 1 and for gene 2
# and returns several intermediate files as well as a report containing the number of chimeric junctions
# of different categories. Note that now the nb for encompassing and overlapping is the same since only
# the internal junction coord are taken into account (due to gtex)
# IMPROVEMENTS: be able to use any position for juncid and for external beg and end

# changed on July 30th 2018:
# - ask for annot.gtf file with at least CDS rows instead of CDS file
# - can deal with any kind of gtf file with gene_id and transcript_id at any place

# usage
#######
# find_orf_preservation_for_chimeras.sh chimeric_junction_file.txt annot.gtf fldgnbt1 fldgnbt2

# example
#########
# cd ~/Blueprint/Gem4/Chimeras/AllCasesPaper/AllMappings/OrfPreservation
# cdsfile=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gen15_CDS.gtf
# time find_orf_preservation_for_chimeras.sh ../distinct_junctions_12runs_collapsed_withmaxbegandend_withexpr_eachrun_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_maxgnsim.txt $cdsfile 20 21 2> find_orf_preservation_for_chimeras.err
# where chimeric junction file looks like this
##############################################
# junc_id beg end ERR180942 ERR180943 ERR180944 ERR180945 ERR180948 ERR180950 ERR180951 ERR186015 ERR230581 ERR232403 ERR232404 ERR244135  gnlist1 gnlist2 gnnamelist1 gnnamelist2 btlist1 btlist2 maxgnsim
# chr6_29976869_+:chr6_30460129_+ 29976847 30460203 0 0 0 0 0 3 0 7 5 6 5 7 ENSG00000204622.6, ENSG00000204592.5, HLA-J, HLA-E, pseudogene, protein_coding, 96.55
# 11623 (22 fields)
# and where gtf file looks like this
####################################
# chr7	HAVANA	CDS	127228553	127228619	.	+	0	gene_id "ENSG00000004059.5"; transcript_id "ENST00000000233.5"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "ARF5"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "ARF5-001"; exon_number 1;  level 2; tag "basic"; tag "CCDS"; ccdsid "CCDS34745.1"; havana_gene "OTTHUMG00000023246.5"; havana_transcript "OTTHUMT00000059567.2";
# 1487 (30 fields)
# 116302 (32 fields)
# 104898 (34 fields)
# 141392 (36 fields)
# 304714 (38 fields)
# 52751 (40 fields)
# 6383 (42 fields)
# 473 (44 fields)
# 30 (46 fields)

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts
# - cannot be run several times in the same dir without erasing previous files


# Input from the user
#####################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] ||  [ ! -n "$14" ]
then
    echo "" >&2
    echo Usage: find_orf_preservation_for_chimeras.sh chimeric_junction_file.txt annot.gtf fldgnbt1 fldgnbt2 >&2
    echo "" >&2
    echo "Example: find_orf_preservation_for_chimeras.sh \\" >&2
    echo "distinct_junctions_12runs_collapsed_withmaxbegandend_withexpr_eachrun_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_maxgnsim.txt \\" >&2
    echo "/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/gen15_CDS.gtf 20 21" >&2
    echo "" >&2
    echo Takes as input:
    echo 1\) a chimeric junction file made from stranded experiments with a header and with junction id in col 1 and external beg and end in col 2 and 3 >&2
    echo 2\) an annotation gtf file including at least CDS rows >&2
    echo 3\) the column number in the junction file, for the biotypes of the genes whose exons overlap part 1 of the junction >&2
    echo 4\) the column number in the junction file, for the biotypes of the genes whose exons overlap part 2 of the junction >&2
    echo and produces several intermediate files as well as statistics in the log file >&2
    echo "" >&2
    exit 1
else
juncfile=$1
annot=$2
fldgnbt1=$3
fldgnbt2=$4
fi

path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

basetmp=${juncfile%.tsv}
base=${base%.txt}

# Programs
###########
MAKEOK=$rootDir/make_gff_ok.awk
GFF2GFF=$rootDir/gff2gff.awk
CUTGFF=$rootDir/cutgff.awk
GFF2GFF=$rootDir/gff2gff.awk
OVERLAP=$rootDir/overlap

# 0) make a cds file ok for overlap and that has an additional key,value for the frame (since not possible to report it otherwise)
##################################################################################################################################
echo I am making a cds file with a format that is understood by overlap >&2
awk '$3=="CDS"' $annot | awk -f $MAKEOK | awk -f $GFF2GFF | awk -v to=20 -f $CUTGFF | awk '{print $0, "frame:", $8}' |  awk -f $GFF2GFF > CDS_cut20_frame.gff
echo done >&2
echo Here is the distribution of the frames in this file >&2
awk '{print $NF}' CDS_cut20_frame.gff | sort | uniq -c >&2
# chr1  HAVANA  CDS     69091   70005   .       +       0       gene_id "ENSG00000186092.4"; transcript_id "ENST00000335137.3"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "OR4F5"; transcript_type "protein_coding"; frame: 0
# 723784 (22 fields)

# 1) identify the chimeric junctions that connect two pcg
#########################################################
echo I am identifying the chimeric junctions that connect two pcg >&2
awk -v fldgnbt1=$fldgnbt1 -v fldgnbt2=$fldgnbt2 'NR==1||($fldgnbt1~/protein_coding/&&$fldgnbt2~/protein_coding/)' $juncfile > distinct_junctions_pcgpcg.txt
echo done >&2
nbtotjunc=`wc -l $juncfile | awk '{print ($1)-1}'`
nbjuncpcgpcg=`wc -l distinct_junctions_pcgpcg.txt | awk '{print ($1)-1}'`
echo There are $nbjuncpcgpcg chimeric junctions that connect two pcg out of the initial $nbtotjunc >&2
echo init: $nbtotjunc > orf_report.txt
echo pcg-pcg: $nbjuncpcgpcg >> orf_report.txt
# junc_id       beg     end     samechrstr      okgxorder       dist    ss1     ss2     gnlist1 gnlist2 gnname1 gnname2 gnbt1   gnbt2   LID16627        LID16628        LID16629  LID16630 LID16631        LID16632        LID16633        LID16634        LID16635        LID16636        LID44497        LID44498        LID44499        LID44594        LID45016  LID45017 LID46598        LID46599        LID8461 LID8462 LID8463 LID8464 LID8686 LID8687 LID8692 LID8701 LID8710 LID8711 LID8963 LID8964 LID8965 LID8966 LID8967 LID8968 LID8969 LID8970
# 469 (50 fields)

# 2) identify the chimeric junctions that connect two pct, ie where the two junction parts 
##########################################################################################
#    overlap a coding exon (CDS) on each side
##############################################
echo I am identifying the chimeric junctions that connect two CDSs >&2
# a. make a gff file with the two parts for each junction and the number of the part of the junction
####################################################################################################
echo "   I am making a gff file with the two parts for each junction and the number of the part of the junction" >&2
# correct but just one base on each side, in order to process Gtex samples that is unstranded and done with 0.6.1 where there was a bug in beg and end
awk 'NR>=2{split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); print a1[1], ".", ".", a1[2], a1[2], ".", a1[3], ".", "junc:", $1, "no:", 1; print a2[1], ".", ".", a2[2], a2[2], ".", a2[3], ".", "junc:", $1, "no:", 2;}' distinct_junctions_pcgpcg.txt | awk -f $GFF2GFF > distinct_junctions_pcgpcg.gff
# small improvement, in order to process Gtex samples that is unstranded and done with 0.6.1 where there was a bug in beg and end but not totally correct
# awk 'NR>=2{split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); print a1[1], ".", ".", (($2<=a1[2]) ? $2 : a1[2]), (($2<=a1[2]) ? a1[2] : $2), ".", a1[3], ".", "junc:", $1, "no:", 1; print a2[1], ".", ".", ((a2[2]<=$3) ? a2[2] : $3), ((a2[2]<=$3) ? $3 : a2[2]), ".", a2[3], ".", "junc:", $1, "no:", 2;}' distinct_junctions_pcgpcg.txt | awk -f $GFF2GFF > distinct_junctions_pcgpcg.gff
# initially done, for example on mouse encode where the beg and end were fine but for - strand gx order and not biological order 
# awk 'NR>=2{split($1,a,":"); split(a[1],a1,"_"); split(a[2],a2,"_"); print a1[1], ".", ".", $2, a1[2], ".", a1[3], ".", "junc:", $1, "no:", 1; print a2[1], ".", ".", a2[2], $3, ".", a2[3], ".", "junc:", $1, "no:", 2;}' distinct_junctions_pcgpcg.txt | awk -f $GFF2GFF > distinct_junctions_pcgpcg.gff
echo "   done" >&2
# chr16 .       .       1823216 1823216 .       +       .       junc: chr16_1823216_+:chr16_1823706_+ no: 1
# 936 (12 fields)
# b. intersect this file with the CDS rows, reporting both the tr list and the frame list, with no nr to keep same order
########################################################################################################################
echo "   I am intersecting this gff file with the CDS file" >&2
$OVERLAP distinct_junctions_pcgpcg.gff CDS_cut20_frame.gff -st 1 -f cds -m 12,22 -v | awk '{gsub(/"/,"",$16); gsub(/;/,"",$16); print $0}' | awk -f $GFF2GFF > distinct_junctions_pcgpcg_trlistwithcdsoverpart_framelist.gff
echo "   done" >&2
# chr1  .       .       1164326 1164326 .       -       .       junc: chr1_1170109_-:chr1_1164326_- no: 2 nb_ov_cds_12: 0 list_cds_12: . nb_ov_cds_22: 0 list_cds_22: .
# chr1    .       .       9833452 9833452 .       -       .       junc: chr1_9931245_-:chr1_9833452_- no: 2 nb_ov_cds_12: 3 list_cds_12: ENST00000361311.4,ENST00000377288.3,ENST00000377298.4, nb_ov_cds_22: 3 list_cds_22: 2,2,2,
# 936 (20 fields)
# c. put everything together in the junction file and report the two trlists and the two framelists as well
###########################################################################################################
echo "   I am putting everything together in the junction file and reporting the two trlists and the two framelists as well" >&2
awk -v fileRef=distinct_junctions_pcgpcg_trlistwithcdsoverpart_framelist.gff 'BEGIN{while (getline < fileRef >0){if($14>0){nb[$10]++; trlist[$10,$12]=$16; framelist[$10,$12]=$20}}} NR==1{print $0, "trlist_part1", "trlist_part2", "framelist_part1", "framelist_part2"}nb[$1]==2{print $0, trlist[$1,1], trlist[$1,2], framelist[$1,1], framelist[$1,2]}' distinct_junctions_pcgpcg.txt > distinct_junctions_pcgpcg_pctpct.txt
echo "   done" >&2
echo done >&2
nbjuncpctpct=`wc -l distinct_junctions_pcgpcg_pctpct.txt | awk '{print ($1)-1}'`
echo There are $nbjuncpctpct chimeric junctions that connect two pct out of the initial $nbtotjunc >&2
echo pct-pct: $nbjuncpctpct >> orf_report.txt
# junc_id       beg     end     samechrstr      okgxorder       dist    ss1     ss2     gnlist1 gnlist2 gnname1 gnname2 gnbt1   gnbt2   LID16627        LID16628        LID16629  LID16630 LID16631        LID16632        LID16633        LID16634        LID16635        LID16636        LID44497        LID44498        LID44499        LID44594        LID45016  LID45017 LID46598        LID46599        LID8461 LID8462 LID8463 LID8464 LID8686 LID8687 LID8692 LID8701 LID8710 LID8711 LID8963 LID8964 LID8965 LID8966 LID8967 LID8968 LID8969 LID8970 trlist_part1 trlist_part2 framelist_part1 framelist_part2
# chr5_39074213_-:chr5_39141243_- 39074263        39141184        1 0 na GT AG ENSG00000164327.8 ENSG00000082074.10 RICTOR FYB protein_coding protein_coding      0       0       0 00       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0       0 03       14      0       0       0       0       0       0 ENST00000296782.5,ENST00000357387.3,ENST00000511516.1, ENST00000351578.6,ENST00000505428.1,ENST00000512982.1,ENST00000515010.1,ENST00000540520.1, 2,2,2, 1,1,1,1,1,
# 158 (54 fields)

# 3) identify the chimeric junctions that have a coding exon (CDS) encompassing them on each side
#################################################################################################
#    and report both the lists of transcripts associated to those coding exons on each side and
###############################################################################################
#    their associated frames
############################
echo I am identifying the chimeric junctions that have a coding exon \(CDS\) encompassing them on each side and reporting >&2
echo both the lists of transcripts associated to those coding exons on each side and their associated frames >&2
# a. intersect the two parts of each junction with the CDS exons but asking the CDS to include the junction part
################################################################################################################
# and report both the transcript list and the frame list, with no redundancy removal for the lists to be in same order
######################################################################################################################
echo "   I am intersecting the two parts of each junction with the CDS exons but asking the CDS to include the junction part" >&2
echo "   and reporting both the transcript list and the frame list, with no redundancy removal for the lists to be in same order" >&2
$OVERLAP distinct_junctions_pcgpcg.gff CDS_cut20_frame.gff -st 1 -i 1 -f cds -m 12,22 -v | awk '{gsub(/"/,"",$16); gsub(/;/,"",$16); print $0}' | awk -f $GFF2GFF > distinct_junctions_pcgpcg_trlistwithcdsincludingpart_cdslist.gff
echo "   done" >&2
# chr1  .       .       1164326 1164326 .       -       .       junc: chr1_1170109_-:chr1_1164326_- no: 2 nb_i1_cds_12: 0 list_cds_12: . nb_i1_cds_22: 0 list_cds_22: .
# 936 (20 fields)
# b. put everything together in the junction file = select the cases that have such CDS on both parts 
#####################################################################################################
# and report for those both the two transcript lists and the two CDS frame lists
################################################################################
echo "   I am putting everything together in the junction file = select the cases that have such CDS on both parts" >&2
echo "   and reporting for those both the two transcript lists and the two CDS frame lists" >&2
awk '{$(NF-3)=""; $(NF-2)=""; $(NF-1)=""; $NF=""; print}' distinct_junctions_pcgpcg_pctpct.txt | awk -v fileRef=distinct_junctions_pcgpcg_trlistwithcdsincludingpart_cdslist.gff 'BEGIN{while (getline < fileRef >0){if($14>0){nb[$10]++; trlist[$10,$12]=$16; framelist[$10,$12]=$20}}} NR==1{print $0, "trlist_part1", "trlist_part2", "framelist_part1", "framelist_part2"}nb[$1]==2{print $0, trlist[$1,1], trlist[$1,2], framelist[$1,1], framelist[$1,2]}' > distinct_junctions_pcgpcg_pctpct_cdsincludingpart.txt
echo "   done" >&2
echo done >&2
# junc_id beg end samechrstr okgxorder dist ss1 ss2 gnlist1 gnlist2 gnname1 gnname2 gnbt1 gnbt2 LID16627 LID16628 LID16629 LID16630 LID16631 LID16632 LID16633 LID16634 LID16635 LID16636 LID44497 LID44498 LID44499 LID44594 LID45016 LID45017 LID46598 LID46599 LID8461 LID8462 LID8463 LID8464 LID8686 LID8687 LID8692 LID8701 LID8710 LID8711 LID8963 LID8964 LID8965 LID8966 LID8967 LID8968 LID8969 LID8970     trlist_part1 trlist_part2 framelist_part1 framelist_part2
# 158 (54 fields)

# 4) identify the chimeric junctions where there exists a transcript t1 on one side and a
########################################################################################## 
#    transcript t2 on the other side, where the CDS exons that encompass the junction parts
###########################################################################################
#    are in the same frame
##########################
echo I am identifying the chimeric junctions where there exists a transcript t1 on one side and a transcript t2 >&2
echo on the other side, where the CDS exons that encompass the junction parts are in the same frame >&2
awk '{found=0; for(i=0; i<=2; i++){for(j=1; j<=2; j++){ok[i,j]=0}} split($(NF-1),a,","); k=1; while(a[k]!=""){ok[a[k],1]=1; k++} split($NF,a,","); k=1; while(a[k]!=""){ok[a[k],2]=1; k++} k=0; while((found==0)&&(k<=2)){if((ok[k,1]==1)&&(ok[k,2]==1)){found=1} k++} if((NR==1)||(found==1)){print}}' distinct_junctions_pcgpcg_pctpct_cdsincludingpart.txt > distinct_junctions_pcgpcg_pctpct_cdsincludingpart_sameframebothparts.txt
echo done >&2
nbjuncwithcds1=`wc -l distinct_junctions_pcgpcg_pctpct_cdsincludingpart.txt | awk '{print ($1)-1}'`
nbjuncorfok1=`wc -l distinct_junctions_pcgpcg_pctpct_cdsincludingpart_sameframebothparts.txt | awk '{print ($1)-1}'`
proporfok1=`echo $nbjuncwithcds1 $nbjuncorfok1 | awk '{print ($2/$1)*100}'`
echo The proportion of junctions with a coding exon encompassing each part of the junction and where the ORF is preserved is $proporfok1 >&2
echo with-encompassing-cds: $nbjuncwithcds1 >> orf_report.txt
echo same-frame1: $nbjuncorfok1 >> orf_report.txt
echo prop1: $proporfok1 >> orf_report.txt
# junc_id beg end samechrstr okgxorder dist ss1 ss2 gnlist1 gnlist2 gnname1 gnname2 gnbt1 gnbt2 LID16627 LID16628 LID16629 LID16630 LID16631 LID16632 LID16633 LID16634 LID16635 LID16636 LID44497 LID44498 LID44499 LID44594 LID45016 LID45017 LID46598 LID46599 LID8461 LID8462 LID8463 LID8464 LID8686 LID8687 LID8692 LID8701 LID8710 LID8711 LID8963 LID8964 LID8965 LID8966 LID8967 LID8968 LID8969 LID8970     trlist_part1 trlist_part2 framelist_part1 framelist_part2
# 78 (54 fields)

# 5) what if I compute the same prop but for the cases that have an overlapping CDS (and not encompassing)?
############################################################################################################
echo I am doing the same but relaxing the encompassing requirement into a simple overlap requirement >&2
awk '{found=0; for(i=0; i<=2; i++){for(j=1; j<=2; j++){ok[i,j]=0}} split($(NF-1),a,","); k=1; while(a[k]!=""){ok[a[k],1]=1; k++} split($NF,a,","); k=1; while(a[k]!=""){ok[a[k],2]=1; k++} k=0; while((found==0)&&(k<=2)){if((ok[k,1]==1)&&(ok[k,2]==1)){found=1} k++} if((NR==1)||(found==1)){print}}' distinct_junctions_pcgpcg_pctpct.txt > distinct_junctions_pcgpcg_pctpct_sameframebothparts.txt
echo done >&2
nbjuncwithcds2=`wc -l distinct_junctions_pcgpcg_pctpct.txt | awk '{print ($1)-1}'`
nbjuncorfok2=`wc -l distinct_junctions_pcgpcg_pctpct_sameframebothparts.txt | awk '{print ($1)-1}'`
proporfok2=`echo $nbjuncwithcds2 $nbjuncorfok2 | awk '{print ($2/$1)*100}'`
echo The proportion of junctions with a coding exon overlapping each part of the junction and where the ORF is preserved is $proporfok2 >&2
echo with-overlapping-cds: $nbjuncwithcds2 >> orf_report.txt
echo same-frame2: $nbjuncorfok2 >> orf_report.txt
echo prop2: $proporfok2 >> orf_report.txt
# junc_id       beg     end     samechrstr      okgxorder       dist    ss1     ss2     gnlist1 gnlist2 gnname1 gnname2 gnbt1   gnbt2   LID16627        LID16628        LID16629  LID16630 LID16631        LID16632        LID16633        LID16634        LID16635        LID16636        LID44497        LID44498        LID44499        LID44594        LID45016  LID45017 LID46598        LID46599        LID8461 LID8462 LID8463 LID8464 LID8686 LID8687 LID8692 LID8701 LID8710 LID8711 LID8963 LID8964 LID8965 LID8966 LID8967 LID8968 LID8969 LID8970 trlist_part1 trlist_part2 framelist_part1 framelist_part2
# 78 (54 fields)

# 6) Gather all the important information in a single file
##########################################################
echo I am gathering all important information \(pcg, cds, orf\) \in a single junction file >&2
awk -v fileRef1=distinct_junctions_pcgpcg.txt -v fileRef2=distinct_junctions_pcgpcg_pctpct_cdsincludingpart.txt -v fileRef3=distinct_junctions_pcgpcg_pctpct_sameframebothparts.txt 'BEGIN{OFS="\t"; while (getline < fileRef1 >0){pcg[$1]=1} while (getline < fileRef2 >0){cds[$1]=1} while (getline < fileRef3 >0){orf[$1]=1}} NR==1{print $0, "pcg", "cds", "orf"} NR>=2{print $0, nn(pcg[$1]), nn(cds[$1]), nn(orf[$1])} function nn(x){return (x==1 ? 1 : 0)}' $juncfile > $base.pcg.cds.orf.tsv
echo done >&2

# 7) clean
##########
echo I am cleaning >&2
rm CDS_cut20_frame.gff
rm distinct_junctions_pcgpcg.gff
rm distinct_junctions_pcgpcg_trlistwithcdsoverpart_framelist.gff
rm distinct_junctions_pcgpcg_trlistwithcdsincludingpart_cdslist.gff
echo done >&2





