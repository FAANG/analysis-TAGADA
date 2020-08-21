#!/bin/bash

# make_chimeric_junction_matrix5.sh
###################################
# this script takes as input the CANDIDATE chimeric junctions obtained with chimpipe version 0.9.4 at least
# or the concatenation of the filtered and the final junctions where the second header has been removed
# and produces a matrix with all the chimeric junctions from different experiments 
# with at least x spanning reads and y paired end reads in one expt or tech rep and no too high similarity between the two connected 
# genes (less than 30 nt sim or less than 80% sim), with their beg and end across all experiments, 
# whether the two parts are on the same chr and strand, whether they are on the expected 
# genomic order, their distance, the list of genes which exons overlap each part of the junction, 
# the real names of those genes, the biotypes of those genes and the number of spanning and discordant pe reads
# supporting the junction in each expt or tech rep
# started in april 2017

# !!! should not be run twice in the same dir at the same time (common files) !!!
# !!! the header might be wrong it would be better to take it from input file !!!

# More precisely, this script takes as input:
#############################################
# 1) a 2 column file containing the chimeric junction files associated to each exp or tech rep
#   = where column 1 is a unique identifier for the expt or tech rep and column 2 is the absolute
#   path to the the candidate file of chimeric junctions in this expt with nbspan, nbdiscpe and all info
# for example this one
# /work/project/fragencode/results/rnaseq/sus_scrofa/liver/pig1/chimpipe/Sus-scrofa-11203-Foie_TTAGGC_L001_Pool1_Foie.9752/chimericJunctions_candidates_Sus-scrofa-11203-Foie_TTAGGC_L001_Pool1_Foie.9752.txt
# juncCoord	type	filtered	reason	nbTotal(spanning+consistent)	nbSpanningReads	nbStaggered	percStaggered	nbMulti	percMulti	nbConsistentPE	nbInconsistentPE	percInconsistentPE	overlapA	overlapB	distExonBoundaryA	distExonBoundaryB	blastAlignLen	blastAlignSim	donorSS	acceptorSS	beg	end	sameChrStr	okGxOrder	dist	gnIdsA	gnIdsB	gnNamesA	gnNamesB	gnTypesA	gnTypesB	juncSpanningReadsIds	consistentPEIds	inconsistentPEIds
# JH118593.1_108648_-:JH118593.1_86156_-	readthrough	1	consistentPE,spanningReads,totalSupport,	2	1	1	100	0	0	1	0	0	100	22.1311	0	95	na	na	GT	AG	108675	86035	1	1	22492	ENSSSCG00000023895	ENSSSCG00000027705	na	na	na	na	ST-J00115:48:H# 12428 (35 fields)
# 12428 (35 fields)

# note: it used to be like this for make_chimeric_junction_matrix4.sh
######################################################################
# juncCoord     totalNbPE       nbSpanningPE    nbStag  percStag        nbMulti percMulti       nbDiscordantPE  nbInconsistentPE        percInconsistentPE      overlapA        overlapB   distExonBoundaryA       distExonBoundaryB       maxBLastAlignLen        blastAlignSim   donorSS acceptorSS      beg     end     sameChrStr      okGxOrder       dist    gnIdA       gnIdB  gnNameA gnNameB gnTypeA gnTypeB juncSpanningReadIds     supportingPairsIds      inconsistentPairsIds
# chrX_100646810_+:chr18_13919745_-       28      28      12      42.8571 28      100     0       0       na      100     100     0       177     321     96.26   GT      AG      100646777  13919690        0       na      na      ENSG00000241343.3,ENSG00000257529.1     ENSG00000243779.1       RPL36A,RP1-164F3.9      RP11-681N23.1   protein_coding,protein_coding      pseudogene      SINATRA_0006:3:109:18510:14395#0/1,SINATRA_0006:1:95:16661:9220#0/1,SINATRA_0006:1:95:19257:7963#0/1,SINATRA_0006:1:78:8979:7945#0/1,SINATRA_0006:2:103:7300:8834#0/1,SINATRA_0006:1:95:6274:6674#0/1,SINATRA_0006:1:95:12119:12408#0/1,SINATRA_0006:1:65:19204:4609#0/1,SINATRA_0006:2:26:15962:4612#0/1,SINATRA_0006:1:78:9741:2208#0/2,SINATRA_0006:3:55:3004:18029#0/1,SINATRA_0006:2:28:17045:15305#0/1,SINATRA_0006:3:28:12150:19552#0/1,SINATRA_0006:1:57:14946:18483#0/1,SINATRA_0006:3:40:17960:10169#0/1,SINATRA_0006:2:78:4492:5160#0/1,SINATRA_0006:1:95:12206:3223#0/1,SINATRA_0006:1:87:6224:19402#0/2,SINATRA_0006:3:117:19626:18481#0/1,SINATRA_0006:2:13:1765:2306#0/1,SINATRA_0006:3:1:10757:14190#0/1,SINATRA_0006:2:69:18762:13853#0/2,SINATRA_0006:1:58:4291:4933#0/1,SINATRA_0006:1:36:4284:12343#0/2,SINATRA_0006:1:95:7884:9420#0/1,SINATRA_0006:3:20:1714:16281#0/1,SINATRA_0006:2:76:17683:16082#0/1,SINATRA_0006:2:72:15990:12773#0/1,      na      na
# 163745 (32 fields)
# 2) a minimum number of spanning reads for a junction to be reported in the matrix (e.g. 2) (default 1)
# 3) a minimum number of discordant pe reads for a junction to be reported in the matrix in the same expt where 2) is achieved (e.g. 1) (default 1)
# and this script will produce as output a matrix like this one where the script has been run
#############################################################################################
# junc_id       beg     end     samechrstr      okgxorder       dist    ss1     ss2     gnlist1 gnlist2 gnname1 gnname2 gnbt1   gnbt2   LID16627        LID16628        LID16629  LID16630 LID16631        LID16632        LID16633        LID16634        LID16635        LID16636        LID44497        LID44498        LID44499        LID44594        LID45016  LID45017 LID46598        LID46599        LID8461 LID8462 LID8463 LID8464 LID8686 LID8687 LID8692 LID8701 LID8710 LID8711 LID8963 LID8964 LID8965 LID8966 LID8967 LID8968 LID8969 LID8970
# chr7_32768131_-:chr7_32766082_- 32768206        32765984        1 1 2049 GT AG ENSG00000237004.2 ENSG00000229358.2 ZNRF2P1 AC018633.4 pseudogene pseudogene     0       0       0 00       0       0       5       0       6       0       0       0       2       0       1       12      12      4       4       3       0       1       11      0       2       0 60       5       0       0       0       0       0       5
# 979 (50 fields)


# usage
#######
# make_chimeric_junction_matrix5.sh file_with_junc_each_exp.txt minspan mindiscpe

# example
#########

# Notes
#######
# - Made for using on a 64 bit linux architecture
# - uses awk scripts


# In case the user does not provide any input file
###################################################
if [ ! -n "$1" ]
then
    echo "" >&2
    echo "Usage: make_chimeric_junction_matrix5.sh file_with_junc_each_exp.txt [minspan] [minpe]" >&2
    echo "" >&2
    echo Example: make_chimeric_junction_matrix5.sh lid_junctionfile.txt 2 >&2
    echo "" >&2
    echo Takes as input:
    echo 1\) a 2 column file containing the candidate chimeric junction files associated to each expt or tech rep >&2
    echo where column 1 is a unique identifier for the expt or tech rep, and columns 2 is the absolute path to the >&2
    echo candidate chimeric junction file \(or proper concatenation of filtered and final junction files\) with the number >&2
    echo of spanning and discordant pe reads and all info from ChimPipe version 0.9.4 at least >&2
    echo 2\) an optional minimum number of spanning reads for a junction to be reported in the matrix \(default 1\) >&2
    echo 3\) an optional minimum number of discordant pe reads in the same expt as 2\) was achieved for a junction >&2
    echo to be reported in the matrix \(default 0\) >&2
    echo and produces a matrix containing the chimeric junctions obtained across all experiments with the min nb >&2
    echo of split and the min nb of discordant pe reads \in at least one experiment \(based on the candidate files\) >&2
    echo "" >&2
    exit 1
else
input=$1
fi

# In case the user does not provide any min stag or min pe we provide default values of 1 and 0 respectively
############################################################################################################
if [ ! -n "$2" ]
then
    minspan=1
    minpe=0
else
    minspan=$2
    if [ ! -n "$3" ]
    then
	minpe=0
    else
	minpe=$3
    fi
fi

# Path to the scripts/programs
##############################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

# Scripts/programs used by this script
######################################
CUTGFF=$rootDir/cutgff.awk
BEGEND=$rootDir/add_maxbegend_to_chimjunc2.awk
GFF2GFF=$rootDir/gff2gff.awk
OVERLAP=$rootDir/overlap

echo "minspan is $minspan and minpe is $minpe" >&2

# 1) In each experiment, select the junctions supported by the min number of spanning reads and min number of discordant paired end reads
##########################################################################################################################################
# input is like this
# juncCoord	type	filtered	reason	nbTotal(spanning+consistent)	nbSpanningReads	nbStaggered	percStaggered	nbMulti	percMulti	nbConsistentPE	nbInconsistentPE	percInconsistentPE	overlapA	overlapB	distExonBoundaryA	distExonBoundaryB	blastAlignLen	blastAlignSim	donorSS	acceptorSS	beg	end	sameChrStr	okGxOrder	dist	gnIdsA	gnIdsB	gnNamesA	gnNamesB	gnTypesA	gnTypesB	juncSpanningReadsIds	consistentPEIds	inconsistentPEIds
# JH118593.1_108648_-:JH118593.1_86156_-	readthrough	1	consistentPE,spanningReads,totalSupport,	2	1	1	100	0	0	1	0	0	100	22.1311	0	95	na	na	GT	AG	108675	86035	1	1	22492	ENSSSCG00000023895	ENSSSCG00000027705	na	na	na	na	ST-J00115:48:H7WJVBBXX:1:1225:7243:34723/2,	ST-J00115:48:H7WJVBBXX:1:1208:30350:38328,	na
# 12428 (35 fields)

# it used to be like this before
# juncCoord     totalNbPE       nbSpanningPE    nbStag  percStag        nbMulti percMulti       nbDiscordantPE  nbInconsistentPE        percInconsistentPE      overlapA        overlapB   distExonBoundaryA       distExonBoundaryB       maxBLastAlignLen        blastAlignSim   donorSS acceptorSS      beg     end     sameChrStr      okGxOrder       dist    gnIdA       gnIdB  gnNameA gnNameB gnTypeA gnTypeB juncSpanningReadIds     supportingPairsIds      inconsistentPairsIds
# chrX_100646810_+:chr18_13919745_-       28      28      12      42.8571 28      100     0       0       na      100     100     0       177     321     96.26   GT      AG      100646777  13919690        0       na      na      ENSG00000241343.3,ENSG00000257529.1     ENSG00000243779.1       RPL36A,RP1-164F3.9      RP11-681N23.1   protein_coding,protein_coding      pseudogene      SINATRA_0006:3:109:18510:14395#0/1,SINATRA_0006:1:95:16661:9220#0/1,SINATRA_0006:1:95:19257:7963#0/1,SINATRA_0006:1:78:8979:7945#0/1,SINATRA_0006:2:103:7300:8834#0/1,SINATRA_0006:1:95:6274:6674#0/1,SINATRA_0006:1:95:12119:12408#0/1,SINATRA_0006:1:65:19204:4609#0/1,SINATRA_0006:2:26:15962:4612#0/1,SINATRA_0006:1:78:9741:2208#0/2,SINATRA_0006:3:55:3004:18029#0/1,SINATRA_0006:2:28:17045:15305#0/1,SINATRA_0006:3:28:12150:19552#0/1,SINATRA_0006:1:57:14946:18483#0/1,SINATRA_0006:3:40:17960:10169#0/1,SINATRA_0006:2:78:4492:5160#0/1,SINATRA_0006:1:95:12206:3223#0/1,SINATRA_0006:1:87:6224:19402#0/2,SINATRA_0006:3:117:19626:18481#0/1,SINATRA_0006:2:13:1765:2306#0/1,SINATRA_0006:3:1:10757:14190#0/1,SINATRA_0006:2:69:18762:13853#0/2,SINATRA_0006:1:58:4291:4933#0/1,SINATRA_0006:1:36:4284:12343#0/2,SINATRA_0006:1:95:7884:9420#0/1,SINATRA_0006:3:20:1714:16281#0/1,SINATRA_0006:2:76:17683:16082#0/1,SINATRA_0006:2:72:15990:12773#0/1,      na      na
# 163745 (32 fields)


# we want to have this information in the output file
#####################################################
# junc_id\tbeg\tend\tsamechrstr\tokgxorder\tdist\ttype\tss1\tss2\tgnlist1\tgnlist2\tgnname1\tgnname2\tgnbt1\tgnbt2
echo I am selecting the junctions supported by at least $minspan spanning reads and $minpe paired end reads \in each experiment or technical replicate >&2
cat $input | while read lid f
do
awk -v minspan=$minspan -v minpe=$minpe 'BEGIN{OFS="\t"}(($1!~/unc/)&&($6>=minspan)&&($11>=minpe)&&(($18=="na")||($18<30)||($19<80))){print $1, $22, $23, $24, $25, $26, $2, $20, $21, $27, $28, $29, $30, $31, $32}' $f > ${f%.txt}\_$minspan\spanning_$minpe\discordantpe_withmaxbegandend.tsv
done 
# chr10_70700950_+:chr10_70741293_+     70700876        70741315        1       1       40343   readthrough GT      AG      ENSG00000107625.6       ENSG00000165732.7       DDX50   DDX21   protein_coding     protein_coding
# 2194 (15 fields)

# 2) Make the junctions obtained in all experiments with their most extreme beg and end across all experiments
##############################################################################################################
echo I am determining the set of junctions found \in all experiments with their most extreme beg and end \in all expts or tech replicates >&2
echo as well as with all other information from the input file that is not expression >&2
cat $input | while read lid f
do
cat ${f%.txt}\_$minspan\spanning_$minpe\discordantpe_withmaxbegandend.tsv
done | awk -v fld1=2 -v fld2=3 -f $BEGEND > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo.tsv
# chr1_174128175_+:chr1_174188263_+     174128113       174188337       1 1 60088 GT AG ENSG00000227373.1 ENSG00000152061.16 RP11-160H22.5 RABGAP1L lincRNA protein_coding 
# 70 (14 fields)

# 3) Add the number of spanning reads and the number of discordant pe reads for each experiment together with a header to the whole matrix
##########################################################################################################################################
echo I am adding the number of spanning reads and the number of discordant pe reads for each experiment or technical replicate >&2
awk 'BEGIN{print "junc_id\tbeg\tend\tsamechrstr\tokgxorder\tdist\ttype\tss1\tss2\tgnlist1\tgnlist2\tgnname1\tgnname2\tgnbt1\tgnbt2"} {print $0}' allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo.tsv > allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo_withexpr_eachexp.tsv
cat $input | while read lid f
do
awk -v lid=$lid -v fileRef=$f 'BEGIN{OFS="\t"; while (getline < fileRef >0){nb[$1]=$6":"$11}} NR==1{print $0, lid}NR>=2{print $0, (nb[$1]!="" ? nb[$1] : 0":"0)}' allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo_withexpr_eachexp.tsv > tmp
mv tmp allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo_withexpr_eachexp.tsv
done
# junc_id       beg     end     samechrstr      okgxorder       dist    ss1     ss2     gnlist1 gnlist2 gnname1 gnname2 gnbt1   gnbt2   LID16627        LID16628        LID16629  LID16630 LID16631        LID16632        LID16633        LID16634        LID16635        LID16636        LID44497        LID44498        LID44499        LID44594        LID45016  LID45017 LID46598        LID46599        LID8461 LID8462 LID8463 LID8464 LID8686 LID8687 LID8692 LID8701 LID8710 LID8711 LID8963 LID8964 LID8965 LID8966 LID8967 LID8968 LID8969 LID8970
# 979 (50 fields)  **** as x:y x represents the number of split (or spanning) reads and y represents the number of pe reads

# 4) Clean
##########
echo I am cleaning >&2
cat $input | while read lid f
do
#rm ${f%.txt}\_$minspan\spanning_$minpe\discordantpe_withmaxbegandend.tsv
done
# rm allexp_distinct_junctions_reliable_ineachexp_withmaxbegandend_allinfo.tsv
echo I am done >&2
