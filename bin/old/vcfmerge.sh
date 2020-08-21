#!/bin/bash

# vcfmerge.sh
#############
# This is to merge vcf.gz files including the type of SV (DEL, INV, ...) but from different individuals
# into one on the basis of reciprocal overlap of 66% of the elements that are merged and while keeping
# the genotype information of all individuals in the final file
# Its basic steps are
# - filter the input vcf.gz files to only have the event type we want
# - concatenate them into a non gzipped vcf file 
# - use intersectBed on this file with itself to have the events that are recip overlapping by 66%
# - make the nodes and edges for the connected component python script from CM and run it to get connected components as a tsv file
# - convert this tsv file into a vcf file
# This script takes as input:
#############################
# - a 2 column tsv file without header with the vcf.gz files of each individual in the order we want to see
#   them in the final output. The first column is the individual identifier and the second one the vcf file
#   these files could contain more than the type of event we are considering but should have the SVTYPE field
#   correctly informed 
# - the type of SV event we are considering (DEL, INV, ...)
# - an id for the merge file that will be output (not including the SV type)
# This script provides as output:
#################################
# - a merged vcf.gz file with the PASS elements of the sv type of interest that reciprocally overlap by 66%
#   and with the genotypes of each individual in the order given in the input tsv file
# - a merged vcf file with the events that have conflicting genotypes for a given individual


# example
# cd /work2/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/test
# input=/work2/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/svim.minimap.DEL.allrun.indiv.vcffiles.tsv
# st=DEL
# outid=svim.minimap.DEL.allrun.trio1
# pgm=/work2/project/fragencode/tools/multi/Scripts/Bash/vcfmerge.sh
# time $pgm $input $st $outid > vcfmerge.out 2> vcfmerge.err
# real	0m6.185s

# input
# offspring	/work2/project/seqoccin/workspace/vcf/persvtype/minimap.svim/trio1_offspring_allruns.minimap.svim.DEL.vcf.gz
# father	/work2/project/seqoccin/workspace/vcf/persvtype/minimap.svim/trio1_father_allruns.minimap.svim.DEL.vcf.gz
# mother	/work2/project/seqoccin/workspace/vcf/persvtype/minimap.svim/trio1_mother_allruns.minimap.svim.DEL.vcf.gz

# output body
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	offspring	father	mother
# 1	165867	indoffspring:idsvim.DEL.74,indmother:idsvim.DEL.200,	N	<DEL>	.	PASS	END=166361;SVLEN=-495;SVTYPE=DEL;POSLIST=165869,165867;ENDLIST=166361,166361;SVLENLIST=-492,-494;GTLIST=0/1,0/1;QUALLIST=7,23;DETLIST=1,0,1;CONFLLIST=0,0,0,	GT	0/1	0/0	0/1

# Check 3 inputs are provided
#############################
if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] 
then
    echo "" >&2
    echo Usage: vcfmerge.sh ind.vcffile.tsv SVtype outputID >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a 2 column tsv file without header with the vcf.gz files of each individual in the order we want to see" >&2
    echo "  them in the final output. The first column is the individual identifier and the second one the vcf file" >&2
    echo "  these files could contain more than the type of event we are considering but should have the SVTYPE field" >&2
    echo "  correctly informed" >&2
    echo "- the type of SV event we are considering (DEL, INV, ...)" >&2
    echo "- an id for the merge file that will be output (not including the SV type)" >&2
    echo "produces as output:" >&2
    echo "- a merged vcf.gz file with the PASS elements of the sv type of interest that reciprocally overlap by 66%" >&2
    echo "  and with the genotypes of each individual in the order given in the input tsv file" >&2
    echo "- a merged vcf file with the events that have conflicting genotypes for a given individual" >&2
    echo "Notes: needs to have awk, bcftools and bedtools in the path or the environment where it is run" >&2
    echo "" >&2
    exit 1
fi


# Check the input file exists
#############################
if [ ! -f "$1" ]
then
    echo "Your input tsv file does not exist" >&2
    exit 1
fi


# Variable assignment
#####################
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
input=$1
st=$2
outid=$3


# Programs
##########
CONVERT=$rootDir/../Awk/conncompSV2vcf_better.awk
CC=$rootDir/../Python/connected_components.py



# Start of the script
#####################
# 1. Filter the input vcf.gz files to only have the kind of events we are interested in and the events that are PASS
####################################################################################################################
echo "I am filtering the input vcf files to only have the kind of events we are interested in and that are PASS" >&2
cat $input | while read ind vcf
do
    b=`basename $vcf`
    b2=${b%.vcf.gz}
    bcftools view $vcf | awk -v st=$st 'BEGIN{OFS="\t"} $1~/#/{print} $1!~/#/&&$7=="PASS"{found=0; split($8,a,";"); k=1; while(a[k]!=""){split(a[k],b,"="); if(b[1]=="SVTYPE"){found=b[2]} k++} if(found!=0){print}}' > $b2.$st.pass.vcf
done
echo "done" >&2
# checked ok files


# 2. Concatenate all those vcf files into one but remembering which individual it was coming from and the id of the event (3rd field)
#####################################################################################################################################
echo "I am concatenating all the resulting vcf files into one and remembering the individual and the id of the events" >&2
fst=`head -1 $input | awk -v st=$st '{n=split($2,a,"/"); split(a[n],b,".vcf.gz"); print b[1]"."st".pass.vcf"}'`
bcftools view -h $fst > $outid.concat.vcf
cat $input | while read ind vcf
do
    bcftools view -H $vcf | awk -v ind=$ind 'BEGIN{OFS="\t"} {$3="ind"ind":id"$3; print}' >> $outid.concat.vcf
done
echo "done" >&2
# checked ok header and nb rows in body


# 3. Make an intersection of this file against itself asking for a reciprocal overlap of 66%
############################################################################################
echo "I am making an intersection of this file against itself asking for a reciprocal overlap of 66% (1 tsv file)" >&2
intersectBed -a $outid.concat.vcf -b $outid.concat.vcf -f 0.66 -r -wao > $outid.concat.vs.itself.tsv
echo "done" >&2
# 1	165869	indoffspring:idsvim.DEL.74	N	<DEL>	7	PASS	SVTYPE=DEL;END=166361;SVLEN=-492;SUPPORT=5;STD_SPAN=3.49;STD_POS=1.43	GT:DP:AD	0/1:10:5,5	1	165869	indoffspring:idsvim.DEL.74	N	<DEL>	7	PASS	SVTYPE=DEL;END=166361;SVLEN=-492;SUPPORT=5;STD_SPAN=3.49;STD_POS=1.43	GT:DP:AD	0/1:10:5,5	493
# 100091 (21 fields)


# 4. Make the node and edge files out of this intersection to be able to run the connected component script (remove redund)
##########################################################################################################################
echo "I am making the node and edge files out of this intersection to be able to run the connected component script (2 tsv files)" >&2
awk '{seen[$3]++; if(seen[$3]==1){print $3}}' $outid.concat.vs.itself.tsv > $outid.nodes.forcc.txt
awk '$3!=$13{seen[$3,$13]++; if(seen[$3,$13]==1){print $3"\t"$13}}' $outid.concat.vs.itself.tsv > $outid.edges.forcc.tsv
echo "done" >&2
# indoffspring:idsvim.DEL.74
# 46153 (1 fields)
# indoffspring:idsvim.DEL.74	indmother:idsvim.DEL.200
# 53938 (2 fields)


# 5. Run the connected component script
#######################################
echo "I am running the connected component (cc) script (1 tsv file)" >&2
python $CC $outid.nodes.forcc.txt $outid.edges.forcc.tsv > $outid.connected_components.tsv
echo "done" >&2
# indoffspring:idsvim.DEL.74	indmother:idsvim.DEL.200
# indoffspring:idsvim.DEL.92	indfather:idsvim.DEL.174	indmother:idsvim.DEL.261
# 12916 (1 fields)
# 6222 (2 fields)
# 6773 (3 fields)
# 71 (4 fields)
# 26 (5 fields)
# 10 (6 fields)

# 6. Convert the cc tsv file into a cc vcf file with the genotype info
#######################################################################
echo "I am converting the cc tsv file into a cc vcf file with the genotype info" >&2
str=`cat $input | awk '{s=(s)($1)(",")} END{print substr(s,1,length(s)-1)}'`
awk -v st=$st -v runid=$outid -v indlist=$str -v fileRef=$outid.concat.vcf -f $CONVERT $outid.connected_components.tsv > $outid.cc.vcf
bcftools view $outid.cc.vcf -O z -o $outid.cc.vcf.gz
echo "done" >&2
# 25928 rows in the body of the main input file
# 90 rows in the body of the conflict file
# CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	offspring	father	mother
# 1	165867	indoffspring:idsvim.DEL.74,indmother:idsvim.DEL.200,	N	<DEL>	.	PASS	END=166361;SVLEN=-495;SVTYPE=DEL;POSLIST=165869,165867;ENDLIST=166361,166361;SVLENLIST=-492,-494;GTLIST=0/1,0/1;QUALLIST=7,23;DETLIST=1,0,1;CONFLLIST=0,0,0,	GT	0/1	0/0	0/1


# 7. Clean
###########
echo "I am cleaning" >&2
# remove the intermediate sv type files
cat $input | while read ind vcf
do
    b=`basename $vcf`
    b2=${b%.vcf.gz}
    rm $b2.$st.pass.vcf
done
rm $outid.cc.vcf
rm $outid.concat.vcf
rm $outid.nodes.forcc.txt
rm $outid.edges.forcc.tsv
echo "done" >&2




