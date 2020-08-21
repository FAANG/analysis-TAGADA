
# TrtoGn_RPKM.sh
# takes as input a gtf file of transcripts produced by the flux capacitor
# and makes as output a gene gff file with an rpkm value associated to each gene
# and computed as the sum of the rpkm of the transcripts belonging to this gene
# as well as the number of reads supporting the gene computed as the sum of the
# reads supporting each transcript of this gene.
# Be careful: makes the following assumptions:
##############################################
# 1. the annotation file has gene_id gene as first key/value pair, and transcript_id transcript as second key/value pair
# 2. the flux output file has transcript_id transcript as first key/value pair

# usage
# ~/bin/TrtoGn_RPKM.sh annot.gtf tr.gtf 

# input
# 1	HAVANA	transcript	12010	13670	.	+	.	transcript_id "ENST00000450305.2"; locus_id "1:11869-14412W"; gene_id "ENSG00000223972.4"; reads 0.000000; length 632; RPKM 0.000000

# output
# 16	HAVANA	gene	68279207	68294961	.	+	.	gene_id "ENSG00000103066.8"; transcript_ids "ENST00000413021.2,ENST00000565744.1,ENST00000219345.5,ENST00000566978.1,ENST00000564827.2,ENST00000566188.1,ENST00000444212.2,ENST00000562966.1,ENST00000568082.1,ENST00000568599.1,ENST00000562449.2,ENST00000565460.1,"; RPKM 0.494646; reads 78;



# In case the user does not provide any annotation file or any flux output file
###############################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] 
then
    echo "" >&2
    echo "Usage:   TrtoGn_RPKM.sh annot.gtf tr.gtf" >&2
    echo "         where annot.gtf is an annotation gtf file and where tr.gtf is a flux output file with rpkm value;" >&2
    echo "         will produce a gene gff file with for each gene an rpkm value and a number of read value, corresponding" >&2
    echo "         to the sum of the rpkm and of the number of reads of the transcripts belonging to this gene." >&2
    echo "" >&2
    exit 1
fi

path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path

annot=$1
tr=$2

annotbase=`basename $annot`
trbase=`basename $tr`

# Programs
##########
GFF2GFF=$rootDir/gff2gff.awk

echo "I am making the file of genes with associated transcripts from the annotation" >&2
awk '$3=="transcript"{split($12,a,"\""); trlist[$10]=(trlist[$10])(a[2])(",");}$3=="gene"{line[$10]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9" "$10}END{for(g in line){print line[g], "transcript_ids", "\""trlist[g]"\"\;"}}' $annot | awk -f $GFF2GFF > ${annotbase%.gtf}.gene.withtrlist.gff

echo "I am making the gene file with rpkm and number of reads" >&2
awk -v fileRef=$tr 'BEGIN{while (getline < fileRef >0){k=9; while(k<=(NF-1)){split($10,a,"\""); if($k=="RPKM"){split($(k+1),b,";"); rpkm[a[2]]=b[1];} if($k=="reads"){split($(k+1),b,";"); reads[a[2]]=b[1];} k+=2}}} {split($12,a,"\""); split(a[2],b,","); s1=0; k=1; while(b[k]!=""){s1+=rpkm[b[k]]; k++} s2=0; k=1; while(b[k]!=""){s2+=reads[b[k]]; k++} print $0, "RPKM", s1"\;", "reads", s2"\;"}' ${annotbase%.gtf}.gene.withtrlist.gff | awk -f $GFF2GFF > ${trbase%.gtf}\_gene_with_rpkm.gff

echo "I am removing unuseful files" >&2
rm ${annotbase%.gtf}.gene.withtrlist.gff

echo "I am done" >&2


