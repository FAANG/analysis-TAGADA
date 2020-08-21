
# TrtoEx_RPKM.sh
# takes as input a gtf file of transcripts produced by the flux capacitor
# and makes as output an exon gff file with an rpkm value associated to each exon
# and computed as the sum of the rpkm of the transcripts sharing this exon

# Be careful: makes the following assumptions:
##############################################
# 1. the annotation file has gene_id gene as first key/value pair, and transcript_id transcript as second key/value pair
# 2. the flux output file has transcript_id transcript as first key/value pair

# usage
# TrtoGn_RPKM.sh annot.gtf tr.gtf 

# input
# 1	HAVANA	transcript	12010	13670	.	+	.	transcript_id "ENST00000450305.2"; locus_id "1:11869-14412W"; gene_id "ENSG00000223972.4"; reads 0.000000; length 632; RPKM 0.000000

# output
# 16	Gencv10	exon	718560	718699	.	+	.	gene_ids ENSG00000140983.7, transcript_ids ENST00000566965.1, RPKM 2.16517; reads 251.931;



# In case the user does not provide any annotation file or any flux output file
###############################################################################
if [ ! -n "$1" ] || [ ! -n "$2" ] 
then
    echo "" >&2
    echo "Usage:   TrtoEx_RPKM.sh annot.gtf tr.gtf" >&2
    echo "         where annot.gtf is an annotation gtf file and where tr.gtf is a flux output file with rpkm value;" >&2
    echo "         will produce an exon gff file with for each exon an rpkm value, correspondin gto the sum of the rpkm" >&2
    echo "         of the transcripts sharing this exon." >&2
    echo "" >&2
    exit 1
fi

# Variables
path="`dirname \"$0\"`" # relative path
rootDir="`( cd \"$path\" && pwd )`" # absolute path
annot=$1
tr=$2
annotbase=`basename $annot`
trbase=`basename $tr`

# Programs
##########
GFF2GFF=$rootDir/../Awk/gff2gff.awk

echo "I am making the file of exons with associated transcripts from the annotation" >&2
awk '$3=="exon"{split($10,a,"\""); split($12,b,"\""); gnlist[$1"?"$4"?"$5"?"$7]=(gnlist[$1"?"$4"?"$5"?"$7])(a[2])(","); trlist[$1"?"$4"?"$5"?"$7]=(trlist[$1"?"$4"?"$5"?"$7])(b[2])(",");}END{for(e in gnlist){split(e,a,"?"); print a[1], "annot", "exon", a[2], a[3], ".", a[4], ".", "gene_ids", "\""gnlist[e]"\"\;", "transcript_ids", "\""trlist[e]"\"\;"}}' $annot | awk -f $GFF2GFF > ${annotbase%.gtf}.distinct.exon.withtrlist.gff 

echo "I am making the exon file with rpkm and number of reads" >&2
awk -v fileRef=$tr 'BEGIN{while (getline < fileRef >0){k=9; while(k<=(NF-1)){split($10,a,"\""); if($k=="RPKM"){split($(k+1),b,";"); rpkm[a[2]]=b[1];} k+=2}}} {split($12,a,"\""); split(a[2],b,","); s1=0; k=1; while(b[k]!=""){s1+=rpkm[b[k]]; k++} print $0, "RPKM", s1"\;"}' ${annotbase%.gtf}.distinct.exon.withtrlist.gff | awk -f $GFF2GFF > ${trbase%.gtf}\_distinct_exon_with_rpkm.gff

echo "I am removing unuseful files" >&2
rm ${annotbase%.gtf}.distinct.exon.withtrlist.gff

echo "I am done" >&2