#!/bin/bash
myname=`basename $0`
FEATURES=( gene exon transcript CDS UTR UTR5 UTR3 utr utr5 utr3 five_prime_utr three_prime_utr intron noncoding non_coding non-coding miRNA snoRNA ncRNA lncRNA tRNA mRNA rRNA protein);
usage() { echo "#############################################################################################################################
                $myname: extract various features from a gtf annotation file
                                    and generate distinct bed files

Usage: $myname <annotation.gtf>

Input: 
<annotation.gtf> : annotation file in gtf format
                   compatible with gff and gff3 format for all features but introns (\"transcript_id\" is required)

Output: 
For each of the following features: ${FEATURES[@]} 
generates two files if the feature is described in the gtf input (in the same directory than the input):

<annotation>.<feature>.distinct.bed :  contains a non-redondant and oriented (stranded) set of intervals from the input.
                                       intervals are sorted (1st field on sort -V) and distinct but can overlap each other.
<annotation>.<feature>.positions.bed : contains a non-overlapping and non-oriented (unstranded) set of sorted intervals
                                       with all the positions that are annotated as feature at least once in the gtf input.
                                       intervals can overlap across feature files (e.g. exon and intron)

Note: introns are generated if exons and transcript_id are found and TSS/TTS are generated if transcripts or mRNA are found.
Warning: needs bedtools, fastalength and gtf2bed - does not operate on .gz files
Warning: needs in the same directory either a tab-separated file called 'chrom.len' with name and size of each chrom, or 
         otherwise a fasta file of the genomic sequence with the same basename than the gtf and the .fa extension 
#############################################################################################################################" 1>&2; exit 1;}

PROMOTERFLANK=( 1000 5000 )
# warning: flanks needs to be a multiple of 1000 (to convert into Kb in filename)

echo "# $myname - `date`"

gtf=$1;

if [[ ! -r "$gtf" ]] ; then
    echo "# ERROR: could not read input annotation \"$gtf\"" ;
    usage ;
fi

if [ ! "`sortBed -h 2>&1 | grep -im1 version`" ] || [ ! "`mergeBed -h 2>&1 | grep -im1 version`" ] ; then
    echo "# ERROR: needs bedtools sort and merge (sortBed & mergeBed) to run" ; 
    usage ;
fi

if [ ! "`fastalength 2>&1 | grep -im1 version`" ] ; then
    echo "# ERROR: needs fastalength to run" ; 
    usage ;
fi

if [ ! "`gtf2bed 2>&1 | grep -im1 usage`" ] ; then
    echo "# ERROR: needs gtf2bed to run" ; 
    usage ;
fi

echo "# Processing $gtf"

# chrom size
chromlen=`dirname $gtf`/chrom.len
fasta=${gtf%.gtf}.fa
if [[ ! -r "$chromlen" ]] ;
    then
    echo "# computing chromosome length..."
    if [[ ! -r "$fasta" ]] ;
    then
	echo "# ERROR: could not find chrom len file $chromlen nor the required fasta file $fasta to generate it"; exit 2;
    fi
    fastalength $fasta | awk '{print $2"\t"$1}' | sort -V > $chromlen ;
fi

# looking for features in the annotation file:
for f in "${FEATURES[@]}" ; do
# bed=${gtf%.gtf}.`echo $f | awk '{print tolower($1)}'`.bed ;
 pos=${gtf%.gtf}.`echo $f | awk '{print tolower($1)}'`.positions.bed ;
 dis=${gtf%.gtf}.`echo $f | awk '{print tolower($1)}'`.distinct.bed ;
 # clean out previous results
 rm -f $pos $dis ;
 # if the annotation feature is found at least once in the gtf file (first quick check with grep)
 if [ "`fgrep -wm1 $f $gtf`" ] ; then
     echo "# looking for $f features to extract intervals and regions..." ;
     cat $gtf | awk -v f=$f '!/^#/ && $3==f {print $1,$4-1,$5,".",".",$7}' OFS="\t" | uniq | sort -Vk1,1 -k2,2n -k3n | uniq > $dis &&
	 cut -f-3 $dis | sortBed | mergeBed | sort -Vk1,1 -k2,2n -k3,3n > $pos ;
     if [[ "`head $pos`" ]] ;
     then
	 echo "# ... done in:"; echo " $dis "; echo " $pos ";
     else
	 echo "# ... no $f found" && rm $dis $pos ; # remove files if empty (feature needs to be in the third field)
     fi
 fi ;
done

# introns
f=intron
pos=${gtf%.gtf}.intron.positions.bed ;
dis=${gtf%.gtf}.intron.distinct.bed ;
# needs exons and transcripts to be defined; also, not required if introns have already been extracted
if [[ -r "${gtf%.gtf}.exon.positions.bed" ]] && [[ "`fgrep -wm1 transcript_id $gtf`" ]] && [[ ! -r "$pos" ]] ; then
    echo "# generating $f intervals and regions..." ;
    idfield=`fgrep -m1 transcript_id $gtf | awk '{for (i=8;i<NF;i++){if ($i=="transcript_id"){print i+1;exit}}print 1}'`
    cat $gtf | awk '!/^#/ && $3=="exon"' | sort -k$idfield,$idfield -Vk1,1 -k4,4n -k5n | 
	awk -v idfield=$idfield '$1==chrom && $idfield==tid && $(idfield-1)=="transcript_id" && ($4-1)>prev {print $1,prev,$4-1,".",".",$7}{chrom=$1;tid=$idfield;prev=$5}' OFS="\t" | 
	uniq | sort -Vk1,1 -k2,2n -k3n | uniq > $dis &&
	cut -f-3 $dis | sortBed | mergeBed | sort -Vk1,1 -k2,2n -k3,3n > $pos ;
    if [[ "`head $pos`" ]] ; then
	echo "#... done in:"; echo " $dis "; echo " $pos ";
    else
	echo "#... FAILED ; skipped.";
	rm $dis $pos ; # remove files if empty (feature needs to be in the third field)
    fi ;
fi ;

# if genes and transcripts were not defined as a feature, try to get them from the gene_id and transcript_id field
dis=${gtf%.gtf}.transcript.distinct.bed;
pos=${gtf%.gtf}.transcript.pos.bed;
if [[ ! -r "$dis" ]] ; then
    echo "# no transcript feature => looking for transcript_id...";
    idfield=`fgrep -m1 transcript_id $gtf | awk '{for (i=8;i<NF;i++){if ($i=="transcript_id"){print i+1;exit}}print 1}'`
    if [[ $idfield -gt 1 ]] ; then
	echo "# ... found in field $idfield ...";
	cat $gtf | awk  -v idfield=$idfield '!/^#/ && $3=="exon" && $(idfield-1)=="transcript_id"' | sort -k$idfield,$idfield -Vk1,1 -k4,4n -k5n | 
	    awk -v idfield=$idfield '($idfield!=tid || $1!=chrom) {if(NR>1 && end>beg)print chrom,beg,end,".",".",strand;chrom=$1;beg=$4-1;end=$5;strand=$7}{chrom=$1;tid=$idfield;end=$5}END{if(end>beg)print chrom,beg,end,".",".",strand}' OFS="\t" |
	    uniq | sort -Vk1,1 -k2,2n -k3n | uniq > $dis &&
	    cut -f-3 $dis | sortBed | mergeBed | sort -Vk1,1 -k2,2n -k3,3n > $pos ;
	if [[ "`head $pos`" ]] ; then
	    echo "# ... done in:"; echo " $dis "; echo " $pos ";
	else
	    echo "#... no valid transcript found" && rm $dis $pos ; # remove files if empty
	fi
    else
	echo "# ... not found";
    fi
fi
# Same for gene: if not in 3rd field, try to find gene_id field
dis=${gtf%.gtf}.gene.distinct.bed;
pos=${gtf%.gtf}.gene.pos.bed;
if [[ ! -r "$dis" ]] ; then
    echo "# no gene feature => looking for gene_id...";
    idfield=`fgrep -m1 gene_id $gtf | awk '{for (i=8;i<NF;i++){if ($i=="gene_id"){print i+1;exit}}print 1}'`
    if [[ $idfield -gt 1 ]] ; then
	echo "# ... found in field $idfield ...";
	cat $gtf | awk  -v idfield=$idfield '!/^#/ && $3=="exon" && $(idfield-1)=="gene_id"' | sort -k$idfield,$idfield -Vk1,1 -k4,4n -k5n | 
	    awk -v idfield=$idfield '($idfield!=tid || $1!=chrom) {if(NR>1 && end>beg)print chrom,beg,end,".",".",strand;chrom=$1;beg=$4-1;end=$5;strand=$7}{chrom=$1;tid=$idfield;end=$5}END{if(end>beg)print chrom,beg,end,".",".",strand}' OFS="\t" |
	    uniq | sort -Vk1,1 -k2,2n -k3n | uniq > $dis &&
	    cut -f-3 $dis | sortBed | mergeBed | sort -Vk1,1 -k2,2n -k3,3n > $pos ;
	if [[ "`head $pos`" ]] ;
	then
	    echo "# ... done in:"; echo " $dis "; echo " $pos ";
	else
	    echo "# ... no valid gene found" && rm $dis $pos ; # remove files if empty
	fi
    else
	echo "# ... not found";
    fi
fi

# promoter= upstream region before each transcript
#t=${gtf%.gtf}.transcript.distinct.bed;
t1=${gtf%.gtf}.transcript.distinct.bed;
t2=${gtf%.gtf}.mrna.distinct.bed;
#if [ -r "$t" ] ;then
if [[ -r "$t1" ]] || [[ -r "$t2" ]] ; then
    for flankn in "${PROMOTERFLANK[@]}" ; do
	flank=${flankn%000}Kb
	dis=${gtf%.gtf}.prom$flank.distinct.bed ;
	pos=${gtf%.gtf}.prom$flank.positions.bed ;
	echo "# generating promoter intervals and regions ($flank bp before each transcript)...";
	(
	    if [ -r "$t1" ] ; then
		cat $t1 ;
	    fi
	    if [ -r "$t2" ] ; then
		cat $t2 ;
	    fi
	    #	cat $t | awk -v f=$chromlen -v flank=$flankn 'BEGIN{while(getline < f > 0)len[$1]=$2}{if($6=="+"){if($2>0)print $1,($2-flank>0?$2-flank:0),$2,".",".",$6}else{if($3<len[$1])print $1,$3,($3+flank<len[$1]?$3+flank:len[$1]),".",".",$6}}' OFS="\t" | uniq | sort -Vk1,1 -k2,2n -k3n | uniq > $dis &&
	) | awk -v f=$chromlen -v flank=$flankn 'BEGIN{while(getline < f > 0)len[$1]=$2}{if($6=="+"){if($2>0)print $1,($2-flank>0?$2-flank:0),$2,".",".",$6}else{if($3<len[$1])print $1,$3,($3+flank<len[$1]?$3+flank:len[$1]),".",".",$6}}' OFS="\t" | uniq | sort -Vk1,1 -k2,2n -k3n | uniq > $dis &&
	    cut -f-3 $dis | sortBed | mergeBed | sort -Vk1,1 -k2,2n -k3,3n > $pos ;
	echo "... done in:"; echo " $dis "; echo " $pos ";
    done
fi

echo "# $myname - `date` - End of processing"
