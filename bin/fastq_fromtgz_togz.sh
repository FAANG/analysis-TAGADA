
# ~sdjebali/bin/fastq_fromtgz_togz.sh

# usage
# fastq_fromtgz_togz.sh input.fastq.tgz output.fastq.gz 
# will produce a file called output.fastq.gz that has the different parts concatenated in the ls order 

if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: fastq_fromtgz_togz.sh input.fastq.tgz output.fastq.gz >&2
    echo "" >&2
    exit 1
fi

input=$1
base=`basename $input` 
output=$2

# extract the files from the archive
tar -xvzf $input > $base.list.txt

# create the fastq.gz
str=`sort -k1,1 $base.list.txt | awk '{s=(s)($1)(" ")}END{print s}'`
cat $str | gzip > $output

# clean
cat $base.list.txt | while read f
do
rm $f
done
rm $base.list.txt