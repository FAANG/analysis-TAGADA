# truncate_and_ens2ucsc_on_vcf.awk
# This is to truncate a vcf file to its first 9th fields (without genotypes)
# and optionally (if ucsc is set to 1) change chromosome names from ensembl
# convention (integers) to ucsc convention (chr<integer>)

# example
# cd ~/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
# pgm=~/fragencode/tools/multi/Scripts/Awk/ens2ucsc_vcf.awk
# time bcftools view 21.vcf.gz | awk -v ucsc=1 -f $pgm | bcftools view -O z -o 21.basic.vcf.gz  
# 12 minutes

BEGIN{OFS="\t"}

{
    if($1=="#CHROM")
    {
	print $1, $2, $3, $4, $5, $6, $7, $8, $9;
    }
    else
    {
	if($0!~/^#/)
	{
	    if(ucsc==1)
	    {
		print "chr"$1, $2, $3, $4, $5, $6, $7, $8, $9;
	    }
	    else
	    {
		print $1, $2, $3, $4, $5, $6, $7, $8, $9;
	    }
	}
	else
	{
	    if(match($0,/(##contig=<ID=)(.*)/,m))
	    {
		if(ucsc==1)
		{
		    print m[1]"chr"m[2];
		}
		else
		{
		    print m[1]m[2]
		}
	    }
	    else
	    {
		print $0;
	    }
	}
    }
}
