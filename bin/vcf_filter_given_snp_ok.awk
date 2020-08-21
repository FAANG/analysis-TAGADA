# vcf_filter_given_snp_ok.awk
# takes as input a vcf file and a fileRef of chr:pos:ref:alt for snps that we should keep, and filters the vcf file accordingly
# and also puts in 3rd field as unique identifier the chr:pos after having removed the non unique ones in terms of chrom:pos

# example
# datadir=/work/project/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
# cd $datadir
# pgm=~/fragencode/tools/multi/Scripts/Awk/vcf_filter_given_snp_ok.awk
# module load bioinfo/bcftools-1.9
# time bcftools view 22.vcf.gz | awk -v fileRef=22.snp.3filters.ok.txt -f $pgm | bcftools view -Oz -o 22.filtered.vcf.gz

# input 22.snp.3filters.ok.txt
# chr22:16050822:G:A
# chr22:16053730:C:A
# 86245 (1 fields)

# output
# 22	16050822	22:16050822	G	A	.	PASS	RefPanelAF=0.137496;AN=6046;AC=685;INFO=0.687915	GT:ADS:DS:GP	0|0:0,0:0:1,0,0	1|0:0.65,0.15:0.8:0.2975,0.605,0.0975	0|0:0.05,0.05:0.1:0.9025,0.095,0.0025	1|0:1,0:1:0,1,0	0|0:0.05,0:0.05:0.95,0.05,0  ...

BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	split($1,a,":");
	nb[a[1]":"a[2]]++;
	ok[$1]=1;
    }
}

($1~/^#/)||(ok[$1":"$2":"$4":"$5]==1){
    if(ok[$1":"$2":"$4":"$5]==1)
    {
	if(nb[$1":"$2]==1)
	{
	    $3=$1":"$2;
	    print;
	}
    }
    else
    {
	if($1~/^#/)
	{
	    print;
	}
    }
}
