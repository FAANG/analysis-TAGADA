# vcf_filterinfo_2_snpok.awk
# from a vcf file with no header (snp file in tsv format)
# and a comma separated list of 8th field subfield attributes corresponding to boolean values for different filters (0 not ok snp, 1 ok snp)
# make a txt file of chrom:pos:ref:alt of the snps to keep after those filters

# example
# datadir=/work/project/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
# cd $datadir
# pgm=~/fragencode/tools/multi/Scripts/Awk/vcf_filterinfo_2_snpok.awk
# module load bioinfo/bcftools-1.9
# time bcftools view -H 22.basic.INFOin0.5-1.ACin60-5086.notstrambig.vcf.gz | awk -v fldlist=INFOin0.5-1,ACin60-5086,notstrambig -f $pgm | awk '{gsub(/chr/,"",$1); print $1}' > 22.snp.3filters.ok.txt
# *** real	8m57.927s

# input
# bcftools view -H 22.basic.INFOin0.5-1.ACin60-5086.notstrambig.vcf.gz
# chr22	16050435	.	T	C	.	PASS	RefPanelAF=0.000323375;AN=6046;AC=8;INFO=0.593662;INFOin0.5-1=1;ACin60-5086=0;notstrambig=1

# output
# 22:16050822:G:A
# 22:16053730:C:A
# 86245 (1 fields)

BEGIN{
    OFS="\t";
    n=split(fldlist,a,",");
    for(i=1; i<=n; i++)
    {
	fld[i]=a[i];
	okfld[a[i]]=1;
	idx[a[i]]=i;
    }
}

($1!~/^#/){
    s=0;
    for(i=1; i<=n; i++)
    {
	ok[i]=0;
    }
    split($8,a,";");
    k=1;
    while(a[k]!="")
    {
	split(a[k],a1,"=");
	if(okfld[a1[1]]==1)
	{
	    ok[idx[a1[1]]]=a1[2];
	}
	k++;
    }
    for(i=1; i<=n; i++)
    {
	s+=ok[i];
    }
    if(s==n)   # in case the snp passed all filters then we put it in the output file
    {
	print $1":"$2":"$4":"$5;
    }
} 
