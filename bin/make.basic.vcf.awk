# make.basic.vcf.awk

# example
# cd ~/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
# module load bioinfo/bcftools-1.9
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/make.basic.vcf.awk
# time bcftools view 22.vcf.gz | awk -f $pgm | bcftools view -O z -o 22.basic.vcf.gz
# real	11m58.904s  *** 9M instead of 370M for 37 millions snps

# input
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SPD-999_3	SPD-997_3 ...
# too big to be looked at

# output
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
# 22	16050435	.	T	C	.	PASS	RefPanelAF=0.000323375;AN=6046;AC=8;INFO=0.593662


BEGIN{
    OFS="\t";
}

$1~/\#/&&$1!="#CHROM"{
    print;
}

$1!~/\#/||$1=="#CHROM"{
    print $1, $2, $3, $4, $5, $6, $7, $8, $9;
}
