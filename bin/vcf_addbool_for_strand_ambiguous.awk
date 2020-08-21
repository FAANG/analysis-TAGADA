# vcf_addbool_for_strand_ambiguous.awk
# takes as input a vcf file and adds in the INFO 8th field a subfield with strnonambig=<bool> where
# bool is a boolean saying whether the alleles are not strand ambiguous
# note that this new subfield description is for the moment not added in the header

# example
# cd ~/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
# pgm=~/fragencode/tools/multi/Scripts/Awk/vcf_addbool_for_strand_ambiguous.awk
# module load bioinfo/bcftools-1.9
# time bcftools view 22.basic.INFOin0.5-1.ACin60-5086.vcf.gz | awk -f $pgm | bcftools view -O z -o 22.basic.INFOin0.5-1.ACin60-5086.allelenonstrambig.vcf.gz
# real	0m2.837s

# input
# chr22	16050435	.	T	C	.	PASS	RefPanelAF=0.000323375;AN=6046;AC=8;INFO=0.593662;INFOin0.5-1=1;ACin60-5086=0

# output
# chr22	16050435	.	T	C	.	PASS	RefPanelAF=0.000323375;AN=6046;AC=8;INFO=0.593662;INFOin0.5-1=1;ACin60-5086=0;notstrambig=1


BEGIN{
    OFS="\t";
}

$1~/\#/{
    print;
}

$1!~/\#/{
    ok=1;
    if(($4=="A"&&$5=="T")||($4=="T"&&$5=="A")||($4=="C"&&$5=="G")||($4=="G"&&$5=="C"))
    {
	ok=0;
    }
    $8=$8";notstrambig="ok;
    print $0;
}
