# vcf_addbool_for_INFO_subfield.awk
# takes as input a vcf file, a subfield called f of the 8th field and a min and max thresholds m and M
# and adds in the 8th field an attribute=value pair where attribute is <f>in<t1>-<t2>
# and value is a boolean saying whether the value in subfield f in the 8th field was in the range m-M
# note that this new subfield description is for the moment not added in the header

# example
# cd ~/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
# pgm=~/fragencode/tools/multi/Scripts/Awk/vcf_addbool_for_INFO_subfield.awk
# module load bioinfo/bcftools-1.9
# time bcftools view 22.basic.vcf.gz | awk -v f=INFO -v m=0.5 -v M=1 -f $pgm | bcftools view -O z -o 22.basic.INFOin0.5-1.vcf.gz

# input
# chr22	16050435	.	T	C	.	PASS	RefPanelAF=0.000323375;AN=6046;AC=8;INFO=0.593662

# output
# chr22	16050435	.	T	C	.	PASS	RefPanelAF=0.000323375;AN=6046;AC=8;INFO=0.593662;INFOin0.5-1=1

BEGIN{
    OFS="\t";
}

$1~/\#/{
    print;
}

$1!~/\#/{
    ok=0;
    found=0;
    split($8,a,";");
    k=1;
    while(found==0&&a[k]!="")
    {
	split(a[k],a1,"=");
	if(a1[1]==f)
	{
	    found=1;
	    if(a1[2]>=m&&a1[2]<=M)
	    {
		ok=1;
	    }
	}
	k++;
    }
    $8=$8";"f"in"m"-"M"="ok;
    print $0;
}
