# vcf_addinfo_make_tsv.awk
# from a vcf file without header, a tsv file and a comma separated list of field number in the tsv file as fldlist
# outputs as a tsv file the vcf file to which have been added separated by tabs the info from the fields
# in question from the input tsv file

# example
# datadir=/work/project/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
# pgm=~/fragencode/tools/multi/Scripts/Awk/vcf_addinfo_make_tsv.awk
# module load bioinfo/bcftools-1.9
# cd $datadir
# time bcftools view -H 22.basic.INFOin0.5-1.ACin60-5086.notstrambig.vcf.gz | awk -v fldlist=9,10,11,12,13,14 -v fileRef=allchr.basic.uniqlist.prom.other.gn.HCmerge.tsv -f $pgm > 22.basic.INFOin0.5-1.ACin60-5086.notstrambig.uniqlist.prom.other.gn.HCmerge.tsv
# real	0m1.488s

# inputs
# allchr.basic.uniqlist.prom.other.gn.HCmerge.tsv
# chr1	13380	rs571093408	C	G	.	PASS	RefPanelAF=7.69941e-05;AN=6046;AC=1;INFO=1	0	NA	0	NA	0	NA

# bcftools view 22.basic.INFOin0.5-1.ACin60-5086.notstrambig.vcf.gz
# chr22	16050435	.	T	C	.	PASS	RefPanelAF=0.000323375;AN=6046;AC=8;INFO=0.593662;INFOin0.5-1=1;ACin60-5086=0;notstrambig=1

# output


     
BEGIN{
    OFS="\t";
    n=split(fldlist,a,",");
    while (getline < fileRef >0)
    {
	for(i=1; i<=n-1; i++)
	{
	    info[$1":"$2":"$3":"$4":"$5]=(info[$1":"$2":"$3":"$4":"$5])($(a[i]))("\t");
	}
	info[$1":"$2":"$3":"$4":"$5]=(info[$1":"$2":"$3":"$4":"$5])($(a[i]));
    }
}

{
    print $0, info[$1":"$2":"$3":"$4":"$5];
}
       
