# conncompSV2vcf.awk
# takes as input a tsv file output by connected_components.py on SVs of type t (DEL, INV or DUP)
# and that has resulted in merged SVs of type $t from several individuals when the SVs had a reciprocal overlap of 66% and
# as fileRef a bcftools merge file out of the events in a set of samples with genotypes in each sample but that can have redundancy due
# to the appearance of the same SV event in x individuals with different attributes
# and output a vcf file from this more laxist merge made by the connected component script and
# with genotypes in each of the individuals
# !!! this is a generalization of the script conncompSV2vcf_trio.awk where we only had 3 samples !!!
# !!! only considers genotypes for which we have reads so not gt:.:. for example where gt=0/0 !!!
# !!! if even with this rule a genotype is not found then homozygous reference !!!
# !!! svlen and svtype are needed in the resulting file and will therefore take the average of the svlen of the individual svs from which we take the genotypes in all the samples (indiv) !!!

# example
# dir=/work/project/seqoccin/svdetection/nanopore/bos_taurus
# cd $dir/merge/trio1.1run.eachindiv
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/conncompSV2vcf.awk 
# time for t in DEL INV DUP
# do
# awk -v t=$t -v fileRef=trio1.1run.eachindiv.$t.vcf -f $pgm trio1.1run.eachindiv.$t.connected.components.tsv > trio1.1run.eachindiv.$t.connected.components.vcf
# done
# real	0m0.074s


# input files
##############
# trio1.1run.eachindiv.DEL.vcf
# ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
# ##FORMAT=<ID=DR,Number=1,Type=Integer,Description="# high-quality reference reads">
# ##FORMAT=<ID=DV,Number=1,Type=Integer,Description="# high-quality variant reads">
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mapping/trio1_offspring_run1.minimap.forsniffles.bam	mapping/trio1_father_run2.minimap.forsniffles.bam	mapping/trio1_mother_run3.minimap.forsniffles.bam
# 1	247245	11;10_1;24_1	N	<DEL>	.	PASS	IMPRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=1;END=503308;SVTYPE=DEL;SVLEN=-256063;STRANDS2=5,4,5,4;RE=9;REF_strand=17,23;PRECISE;STD_quant_start=0.948683;STD_quant_stop=1.70294;Kurtosis_quant_start=6.74573;Kurtosis_quant_stop=-0.356994;SUPTYPE=SR;STRANDS=+-;AF=0.340909	GT:DR:DV	0/0:40:9	0/0:18:3	0/1:29:15

# trio1.1run.eachindiv.DEL.connected.components.tsv
# 10:10010297:DEL
# 10:100838356:DEL	10:100838361:DEL

# output file
##############
# ##fileformat=VCFv4.1
# ##FILTER=<ID=PASS,Description="All filters passed">
# 2222 (1 fields)
# 2 (2 fields)
# 3 (3 fields)
# 6 (4 fields)
# 1 (5 fields)
# 2 (6 fields)
# 4 (7 fields)
# 2 (8 fields)
# 4 (9 fields)
# 1 (10 fields)
# 1 (11 fields)
# 1116 (12 fields)
# 1 (14 fields)
# 1 (15 fields)
# 1 (16 fields)

BEGIN{
    if(t=="")
    {
	t="DEL";
    }
    OFS="\t";
    while (getline < fileRef >0)
    {
	if($1~/#/)
	{
	    print;
	}
	else
	{
	    # remember the length of each SV
	    split($8,a,";");
	    m=1;
	    fnd=0;
	    while(fnd==0&&a[m]!="")
	    {
		split(a[m],b,"=");
		if(b[1]=="SVLEN")
		{
		    len[$1":"$2":"t]=b[2];
		}
		m++;
	    }
	    # make the 9th field string of the final vcf file
	    s2=$9;
	    # for each SV of type t and each sample (indiv) remember the genotype of this SV in this sample (indiv) 
	    for(i=10; i<=NF; i++)
	    {
		split($i,a,":");
		if(a[2]!="."||a[3]!=".")
		{
		    gt[i-10+1,$1":"$2":"t]=$i;
		}
	    }
	    n=NF-10+1;
	}
    }
}

{
    s3="";
    l=0;
    # for each of the samples (indiv), finds the first SV of the connected component that has reads for the ref and alt alelles and remember it
    for(i=1; i<=n; i++)
    {
	k[i]=1;
	found[i]=0;
	sv[i]="";
	while(found[i]==0&&$(k[i])!="")
	{
	    split(gt[i,$(k[i])],a,":");
	    if(gt[i,$(k[i])]!=""&&(a[2]!="."||a[3]!="."))
	    {
		found[i]=1;
		sv[i]=$(k[i]);
		l+=len[sv[i]];
	    }
	    k[i]++;
	}
    }
    split($1,a,":");
    tot=0;
    for(i=1; i<=n-1; i++)
    {
	s3=(s3)(found[i]==1 ? gt[i,sv[i]] : "0/0:.:.")("\t");
	tot+=(found[i]);
    }
    s3=(s3)(found[i]==1 ? gt[i,sv[i]] : "0/0:.:.");
    tot+=(found[i]);
    s1="SVLEN="l/tot";SVTYPE="t;
    print a[1], a[2], "conncomp"NR, "N", "<"t">", ".", "PASS", s1, (s2!="" ? s2 : "GT:DR:DV"), s3;
}

