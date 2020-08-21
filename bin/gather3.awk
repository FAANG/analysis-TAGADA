# gather3.awk
# takes a set of merged projected peaks on human from 4 species that all have a max of 1 atac-seq peak from each species and 8 files
# of header of normalized peak activities and normalized peak activities in the 4 species (by alphabetical order) and makes a single
# matrix with human hits with 1 and only 1 hit from each species in rows and with atac-seq samples in columns and with n-1 rows in the
# header compared to the body. It also takes care of inverting the tissue and the animal in the lids from the 4 header files

# example
# cd /work/project/fragencode/results/atacseq/multi/tissue.peaks.merged/project.peaks.to.human/hierarch.clust/4species
# header1=/work/project/fragencode/results/atacseq/bos_taurus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed.readme.idx
# peak1=/work/project/fragencode/results/atacseq/bos_taurus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# header2=/work/project/fragencode/results/atacseq/capra_hircus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed.readme.idx
# peak2=/work/project/fragencode/results/atacseq/capra_hircus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# header3=/work/project/fragencode/results/atacseq/gallus_gallus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed.readme.idx
# peak3=/work/project/fragencode/results/atacseq/gallus_gallus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# header4=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed.readme.idx
# peak4=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/gather3.awk
# time zcat ../../4species.besthit.merged.max1eachspecies.tsv.gz | awk -v fileRef1=$header1 -v fileRef2=$peak1 -v fileRef3=$header2 -v fileRef4=$peak2 -v fileRef5=$header3 -v fileRef6=$peak3 -v fileRef7=$header4 -v fileRef8=$peak4 -f $pgm > 4species.orth.peaks.access.each.sample.tsv

# inputs
# chr	beg	end	id	peaknb	str	cattle3.cd4	cattle2.cd4	cattle4.cd4	cattle1.cd4	cattle3.cd8	cattle2.cd8	cattle4.cd8	cattle1.cd8
# 1 (14 fields)
# 1	5122	7816	.	3	+	163.63	175.85	148.53	180.21	165.89	148.8	150.05	133.8
# 104985 (14 fields)
# chr1	1574444	1574688	bos_taurus:16:52191332-52192318:798.697:24.4422:80.4979:intergenic,capra_hircus:16:49552453-49553534:286.484:22.2942:80.9129:tss,gallus_gallus:21:2112843-2113798:203.61:17.5916:75:tss,sus_scrofa:6:63754391-63755728:657.248:18.2498:79.918:tss,	4

# output
# cd4.cattle3	cd4.cattle2	cd4.cattle4	cd4.cattle1	cd8.cattle3	cd8.cattle2	cd8.cattle4	cd8.cattle1	liver.goat2	liver.goat3	cd4.goat3	cd4.goat1	cd4.goat4	cd8.goat3	cd8.goat2	cd8.goat1	cd8.goat4	liver.goat1	liver.goat4	liver.chicken2	liver.chicken3	cd4.chicken3	cd4.chicken2	cd4.chicken1	cd8.chicken2	liver.chicken1	liver.chicken4	cd4.pig1	cd4.pig2	cd4.pig3	cd4.pig4	cd8.pig1	cd8.pig2	cd8.pig3	cd8.pig4	liver.pig1	liver.pig2	liver.pig3
# chr1:1574444-1574688	720.44	693.8	826.49	730.25	925.12	810.17	811.66	871.65	236.03	265.46	335.01	347.2	313.56	295.79	333.76	245.83	301.12	209.69	267.87	227.29	160.87214.67	237.82	213.31	212.33	198.39	164.2	702.78	678.03	719	750.84	760.45	667.04	689.28	745.27	599.83	512.87	404.34
# 1 (38 fields)
# 1083 (39 fields)

# NOTE:
# works with raw counts but need to postprocess since in the readme.idx of those files the separator is -


BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    split($i,a,".");
	    header=(header)(a[2]"."a[1])("\t");
	}
    }
    while (getline < fileRef2 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    val["bos_taurus:"$1":"$2"-"$3]=(val["bos_taurus:"$1":"$2"-"$3])($i)("\t");
	}
    }
    while (getline < fileRef3 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    split($i,a,".");
	    header=(header)(a[2]"."a[1])("\t");
	}
    }
    while (getline < fileRef4 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    val["capra_hircus:"$1":"$2"-"$3]=(val["capra_hircus:"$1":"$2"-"$3])($i)("\t");
	}
    }
    while (getline < fileRef5 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    split($i,a,".");
	    header=(header)(a[2]"."a[1])("\t");
	}
    }
    while (getline < fileRef6 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    val["gallus_gallus:"$1":"$2"-"$3]=(val["gallus_gallus:"$1":"$2"-"$3])($i)("\t");
	}
    }
    while (getline < fileRef7 >0)
    {
	for(i=7; i<=NF-1; i++)
	{
	    split($i,a,".");
	    header=(header)(a[2]"."a[1])("\t");
	}
	split($i,a,".");
	header=(header)(a[2]"."a[1]);
    }
    while (getline < fileRef8 >0)
    {
	for(i=7; i<=NF-1; i++)
	{
	    val["sus_scrofa:"$1":"$2"-"$3]=(val["sus_scrofa:"$1":"$2"-"$3])($i)("\t");
	}
	val["sus_scrofa:"$1":"$2"-"$3]=(val["sus_scrofa:"$1":"$2"-"$3])($i);
    }
    print header;
}

# chr1	1574444	1574688	bos_taurus:16:52191332-52192318:798.697:24.4422:80.4979:intergenic,capra_hircus:16:49552453-49553534:286.484:22.2942:80.9129:tss,gallus_gallus:21:2112843-2113798:203.61:17.5916:75:tss,sus_scrofa:6:63754391-63755728:657.248:18.2498:79.918:tss,	4
$5==4{
    split($4,a,",");
    k=1;
    while(a[k]!="")
    {
	split(a[k],b,":");
	p=b[1]":"b[2]":"b[3];
	coord[b[1]]=p;
	k++;
    }
    print $1":"$2"-"$3, (val[coord["bos_taurus"]])(val[coord["capra_hircus"]])(val[coord["gallus_gallus"]])(val[coord["sus_scrofa"]]);
}
