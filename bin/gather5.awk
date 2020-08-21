# gather5.awk
# takes a set of merged projected peaks on human from 4 species that all have a max of 1 atac-seq peak from each species and 6 files
# of header of normalized peak activities and normalized peak activities in bt, hc and ss (by alphabetical order) and makes a single
# matrix with human hits with 1 and only 1 hit from bt, hc and ss 
# in rows and with atac-seq samples in columns and with n-1 rows in the
# header compared to the body. It also takes care of inverting the tissue and the animal in the lids from the 3 header files
# note : I have to do it for the 4 triplets of species but I cannot really generalize it since
# for the last species we have to avoid writing a last tab in the row and the last species is not always the same ...

# example
# cd /work/project/fragencode/results/atacseq/multi/tissue.peaks.merged/project.peaks.to.human/hierarch.clust/3species
# header1=/work/project/fragencode/results/atacseq/bos_taurus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed.readme.idx
# peak1=/work/project/fragencode/results/atacseq/bos_taurus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# header2=/work/project/fragencode/results/atacseq/capra_hircus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed.readme.idx
# peak2=/work/project/fragencode/results/atacseq/capra_hircus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# header3=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed.readme.idx
# peak3=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/gather5.awk
# time cat ../../4species.besthit.merged.max1eachspecies.tsv | awk -v fileRef1=$header1 -v fileRef2=$peak1 -v fileRef3=$header2 -v fileRef4=$peak2 -v fileRef5=$header3 -v fileRef6=$peak3 -f $pgm > 3species.bthcss.orth.peaks.access.each.sample.tsv

# inputs
# chr	beg	end	id	peaknb	str	cattle3.cd4	cattle2.cd4	cattle4.cd4	cattle1.cd4	cattle3.cd8	cattle2.cd8	cattle4.cd8	cattle1.cd8
# 1 (14 fields)
# 1	5122	7816	.	3	+	163.63	175.85	148.53	180.21	165.89	148.8	150.05	133.8
# 104985 (14 fields)
# chr1	1574444	1574688	bos_taurus:16:52191332-52192318:798.697:24.4422:80.4979:intergenic,capra_hircus:16:49552453-49553534:286.484:22.2942:80.9129:tss,gallus_gallus:21:2112843-2113798:203.61:17.5916:75:tss,sus_scrofa:6:63754391-63755728:657.248:18.2498:79.918:tss,	4

# output
# cd4.cattle3	cd4.cattle2	cd4.cattle4	cd4.cattle1	cd8.cattle3	cd8.cattle2	cd8.cattle4	cd8.cattle1	liver.goat2	liver.goat3	cd4.goat3	cd4.goat1	cd4.goat4	cd8.goat3	cd8.goat2	cd8.goat1	cd8.goat4	liver.goat1	liver.goat4	liver.chicken2	liver.chicken3	cd4.chicken3	cd4.chicken2	cd4.chicken1	cd8.chicken2	liver.chicken1	liver.chicken4
# chr1:1471663-1472133	823.84	846.55	767.59	809.92	807.51	820.4	749.35	793.69	247.64	255.2	303.78	304.95	310.56	234.24	251.13231.91	277.41	236.79	118.24	73.07	105.46	215.97	183.26	200.54	165.55	77	112.31
# 1 (27 fields)
# 1362 (28 fields)

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
	for(i=7; i<=NF-1; i++)
	{
	    split($i,a,".");
	    header=(header)(a[2]"."a[1])("\t");
	}
	split($i,a,".");
	header=(header)(a[2]"."a[1]);
    }
    while (getline < fileRef6 >0)
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
$5>=3&&$4~/bos_taurus/&&$4~/capra_hircus/&&$4~/sus_scrofa/{
    split($4,a,",");
    k=1;
    while(a[k]!="")
    {
	split(a[k],b,":");
	p=b[1]":"b[2]":"b[3];
	coord[b[1]]=p;
	k++;
    }
    print $1":"$2"-"$3, (val[coord["bos_taurus"]])(val[coord["capra_hircus"]])(val[coord["sus_scrofa"]]);
}
