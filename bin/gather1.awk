# gather1.awk
# this is to gather 3 pieces of information about atac-seq peaks on each of 4 species present in input file
# - information about accessibility level for the 4 species (file1-4)
# - information about % of size aligned and % similarity for the 4 species (file5-8)
# - information about a gx simple class for the 4 species (file9-12)

# example
# dir=/work/project/fragencode/results/atacseq
# dir2=$dir/multi/tissue.peaks.merged/project.peaks.to.human
# cd $dir2
# file1=$dir/bos_taurus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# file2=$dir/capra_hircus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# file3=$dir/gallus_gallus/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# file4=$dir/sus_scrofa/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
# file5=$dir2/bos_taurus/maf-convert.besthit.psl
# file6=$dir2/capra_hircus/maf-convert.besthit.psl
# file7=$dir2/gallus_gallus/maf-convert.besthit.psl
# file8=$dir2/sus_scrofa/maf-convert.besthit.psl
# file9=$dir/bos_taurus/tissue.peaks.merged/mergedpeaks_simpleclass.tsv 
# file10=$dir/capra_hircus/tissue.peaks.merged/mergedpeaks_simpleclass.tsv 
# file11=$dir/gallus_gallus/tissue.peaks.merged/mergedpeaks_simpleclass.tsv 
# file12=$dir/sus_scrofa/tissue.peaks.merged/mergedpeaks_simpleclass.tsv 
# awk -v fileRef1=$file1 -v fileRef2=$file2 -v fileRef3=$file3 -v fileRef4=$file4 -v fileRef5=$file5 -v fileRef6=$file6 -v fileRef7=$file7 -v fileRef8=$file8 -v fileRef9=$file9 -v fileRef10=$file10 -v fileRef11=$file11 -v fileRef12=$file12 -f gather1.awk 4species.besthit.merged.max1eachspecies.min.tsv > 4species.besthit.merged.max1eachspecies.tsv

# input
# chr1	180083	180276	gallus_gallus:AADN04007855.1:5808-8532,sus_scrofa:6:63179987-63181712	2
# chr1	180862	180935	gallus_gallus:2:15686887-15687279	1
# 212021 (5 fields)

# file1-4 are of this type
# 1	5122	7816	.	3	+	163.63	175.85	148.53	180.21	165.89	148.8	150.05	133.8
# 1	13519	15925	.	2	+	223.11	209.99	196.65	164.71	209.8	184.57	177.52	159.08
# 104985 (14 fields)

# file5-8 are of this type
# 70	19	0	0	0	0	1	4	+	10:100097024-100097471	447	34	123	chr14	107043718	87336296	87336389	2	43,46,34,77,	87336296,87336343,
# 49	12	0	0	0	0	0	0	-	10:100107305-100107831	526	0	61	chr10	133797422	21343084	21343145	1	61,	465,	21343084,
# 76253 (21 fields)

# file 9-12 are of this type
# chr	beg	end	class
# 1	5122	7816	intergenic
# 104986 (4 fields)

# output
# chr1	180083	180276	gallus_gallus:AADN04007855.1:5808-8532:1492.32:3.63436:76.7677:intergenic,sus_scrofa:6:63179987-63181712:1295.48:8.98551:81.2903:intron,	2
# chr1	180862	180935	gallus_gallus:2:15686887-15687279:206.91:18.1122:84.507:intergenic,	1
# 212021 (5 fields)


BEGIN{
    OFS="\t";
    sp[1]="bos_taurus";
    sp[2]="capra_hircus";
    sp[3]="gallus_gallus";
    sp[4]="sus_scrofa";
    while (getline < fileRef1 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    n=NF;
	    acc[sp[1]":"$1":"$2"-"$3]+=$i;
	}
	acc[sp[1]":"$1":"$2"-"$3]=acc[sp[1]":"$1":"$2"-"$3]/(n-7+1);
    }
    while (getline < fileRef2 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    n=NF;
	    acc[sp[2]":"$1":"$2"-"$3]+=$i;
	}
	acc[sp[2]":"$1":"$2"-"$3]=acc[sp[2]":"$1":"$2"-"$3]/(n-7+1);
    }
    while (getline < fileRef3 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    n=NF;
	    acc[sp[3]":"$1":"$2"-"$3]+=$i;
	}
	acc[sp[3]":"$1":"$2"-"$3]=acc[sp[3]":"$1":"$2"-"$3]/(n-7+1);
    }
    while (getline < fileRef4 >0)
    {
	for(i=7; i<=NF; i++)
	{
	    n=NF;
	    acc[sp[4]":"$1":"$2"-"$3]+=$i;
	}
	acc[sp[4]":"$1":"$2"-"$3]=acc[sp[4]":"$1":"$2"-"$3]/(n-7+1);
    }
    while (getline < fileRef5 >0)
    {
	palign[sp[1]":"$10]=($1+$2)*100/$11;
	psim[sp[1]":"$10]=$1*100/($1+$2);
    }
    while (getline < fileRef6 >0)
    {
	palign[sp[2]":"$10]=($1+$2)*100/$11;
	psim[sp[2]":"$10]=$1*100/($1+$2);
    }
    while (getline < fileRef7 >0)
    {
	palign[sp[3]":"$10]=($1+$2)*100/$11;
	psim[sp[3]":"$10]=$1*100/($1+$2);
    }
    while (getline < fileRef8 >0)
    {
	palign[sp[4]":"$10]=($1+$2)*100/$11;
	psim[sp[4]":"$10]=$1*100/($1+$2);
    }
    while (getline < fileRef9 >0)
    {
	class[sp[1]":"$1":"$2"-"$3]=$4;
    }
    while (getline < fileRef10 >0)
    {
	class[sp[2]":"$1":"$2"-"$3]=$4;
    }
    while (getline < fileRef11 >0)
    {
	class[sp[3]":"$1":"$2"-"$3]=$4;
    }
    while (getline < fileRef12 >0)
    {
	class[sp[4]":"$1":"$2"-"$3]=$4;
    }
}

{
    s="";
    split($4,a,",");
    k=1;
    while(a[k]!="")
    {
	s=(s)(a[k]":"acc[a[k]]":"palign[a[k]]":"psim[a[k]]":"class[a[k]])(",");
	k++;
    }
    $4=s;
    print;
}

