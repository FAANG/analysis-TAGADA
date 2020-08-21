# prom.elt2.5reg.snp.allinfo.awk
# takes as input
# - fileRef1 = a bed4 file of regions with snp coord list associated (overlapping in general)
# - fileRef2 = a tsv snp file output by plink gwas (done with dosage as info in general)
# - standard input = a tsv file with header that has snpid, celltype, prom.elt2.region.gxorder, same extended by 500kb on each side
#   and then the coordinates of 5 subregions which are
#   * 500kb.before.region
#   * first.gx.order.elt.region
#   * prom.elt2.inside.region
#   * second.gx.order.elt.region
#   * 500kb.after.region
#   and then for those 5 subregions their
#   * sizes
#   * nbsnps
#   * snp densities
# and outputs
# - a tsv file with header that has
#   * snpid
#   * celltype
#   * region number
#   * reg.coord
#   * reg.size
#   * reg.nbsnp
#   * reg.denssnp
#   * snp.coords
#   * snp.freqs
#   * snp.info
#   * snp.ORs
#   * snp.SEs
#   * snp.pvals
# so 13 fiels in total

# example
# cd ~/gwas/parkinson.lancet.ipdgc.2019
# snps=~/gwas/parkinson_france/allchr.MP_status.assoc.dosage
# pgm=~/tools/multi/Scripts/Awk/prom.elt2.5reg.snp.allinfo.awk
# time awk -v fileRef1=2cases.promelt2500kb.5subregions.uniq.withpostqcsnps.bed -v fileRef2=$snps -f $pgm 2cases.snpid.promelt2500kb.5subregions.5reg.size.postqcsnp.nb.density.tsv > 2cases.snpid.promelt2500kb.5subregions.5reg.size.postqcsnp.nb.density.coord.freq.dosage.or.se.pval.lists.tsv

# inputs
# 2cases.promelt2500kb.5subregions.uniq.withpostqcsnps.bed
# 5 133858641 133865903 5:133858719,5:133859609,5:133860101,5:133860435,5:133860901,5:133861426,5:133861663,5:133861716,5:133861756,5:133864436,5:133864599,5:133865452,
# 14 (4 fields)
# $snps
# SNP	A1	A2	FRQ	INFO	OR	SE	P
#      1:49298   T   C  0.6058  0.7511  0.8914  0.0638 0.07169
# 5854882 (8 fields)
# 2cases.snpid.promelt2500kb.5subregions.5reg.size.postqcsnp.nb.density.tsv
# snpid	ct	prom.elt2.region.gxorder	prom.elt2.region.gxorder.ext500kbeachside	500kb.before.region	first.gx.order.elt.region	prom.elt2.inside.region	second.gx.order.elt.region	500kb.after.region	reg1.size	reg2.size	reg3.size	reg4.size	reg5.size	reg1.nbsnp	reg2.nbsnp	reg3.nbsnp	reg4.nbsnp	reg5.nbsnp	reg1.denssnp	reg2.denssnp	reg3.denssnp	reg4.denssnp	reg5.denssnp
# 5:134199105	HCmerge	5:133858641-134199358	5:133358641-134699358	5:133358641-133858641	5:133858641-133865903	5:133865903-134192251	5:134192251-134199358	5:134199358-134699358	500000	7262	326348	7107	500000	810	12	358	5	896	1.62	1.65244	1.09699	0.703532	1.792
# 5 (24 fields)  *** from 5 to 9 for reg coord, from 10 to 14 for sizes, from 15 to 19 for nbsnps, from 20 to 24 for densities

# output
# 21 (13 fields)


BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	snplist[$1":"$2"-"$3]=$4;
    }
    while (getline < fileRef2 >0)
    {
	n++;
	if(n>=2)
	{
	    for(j=1; j<=5; j++)
	    {
		info[$1,j]=$(5+j-2);
	    }
	}
    }
}

NR==1{
    print $1"\t"$2"\treg.number\treg.coord\treg.size\treg.nbsnp\treg.denssnp\tsnp.coords\tsnp.freqs\tsnp.infos\tsnp.ORs\tsnp.SEs\tsnp.pvals";
}

NR>=2{
    # for each of the 5 subregions, print the snpid and the cell type
    for(i=1; i<=5; i++)
    {
	# put the snp list into a variable
	l[i]=snplist[$(5+i-1)]
	
        # remember for subreg i, the snpid, the cell type, the reg number, the reg coord, the reg size, the reg nbsnp, the reg snp density and the snp coord
	s[i]=$1"\t"$2"\t"(i)("\t")($(5+i-1))("\t")($(10+i-1))("\t")($(15+i-1))("\t")($(20+i-1))("\t")(l[i])("\t");
	
	# split the comma separated list of snp coord for region i
	split(l[i],a,",");

	# fill in the string of subbreg i with 5 infos about it which are the 5 infos about each of its snps (and therefore adds 5 comma separated lists)
	for(j=1; j<=4; j++)
	{
	    k=1;
	    while(a[k]!="")
	    {
		s[i]=(s[i])(info[a[k],j])(",");
		k++;
	    }
	    s[i]=(s[i])("\t");
	}
	k=1;
	while(a[k]!="")
	{
	    s[i]=(s[i])(info[a[k],j])(",");
	    k++;
	}
	print s[i];
    }
}
