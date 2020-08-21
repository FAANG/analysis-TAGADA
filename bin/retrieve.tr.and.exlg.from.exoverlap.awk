# retrieve.tr.and.exlg.from.exoverlap.awk
# takes as input
# - fileRef which is a complete exon gff file (gene_id, transcript_id) to which we compared the input exon file to using overlap (2nd file is fileRef, 1st file is input file)
# - input file is a complete exon gff file (gene_id, transcript_id) for which we want to know for each of its transcripts, the cumulative exon length that is overlapped with the transcripts of fileRef
# provides as output a with header tsv file that has for each transcript from the input file:
#############################################################################################
# - transcript id
# - gene id
# - number of exons
# - cumulative exon length
# - list of transcript ids from fileRef file whose exons overlap its exons (comma separated list)
# - list of their gene ids (comma separated list)
# - number of such transcripts
# - number of exons of those (comma separated list)
# - cumulative exon length of those (comma separated list)
# - number of exons in the intersection (comma separated list) 
# - cumulative exon length in the intersection (comma separated list)

# example
#########
# cd /work/project/fragencode/workspace/sdjebali/comparativegx/projections/genes_diff_species/livestock_to_human/sus_scrofa
# annot=/work/project/fragencode/data/species/homo_sapiens/GRCh38.90/homo_sapiens_ucsc_exons.gff
# awk -v fld=16 -v fileRef=$annot -f retrieve.tr.and.exlg.from.exoverlap.awk fragtr_to_hg38.besthit.exon.over.ens90.gff > fragtr_to_hg38.besthit.trid.gnid.nbex.lgex.overtr.idlist.gnidlist.nbex.lgex.interex.nb.lg.tsv
# trid	gnid	nbex	exlg	overtrlist	overgnid	nbtr	nbexlist	exlglist	internbexlist	interexlglist
# TCONS_00087727.1	XLOC_025169	9	233	NA	NA	0	37	15416	NA	NA
# TCONS_00062595.1	XLOC_018093	66	4282	ENST00000620394,ENST00000614484,ENST00000618262,ENST00000614370,ENST00000618411,ENST00000611405,ENST00000431645,ENST00000423692,ENST00000641476,ENST00000474896,ENST00000460194,ENST00000448764,ENST00000449389,ENST00000466432,ENST00000487311,ENST00000469602,ENST00000478693,ENST00000445337,ENST00000477529,	ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087,ENSG00000087087	19	2,20,20,5,20,20,4,4,3,3,3,12,5,2,2,2,3,3,2	582,2964,2955,556,2893,2905,397,727,992,726,575,1607,448,378,563,376,655,446,540	9,31,31,9,30,29,3,8,7,7,4,16,6,2,2,2,4,6,7	353,2943,2946,477,2884,2884,333,726,684,680,357,1591,435,134,328,331,530,431,405
# 73506 (11 fields) *** real	1m20.942s  *** we see the genes in the list are usually the same ones  *** file is 31M
# Note:
# supposes that this command was run before and that the output of this command is the one provided as input of the present script (same $annot file for both)
# overlap=/work/project/fragencode/tools/multi/Scripts/bin/overlap
# annot=/work/project/fragencode/data/species/homo_sapiens/GRCh38.90/homo_sapiens_ucsc_exons.gff
# $overlap fragtr_to_hg38.besthit.exon.gff $annot -m -1 -nr -st 1 -o fragtr_to_hg38.besthit.exon.over.ens90.gff

# fileRef file is like this
# chr1	ensembl	exon	100000637	100000739	.	-	.	gene_id "ENSG00000202259"; transcript_id "ENST00000365389";
# chr1	ensembl	exon	100148449	100149097	.	-	.	gene_id "ENSG00000122477"; transcript_id "ENST00000342895";
# 1199596 (12 fields)

# input file is like this
# chr1	p2g	exon	925708	925713	.	+	.	 gene_id "XLOC_024688"; transcript_id "TCONS_00085879.1"; nb_ov_feat2: 0 list_feat2: .
# chr1	p2g	exon	925715	925771	.	+	.	 gene_id "XLOC_024688"; transcript_id "TCONS_00085879.1"; nb_ov_feat2: 2 list_feat2: chr1_925738_925800_+,chr1_925741_925800_+,
# 5056308 (16 fields) 
# for example exon chr1_925738_925800_+ belongs to ENST00000342066 so ENST00000342066 should be in the tr list of TCONS_00085879.1 and indeed
# grep TCONS_00085879.1 fragtr_to_hg38.besthit.trid.gnid.nbex.lgex.overtr.idlist.gnidlist.nbex.lgex.interex.nb.lg.tsv
# TCONS_00085879.1	XLOC_024688	45	2724	ENST00000342066,ENST00000616016,ENST00000616125,ENST00000617307,ENST00000618181,ENST00000618323,ENST00000618779,ENST00000620200,ENST00000622503,ENST00000420190,ENST00000437963,ENST00000341065,ENST00000455979,ENST00000474461,ENST00000478729,ENST00000466827,ENST00000464948,	ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634,ENSG00000187634	17	14,13,12,13,11,12,13,9,14,7,5,12,7,4,3,2,2	2551,2391,2230,2314,2179,2159,2368,1874,2557,1578,387,2191,1731,862,318,542,657	38,37,36,37,35,35,36,32,38,6,4,35,17,6,2,5,5	2451,2291,2175,2262,2124,2110,2319,1825,2454,557,347,2162,1567,794,201,439,499	


BEGIN{
    # prints the header
    OFS="\t";
    print "trid", "gnid", "nbex", "exlg", "overtrlist", "overgnid", "nbtr", "nbexlist", "exlglist", "internbexlist", "interexlglist";
    # defines the field of the input file where we will be looking for the list of exon coordinates from fileRef that overlap the exon of the input file
    # by default it will be the 16th since if the initial exon file only had gene_id and transcript_id and was used as 1st file of overlap with -m -1 -nr -st 1 it will define a new 16 column file
    if(fld=="")
    {
	fld=16;
    }
    while (getline < fileRef >0)
	# for each transcript from the first file (overlap file2 in the comparison to the main input file) get its gene id, its number of exons and its cumulative exon length
	# also get for each exon coordinates, its transcript list (we suppose there is no redundancy in the fileRef file)
    {
	split($12,a,"\"");
	split($10,b,"\"");
	gnid1[a[2]]=b[2];
	nbex1[a[2]]++;
	cumulex1[a[2]]+=$5-$4+1;
	trlist1[$1"_"$4"_"$5"_"$7]=(trlist1[$1"_"$4"_"$5"_"$7])(a[2])(",");
    }
}

{
    # for each transcript present in the input file (= the one for which we want to know for each of its transcripts, the cumulative exon length intersected by each tr of the fileRef file)
    # get is number of exons and its cumulative exon length
    split($12,a,"\"");
    split($10,b,"\"");
    gnid2[a[2]]=b[2];
    nbex2[a[2]]++;
    cumulex2[a[2]]+=$5-$4+1;

    # the list of exons from fileRef that intersects the current exon
    split($fld,b,",");
    k=1;
    # for each exon from fileRef file that intersects the current exon
    while(b[k]!="")
    {
	# split its coordinates to see its beg and end
	split(b[k],c,"_");
	# its corresponding transcript list
	split(trlist1[b[k]],d,",");
	l=1;
	# for each transcript from fileRef whose exon overlaps the current exon, we need to know the number of exons and cumulative length in the intersection
	while(d[l]!="")
	{
	    # the first time we see a transcript associated to the current transcript, we rememeber it as being associated to the current transcript
	    # and we also compute the nb ex and cumul ex in the intersection
	    nbassoc[a[2],d[l]]++;
	    if(nbassoc[a[2],d[l]]==1)
	    {
		nbtr[a[2]]++;
		assoctr[a[2],nbtr[a[2]]]=d[l];
		trlist2[a[2]]=(trlist2[a[2]])(d[l])(",")
	    }
	    internbex[a[2],d[l]]++;
	    interexlg[a[2],d[l]]+=(min($5,c[3]))-(max($4,c[2]))+1;
	    l++;
	}
	k++;
    }
}

END{
    for(t in nbex2)
    {
	sgnid="";
	snbex="";
	slgex="";
	snbexinter="";
	slgexinter="";
	for(i=1; i<=nbtr[t]-1; i++)
	{
	    sgnid=(sgnid)(gnid1[assoctr[t,i]])(",");
	    snbex=(snbex)(nbex1[assoctr[t,i]])(",");
	    slgex=(slgex)(cumulex1[assoctr[t,i]])(",");
	    snbexinter=(snbexinter)(internbex[t,assoctr[t,i]])(",");
	    slgexinter=(slgexinter)(interexlg[t,assoctr[t,i]])(",");
	}
	print t, gnid2[t], nbex2[t], cumulex2[t], nn(trlist2[t]), nn((sgnid)(gnid1[assoctr[t,i]])), (nbtr[t]!="" ? nbtr[t] : 0), nn((snbex)(nbex1[assoctr[t,i]])), nn((slgex)(cumulex1[assoctr[t,i]])), nn((snbexinter)(internbex[t,assoctr[t,i]])), nn((slgexinter)(interexlg[t,assoctr[t,i]]));
    }
}

function min(x,y){
    return (x<=y ? x : y);
}

function max(x,y){
    return (x>=y ? x : y);
}

function nn(x){
    return (x!="" ? x : "NA")
}
