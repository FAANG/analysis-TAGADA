# conncompSV2vcf_better.awk
# takes as input:
# - a parameter st for the type of SV we are considering in the conversion (only present in the two next input files)
# - a string runid that is a unique identifier of the set we are considering
# - a string indlist that is a comma separated list of individuals in the order we want to display them in the output files
#   and that are also used to index $3 of the vcf file and the objects of the tsv file (indfather:id15703 for example
#   and s2 will be offspring,father,mother)  *** absolutely compulsory
# - a string infolist with the pieces of information from the INFO $8 field in the initial concat vcf that we want to keep
# - tsv file of the connected components (cc) made from the SVs of this type from several individuals
# - a fileRef file that is the concatenation of all the SVs of the different individuals as a vcf file (non gzipped)
# and outputs
# - a main vcf file with the SVs corresponding to the connected components (cc) and with genotypes in each of the individuals
#   when there is consistency between all the genotypes of a given individual for the cc
# - an aux vcf file with the conflicting cc SVs in a file named "conflicting.cc.of"runid".vcf"
# It has to be noted that the SV identifiers present in the cc tsv file is the one present in $3 of the vcf file
# Notes:
# - if there are different SVs from the same individual inside a cc and their genotypes do not agree flag as conflicting
# - if for a cc there is no SV from a given individual then put homozygous reference for this individual
# - svlen and svtype are needed in the resulting file
# - $3 in the concatenation file is supposed to be the unique identifier of each individual SV
# - only PASS variants are considered and therefore PASS will be provided in the outputs
# - only the GT information will be provided in $9 and we will suppose it is given as first info in the vcf input file (1st info before first :)

# Important
###########
# !!! it is essential that the vcf file passed as input and the ss tsv file have as identifiers of the objects ind$ind:id$id !!!

# example
# cd /work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds
# map=minimap
# sc=sniffles
# st=DEL
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/conncompSV2vcf_better.awk
# time awk -v st=$st -v runid=$sc.$map.$st.allrun -v indlist=offspring,father,mother -v infolist=RE -v fileRef=$sc.$map.$st.allrun.concat.sorted.vcf -f $pgm $sc.$map.$st.allrun.connected_components.tsv > $sc.$map.$st.allrun.cc.vcf


# input files
##############
# 1. the connected component file 
#################################
# /work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/sniffles.minimap.DEL.allrun.connected_components.tsv
# indfather:id15703	indmother:id23010
# indoffspring:id31365	indmother:id56224_1
# 1021 (1 fields)
# 347 (2 fields)
# 229 (3 fields)
# 31 (4 fields)
# 18 (5 fields)
# 7 (6 fields)
# 7 (7 fields)
# 6 (8 fields)
# 8 (9 fields)
# 1 (10 fields)
# 2 (13 fields)
# 2 (14 fields)
# 1 (15 fields)
# 2 (20 fields)
# 1 (58 fields)
# 2. the concatenation of all the individual SVs of all the individuals
#######################################################################
# bcftools view -H /work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/sniffles.minimap.DEL.allrun.concat.sorted.vcf | head
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	mapping/trio1_offspring_run1.minimap.forsniffles.bam
# 1	158099	indmother:id26_6	N	<DEL>	.	PASS	PRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=1;END=158141;STD_quant_start=6.95222;STD_quant_stop=5.68624;Kurtosis_quant_start=-0.041094;Kurtosis_quant_stop=-0.826443;SVTYPE=DEL;SUPTYPE=AL;SVLEN=-42;STRANDS=+-;STRANDS2=0,3,0,3;RE=3;REF_strand=0,0;AF=1	GT:DR:DV	1/1:0:3

# output files (1 last of header and after)
##############
# main one
# 10	100784298	indfather:id15703,indmother:id23010,	N	<DEL>	.	PASS	END=100793086;SVLEN=-8789;SVTYPE=DEL;POSLIST=100785684,100784298,;ENDLIST=100791904,100793086,;SVLENLIST=-6220,-8788,;GTLIST=0/1,0/1,;DETLIST=0,1,1,;CONFLLIST=0,0,0,	GT	0/0	0/1	0/1
# NKLS02001553.1	29450	indoffspring:id31365,indmother:id56224_1,	N	<DEL>	.	PASS	END=32085;SVLEN=-2636;SVTYPE=DEL;POSLIST=29450,29458,;ENDLIST=32085,32059,;SVLENLIST=-2635,-2601,;GTLIST=1/1,1/1,;DETLIST=1,0,1,;CONFLLIST=0,0,0,	GT	1/1	0/0	1/1
# 10	100786823	indoffspring:id13339,	N	<DEL>	.	PASS	END=100790628;SVLEN=-3806;SVTYPE=DEL;POSLIST=100786823,;ENDLIST=100790628,;SVLENLIST=-3805,;GTLIST=1/1,;DETLIST=1,0,0,;CONFLLIST=0,0,0,	GT	1/1	0/0	0/0

# conflicting rows
# 17	68056561	indoffspring:id20114,indmother:id34816_1,indfather:id23909,indmother:id34817,	N	<DEL>	.	PASS	END=68217030;SVLEN=-160470;SVTYPE=DEL;POSLIST=68057710,68057711,68057710,68056561,;ENDLIST=68216492,68216491,68216492,68217030,;SVLENLIST=-158782,-158780,-158782,-160469,;GTLIST=1/1,1/1,1/1,0/0,;DETLIST=1,1,1,;CONFLLIST=1,1,1,	GT	1/1	1/1	1/1	
# 17	70993072	indfather:id23972_3,indfather:id24017_0,indfather:id24020_0,indfather:id24022,	N	<DEL>	.	PASS	END=71000580;SVLEN=-7509;SVTYPE=DEL;POSLIST=70995250,70993072,70994054,70994271,;ENDLIST=70999577,70999072,70999655,71000580,;SVLENLIST=-4327,-6000,-5601,-6309,;GTLIST=1/1,0/0,0/1,0/0,;DETLIST=0,1,0,;CONFLLIST=0,1,0,	GT	0/0	1/1	0/0	
# 83 (12 fields)  *** indeed we see the conflicts for those rows



BEGIN{
    if(st=="")
    {
	st="DEL";
    }
    if(runid=="")
    {
	runid="runid";
    }
    # infolist is the list of informations from the info no 8 field from vcf input file that we want to remember
    o=split(infolist,c,",");
    n=split(indlist,a,",");
    for(k=1; k<=n-1; k++)
    {
	ind[k]=a[k];
	str0=(str0)(a[k])("\t");
    }
    ind[k]=a[k];
    str0=(str0)(a[k]);
    OFS="\t";
    # print the exact same header as in the vcf input file (apart from last line where we put the indiv names) and otherwise remember info about individual SVs 
    while (getline < fileRef >0)  
    {
	if($1~/#/)
	{
	    if($1!="#CHROM")
	    {
		print;
		print > "conflicting.cc.of."(runid)".vcf";
	    }
	    else
	    {
		print "##INFO=<ID=POSLIST,Number=1,Type=String,Description=\"Comma-separated list of begining positions of the SVs included in the connected component (cc)\">";
		print "##INFO=<ID=ENDLIST,Number=1,Type=String,Description=\"Comma-separated list of end positions of the SVs included in the connected component (cc)\">";
		print "##INFO=<ID=SVLENLIST,Number=1,Type=String,Description=\"Comma-separated list of lengths of the SVs included in the connected component (cc)\">";
		print "##INFO=<ID=GTLIST,Number=1,Type=String,Description=\"Comma-separated list of genotypes of the SVs included in the connected component (cc)\">";
		print "##INFO=<ID=QUALLIST,Number=1,Type=String,Description=\"Comma-separated list of qualities of the SVs included in the connected component (cc)\">";
		for(i=1; i<=o; i++)
		{
		    print "##INFO=<ID="c[i]"LIST,Number=1,Type=String,Description=\"Comma-separated list of "c[i]"s of the SVs included in the connected component (cc)\">";
		}
		print "##INFO=<ID=DETLIST,Number=1,Type=String,Description=\"Comma-separated list of booleans indicating for each individual whether the SV of the cc was detected\">";
		print "##INFO=<ID=CONFLLIST,Number=1,Type=String,Description=\"Comma-separated list of booleans indicating for each individual whether there were conflicting genotypes\">";
		print $1, $2, $3, $4, $5, $6, $7, $8, $9, str0;

                # and same in the conflict file
		print "##INFO=<ID=POSLIST,Number=1,Type=String,Description=\"Comma-separated list of begining positions of the SVs included in the connected component (cc)\">" > "conflicting.cc.of."(runid)".vcf";
		print "##INFO=<ID=ENDLIST,Number=1,Type=String,Description=\"Comma-separated list of end positions of the SVs included in the connected component (cc)\">" > "conflicting.cc.of."(runid)".vcf";
		print "##INFO=<ID=SVLENLIST,Number=1,Type=String,Description=\"Comma-separated list of lengths of the SVs included in the connected component (cc)\">" > "conflicting.cc.of."(runid)".vcf";
		print "##INFO=<ID=GTLIST,Number=1,Type=String,Description=\"Comma-separated list of genotypes of the SVs included in the connected component (cc)\">" > "conflicting.cc.of."(runid)".vcf";
		print "##INFO=<ID=QUALLIST,Number=1,Type=String,Description=\"Comma-separated list of qualities of the SVs included in the connected component (cc)\">" > "conflicting.cc.of."(runid)".vcf";
		for(i=1; i<=o; i++)
		{
		    print "##INFO=<ID="c[i]"LIST,Number=1,Type=String,Description=\"Comma-separated list of "c[i]"s of the SVs included in the connected component (cc)\">" > "conflicting.cc.of."(runid)".vcf";
		}
		print "##INFO=<ID=DETLIST,Number=1,Type=String,Description=\"Comma-separated list of booleans indicating for each individual whether the SV of the cc was detected\">" > "conflicting.cc.of."(runid)".vcf";
		print "##INFO=<ID=CONFLLIST,Number=1,Type=String,Description=\"Comma-separated list of booleans indicating for each individual whether there were conflicting genotypes\">" > "conflicting.cc.of."(runid)".vcf";
		print $1, $2, $3, $4, $5, $6, $7, $8, $9, str0 > "conflicting.cc.of."(runid)".vcf";
	    }
	}
	else
	{
	    # remember the chr, the beg, the quality, the end, the length and all the info from infolist that were asked for, of each individual SV
	    gchr[$3]=$1;
	    gbeg[$3]=$2;
	    qual[$3]=$6;
	    split($8,a,";");
	    m=1;
	    while(a[m]!="")
	    {
		split(a[m],b,"=");
		if(b[1]=="END")
		{
		    gend[$3]=b[2];
		}
		if(b[1]=="SVLEN")
		{
		    len[$3]=b[2];
		}
		# remember all the info asked for in infolist if they exist in the input vcf file (their corresponding header rows will be added in all cases above)
		for(i=1; i<=o; i++)
		{
		    if(b[1]==c[i])
		    {
			info[$3,i]=b[2];
		    }
		}
		m++;
	    }
	    # for each individual SV from the vcf file remember its genotype (here no 0/0:./. since no bcftools merge done before so not test needed for that string)
	    # we suppose the 1st info before 1st : in the 10th field is the genotype
	    split($10,a,":");
	    gt[$3]=a[1];

 	}
    }
}

# indfather:id15703	indmother:id23010
# indoffspring:id31365	indmother:id56224_1
# 1	158099	indmother:id26_6	N	<DEL>	.	PASS	PRECISE;SVMETHOD=Snifflesv1.0.11;CHR2=1;END=158141;STD_quant_start=6.95222;STD_quant_stop=5.68624;Kurtosis_quant_start=-0.041094;Kurtosis_quant_stop=-0.826443;SVTYPE=DEL;SUPTYPE=AL;SVLEN=-42;STRANDS=+-;STRANDS2=0,3,0,3;RE=3;REF_strand=0,0;AF=1	GT:DR:DV	1/1:0:3
{
    # for each connected component, find the genotype of each indiv and if the indiv is not present put 0/0:./.
    # if the indiv has several SV in the cc them output if only if all genotypes are the same
    # and put in an aux file the conflicting connected components (ie with different genotypes for a given individual)
    for(k=1; k<=n; k++)
    {
	det[ind[k]]=0;
    }
    ok=1;
    k=1;
    id="";  # this is the id of the connected component, made of the ids of its individual SVs
    gb=1000000000;
    ge=-1;
    gblist="";
    gelist="";
    lenlist="";
    gtlist="";
    qualist="";
    detlist="";
    conflist="";
    str2="";    # string for the genotypes of the indiv, from $10 to $10+n-1
    addinfostr=""; # complete string for the additional info frominfolist
    # strings for each info from infolist to remember from each SV the cc is made of
    for(i=1; i<=o; i++)
    {
	str11[i]=((c[i])"LIST=");
    }
    # Go along the complete set of SVs of the cc since all the begs, ends and svlen are needed
    # Get the begining of the region as the smallest begining of all the SVs to put in column 2 later on (gb)
    # Get the end of the region as the biggest end of all the SVs to put as info in the INFO field
    # Make the INFO string (col8) with the global end and svlen obtained from global beg and end (svlen is negative for DEL),
    # and also with the list of begs, list of ends and list of lengths
    # and also with the list of qualities provided in the $6 field
    # and also with the list of 
    # and also with the list of 3 booleans specifying whether the SV was present in each of the 3 individuals
    # and also with the list of 3 booleans specifying whether there was any conflict in any of the 3 individuals
    # n2 is the number of SVs in the connected component we are looking at and which to make a vcf row from
    n2=split($0,a,"\t");
    while(k<=n2)
    {
	id=(id)(a[k])(",");
	split(a[k],b,":");
	split(b[1],b1,"ind");
	# remember the SV was detected in the current individual of the current cc
	det[b1[2]]=1;
	gblist=(gblist)(gbeg[a[k]])(",");
	gelist=(gelist)(gend[a[k]])(",");
	lenlist=(lenlist)(len[a[k]])(",");  # should also be end-beg+1 and *1 for del but better to take it from info field (see above)
	gtlist=(gtlist)(gt[a[k]])(",");
	qualist=(qualist)(qual[a[k]])(",");
	if(gbeg[a[k]]<gb)
	{
	    gb=gbeg[a[k]];
	}
	if(gend[a[k]]>ge)
	{
	    ge=gend[a[k]];
	}
	# if there is one conflict for the row and a given individual, regarding genotypes, flag as conflictual
	if(g[b1[2],NR]!=""&&g[b1[2],NR]!=gt[a[k]])
	{
	    ok=0;
	}
	else
	{
	    g[b1[2],NR]=gt[a[k]];
	}
	for(i=1; i<=o; i++)
	{
	    str11[i]=(str11[i])(info[a[k],i])(",");
	}
	k++;
    }
    for(i=1; i<=o; i++)
    {
	addinfostr=(addinfostr)(str11[i])(";");
    }
    # Make the string of genotypes of all the indivs in the order given by indlist and make the DETLIST and CONFLLIST at the same time
    for(k=1; k<=n; k++)
    {
	str2=(str2)(g[ind[k],NR]!="" ? g[ind[k],NR] : "0/0")("\t");
	detlist=(detlist)(det[ind[k]])(",");
	conflist=(conflist)((g[ind[k],NR]!=""&&ok==0) ? 1 : 0)(",");
    }
    str1="END="ge";SVLEN="(st=="DEL" ? -1 : 1)*(ge-gb+1)";SVTYPE="st";POSLIST="gblist";ENDLIST="gelist";SVLENLIST="lenlist";GTLIST="gtlist";QUALLIST="(qualist)";"(addinfostr)"DETLIST="detlist";CONFLLIST="conflist;
    gsub(/\,\;/,";",str1);
    if(ok==1)
    {
	print gchr[$1], gb, id, "N", "<"st">", ".", "PASS", str1, "GT", str2;
    }
    else
    {
	print gchr[$1], gb, id, "N", "<"st">", ".", "PASS", str1, "GT", str2 > "conflicting.cc.of."(runid)".vcf";
    }
}

