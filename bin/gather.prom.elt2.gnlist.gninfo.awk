# gather.prom.elt2.gnlist.gninfo.awk
# takes as input a tsv file with snps and info of overlapping prom frag (from prom capture hic)
# as well as
# - in fileRef1: the same with overlapping elt2 frag (from prom capture hic)
# - in fileRef2: the same with number and list of genes of prom frag (from prom capture hic)
# - in fileRef3: prom-elt connections with prom frag first, elt2 frag then and at the end the list of genes whose tss+-500bp overlap
#   the prom frag that is provided first
# and returns a tsv file corresponding to
# - the input file to which we add the number and gnlist corresponding to prom frag overlapping the snp
# - the elt2 frag and the corresponding number and gnlist corresponding to the prom frag is is connected to (otherwise NA)
# - the number and list of genes of prom frag overlapping it

# example
# cd /work/project/fragencode/workspace/sdjebali/irsd/data/parkinson/lancet.ipdgc.2019
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/gather.prom.elt2.gnlist.gninfo.awk
# ct=HCmerge
# conn=/work/project/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019/$ct/score2/pall/$ct.pall.score2.gninfo.ens.bedpe
# time awk -v fileRef1=PD_lancet_tableS2_known_loci.over.$ct.elt2frag.tsv -v fileRef2=PD_lancet_tableS2_known_loci.over.$ct.gnofpromfrag.tsv -v fileRef3=$conn -f $pgm PD_lancet_tableS2_known_loci.over.$ct.promfrag.tsv > PD_lancet_tableS2_known_loci.over.$ct.promfrag.elt2frag.withgenes.gnofpromlist.tsv


# inputs
# fileRef1
# 1	154898184	154898185	0:PMVK:NA:C-G	0.0112	.	NA
# 107 (7 fields)  *** elt2 frag
# fileRef2
# 1	154898184	154898185	0:PMVK:NA:C-G	0.0112	.	1	1:154897209-154909467:ENSG00000163344.5:PMVK:1:154891519:154930841,:-?
# 107 (8 fields)  *** gn of prom list
# fileRef3
# 1	1183838	1209657	1	1491937	1496017	1:1183838-1209657,1:1491937-1496017	2.18497234006452	.	.	po	22511	4859	6	2972302.0085	2.9873	1.0369	1.9571	1	1:1189288:1209265:-:ENSG00000160087.16:UBE2J2,
# 64514 (21 fields)  *** connections with at the end the list of genes corresponding to the prom frag and therefore to the elt2 frag as well
# standard input
# 1	154898184	154898185	0:PMVK:NA:C-G	0.0112	.	1:154891519-154930841
# 107 (7 fields)  *** prom frag

# output
# 1	154898184	154898185	0:PMVK:NA:C-G	0.0112	.	NA	0	NA	NA	0	NA	0	NA
# 107 (14 fields)


BEGIN{
    OFS="\t";
    nbgn1["NA"]=0;
    nbgn2["NA"]=0;
    gnlist1["NA"]="NA";
    gnlist2["NA"]="NA";
    
    # read elt2 frag overlap file
    while (getline < fileRef1 >0)
    {
	elt2frag[$1":"$2"-"$3]=$7;
    }
    # read gnlist of prom frag overlap file
    while (getline < fileRef2 >0)
    {
	if($7==0)
	{
	    nbgn3[$1":"$2"-"$3]=$7;
	    gnlist3[$1":"$2"-"$3]=$8;
	}
	else
	{
	    split($8,a,"?");
	    k=1;
	    while(a[k]!="")
	    {
		split(a[k],a1,":");
		nbgn3[$1":"$2"-"$3]++;
		gnlist3[$1":"$2"-"$3]=(gnlist3[$1":"$2"-"$3])(a1[4])(",");
		k++;
	    }
	}
    }
    # read prom-elt2 connection file
    # 1	1183838	1209657	1	1491937	1496017	1:1183838-1209657,1:1491937-1496017	2.18497234006452	.	.	po	22511	4859	6	2972302.0085	2.9873	1.0369	1.9571	1	1:1189288:1209265:-:ENSG00000160087.16:UBE2J2,
    while (getline < fileRef3 >0)
    {
	if($20>0)
	{
	    split($21,a,",");
	    k=1;
	    while(a[k]!="")
	    {
		split(a[k],a1,":");
		seen[$1":"$2"-"$3,a1[6]]++;
		if(seen[$1":"$2"-"$3,a1[6]]==1)
		{
		    nbgn1[$1":"$2"-"$3]++;
		    gnlist1[$1":"$2"-"$3]=(gnlist1[$1":"$2"-"$3])(a1[6])(",");
		}
		seen[$4":"$5"-"$6,a1[6]]++;
		if(seen[$4":"$5"-"$6,a1[6]]==1)
		{
		    nbgn2[$4":"$5"-"$6]++;
		    gnlist2[$4":"$5"-"$6]=(gnlist2[$4":"$5"-"$6])(a1[6])(",");
		    # print $4":"$5"-"$6, nbgn2[$4":"$5"-"$6], gnlist2[$4":"$5"-"$6];
		}
		k++;
	    }
	}
    }
}

# read prom frag overlap file
{
    elt2=elt2frag[$1":"$2"-"$3];
    print $0, nbgn1[$7], gnlist1[$7], elt2, (nbgn2[elt2]!="" ? nbgn2[elt2] : 0), (gnlist2[elt2]!="" ? gnlist2[elt2] : "NA"), nbgn3[$1":"$2"-"$3], gnlist3[$1":"$2"-"$3];
}
