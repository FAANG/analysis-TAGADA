# decorate.conncomp.from.promelt2.awk
# from:
#######
# - an input tsv file with a connected component of prom-elt2 objects at each row (nb of columns represent nb of prom-elt2 object in each conn comp)
# - a fileRef of bedpe file of prom-elt2 objects with their score in column 8 and their id as promcoord,elt2coord in column 7 (same as in input tsv file)
# and supposing that some elt2 could be equal to some prom (but not overlapping without being equal)
# outputs a tsv file without header that has 10 fields:
#######################################################
# - a comma separated list of prom
# - their nb 
# - a comma separated list of elt2
# - their nb
# - the score list (of all prom-elt2 objects in the conn comp)
# - the cumulative length of prom (supposed not to overlap)
# - the cumulative length of elt2 (same)
# - the cumulative length of prom and elt2, computed as the sum of both to which the cumul length of elt2 that are totally equal to prom is removed
# - the min distance between a prom and an elt2 in such a conn comp
# - the max distance between a prom and an elt2 in such a conn comp
# - a comma separated list of elt2 that are totally equal to prom (NA if empty)

# example
##########
# dir=/work/project/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/decorate.conncomp.from.promelt2.awk
# ct=HCmerge
# cd $dir/$ct/score2/pall
# time awk -v fileRef=$ct.pall.score2.gninfo.ens.bedpe -f $pgm $ct.pall.score2.ens.conn.connected.components.tsv > $ct.pall.score2.ens.allpromalletl2obj.prom.elt2.list.nb.sclist.prom.elt2.cumullength.sum.dist.min.max.elt2equprom.tsv
# real	0m2.054s

# fileRef
#########
# 1	1183838	1209657	1	1491937	1496017	1:1183838:1209657,1:1491937:1496017	2.18497234006452	.	.	po	22511	4859	6	2972302.0085	2.9873	1.0369	1.9571	1	1:1189288:1209265:-:ENSG00000160087.16:UBE2J2,	282280
# 64514 (22 fields)
# input file
#############
# 1:1183838:1209657,1:1491937:1496017	1:1183838:1209657,1:1647936:1674371	1:1331882:1357608,1:1435680:1451318	1:1331882:1357608,1:1491937:1496017	1:1331882:1357608,1:1757202:1764600	1:1379293:1385838,1:1647936:1674371	1:1379293:1385838,1:1814692:1819360	1:1435680:1451318,1:1331882:1357608	1:1435680:1451318,1:1496018:1502827	1:1435680:1451318,1:1647936:1674371	1:1435680:1451318,1:1752561:1757201	1:1502828:1517794,1:1647936:1674371	1:1502828:1517794,1:1814692:1819360	1:1558274:1582294,1:1647936:1674371	1:1584751:1621263,1:1647936:1674371	1:1584751:1621263,1:1814692:1819360	1:1647936:1674371,1:1183838:1209657	1:1647936:1674371,1:1379293:1385838	1:1647936:1674371,1:1435680:1451318	1:1647936:1674371,1:1502828:1517794	1:1647936:1674371,1:1539937:1549357	1:1647936:1674371,1:1558274:1582294	1:1647936:1674371,1:1584751:1621263	1:1674372:1684389,1:1621264:1634847	1:1674372:1684389,1:1814692:1819360	1:1704140:1729970,1:1752561:1757201	1:1704140:1729970,1:1814692:1819360	1:1704140:1729970,1:1896933:1910676	1:1819361:1831600,1:1739815:1747632	1:1819361:1831600,1:1752561:1757201
# ...
# tsv file with different nb of columns for each row, representing the different prom-elt2 objects that are within a connected component

# output file
##############
# 1:1183838:1209657,1:1331882:1357608,1:1379293:1385838,1:1435680:1451318,1:1502828:1517794,1:1558274:1582294,1:1584751:1621263,1:1647936:1674371,1:1674372:1684389,1:1704140:1729970,1:1819361:1831600,	11	1:1491937:1496017,1:1647936:1674371,1:1435680:1451318,1:1757202:1764600,1:1814692:1819360,1:1331882:1357608,1:1496018:1502827,1:1752561:1757201,1:1183838:1209657,1:1379293:1385838,1:1502828:1517794,1:1539937:1549357,1:1558274:1582294,1:1584751:1621263,1:1621264:1634847,1:1896933:1910676,1:1739815:1747632,	17	2.18497234006452,2.18346841284655,2.50607490345802,2.45370189969025,2.36165800236146,3.32816882381377,3.14906502956605,2.5106741974425,2.29507495792811,2.35460879825405,6.00647209644406,3.33034565728985,4.68596905408866,2.08116070475416,2.41291622150509,2.1717770163861,2.18026919137198,3.32463272367633,2.35083687076021,3.32164185204039,2.46013680123846,2.07802431392853,2.41005829485489,92.3459964117904,4.00235507786993,2.43148040906616,12.4760314900563,2.11419909995048,2.49565164085615,3.86588299256036,	223747	247819	295905	1	687276	1:1183838:1209657,1:1331882:1357608,1:1379293:1385838,1:1435680:1451318,1:1502828:1517794,1:1558274:1582294,1:1584751:1621263,1:1647936:1674371,
# 1:1621264:1634847,	1	1:1674372:1684389,	1	92.2880184499934,	13583	10017	23600	39525	39525	NA
# 4640 (11 fields)


BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	sc[$7]=$8;
    }
}

{
    # comma separated list of unique prom and its length
    s1="";
    n1=0;
    # comma separated list of unique elt2 and its length
    s2="";
    n2=0;
    # comma separated list of scores of prom-elt2
    s3="";
    # comma separated list of elt2 coord that are exactly equal to some prom
    s4="";
    # cumul length of prom
    sum1=0;
    # cumul length of elt2
    sum2=0;
    # cumul length of elt2 that are exactly equal to prom 
    sum21=0;
    # min and max dist between a prom and an elt2 that do not overlap
    m=300000000;
    M=0;
    # put each prom-elt2 of the conn comp into an array a
    n=split($0,a,"\t");
    # loop on the prom-elt2 present in the connected component but leaving the last one out to have comma separated lists that do not end with <,>
    for(k=1; k<=n; k++)
    {
	# put the prom and the elt2 coord in an array b
	split(a[k],b,",");
	# split the coord of prom and elt2
	split(b[1],b1,":");
	split(b[2],b2,":");
	# since some prom and some elt2 could appear several times in the same row, make a hashtable for each of them and use each only once
	seen1[NR,b[1]]++;
	seen2[NR,b[2]]++;
	if(seen1[NR,b[1]]==1)
	{
	    n1++;
	    s1=(s1)(b[1])(",");
	    sum1+=(b1[3]-b1[2]);
	}
	if(seen2[NR,b[2]]==1)
	{
	    n2++;
	    s2=(s2)(b[2])(",");
	    sum2+=(b2[3]-b2[2]);
	}
	s3=(s3)(sc[a[k]])(",");	
    }
    # in order to compute the min and max dist between a prom and an elt2 and also the list of elt2 that are exactly equal to a prom
    # make a double loop over the unique prom and the unique elt2 and compute those
    split(s1,a,",");
    split(s2,b,",");
    k=1;
    while(a[k]!="")
    {
	# put current prom coord in array p
	split(a[k],p,":");
	l=1;
	while(b[l]!="")
	{
	    # put current elt2 coord in array e
	    split(b[l],e,":");
	    # in case an elt2 is exactly equal to a prom, update the s4 string and sum21 to be able to compute the cumul length of prom and elt2
	    if(((p[1]==e[1])&&(p[2]==e[2])&&(p[3]==e[3]))==1)
	    {
		s4=(s4)(b[l])(",");
		sum21+=(e[3]-e[2]);
	    }
	    else
	    {
		if(p[3]<e[2])
		{
		    if((e[2]-p[3])<m)
		    {
			m=e[2]-p[3];
		    }
		    if((e[2]-p[3])>M)
		    {
			M=e[2]-p[3];
		    }
		}
		else
		{
		    if(e[3]<p[2])
		    {
			if((p[2]-e[3])<m)
			{
			    m=p[2]-e[3];
			}
			if((p[2]-e[3])>M)
			{
			    M=p[2]-e[3];
			}
		    }
		}
	    }
	    l++;
	}
	k++;
    }
    print s1, n1, s2, n2, s3, sum1, sum2, sum1+sum2-sum21, m, M, (s4!="" ? s4 : "NA");
}


 
