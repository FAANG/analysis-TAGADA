# add_maxbegend_to_chimjunc2.awk
# improvement of add_maxbegend_to_chimjunc.awk since it enables to write all the information associated to a junction
# and that is located from field no $(fld2+1) on, and which are samechrstr, okgxorder, dist, (don, acc), gnlist1, gnlist2,
# gnnamelist1, gnnamelist2, gnbtlist1, gnbtlist2

# example
# cd ~sdjebali/ENCODE_AWG/Analyses/Mouse_Human/Chimeras/Human
# cat LID16627/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt LID16628/Chimsplice/distinct_junctions_nbstaggered_nbtotalsplimappings_withmaxbegandend_samechrstr_okgxorder_dist_gnlist1_gnlist2_gnname1_gnname2_bt1_bt2_from_split_mappings_part1overA_part2overB_only_A_B_indiffgn_and_inonegn.txt | awk '$2>=2' | awk -v fld1=4 -v fld2=5 -f ~sdjebali/Awk/add_maxbegend_to_chimjunc2.awk > junctions_of_LID16627_and_LID16628_with_maxbegandend.txt

# input
#######
# chr10_102114354_+:chr17_20688573_+ 2 2 102114321 20688619 0 NA NA ENSG00000099194.4, ENSG00000226549.2, SCD, SCDP1, protein_coding, pseudogene,
# chr10_102114354_+:chr17_20688573_+ 3 3 102114306 20688619 0 NA NA ENSG00000099194.4, ENSG00000226549.2, SCD, SCDP1, protein_coding, pseudogene,
# chr10_7844267_+:chr14_54458042_- 5 15 7844245 54457982 0 NA NA ENSG00000165629.14, ENSG00000224004.1, ATP5C1, ATP5C1P1, protein_coding, pseudogene,
# chr10_7844267_+:chr14_54458042_- 3 9 7844250 54457982 0 NA NA ENSG00000165629.14, ENSG00000224004.1, ATP5C1, ATP5C1P1, protein_coding, pseudogene,
# chr11_122931879_-:chr1_247393396_+ 3 21 122931925 247393432 0 NA NA ENSG00000109971.8, ENSG00000214144.2, HSPA8, RP11-488L18.1, protein_coding, pseudogene,
# chr11_122931879_-:chr1_247393396_+ 3 21 122931928 247393432 0 NA NA ENSG00000109971.8, ENSG00000214144.2, HSPA8, RP11-488L18.1, protein_coding, pseudogene,
# chr10_103799770_-:chr10_103588956_- 20 60 103799835 103588887 1 0 NA ENSG00000120029.7, ENSG00000120049.12, C10orf76, KCNIP2, protein_coding, protein_coding,
# chr10_103799770_-:chr10_103588956_- 15 45 103799831 103588891 1 0 NA ENSG00000120029.7, ENSG00000120049.12, C10orf76, KCNIP2, protein_coding, protein_coding,

# output
########
# a n = 3+m column file with junction id first, then beg and end, and then m fields which are the ones from fld2+1 to the end in the input matrix

BEGIN{OFS="\t"}

{
    if(info[$1]=="")
    {
	for(i=(fld2+1); i<=NF-1; i++)
	{
	    info[$1]=(info[$1])($i)("\t");
	}
	info[$1]=(info[$1])($i);
    }
    split($1,a,":"); 
    split(a[1],a1,"_"); 
    split(a[2],a2,"_"); 
    if(a1[3]=="+")
    {
	if((beg[$1]=="")||($fld1<beg[$1]))
	{
	    beg[$1]=$fld1;
	}

	if(a2[3]=="+")    # +/+
	{
	    if((end[$1]=="")||($fld2>end[$1]))
	    {
		end[$1]=$fld2;
	    }
	}
	else              # +/-
	{
	    if((end[$1]=="")||($fld2<end[$1]))
	    {
		end[$1]=$fld2;
	    }
	}
    }
    else  
    {
	if((beg[$1]=="")||($fld1>beg[$1]))
	    {
		beg[$1]=$fld1;
	    }
	
	if(a2[3]=="+")    # -/+
	{
	    if((end[$1]=="")||($fld2>end[$1]))
	    {
		end[$1]=$fld2;
	    }
	}
	else              # -/-  *** most complicated case since here the beg of the junction has to do with the second part of the junction, and the end of the junction with the first
	{
	    if((end[$1]=="")||($fld2<end[$1]))
	    {
		end[$1]=$fld2;
	    }
	}
    }
}

END{
    for(j in beg)
    {
	print j, beg[j], end[j], info[j];
    }
}
