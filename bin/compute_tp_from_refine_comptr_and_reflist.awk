# compute_tp_from_refine_comptr_and_reflist.awk
# takes as input an output from refine_comptr_output.sh (tsv file with 7 fields and 3 kinds of classes in field no 2, 4 and 6
# and a label and a fileRef file as parameter
# this last file is a tsv file of reference transcripts (transcript id, gene id) on which we want to make the assessment
# = compute the number of predicted transcripts from the input that are matching those reference transcripts (those should be
# a subset of the tr used as annotation input to refine_comptr_output.sh)
# with different criteria, from the least to the most stringent (number and prop of total predicted):
# 1. matching any transcript from the annotation used as input to refine_comptr_output.sh (Exact or Extension or Inclusion or Overlap)
# 2. same but with a reference transcript from fileRef
# 3. matching any transcript from the annotation used as input to refine_comptr_output.sh but with Exact or Extension or Inclusion class
# 4. same but with a reference transcript from fileRef
# 5. matching any transcript from the annotation used as input to refine_comptr_output.sh but with Exact or Extension class
# 6. same but with a reference transcript from fileRef
# 7. matching any transcript from the annotation used as input to refine_comptr_output.sh but with Exact class
# 8. same but with a reference transcript from fileRef
# it also adds at the end (maybe later add the same but with the annot tr not just the ref ones):
# 9. number and proportion of reference transcripts from fileRef matched with Exact or Extension or Inclusion or Overlap class
# 10. number and proportion of reference transcripts from fileRef matched with Exact or Extension or Inclusion class
# 11. number and proportion of reference transcripts from fileRef matched with Exact or Extension class
# 12. number and proportion of reference transcripts from fileRef matched with Exact class

# example
# dir=/work/project/fragencode/workspace/sdjebali/trreconstruct/Livestock/Spring2018_Oana_Mouse
# cd $dir/star.cufflinks.cuffmerge
# ref=$dir/gencvM16.complete/gencode.vM16.reflnc.trid.gnid.tsv
# awk -v label=star.cufflinks.cuffmerge -v fileRef=$ref -f $dir/compute_tp_from_refine_comptr_and_reflist.awk lncRNA_exon_comp_refinedclass_nbex_intermclass.tsv > lncRNA.TP.info.tsv


# kind of input
# trid	comptrclass	annottrlist	refinedtrclass	nbex	interm_class	interm_gnlist
# tTCONS_00000082	Intergenic_or_antisense	.	intergenic	n3	novel	.
# tMSTRG.1224.1	Overlap	ENSMUST00000193042.5,ENSMUST00000046110.15,ENSMUST00000195311.5,ENSMUST00000194217.1,ENSMUST00000192821.1,ENSMUST00000194658.1,	extension	n2	novel	.
# 26060 (7 fields)

# kind of output
# label	nbtot	match.4classes.any.nb	match.4classes.any.pcent	match.4classes.ref.nb	match.4classes.ref.pcent	match.3classes.any.nb	match.3classes.any.pcent	match.3classes.ref.nb	match.3classes.ref.pcent	match.2classes.any.nb	match.2classes.any.pcent	match.2classes.ref.nb	match.2classes.ref.pcent	match.1class.any.nb	match.1class.any.pcent	match.1class.ref.nb	match.1class.ref.pcent	nbtotref	matchedbypred.4classes.nb	matchedbypred.4classes.pcent	matchedbypred.3classes.nb	matchedbypred.3classes.pcent	matchedbypred.2classes.nb	matchedbypred.2classes.pcent	matchedbypred.1class.nb	matchedbypred.1class.pcent
# star.cufflinks.cuffmerge	26059	14748	56.5947	1787	6.85752	7209	27.6641	1012	3.8835	6202	23.7998	828	3.17741	1223	4.6932	196	0.752139	2425	1003	41.3608	467	19.2577	350	14.433	120	4.94845
# 2 (27 fields)


BEGIN{
    OFS="\t";
    if(label=="")
    {
	label="my_label";
    }
    # remember the reference transcript ids
    while (getline < fileRef >0)
    {
	ok[$1]=1;
	nbtotref++;
    }
}


NR>=2{
    # for each predicted transcripts
    # increment the counter
    ntot++;

    # for the ones that are of the 4 matching classes (Overlap or Inclusion or Extension or Exact)
    # increment the 1st counter and find out whether one of the associated transcript is from the reference
    # and in this case increment the 1st counter ok
    if($3!=".")
    {
	nok1++;
	split($3,a,",");
	found=0;
	k=1;
	while(a[k]!="")
	{
	    if(ok[a[k]]==1)
	    {
		found=1;
		nb[a[k]]++;
		if(nb[a[k]]==1)
		{
		    nbref1++;
		    if($2=="Inclusion")
		    {
			nbref2++;
		    }
		    else
		    {
			if($2=="Extension")
			{
			    nbref2++;
			    nbref3++;
			}
			else
			{
			    if($2=="Exact")
			    {
				nbref2++;
				nbref3++;
				nbref4++
			    }
			}
		    }
		}
	    }
	    k++;
	}
	if(found==1)
	{
	    nok1ref++;
	}	

	# for the ones that are of the class Inclusion increment the second counter and the second counter ok if there is a match from the ref
	if($2=="Inclusion")
	{
	    nok2++;
	    if(found==1)
	    {
		nok2ref++;
	    }
	}
	else
	{
	    # for the ones that are of the class Extension increment the third counter and the third counter ok if there is a match from the ref
	    if($2=="Extension")
	    {
		nok2++;
		nok3++;
		if(found==1)
		{
		    nok2ref++;
		    nok3ref++;
		}
	    }
	    else
	    {
		# for the ones that are of the class Exact increment the fourth counter and the fourth counter ok if there is a match from the ref
		if($2=="Exact")
		{
		    nok2++;
		    nok3++;
		    nok4++;
		    if(found==1)
		    {
			nok2ref++;
			nok3ref++;
			nok4ref++;
		    }
		}
	    }
	    
	}
    } 
}

END{
    print "label", "nbtot", "match.4classes.any.nb", "match.4classes.any.pcent", "match.4classes.ref.nb", "match.4classes.ref.pcent", "match.3classes.any.nb", "match.3classes.any.pcent", "match.3classes.ref.nb", "match.3classes.ref.pcent", "match.2classes.any.nb", "match.2classes.any.pcent", "match.2classes.ref.nb", "match.2classes.ref.pcent", "match.1class.any.nb", "match.1class.any.pcent", "match.1class.ref.nb", "match.1class.ref.pcent", "nbtotref", "matchedbypred.4classes.nb", "matchedbypred.4classes.pcent", "matchedbypred.3classes.nb", "matchedbypred.3classes.pcent", "matchedbypred.2classes.nb", "matchedbypred.2classes.pcent", "matchedbypred.1class.nb", "matchedbypred.1class.pcent";
    print label, ntot, nok1, nok1/ntot*100, nok1ref, nok1ref/ntot*100, nok2, nok2/ntot*100, nok2ref, nok2ref/ntot*100, nok3, nok3/ntot*100, nok3ref, nok3ref/ntot*100, nok4, nok4/ntot*100, nok4ref, nok4ref/ntot*100, nbtotref, nbref1, nbref1/nbtotref*100, nbref2, nbref2/nbtotref*100, nbref3, nbref3/nbtotref*100, nbref4, nbref4/nbtotref*100;
}
