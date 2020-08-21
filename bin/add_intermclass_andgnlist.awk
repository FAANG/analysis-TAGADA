# add_intermclass_andgnlist.awk
# this script takes as input
# - a complete file of predicted transcripts in gff format (includes exons and transcripts)
# - an annotated transcript gff file in gff format
# - a 2 column tsv file with for each predicted spliced transcript with an intron common with an annotated transcript, its id and the list of these transcript genes
# - a 5 column tsv file that was the previous output of refine_comptr_output.sh which has
#   predicted transcript id with t in front of it, its comptr class, its associated annotated transcript list, its refined class and its nb ex
# - a 3 column tsv file with for each monoexonic predicted transcript its new class and the list of associated reference monoexonic transcripts
# on March 27th 2018 I made the following changes
# - changed the class names (known instead of annot and alternative instead of novel_of_annot)
# - provide another file as input which is a 3 column tsv file with for each monoex pred tr has the new class and the associated reference transcripts

# and outputs the last file but with 2 additional fields:
##########################################################
# - its intermediate class which is one of those 4:
# 1) known = predicted transcript that have exactly the same exon structure as a reference tr but that do not
#    extend the annotated transcript on either side (for spliced transcripts I can use my Exact class but add
#    the constraint that there is no extension and for monoexonic transcripts I have to check whether there is any
#    monoexonic annotated transcript with exactly the same coordinates as it)
# 2) extension = predicted transcript that have exactly the same exon structure as a reference tr but that is
#    not in 1) (for spliced transcripts I can use my Exact class but check they are not in class 1) but for monoexonic
#    transcripts I have to check if it overlaps an annotated monoexonic transcript and that it is not in class 1)
# 3) alternative = new isoform or variant of reference genes = predicted spliced transcripts that are not in 1) or 2)
#    but that have at least one common intron with a reference annotated tr (I have to compute it from scratch by
#    asking for one common intron same strand as the reference and then ask that the transcript is neither
#    in class 1 nor in class 2                                                 
# 4) novel = new transcripts = the ones not in 1) or 2) or 3) (for monoex they can only be from class 1, 2 or 4, not 3)

# note: to be considered annot we used to have this condition before
# if((abs(gbeg[$1]-gbeg2[a[k]])<=10)&&(abs(gend[$1]-gend2[a[k]])<=10))
# but in fact it was wrong because then a transcript more than 10bp inward
# was considered an extension

# example
#########
# cd /work/project/fragencode/results/rnaseq/sus_scrofa/assembled
# time awk -v fileRef1=sus_scrofa_cuff_tpm0.1_2sample_complete.gff -v fileRef2=annot_tr.gff -v fileRef3=sus_scrofa_cuff_tpm0.1_2sample_transcriptid_gnlistwithcommonintron.tsv -v fileRef4=sus_scrofa_cuff_tpm0.1_2sample_monoextr_id_class_reftrlist.tsv -f $ADDCLASS sus_scrofa_cuff_tpm0.1_2sample_comp_refinedclass_nbex.tsv > sus_scrofa_cuff_tpm0.1_2sample_comp_refinedclass_nbex_intermclass.tsv

# input files
#############
# sus_scrofa_cuff_tpm0.1_2sample_complete.gff
# 1	Cufflinks	exon	561	961	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001";
# 1	Cufflinks	exon	2371	2465	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001";
# 1652756 (12 fields)

# annot_tr.gff
# AEMK02000137.1	ensembl	transcript	18908	19849	.	+	.	gene_id "ENSSSCG00000035892"; transcript_id "ENSSSCT00000066138";
# 5	ensembl	transcript	64072725	64105172	.	+	.	gene_id "ENSSSCG00000000697"; transcript_id "ENSSSCT00000057514";
# 49448 (12 fields)

# sus_scrofa_cuff_tpm0.1_2sample_transcriptid_gnlistwithcommonintron.tsv
# TCONS_00024222	ENSSSCG00000028855,
# TCONS_00083097	ENSSSCG00000002680,
# 61261 (2 fields)

# sus_scrofa_cuff_tpm0.1_2sample_comp_refinedclass_nbex.tsv
# trid	comptrclass	annottrlist	refinedtrclass	nbex
# tTCONS_00000001	Intergenic_or_antisense	.	antisense	n3
# 77541 (5 fields)

# sus_scrofa_cuff_tpm0.1_2sample_monoextr_id_class_reftrlist.tsv
# TCONS_00010286	novel	.
# TCONS_00068266	novel	.
# 9519 (3 fields)


# output file
#############
# trid	comptrclass	annottrlist	refinedtrclass	nbex	interm_class	interm_gnlist
# tTCONS_00000001	Intergenic_or_antisense	.	antisense	n3	novel	.
# 77541 (7 fields)
# and sum stats
# awk 'NR>=2{print $6}' sus_scrofa_cuff_tpm0.1_2sample_comp_refinedclass_nbex_intermclass.tsv | sort | uniq -c | sort -k1,1nr
# 35062 alternative
# 23921 known
# 15855 novel
#  2702 extension

BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)  # file of predicted transcripts with at least transcript rows in gff format
    {
	if($3=="transcript")    # for each predicted transcript, store its gene id as well as its beg and end (beg<end)
	{
	    split($10,a,"\"");
	    split($12,b,"\"");
	    gn["t"b[2]]=a[2];
	    gbeg["t"b[2]]=$4;
	    gend["t"b[2]]=$5;
	}
    }
    while (getline < fileRef2 >0)  # file of annotated transcripts with gene_id and transcript_id as 1st subfields in 9th field of a gff formated file
	                           # for each annotated transcript, store its gene id as well as its beg and end (beg<end)
    {
	split($10,a,"\"");
	split($12,b,"\"");
	gn2[b[2]]=a[2];
	gbeg2[b[2]]=$4;
	gend2[b[2]]=$5;
    }
    while (getline < fileRef3 >0)        # tsv file of predicted spliced transcripts with the nr list of annotated spliced genes with an intron common with them
                                         # for each predicted transcript, store the nr list of annotated genes with an intron common with it
    {
	commonintron_gnlist["t"$1]=$2;
    }
    while (getline < fileRef4 >0)        # tsv file of predicted monoexonic transcripts with it new class and the list of annotated monoex transcripts corresponding to it
                                         # for each predicted monoexonic transcript, store its new class and the list of annotated monoex transcripts corresponding to it
    {
	monoex_newclass["t"$1]=$2;
	monoex_reftrlist["t"$1]=$3;
    }
}

NR==1{
    # print the header with 2 additional labels
    print $0, "interm_class", "interm_gnlist";   
}

NR>=2{
    # for each predicted transcript, find its new class together with the corresponding gene list
    class="";
    gnlist="";
    if($2=="Exact")   # if the predicted transcript is from the comptr exact class (meaning it is spliced)
    {
	split($3,a,",");  
	k=1;
	while(a[k]!="")   # for each annotated transcript corresponding to this predicted transcript
	{
	    if((gbeg[$1]<gbeg2[a[k]])||(gend[$1]>gend2[a[k]]))   # check if the predicted transcript extends the annotated transcript even by a single bp on either side
		                                                 # in this case it is of the class extension but only associated to this transcript
	    {
		class="extension";   # remember it is of class "extension" but just associated to annotated transcript number k
		ok[$1,k]=1;          # and remember the number of the annotated transcript associated to it
	    }
	    k++;
	}
	if(class=="extension")
	{
	    for(i=1; i<=k; i++)
	    {
		if(ok[$1,i]==1)
		{
		    if(gn2[a[i]]!="")
			gnlist=(gnlist)(gn2[a[i]])(",");   # we only remember the genes that correspond to the annotated transcripts that are extended by the current predicted tr
		                                       # note that this list could be redundant
		}
	    }
	}
	else
	{
	    class="known";
	    for(i=1; i<=k; i++)
	    {
		if(gn2[a[i]]!="")
		    gnlist=(gnlist)(gn2[a[i]])(",");   # in case it was never an extension than it is an annotated one and it takes the genes of all the transcripts to which it is similar
	    }
	}
    }
    else   # if the predicted transcript is not from the comptr first Exact class 
    {
	if(commonintron_gnlist[$1]!="") # if it has an intron in common with an annotated gene (meaning it is spliced)
	{
	    class="alternative";
	    gnlist=commonintron_gnlist[$1];
	}
	else
	{
	    if($2=="Monoexonic")     # if it is monoexonic then we look whether it is of the new class known or the new class extension and in this case
		                     # assign this class to it and find the annot gene list from the annot tr list
	    {
		if(monoex_newclass[$1]!="")
		{
		    class=monoex_newclass[$1];
		    if(class=="novel")
		    {
			gnlist="."
		    }
		    else
		    {
			split(monoex_reftrlist[$1],a,",");
			k=1;
			while(a[k]!="")   # for each annotated transcript corresponding to this predicted transcript find its corresponding gene
			{
			    gnlist=(gnlist)(gn2[a[k]])(",");
			    k++;
			}
		    }
		}
	    }
	    else
	    {
		class="novel";
		gnlist=".";
	    }
	}
    }
    print $0, class, gnlist;
}


function abs(x){
    return (x>=0 ? x : -1*x);
}
