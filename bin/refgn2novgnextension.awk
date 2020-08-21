# refgn2novgnextension.awk
# takes as input three files (3 first as fileRef and last one as input as well)
# - a tsv file from refine_comptr_output.sh with for each predicted transcript its intermediate class with respect
#   to the reference transcript and the associated ref gene list at last 2 columns
# - a complete annot gff2 file for the predicted transcripts
# - a complete annot gff2 file for the reference transcripts
# and outputs a file of annotated genes with:
# - annot gene id
# - list of predicted transcripts that are extension for the annot gene
# - list of predicted genes for those
# - list of booleans saying whether each predicted transcript extends the annot gene in 5'
# - list of booleans saying whether each predicted transcript extends the annot gene in 3'
# - list of predicted transcripts that are alternative for the annot gene
# - list of predicted genes for those
# - list of booleans saying whether each predicted transcript extends the annot gene in 5'
# - list of booleans saying whether each predicted transcript extends the annot gene in 3'

# Example
# resdir=/work/project/fragencode/results/rnaseq
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/refgn2novgnextension.awk
# sp=sus_scrofa
# cd $resdir/$sp/assembled
# awk -v fileRef1=$sp\_cuff_tpm0.1_2sample_comp_refinedclass_nbex_intermclass.tsv -v fileRef2=tmppred.gff -v fileRef3=tmpref.gff -f $pgm tmpref.gff > anngn.fragtrlist.fraggnlist.ext5p.ext3p.ext.alt.tsv

# fileRef1 = $sp\_cuff_tpm0.1_2sample_comp_refinedclass_nbex_intermclass.ts
# trid	comptrclass	annottrlist	refinedtrclass	nbex	interm_class	interm_gnlist
# tTCONS_00000001	Intergenic_or_antisense	.	antisense	n3	novel	.
# 77541 (7 fields)  *** in the last column and appart from the novel class, the list of annot gene is given

# fileRef2 = tmppred.gff
# 1	Cufflinks	exon	561	961	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001";
# 1	Cufflinks	exon	2371	2465	.	+	.	gene_id "XLOC_000001"; transcript_id "TCONS_00000001";
# 1652756 (12 fields)

# fileRef3 and input = tmpref.gff
# 5	ensembl	exon	3227814	3228279	.	-	.	gene_id "ENSSSCG00000000002"; transcript_id "ENSSSCT00000000003"; gene_version "3"; transcript_version "3"; exon_number "12"; gene_name "GTSE1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "GTSE1-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSSSCE00000000026"; exon_version "3";
# 5	ensembl	exon	3228687	3228887	.	-	.	gene_id "ENSSSCG00000000002"; transcript_id "ENSSSCT00000000003"; gene_version "3"; transcript_version "3"; exon_number "11"; gene_name "GTSE1"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_name "GTSE1-201"; transcript_source "ensembl"; transcript_biotype "protein_coding"; exon_id "ENSSSCE00000000025"; exon_version "1";
# 75328 (12 fields)
# 58651 (26 fields)
# 498784 (30 fields)
# 466736 (34 fields)
# 99 (36 fields)

# output
# ENSSSCG00000009432	TCONS_00040185,TCONS_00040191,	XLOC_011491,XLOC_011491,	0,0,	0,0,	TCONS_00040190,TCONS_00040192,	XLOC_011491,XLOC_011491,	0,0,	0,0,	
# ENSSSCG00000004074	.	.	.	.	.	.	.	.	
# 25880 (9 fields)  

BEGIN{
    OFS="\t";
    c[1]="extension";
    c[2]="alternative";
    # remember for each reference gene, the list of predicted transcripts that extend or are alternative wrt its transcripts
    while (getline < fileRef1 >0)
    {
	split($1,a,"t");
	if($6=="extension"||$6=="alternative")
	{
	    split($7,b,",");
	    k=1;
	    while(b[k]!="")
	    {
		trlist[b[k],$6]=(trlist[b[k],$6])(a[2])(",");
		k++;
	    }
	}
    }
    # remember the coordinates of the predicted transcripts and the list of transcripts of each gene
    while (getline < fileRef2 >0)
    {
	if($3=="transcript")
	{
	    split($10,a,"\"");
	    split($12,b,"\"");
	    gnpred[b[2]]=a[2];
	    strpred[b[2]]=$7;
	    gbegpred[b[2]]=$4;
	    gendpred[b[2]]=$5;
	}
    }
    # remember the coordinates of the reference genes
    while (getline < fileRef3 >0)
    {
	if($3=="gene")
	{
	    split($10,a,"\"");
	    strref[a[2]]=$7;
	    gbegref[a[2]]=$4;
	    gendref[a[2]]=$5;
	}
    }
}

$3=="gene"{
    # start the string to print with the reference gene id
    split($10,a,"\"");
    s=a[2]"\t";
    
    # for each of the two classes (extension and alternative), increment the string to be printted (s) with 4 substrings (commas separated lists) as stated above
    for(i=1; i<=2; i++)
    {
	# s1 is the list of predicted transcripts that are of the class c[i] wrt the current reference gene
	s1=trlist[a[2],c[i]];
	s2="";
	s3="";
	s4="";
	split(s1,b,",");
	k=1;
	while(b[k]!="")
	{
	    # Remember the predicted genes of those predicted transcripts that are of the class c[i] wrt the current reference gene
	    s2=(s2)(gnpred[b[k]])(",");
	    # Make lists of booleans saying whether those predicted transcripts extend the ref gene in 5' and in 3'
	    if(strref[a[2]]=="+")
	    {
		s3=(s3)((gbegpred[b[k]]<gbegref[a[2]]) ? "1" : "0")(",");
		s4=(s4)((gendpred[b[k]]>gendref[a[2]]) ? "1" : "0")(",");
	    }
	    else
	    {
		s3=(s3)((gendpred[b[k]]>gendref[a[2]]) ? "1" : "0")(",");
		s4=(s4)((gbegpred[b[k]]<gbegref[a[2]]) ? "1" : "0")(",");
	    }
	    k++;
	}
	s=(s)(nn(s1))("\t")(nn(s2))("\t")(nn(s3))("\t")(nn(s4))("\t");
    }
    print s;
} 

function nn(x)
{
    return (x!="" ? x : ".");
}
