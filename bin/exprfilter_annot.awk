# exprfilter_annot.awk

# takes as
# - fileRef a transcript expression matrix (tsv file wuith header) with transcripts in rows and samples in columns
#   that has transcript id, some information and then the TPMs in a set of samples from field number fstexpr (see below)
# - input file a gtf or gff2 (with at least exons) file, with transcript_id in the 9th field
# - mexpr: a min threshold for expression (0.1 by default)
# - msamp: a min number of samples where we want to this the min expression to consider the transcript in the annotation (2 by default)
# - fstexpr: the first field with TPM information in the transcript expression matrix (3 by default)

# and outputs
# -  a gff2 exon file where only transcripts with expression above but that needs to pass to gffok.awk to be a real gff file wrt tabs and spaces

# example 
# cd /work/project/fragencode/workspace/geneswitch/analyses/qc_and_first_results/elements
# expr=/work/project/fragencode/workspace/geneswitch/results/counts/transcripts_TPM.tsv
# ref=/work/project/fragencode/workspace/geneswitch/pipelines/rnaseq/tests/sus_scrofa/allbig/data/species/sus_scrofa.gtf
# exprfilter=/work/project/fragencode/tools/multi/Scripts/Awk/exprfilter_annot.awk
# gffok=/work/project/fragencode/tools/multi/Scripts/Awk/make_gff_ok.awk
# time awk -v mexpr=0.1 -v msamp=2 -v fstexpr=3 -v fileRef=$expr -f $exprfilter $ref | awk -f $gffok > ref.annot.tpm0.1.2samples.exons.gff


BEGIN{
    OFS="\t";
    if(mexpr=="")
    {
	mexpr=0.1;
    }
    if(msamp=="")
    {
	msamp=2;
    }
    if(fstexpr=="")
    {
	fstexpr=3;
    }
    while (getline < fileRef >0)
    {
	ok=0;
	k=fstexpr;
	while(ok<msamp&&k<=NF)
	{
	    if($k>=mexpr)
	    {
		ok++;
	    }
	    k++;
	}
	if(ok==msamp)
	{
	    keep[$1]=1;
	}
    }
}

$3=="exon"{
    fine=0;
    split($0,a,"\t");
    split(a[9],a1,"; ");
    k=1;
    while(fine==0&&a1[k]!="")
    {
	split(a1[k],a2," ");
	split(a2[2],a3,"\"");
	if(a2[1]=="transcript_id"&&keep[a3[2]]==1)
	{
	    fine=1;
	    print;
	}
	k++;
    }
}
