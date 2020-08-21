# pair_dnasetss1kb_to_trexpr.awk
# takes as input fld1 and fld2 integer variables and 4 files
# - fileRef1 is an rnaseq metadata file with samples to consider as rows where field fld1 is not null, labExpId in 1st field and sample name in field fld2
# - fileRef2 is a dnase-seq metadata file with samples to consider as rows where field fld1 is not null, labExpId in 1st field and sample name in field fld2
# - fileRef3 is a matrix with n-1 columns in the header and n columns in the rest and with transcripts or genes or tss in rows and rnaseq labExpIds in columns
# - an input file that is a matrix with n-1 columns in the header and n columns in the rest and with tss-1kb openness of the transcripts or genes in rows and dnase labExpIds in columns
#   this file has in 1st column a ":" separated string with the coordinates of the tss-1kb and the corresponding transcript or gene or tss
# and provides as output:
# - in the standard output all the pairs of tr and tss-1kb windows for which a tsv file will be made in Aux
# - in the Aux directory that is supposed to already exist, as many files as there are rows in the input file with transcripts or genes or tss in the input files,
#   these files are named after the 1st column of the input file and with tsv extension, without headers and with sample names in rows and rnaseq and dnase-seq values in columns
# !!! note: the input can contain less transcripts than fileRef3 !!!

# cd /work/project/dynagen/sdjebali/enhancer.gene/exploratory.analysis/corr.atac.prom.rna/10remc
# mkdir -p Aux
# time awk -v fld1=9 -v fld2=7 -v fileRef1=rnaseq/15.remc.metadata.tsv -v fileRef2=dnaseseq/15.remc.metadata.tsv -v fileRef3=rnaseq/pcgtrid_with_TPM_10samples.tsv -f ../pair_dnasetss1kb_to_trexpr.awk dnaseseq/multiBamCov.10samp.libsizenorm.tsv > tss1kb.trpairs.txt
# ls Aux/ | wc -l
# 58955

# rnaseq/15.remc.metadata.tsv
# labExpId	long.tiss	short.tiss	long.indiv	short.indiv	geo	labExpId1	labExpId2	frag.config	trfile.id	gnfile.id
# ENCSR592EZK	Fetal_Lung_Left	lungL	23964	h64	GSM1101687	h64-lungL	18216	0	ENCFF238LVW	ENCFF735RLM
# 16 (11 fields)

# dnaseseq/15.remc.metadata.tsv
# labExpId	long.tiss	short.tiss	long.indiv	short.indiv	geo	labExpId1	labExpId2	frag.config	bamfile.id
# ENCSR357RED	Fetal_Muscle_Back	muscleB	24111	h11	GSM817182	h11-muscleB	18468	0	ENCFF951KSP
# 16 (10 fields)

# rnaseq/pcgtrid_with_TPM_10samples.tsv
# ENCSR074APH	ENCSR560MDQ	ENCSR044JAQ	ENCSR317LMH	ENCSR222IGR	ENCSR620ZNQ	ENCSR332MTG	ENCSR572FXC	ENCSR499NEL	ENCSR176WMG
# ENST00000373020.4	87.48	60.78	85.35	19.19	72.21	20.77	0.00	75.08	61.38	90.35
# 1 (10 fields)
# 81814 (11 fields)

# dnaseseq/multiBamCov.10samp.libsizenorm.tsv
# ENCSR445XYW	ENCSR949JYS	ENCSR846CTA	ENCSR639XWF	ENCSR468OVW	ENCSR299INS	ENCSR964VVW	ENCSR517NHP	ENCSR953WJO	ENCSR753QKD
# chr1:68090:69091:+:ENST00000335137.3	0	1.56937	1.28517	1.02958	1.19147	0	0	0	0	0
# 1 (10 fields)
# 81814 (11 fields)

# Aux/chr10:100028006:100029007:-:ENST00000260702.3.tsv
# h64-lungR	4.34	74.5116
# h05-lungR	4.27	83.7261
# 10 (3 fields)


BEGIN{
    OFS="\t";
    while (getline < fileRef1 >0)
    {
	if($fld1!=0)
	{
	    # count the samples with n 
	    n++;
	    # store the sample name in rnasamp and remember the index of this sample in idx
	    rnasamp[n]=$fld2;
	    rnasamp2[$1]=$fld2;
	}
    }

    while (getline < fileRef2 >0)
    {
	if($fld1!=0)
	{
	    # record the sample name for each dnase-seq labExpId to be able to make the correspondance between rnaseq lid and dnaseseq lid later on
	    dnasesamp2[$1]=$fld2;
	}
    }

    # read the rna matrix
    while (getline < fileRef3 >0)
    {
	m++;
	# from the header record the sample name corresponding to each column of the file body
	if(m==1)
	{
	    for(i=1; i<=NF; i++)
	    {
		rnamatsampname[i+1]=rnasamp2[$i];
	    }
	}

	# from the rna matrix body record the expression values of each transcript in each sample
	else
	{
	    if(m>=2)
	    {
		for(i=2; i<=NF; i++)
		{
		    rnaval[$1,rnamatsampname[i]]=$i;
		}
	    }
	}
    }
}

# go over the header of the dnase matrix and record the sample name corresponding to each column of the file body
NR==1{
    for(i=1; i<=NF; i++)
    {
	dnasematsampname[i+1]=dnasesamp2[$i];
	dnaseidx[dnasesamp2[$i]]=i+1;
    }
}

# go over the body of the dnase matrix and for each dnase/tr pair for which tr has expr more than 0.1 in at least 2 samples
# print a file in Aux labelled after the 1st column value and ending in tsv and without header that has sample names as rows
# and rna and dnase values as columns
NR>=2{
    split($1,a,":");
    s="";
    for(i=2; i<=n-1; i++)
    {
	s=(s)(rnasamp[i])("\t")(rnaval[a[5],rnasamp[i]])("\t")($(dnaseidx[rnasamp[i]]))("\n");
    }
    s=(s)(rnasamp[i])("\t")(rnaval[a[5],rnasamp[i]])("\t")($(dnaseidx[rnasamp[i]]));
    print s > "Aux/"$1".tsv";
    print $1;
}
