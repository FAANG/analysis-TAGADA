# norm.lib.size.and.select.samples.awk
# take as inputs
# - fldno = field number in fileRef corresponding to a string which when non equal to "0" says that labExpId needs to be kept
# - fileRef= metadata file with header and having labExpId and a string which when non equal to "0" says that labExpId needs to be kept in column fldno
# - a matrix tsv file with a header with n-1 column than the rest with read counts for each object/row in each sample/column and with header
#   being a list of labExpId
# and returns as output
# the matrix tsv file with only the samples to be kept and where the counts have been normalized by library size and multiplied
# by the average read number in all libraries to keep (in order to have a number that looks like a read count)

# Example
#########
# cd /work/project/dynagen/sdjebali/enhancer.gene/exploratory.analysis/corr.atac.prom.rna/15remc/dnaseseq
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/norm.lib.size.and.select.samples.awk
# awk -v fldno=9 -v fileRef=15.remc.metadata.tsv -f $pgm multiBamCov.tsv > multiBamCov.10samp.libsizenorm.tsv

# 15.remc.metadata.tsv
# labExpId	long.tiss	short.tiss	long.indiv	short.indiv	geo	labExpId1	labExpId2	frag.config
# ENCSR357RED	Fetal_Muscle_Back	muscleB	24111	h11	GSM817182	h11-muscleB	18468	0
# 16 (9 fields)

# multiBamCov.tsv
# ENCSR445XYW	ENCSR949JYS	ENCSR846CTA	ENCSR639XWF	ENCSR468OVW	ENCSR888GBS	ENCSR299INS	ENCSR033STL	ENCSR072NBR	ENCSR964VVW	ENCSR357RED	ENCSR517NHP	ENCSR006XED	ENCSR953WJO	ENCSR753QKD
# chr1:68090:69091:+:ENST00000335137.3	0	1	1	1	1	0	0	0	0	0	0	0	0	0	0
# 1 (15 fields)
# 81814 (16 fields)

# multiBamCov.10samp.libsizenorm.tsv
# ENCSR445XYW	ENCSR949JYS	ENCSR846CTA	ENCSR639XWF	ENCSR468OVW	ENCSR299INS	ENCSR964VVW	ENCSR517NHP	ENCSR953WJO	ENCSR753QKD
# chr1:68090:69091:+:ENST00000335137.3	0	1.56937	1.28517	1.02958	1.19147	0	0	0	0	0	
# 1 (10 fields)
# 81814 (11 fields)

BEGIN{
    # reads the metadafile and remember the labExpId to keep and how many they are
    OFS="\t";
    while (getline < fileRef >0)
    {
	n++;
	if(n==1)
	{
	    found=0;
	    i=1;
	    while(found==0&&i<NF)
	    {
		if($i=="labExpId")
		{
		    found=1;
		    idxlid=i;
		}
		i++;
	    }
	}
	else
	{
	    if($fldno!=0)
	    {
		ok[$idxlid]=1;
		nbsamp++;
	    }
	}
    }
}

NR==1{
    # from the header of the matrix only keep the labExpIds we want
    for(i=1; i<=NF-1; i++)
    {
	if(ok[$i]==1)
	{
	    okidx[i+1]=1;
	    s=(s)($i)("\t");
	}
    }
    if(ok[$i]==1)
    {
	okidx[i+1]=1;
	s=(s)($i);
    }
    # note that if the last one is not to be kept we will have 1 more tab at the end, to solve later on
    print s;
}

NR>=2{
    # j is a counter for the number of rows, each count for each labExpId to keep and each row will be stored in count[i,j]
    j++;
    elt[j]=$1;
    for(i=2; i<=NF; i++)
    {
	if(okidx[i]==1)
	{
	    count[i,j]=$i;
	    # total count of all reads in the matrix for the labExpIds to keep
	    tot+=$i;
	    # total count of all reads in library i
	    N[i]+=$i;
	}
    }
}

END{
    # average count per library
    avlib=tot/nbsamp;
    # k is a counter for the rows
    for(k=1; k<=j; k++)
    {
	s=elt[k]"\t";
	# c is a counter for the columns
	for(c=2; c<=i-2; c++)
	{
	    if(okidx[c]==1)
	    {
		s=(s)(count[c,k]*avlib/N[c])("\t");
	    }
	}
	if(okidx[c]==1)
	{
	    s=(s)(count[c,k]*avlib/N[c]);
	}
	# note that if the last one is not to be kept we will have 1 more tab at the end, to solve later on
	print s;
    }
} 
