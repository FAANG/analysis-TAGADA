# make_matrix_from_redund_methyl_bed.awk

# example
# cd ~/dynagen/sdjebali/enhancer.gene/methylation.rnaseq/chicken
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/make_matrix_from_redund_methyl_bed.aw
# awk -v fileRef=sample.metadata.tsv -f $pgm unite.bed > cpgid.methylpcent.7samples.tsv
# on cluster because needs more than 8G of RAM

# input: unite.bed
# chr	start	end	id	score	strand	coverage	numCs	numTs
# 1	1583	1585	7127emb08	.	+	8	0	8
# 41880700 (9 fields)  *** 7 rows for each cpg (since 7 embryos)

# output: cpgid.methylpcent.7samples.tsv
# 	lid.7127.F23	lid.7127.F8	lid.7127.M15	lid.7127.M20	lid.8303.F12	lid.8303.F22	lid.8303.M14
# 1:69121798:69121800:+	100	100	72.7273	44.4444	100	84.6154	66.6667
# 1 (7 fields)
# 5982957 (8 fields) *** note that 5982957*7+1 = 41880700


BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	n++;
	if(n>=2)
	{
	    lid[$2]=$1;
	    samp[n]=$1;
	}
    }
    for(i=1; i<=n-1; i++)
    {
	s=(s)(samp[i])("\t");
    }
    print (s)(samp[i]);
}

NR>=2{
    cpg[$1":"$2":"$3":"$6]=1;
    val[$1":"$2":"$3":"$6,lid[$4]]=$8/$7*100;
}

END{
    for(c in cpg)
    {
	s=c"\t";
	for(i=2; i<=n-1; i++)
	{
	    s=(s)(val[c,samp[i]])("\t");
	}
	print (s)(val[c,samp[i]])
    }
}
