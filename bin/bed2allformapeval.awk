# bed2allformapeval.awk
# from a the min overlap ratio and a bed file of primary longest mapped reads that were simulated and for which we know where they should align
# (indicated in field 4 of the bed file), outputs a file that has everything needed to make a mapeval output like file except the ordering of the Q rows, the U row
# and the 2 last cumulative numbers (for which ordering is needed):
# - 1 Q line for each mapping quality in the input file that has after Q
#   * the quality threshold
#   * the number of reads mapped with mapping quality equal to or greater than the threshold
#   * the number of wrong mappings

# Takes as input an r variable indicating the minimum intersection over union bases to consider the mapping as correct
# and a bed file of mapped reads that has
# - in name field 4: the read id in the form <S1_readno>!<mapseqid>!<gbeg>!<gend>!<strand> where the coorinates represent where the read should map according to simulation
# - in score field 5: the mapping quality of the read
# and provides as output a tsv file with as many Q lines as mapq in the bed file 
# each Q line is like this
# Q       60      32478   0       0.000000000     32478
# Q       22      16      1       0.000030775     32494
# and has successively 
# - Q
# - quality threshold
# - number of reads mapped with mapping quality equal to or greater than the threshold
# - number of reads that are not correctly mapped among those (according to overlap amount and r threshold)

# Example:
##########
# cd /work/project/dynagen/sdjebali/SeqOccIn/snakemake/pipeline/complete/minimap2
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/bed2allformapeval.awk
# awk -v r=0.1 -f $pgm simulated_reads.genome.aln.primary.longest.sorted.bed | sort -k2,2nr > simulated_reads.genome.aln.primary.longest.all4mapeval.out

# NOTE: in case we want a real mapeval output file need to add another awk command after the sort
# awk -v r=0.1 -f $pgm simulated_reads.genome.aln.primary.longest.sorted.bed | sort -k2,2nr | awk -v tot=10000 'BEGIN{OFS="\t"} {cummap+=$3; cumerr+=$4; print $0, cumerr/cummap, cummap} END{print "U", tot-cummap}' > simulated_reads.genome.aln.primary.longest.homemapeval.out
															 
# input
# scaffold19	25465520	25484832	S1_5659!scaffold19!25465520!25484863!+	60	+
# scaffold12	13884280	13885626	S1_2829!scaffold12!13882871!13884280!-	60	-
# 8689 (6 fields)

# output
# Q	60	8593	0
# Q	59	2	0
# 43 (4 fields)

{
    ntot[$5]++;
    split($4,a,"!");
    if($1==a[2]&&$6==a[5])   # if the chr and the strand are fine then look at overlap and if smaller than r then put as wrong for this mapq
    {
	inter=min($3,a[4]);
	union=max($2,a[3]);
	if((inter/union)<r)
	{
	    nwrong[$5]++;
	}
    }
    else  # otherwise the mapping is incorrect for this mapq
    {
	nwrong[$5]++;
    }
}
END{
    OFS="\t";
    for(q in ntot)
    {
	print "Q", q, ntot[q], (nwrong[q]!="" ? nwrong[q] : 0);
    }
}


function min(x,y)
{
   return (x<=y ? x : y)
}

function max(x,y)
{
    return (x>=y ? x : y)
}

