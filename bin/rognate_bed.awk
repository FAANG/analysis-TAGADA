# rognate_bed.awk
# rognate end features going above the end of the chr, and completely discard if it becomes smaller than the beg

# usage
# awk -v fileRef=chrom.len -f rognate_bed.awk features.bed > features_rognate.bed

# example
# chrom=/work/project/fragencode/data/species/gallus_gallus/Gallus_gallus-5.0.87/chrom.len
# infile=/work/project/fragencode/results/atacseq/gallus_gallus/cd4/chicken1/peaks/ATAC66_atacseq_combined.q10.noMT.nodup.bedgraph
# awk -v fileRef=$chrom -f rognate_bed.awk $infile > ${infile%.bedgraph}_ok.bedgraph
# 1	196202544
# 23475 (2 fields)
# 1	139	287	0.0350715

BEGIN{
    OFS="\t";
    while (getline < fileRef >0)
    {
	gend[$1]=$2
    }
}

{
    $3=($3<=gend[$1] ? $3 : gend[$1]);
    if($2<$3)
	print $0;
}
