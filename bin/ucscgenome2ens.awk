# ucscgenome2ens.awk

# usage
# awk -f ucscgenome2ens.awk ucsc_genome.fa > ens_genome.fa

# example
# cd /work/project/fragencode/data/species/capra_hircus/CHIR_1.0.102/softmasked
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/ucscgenome2ens.awk
# time awk -f $pgm chi_ref_ASM170441v1_allchr.mfa > ~/dynagen/sdjebali/comparativegx/data/species/capra_hircus/capra_hircus.fa
# real	0m48.576s

# input
# >chr1 gi|1060061422|ref|NC_030808.1| Capra hircus breed San Clemente chromosome 1, ASM170441v1, whole genome shotgun sequence
# GAAGAGCCAGGAGAACTATGAgatgtcaagggaacatttcatgcaagatggccattaaaggacagaaaca

BEGIN{
    OFS="\t";
}

{
    if($1~/^>chr[0-9]+$/||$1==">chrX"||$1==">chrY"||$1==">chrM"||$1==">chrMT"||$1==">chrW"||$1==">chrZ")
    {
	new=substr($1,5,length($1)-4);
	if($1!=">chrM")
	{
	    print ">"new;
	}
	else
	    print ">MT";
    }
    else
    {
	print;
    }
}
