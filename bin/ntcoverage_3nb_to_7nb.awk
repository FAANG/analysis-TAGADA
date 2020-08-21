# ~/Awk/ntcoverage_3nb_to_7nb.awk
# this script takes as input a coverage file like the one produced by ~/bin/compute_coverage_from_elements_gal_nocluster_better.sh
# = a file with for a given RNAseq experiment or set of experiments, the nt covered by its contigs
# the exonic nt covered by its contigs and the genic nt covered by its contigs
# and produces as output the 7 proportions that a figure like figs4 of the rna paper uses
# = a file with 7 proportions:
# - % of genome covered 
# - % of rnaseq nt in exons
# - % of exonic nt covered
# - % of rnaseq nt in introns
# - % of intronic nt covered
# - % of rnaseq nt in intergenic
# - % of intergenic nt covered
# it also requires nttotalnogap from the genome, ntinprojex and ntinprojgn  from the annotation
# as well as label1 and label2 for the experiment or set of experiments. when it is a set of experiments
# then label1 must be all and label2 must be "cumulative", when it is an experiment then label1 is usually
# the labexpid and label2 must be "individual"

# example
# cd ~sdjebali/Bean/Contigs
# awk -v nttotalnogap=$nttotalnogap -v ntinprojex=$ntinprojex -v ntinprojgn=$ntinprojgn -v label1=all -v label2=cumulative -f ~/Awk/ntcoverage_3nb_to_7nb.awk ~/Bean/Contigs/Coverages/32exp_03212013_contigs_coverage.txt >> basename_indiv_and_proj_notshort_notnucsubcomp_7propforboxplots_nonnull.txt

{
    print "total", label1, $1, $1/nttotalnogap*100, label2, "0";
    print "exonic", label1, $2, ($1!=0 ? $2/$1*100 : 0), label2, "2";
    print "exonic", label1, $2, $2/ntinprojex*100, label2, "1";
    print "intronic", label1, ($3-$2), ($1!=0 ? ($3-$2)/$1*100 : 0), label2, "2"; 
    print "intronic", label1, ($3-$2), ($3-$2)/(ntinprojgn-ntinprojex)*100,  label2, "1"; 
    print "intergenic", label1, ($1-$3), ($1!=0 ? ($1-$3)/$1*100 : 0), label2, "2"; 
    print "intergenic", label1, ($1-$3), ($1-$3)/(nttotalnogap-ntinprojgn)*100, label2, "1";
}