# atac.meta.2.run.sheet.awk
# From a metadata tsv file of atac-seq combined fastq files per biorep, make a run sheet for the faang submission
# An improvement would be to have a single awk script for atacseq and hic but for this we would need to
# - pass the field no for R1 and R2 or get it from metadata header (9 for atac, 7 for hic)
# - pass the type of attribute we have for R1 and R2 (1 and 2 for atac, R1 and R2 for hic)
# - pass the way the files are ending (_R1.fastq.gz for atac, .R1.fastq.gz for hic)
# - pass the field no for animal init id or get it from metadata  header (6 for both here but only by chance)
# - pass the field no for species or get it from metadata  header (2 for both but again by chance)
# - pass the field no for tissues or get it from metadata  header (4 for atac and 3 for hic)
# - pass the kind of info from the metadata file that has to be connected to lib insert size (atacsample ($8 in meta) for atac, animal_old ($4 in meta) for hic)


# example
# cd /work/project/fragencode/workspace/sdjebali/fragencode/data_submission/atacseq
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/atac.meta.2.run.sheet.awk 
# meta=/work/project/fragencode/data/metadata/atacseq/atacseq_combinedreadfile_metadata.tsv
# awk -v fileRef=atacseq.fastq.md5sum.tsv -f $pgm $meta > atacseq.run.sheet.tsv

# atacseq.fastq.md5sum.tsv
# 2667d5177aba5dc9952253f2e4769a96  /work/project/fragencode/data/reads/atacseq/bos_taurus/cd4/cattle3/ATAC36_atacseq_combined_R1.fastq.gz
# 1ba18427949573fe4efc776c449763e3  /work/project/fragencode/data/reads/atacseq/bos_taurus/cd4/cattle3/ATAC36_atacseq_combined_R2.fastq.gz
# 76 (2 fields)

# $meta
# /work/project/fragencode/data/metadata/atacseq/atacseq_combinedreadfile_metadata.tsv
# path	species	olabExpId	tissue	sex	animal	animal_short	atacsample	read_no	labExpId
# /work/project/fragencode/data/reads/atacseq/bos_taurus/cd4/cattle3/ATAC36_atacseq_combined_R1.fastq.gz	bos_taurus	bostauruscd4FR6125612130	cd4	female	FR6125612130	cattle3	ATAC36	1	cattle3.cd4
# 77 (10 fields)

# output file should look like this
# alias	run_center	run_date	EXPERIMENT_REF	filename	filetype	checksum_method	checksum	filename_pair	filetype_pair	checksum_method_pair	checksum_pair
# ATAC85_atacseq_combined	INRA_Toulouse	2017-09-27T00:00:00	INRA_FRAGENCODE_20170927_ATACSEQ_GAL_1346_LIVER	ATAC85_atacseq_combined_R1.fastq.gz	fastq	MD5	5d691705080ba80bfb8b50a94032a7f9	ATAC85_atacseq_combined_R2.fastq.gz	fastq	MD5	dec5821800ef3fa3626e7e605f6af376
# 39 (12 fields)


BEGIN{
    OFS="\t";
    # correspondance between long species and short species names
    corr["bos_taurus"]="BOS";
    corr["capra_hircus"]="CAP";
    corr["gallus_gallus"]="GAL";
    corr["sus_scrofa"]="SUS";
    # correspondance between lowercase and uppercase tissue names
    corr["cd4"]="CD4";
    corr["cd8"]="CD8";
    corr["liver"]="LIVER";

    # read the md5sum file in order to associate an md5sum to each fastq file 
    while (getline < fileRef >0)
    {
	n=split($2,a,"/");
	md5[a[n]]=$1;
    }

    # print the header of the run file we want
    print "alias", "run_center", "run_date", "EXPERIMENT_REF", "filename", "filetype", "checksum_method", "checksum", "filename_pair", "filetype_pair", "checksum_method_pair", "checksum_pair";
}

# Reads the metadata file of combined reads (1 and 2) for each biorep
# record all info for /1
$9==1{
    n=split($1,a,"/");
    split(a[n],b,"_R1.fastq.gz");
    nb[b[1]]++;
    # the initial animal id
    n2=split($6,c,"");
    s="";
    for(i=n2-3; i<=n2; i++)
    {
	s=(s)(c[i]);
    }
    info[b[1],$9]="INRA_Toulouse\t2017-09-27T00:00:00\tINRA_FRAGENCODE_20170927_ATACSEQ_"corr[$2]"_"s"_"corr[$4]"\t"a[n]"\tfastq\tMD5\t"md5[a[n]];
}

# record all complementary info for /2
$9==2{
    n=split($1,a,"/");
    split(a[n],b,"_R2.fastq.gz");
    nb[b[1]]++;
    info[b[1],$9]=a[n]"\tfastq\tMD5\t"md5[a[n]];
}

# At the end write all the info about each biorep (/1 and /2 on the same row)
END{
    for(samp in nb)
    {
	print samp"\t"info[samp,1]"\t"info[samp,2];
    }
}
