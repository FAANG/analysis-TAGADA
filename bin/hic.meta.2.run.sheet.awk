# hic.meta.2.run.sheet.awk
# From a metadata tsv file of hic combined fastq files per biorep, make a run sheet for the faang submission
# In order to have a single script for atacseq and hic we would need to
# - pass the field no for R1 and R2 or get it from metadata header (9 for atac, 7 for hic)
# - pass the type of attribute we have for R1 and R2 (1 and 2 for atac, R1 and R2 for hic)
# - pass the way the files are ending (_R1.fastq.gz for atac, .R1.fastq.gz for hic)
# - pass the field no for animal init id or get it from metadata  header (6 for both here but only by chance)
# - pass the field no for species or get it from metadata  header (2 for both but again by chance)
# - pass the field no for tissues or get it from metadata  header (4 for atac and 3 for hic)

# example
# cd /work/project/fragencode/data/ebi_submission/hic
# md5=/work/project/fragencode/data/ebi_submission/hic/sequences/hic.md5sum.file.tsv
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/hic.meta.2.run.sheet.awk
# meta=/work/project/fragencode/data/ebi_submission/hic/sequences/hic_fastq_biorep_metadata.tsv
# awk -v fileRef=$md5 -f $pgm $meta > hic.run.sheet.tsv

# hic.md5sum.file.tsv
# e8c283e2dc8e86e3fcd0ca371411a760  /work/project/fragencode/data/ebi_submission/hic/sequences/hic.capra_hircus.goat1.R1.fastq.gz
# 3203322feb487f3af631a4f742ce0838  /work/project/fragencode/data/ebi_submission/hic/sequences/hic.capra_hircus.goat1.R2.fastq.gz
# 24 (2 fields)

# $meta
# path	species	tissue	animal_old	animal_short	animal	readno	Species
# /work/project/fragencode/data/ebi_submission/hic/sequences/hic.capra_hircus.goat2.R1.fastq.gz	capra_hircus	liver	Goat1	goat2	FR250020007415157	R1	Goat
# 25 (8 fields)
# for atacseq we had this
# /work/project/fragencode/data/metadata/atacseq/atacseq_combinedreadfile_metadata.tsv
# path	species	olabExpId	tissue	sex	animal	animal_short	atacsample	read_no	labExpId
# /work/project/fragencode/data/reads/atacseq/bos_taurus/cd4/cattle3/ATAC36_atacseq_combined_R1.fastq.gz	bos_taurus	bostauruscd4FR6125612130	cd4	female	FR6125612130	cattle3	ATAC36	1	cattle3.cd4
# 77 (10 fields)

# output file hic.run.sheet.tsv
# alias	run_center	run_date	EXPERIMENT_REF	filename	filetype	checksum_method	checksum	filename_pair	filetype_pair	checksum_method_pair	checksum_pair
# hic.gallus_gallus.chicken1	INRA_Toulouse	2016-07-04T00:00:00	INRA_FRAGENCODE_20160704_HIC_GAL_2756_LIVER	hic.gallus_gallus.chicken1.R1.fastq.gz	fastq	MD5	aea0bd5f895cad0f2bc8deb85b12d7f8	hic.gallus_gallus.chicken1.R2.fastq.gz	fastq	MD5	fe3bc7d58d0f98d0104f8420c4aac436
# 13 (12 fields)


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
$7=="R1"{
    n=split($1,a,"/");
    split(a[n],b,".R1.fastq.gz");
    nb[b[1]]++;
    # the initial animal id
    n2=split($6,c,"");
    s="";
    for(i=n2-3; i<=n2; i++)
    {
	s=(s)(c[i]);
    }
    info[b[1],$7]="INRA_Toulouse\t2016-07-04T00:00:00\tINRA_FRAGENCODE_20160704_HIC_"corr[$2]"_"s"_"corr[$3]"\t"a[n]"\tfastq\tMD5\t"md5[a[n]];
}

# record all complementary info for /2
$7=="R2"{
    n=split($1,a,"/");
    split(a[n],b,".R2.fastq.gz");
    nb[b[1]]++;
    info[b[1],$7]=a[n]"\tfastq\tMD5\t"md5[a[n]];
}

# At the end write all the info about each biorep (/1 and /2 on the same row)
END{
    for(samp in nb)
    {
	print samp"\t"info[samp,"R1"]"\t"info[samp,"R2"];
    }
}
