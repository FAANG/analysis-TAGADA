# rnaseq.meta.2.run.sheet.awk
# From a metadata tsv file of rnaseq combined fastq files per biorep and an md5sum file, make a run sheet for the faang submission
# !!! changed on June 25th 2018 not to have the same experiment alias as in SF's first submission and for this adding _2 !!!

# In order to have a single script for atacseq, hic and rnaseq we would need to
# - pass the field no for R1 and R2 or get it from metadata header (9 for atac, 7 for hic and rnaseq)
# - pass the type of attribute we have for R1 and R2 (1 and 2 for atac and rnaseq, R1 and R2 for hic)
# - pass the way the files are ending (_R1.fastq.gz for atac, .R1.fastq.gz for hic and rnaseq)
# - pass the field no for animal init id or get it from metadata  header (6 for atacseq and hic, 5 for rnaseq)
# - pass the field no for species or get it from metadata  header (2 for the 3 techniques but again by chance)
# - pass the field no for tissues or get it from metadata  header (4 for atacseq and rnaseq and 3 for hic)
# or at least have the same keys in the 3 metadata files for those attributes

# Here we need this information in the output
#############################################
# - alias = basename of the file before 1/2.fastq.gz = easy to get from the file (need to do but will do similar to hic awk)
# - run center = INRA_Toulouse = easy (as is)
# - run date = Liver: 2016-07-08T00:00:00 and CD: 2016-08-12T00:00:00
# - EXPERIMENT_REF = must correspond to the library name in the experiment_ena sheet or to the experiment_alias in the experiment_faang sheet
#   in fact SF already put a label like this one INRA_FRAGENCODE_20170716_RNASEQ_BOS_2130_CD4 so I will do the same here
#   this one seems the most difficult to get but in fact it appears that SF has already made a correspondence between
#   SAMPLE_DESCRIPTOR and EXPERIMENT_alias such as between SAMEA1088320	and INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4
#   so if I manage to get a correspondance between his short animal name such as 1346 and our nickname such as pig1 then it is fine (see hic awk)
#   !!! however here need to add _2 after the date so as not to have any conflict with SF's previous submission !!!
# - filename = the one I have for each biorep = easy (need to do but will do similar to hic awk)
# - filetype = fastq = easy (as is)
# - checksum_method = MD5 = easy (as is)
# - checksum = to compute with md5sum (gnu core_utils 8.24) = easy (done above)
# - filename_pair = the reciprocal one (if /1 then /2, if /2 then /1) = easy (need to do but will do similar to hic awk)
# - filetype_pair = fastq = easy (as is)
# - checksum_method_pair = MD5 (as is)
# - checksum_pair = to compute with md5sum (gnu core_utils 8.24) = easy (done above)

# example
# cd /work/project/fragencode/data/ebi_submission/rnaseq
# md5=/work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq.md5sum.file.tsv
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/rnaseq.meta.2.run.sheet.awk
# meta=/work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq_fastq_biorep_metadata.tsv
# awk -v fileRef=$md5 -f $pgm $meta > rnaseq.run.sheet.tsv

# rnaseq.md5sum.file.tsv
# ddf152da4c0419c51651cde2cb1bdf2a  /work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq.bos_taurus.cd4.cattle1.R1.fastq.gz
# b14a7e376eccda17e43ca5e875611e98  /work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq.bos_taurus.cd4.cattle1.R2.fastq.gz
# 82 (2 fields)

# $meta
# path	species	labExpId	tissue	animal_old	animal_short	readno
# /work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq.bos_taurus.cd4.cattle1.R1.fastq.gz	bos_taurus	cd4.cattle1	cd4	FR4934530986	cattle1	1
# 83 (7 fields)
# for hic we had this
# path	species	tissue	animal_old	animal_short	animal	readno	Species
# /work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq.capra_hircus.goat2.R1.fastq.gz	capra_hircus	liver	Goat1	goat2	FR250020007415157	R1	Goat
# 25 (8 fields)

# output file rnaseq.run.sheet.tsv
# alias	run_center	run_date	EXPERIMENT_REF	filename	filetype	checksum_method	checksum	filename_pair	filetype_pair	checksum_method_pair	checksum_pair
# rnaseq.bos_taurus.cd4.cattle1	INRA_Toulouse	2016-08-12T00:00:00	INRA_FRAGENCODE_20170716_2_RNASEQ_BOS_0986_CD4	rnaseq.bos_taurus.cd4.cattle1.R1.fastq.gz	fastq	MD5	ddf152da4c0419c51651cde2cb1bdf2a	rnaseq.bos_taurus.cd4.cattle1.R2.fastq.gz	fastq	MD5	b14a7e376eccda17e43ca5e875611e98
# 42 (12 fields)

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
$7==1{
    n=split($1,a,"/");
    split(a[n],b,".R1.fastq.gz");
    nb[b[1]]++;
    # the initial animal id
    n2=split($5,c,"");
    s="";
    for(i=n2-3; i<=n2; i++)
    {
	s=(s)(c[i]);
    }
    info[b[1],$7]="INRA_Toulouse\t"($4=="liver" ? "2016-07-08T00:00:00" : "2016-08-12T00:00:00")"\tINRA_FRAGENCODE_20170716_2_RNASEQ_"corr[$2]"_"s"_"corr[$4]"\t"a[n]"\tfastq\tMD5\t"md5[a[n]];
}

# record all complementary info for /2
$7==2{
    n=split($1,a,"/");
    split(a[n],b,".R2.fastq.gz");
    nb[b[1]]++;
    info[b[1],$7]=a[n]"\tfastq\tMD5\t"md5[a[n]];
}

# At the end write all the info about each biorep (/1 and /2 on the same row)
END{
    for(samp in nb)
    {
	print samp"\t"info[samp,1]"\t"info[samp,2];
    }
}
