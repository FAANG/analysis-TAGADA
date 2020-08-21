# atac.meta.2.atacseq.sheet.awk

# example
# cd /work/project/fragencode/workspace/sdjebali/fragencode/data_submission/atacseq
# rna=/work/project/fragencode/data/ebi_submission/rnaseq/by.SF/experiment_ena.tsv
# meta=/work/project/fragencode/data/metadata/atacseq/atacseq_combinedreadfile_metadata.tsv
# awk -v fileRef=$rna -f atac.meta.2.atacseq.sheet.awk $meta > atacseq.atacseq.tsv

# $rna
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	TITLE	STUDY_REF	DESIGN_DESCRIPTION	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	LIBRARY_LAYOUTNOMINAL_LENGTH	NOMINAL_SDEV	LIBRARY_CONSTRUCTION_PROTOCOL	PLATFORM	INSTRUMENT_MODEL
# SAMEA1088320	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	Pig transcriptome profiling in cd4 cells by the FAANG pilot project FR-AgENCODE	INRA_FRAGENCODE_20180228_RNASEQ	Paired-end RNA-Seq of polyA+ transcripts from cd4 cells on an Illumina HiSeq 3000	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	RNA-Seq	TRANSCRIPTOMIC	Oligo-dT	PAIRED	515	not available TruSeq Stranded mRNA Sample Preparation kit V2 Guide (with a reduction of the fragmentation time from 8 to 2 minutes)	ILLUMINA	Illumina HiSeq 3000
# 1 (15 fields)
# 42 (60 fields)

# $meta
# path	species	olabExpId	tissue	sex	animal	animal_short	atacsample	read_no	labExpId
# /work/project/fragencode/data/reads/atacseq/bos_taurus/cd4/cattle3/ATAC36_atacseq_combined_R1.fastq.gz	bos_taurus	bostauruscd4FR6125612130	cd4	female	FR6125612130	cattle3	ATAC36	1	cattle3.cd4
# 77 (10 fields)

# output
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	experiment target	transposase protocol
# SAMEA1088339	INRA_FRAGENCODE_20170927_ATACSEQ_BOS_2130_CD4	open_chromatin_region	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_ATAC-seq_AG_v1_20160805.pdf
# 38 (4 fields)
# 1 (6 fields) *** all rows have 4 fields separated by tabs


# - SAMPLE_DESCRIPTOR (SAMEA1088320): OK I have from run sheet and can make the same way
# - EXPERIMENT_alias (INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4): OK I have from run sheet and can make the same way
# - experiment target = open_chromatin_region
# - transposase protocol = ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_ATAC-seq_AG_v1_20160805.pdf


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
    
    # for 3 chickens, 1 cattle and 1 pig tcell atacseq experiments, provide correspondance between metadata and biosample id manually
    # since biosample metadata is badly organized (2 files, one for liver, one for tcells, not same format and more in the file)
    # an added complication is that tcell samples are not directly linked to the animal sample but via spleen/splenocyte samples
    # for chicken and via blood sample for cattle/pig
    biosamp["GAL_2657_CD4"]="SAMEA1088331";
    biosamp["GAL_2657_CD8"]="SAMEA1088311";
    biosamp["GAL_2756_CD4"]="SAMEA1088349";
    biosamp["BOS_0986_CD8"]="SAMEA1088340";
    biosamp["SUS_1203_CD4"]="SAMEA1088372";
    
    # read the rnaseq ena file in order to associate a sample to 
    while (getline < fileRef >0)
    {
	split($2,a,"_");
	biosamp[a[5]"_"a[6]"_"a[7]]=$1;
    }

    # print the header of the faang experiment file we want
    print "SAMPLE_DESCRIPTOR", "EXPERIMENT_alias", "experiment target", "transposase protocol";
}

# Read the metadata file of fastq files but just read the /1 rows since exact same info in /2 and we need it once
$9==1{
    n=split($1,a,"/");
    split(a[n],b,"_R1.fastq.gz");
    # the initial animal id
    n2=split($6,c,"");
    s="";
    for(i=n2-3; i<=n2; i++)
    {
	s=(s)(c[i]);
    }
    auxlib=corr[$2]"_"s"_"corr[$4];
    libid="INRA_FRAGENCODE_20170927_ATACSEQ_"auxlib;
    print biosamp[auxlib], libid, "open_chromatin_region", "ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_ATAC-seq_AG_v1_20160805.pdf";
}
