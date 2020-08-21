# hic.meta.2.ena.sheet.awk
# From 3 information files make a tsv file corresponding to the submission form experiment_ena sheet for hic data
# Input files are
# - $rna = the experiment_ena file used for rnaseq in order to get the biosample id (sample_descriptor)
# - $insert = the library insert size tsv file
# - $meta = the metadata file for the biorep
# An improvement would be to have a single such script for hic and atacseq but for this we would need to
# - pass the field no for R1 and R2 or get it from metadata header (9 for atac, 7 for hic)
# - pass the type of attribute we have for R1 and R2 (1 and 2 for atac, R1 and R2 for hic)
# - pass the way the files are ending (_R1.fastq.gz for atac, .R1.fastq.gz for hic)
# - pass the field no for animal init id or get it from metadata  header (6 for both here but only by chance)
# - pass the field no for species or get it from metadata  header (2 for both but again by chance)
# - pass the field no for tissues or get it from metadata  header (4 for atac and 3 for hic)
# - pass the kind of info from the metadata file that has to be connected to lib insert size (atacsample ($8 in meta) for atac, animal_old ($4 in meta) for hic)


# example
# cd /work/project/fragencode/data/ebi_submission/hic
# rna=/work/project/fragencode/data/ebi_submission/rnaseq/by.SF/experiment_ena.tsv
# insert=/work/project/fragencode/data/metadata/hic/animal_old.lib_insert_size.tsv
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/hic.meta.2.ena.sheet.awk 
# meta=/work/project/fragencode/data/ebi_submission/hic/sequences/hic_fastq_biorep_metadata.tsv
# awk -v fileRef1=$rna -v fileRef2=$insert -f $pgm $meta > hic.experiment_ena.tsv

# $rna
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	TITLE	STUDY_REF	DESIGN_DESCRIPTION	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	LIBRARY_LAYOUTNOMINAL_LENGTH	NOMINAL_SDEV	LIBRARY_CONSTRUCTION_PROTOCOL	PLATFORM	INSTRUMENT_MODEL
# SAMEA1088320	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	Pig transcriptome profiling in cd4 cells by the FAANG pilot project FR-AgENCODE	INRA_FRAGENCODE_20180228_RNASEQ	Paired-end RNA-Seq of polyA+ transcripts from cd4 cells on an Illumina HiSeq 3000	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	RNA-Seq	TRANSCRIPTOMIC	Oligo-dT	PAIRED	515	not available TruSeq Stranded mRNA Sample Preparation kit V2 Guide (with a reduction of the fragmentation time from 8 to 2 minutes)	ILLUMINA	Illumina HiSeq 3000
# 1 (15 fields)
# 42 (60 fields)

# $insert
# Pig1	403
# Pig2	474
# 12 (2 fields)

# $meta
# path	species	tissue	animal_old	animal_short	animal	readno	Species
# /work/project/fragencode/data/ebi_submission/hic/sequences/hic.capra_hircus.goat2.R1.fastq.gz	capra_hircus	liver	Goat1	goat2	FR250020007415157	R1	Goat
# 25 (8 fields)
# for atacseq we had this
# path	species	olabExpId	tissue	sex	animal	animal_short	atacsample	read_no	labExpId
# /work/project/fragencode/data/reads/atacseq/bos_taurus/cd4/cattle3/ATAC36_atacseq_combined_R1.fastq.gz	bos_taurus	bostauruscd4FR6125612130	cd4	female	FR6125612130	cattle3	ATAC36	1	cattle3.cd4
# 77 (10 fields)

# output
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	TITLE	STUDY_REF	DESIGN_DESCRIPTION	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	LIBRARY_LAYOUTNOMINAL_LENGTH	NOMINAL_SDEV	LIBRARY_CONSTRUCTION_PROTOCOL	PLATFORM	INSTRUMENT_MODEL
# SAMEA1088339	INRA_FRAGENCODE_20170927_ATACSEQ_BOS_2130_CD4	Cattle chromatin accessibility profiling in cd4 cells by the FAANG pilot project FR-AgENCODE	INRA_FRAGENCODE_20180228_ATACSEQ	Paired-end ATAC-Seq of chromatin from cd4 cells on an Illumina HiSeq 3000	INRA_FRAGENCODE_20170927_ATACSEQ_BOS_2130_CD4	ATAC-Seq	GENOMIC	not applicable	PAIRED	265	not available	ATAC-seq from Buenrostro, 2014 protocol adapted for fresh tissue and primary cells	ILLUMINA	Illumina HiSeq 3000
# 1 (15 fields)
# 38 (53 fields)

# - SAMPLE_DESCRIPTOR = this is the biosample id, so take it from rnaseq as the one associated to a given combination of animal and tissue (in second column) (ex SAMEA1088320)
# - EXPERIMENT_alias = this is the experiment alias I made up and put in the run sheet (eg INRA_FRAGENCODE_20160704_HIC_SUS_1346_LIVER)
# - TITLE = do same as rnaseq, example Pig in situ HiC chromosome conformation capture experiments in the liver tissue by the FAANG pilot project FR-AgENCODE (change the species each time)
# - STUDY_REF = do same as rnaseq and take from study_alias in study sheet so here will be INRA_FRAGENCODE_20180529_HIC
# - DESIGN_DESCRIPTION = do same as rnaseq , ex Paired-end in situ HiC for chromosome conformation capture in the liver tissue on an Illumina HiSeq 3000
# - LIBRARY_NAME = same as the experiment_alias in field no 2
# - LIBRARY_STRATEGY = Hi-C
# - LIBRARY_SOURCE = GENOMIC
# - LIBRARY_SELECTION = not applicable ?
# - LIBRARY_LAYOUT = PAIRED
# - NOMINAL_LENGTH = take from tsv file I made from NG6
# - NOMINAL_SDEV = not available
# - LIBRARY_CONSTRUCTION_PROTOCOL = look in the paper we submitted
# - PLATFORM = ILLUMINA
# - INSTRUMENT_MODEL = Illumina HiSeq 3000

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
    # correspondance between long species and nice species names
    corr2["bos_taurus"]="Cattle";
    corr2["capra_hircus"]="Goat";
    corr2["gallus_gallus"]="Chicken";
    corr2["sus_scrofa"]="Pig";

    # for 3 chicken, 1 cattle and 1 pig tcell atacseq experiments, provide correspondance between metadata and biosample id manually
    # since biosample metadata is badly organized (2 files, one for liver, one for tcells, not same format and more in the file)
    # an added complication is that tcell samples are not directly linked to the animal sample but via spleen/splenocyte samples
    # for chicken and via blood sample for cattle/pig
    biosamp["GAL_2657_CD4"]="SAMEA1088331";
    biosamp["GAL_2657_CD8"]="SAMEA1088311";
    biosamp["GAL_2756_CD4"]="SAMEA1088349";
    biosamp["BOS_0986_CD8"]="SAMEA1088340";
    biosamp["SUS_1203_CD4"]="SAMEA1088372";
    
    # read the rnaseq ena file in order to associate a sample to 
    while (getline < fileRef1 >0)
    {
	split($2,a,"_");
	biosamp[a[5]"_"a[6]"_"a[7]]=$1;
    }
    # read the insert size file in order to associate a size to each library
    while (getline < fileRef2 >0)
    {
	size[$1]=$2;
    }

    # print the header of the ena file we want
    print "SAMPLE_DESCRIPTOR", "EXPERIMENT_alias", "TITLE", "STUDY_REF", "DESIGN_DESCRIPTION",	"LIBRARY_NAME", "LIBRARY_STRATEGY", "LIBRARY_SOURCE", "LIBRARY_SELECTION", "LIBRARY_LAYOUT", "NOMINAL_LENGTH", "NOMINAL_SDEV", "LIBRARY_CONSTRUCTION_PROTOCOL", "PLATFORM", "INSTRUMENT_MODEL";
}

# Read the metadata file of fastq files but just read the /1 rows since exact same info in /2 and we need it once
$7=="R1"{
    n=split($1,a,"/");
    split(a[n],b,".R1.fastq.gz");
    # the initial animal id
    n2=split($6,c,"");
    s="";
    for(i=n2-3; i<=n2; i++)
    {
	s=(s)(c[i]);
    }
    auxlib=corr[$2]"_"s"_"corr[$3];
    libid="INRA_FRAGENCODE_20160704_HIC_"auxlib;
    print biosamp[auxlib], libid, corr2[$2]" in situ HiC chromosome conformation capture experiments in the "$3" tissue by the FAANG pilot project FR-AgENCODE", "INRA_FRAGENCODE_20180529_HIC", "Paired-end in situ HiC for chromosome conformation capture in the "$3" tissue on an Illumina HiSeq 3000", libid, "Hi-C", "GENOMIC", "unspecified", "PAIRED", size[$4], "not available", "in situ Hi-C from Rao, 2014 protocol adapted for frozen tissue", "ILLUMINA", "Illumina HiSeq 3000";
}
