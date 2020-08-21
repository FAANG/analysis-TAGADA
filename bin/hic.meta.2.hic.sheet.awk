# hic.meta.2.hic.sheet.awk
# From 2 information files make a tsv file corresponding to the submission form hic sheet for hic data
# Input files are
# - $rna = the experiment_ena file used for rnaseq in order to get the biosample id (sample_descriptor)
# - $meta = the metadata file for the biorep

# example
# cd /work/project/fragencode/data/ebi_submission/hic
# rna=/work/project/fragencode/data/ebi_submission/rnaseq/by.SF/experiment_ena.tsv
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/hic.meta.2.hic.sheet.awk
# meta=/work/project/fragencode/data/metadata/atacseq/atacseq_combinedreadfile_metadata.tsv
# awk -v fileRef=$rna -f $pgm $meta > hic.hic.tsv

# $rna
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	TITLE	STUDY_REF	DESIGN_DESCRIPTION	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	LIBRARY_LAYOUTNOMINAL_LENGTH	NOMINAL_SDEV	LIBRARY_CONSTRUCTION_PROTOCOL	PLATFORM	INSTRUMENT_MODEL
# SAMEA1088320	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	Pig transcriptome profiling in cd4 cells by the FAANG pilot project FR-AgENCODE	INRA_FRAGENCODE_20180228_RNASEQ	Paired-end RNA-Seq of polyA+ transcripts from cd4 cells on an Illumina HiSeq 3000	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	RNA-Seq	TRANSCRIPTOMIC	Oligo-dT	PAIRED	515	not available TruSeq Stranded mRNA Sample Preparation kit V2 Guide (with a reduction of the fragmentation time from 8 to 2 minutes)	ILLUMINA	Illumina HiSeq 3000
# 1 (15 fields)
# 42 (60 fields)

# $meta
# path	species	tissue	animal_old	animal_short	animal	readno	Species
# /work/project/fragencode/data/ebi_submission/hic/sequences/hic.capra_hircus.goat2.R1.fastq.gz	capra_hircus	liver	Goat1	goat2	FR250020007415157	R1	Goat
# 25 (8 fields)

# output


# - SAMPLE_DESCRIPTOR (SAMEA1088320) = OK I have from run sheet and can make the same way
# - EXPERIMENT_alias (INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4) = OK I have from run sheet and can make the same way
# - experiment target = chromatin
# - restriction enzyme = HindIII
# - restriction site = AAGCTAGCTT


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
    print "SAMPLE_DESCRIPTOR", "EXPERIMENT_alias", "experiment target",	"restriction enzyme", "restriction site";
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
    print biosamp[auxlib], libid, "chromatin", "HindIII", "AAGCTAGCTT";
}
