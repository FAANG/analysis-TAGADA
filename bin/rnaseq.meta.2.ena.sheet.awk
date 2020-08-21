# rnaseq.meta.2.ena.sheet.awk
# From 2 information files make a tsv file corresponding to the submission form experiment_ena sheet for rnaseq data
# Input files are
# - $rna = the experiment_ena file used for rnaseq in order to get the biosample id (sample_descriptor)
# - $meta = the metadata file for the biorep
# here the insert size is obtained from the same file as the biosample id so it will be difficult to make a single awk script
# of this kind for rnaseq and the other two techniques (atacseq and hic) unless we make a specific file for rnaseq insert sizes
# and change the present script ...
# !!! changed on June 25th 2018 not to have the same experiment alias as in SF's first submission and for this adding _2 !!!

# Here we need this information in the output
#############################################
# - SAMPLE_DESCRIPTOR = this is the biosample id, so take it from rnaseq as the one associated to a given combination of animal and tissue (in second column) (ex SAMEA1088320)
# - EXPERIMENT_alias = this is the experiment alias I made up and put in the run sheet (eg INRA_FRAGENCODE_20170716_2_RNASEQ_BOS_2361_CD4) !!! be careful to put _2 to differ from SF !!!
# - TITLE = Pig transcriptome profiling in cd4 cells by the FAANG pilot project FR-AgENCODE (change the tissue and species each time)
# - STUDY_REF = take from study sheet = INRA_FRAGENCODE_20180620_RNASEQ
# - DESIGN_DESCRIPTION = take from old SF = Paired-end RNA-Seq of polyA+ transcripts from cd4 cells on an Illumina HiSeq 3000 (change tissue each time)
# - LIBRARY_NAME = same as the experiment_alias in field no 2
# - LIBRARY_STRATEGY = RNA-Seq
# - LIBRARY_SOURCE = TRANSCRIPTOMIC
# - LIBRARY_SELECTION = Oligo-dT
# - LIBRARY_LAYOUT = PAIRED
# - NOMINAL_LENGTH = take from SF's old file
# - NOMINAL_SDEV = not available
# - LIBRARY_CONSTRUCTION_PROTOCOL = TruSeq Stranded mRNA Sample Preparation kit V2 Guide (with a reduction of the fragmentation time from 8 to 2 minutes)
# - PLATFORM = ILLUMINA
# - INSTRUMENT_MODEL = Illumina HiSeq 3000

# example
# cd /work/project/fragencode/data/ebi_submission/rnaseq
# rna=/work/project/fragencode/data/ebi_submission/rnaseq/by.SF/experiment_ena.tsv
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/rnaseq.meta.2.ena.sheet.awk 
# meta=/work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq_fastq_biorep_metadata.tsv
# awk -v fileRef=$rna -f $pgm $meta > rnaseq.experiment_ena.tsv

# $rna
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	TITLE	STUDY_REF	DESIGN_DESCRIPTION	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	LIBRARY_LAYOUTNOMINAL_LENGTH	NOMINAL_SDEV	LIBRARY_CONSTRUCTION_PROTOCOL	PLATFORM	INSTRUMENT_MODEL
# SAMEA1088320	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	Pig transcriptome profiling in cd4 cells by the FAANG pilot project FR-AgENCODE	INRA_FRAGENCODE_20180228_RNASEQ	Paired-end RNA-Seq of polyA+ transcripts from cd4 cells on an Illumina HiSeq 3000	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	RNA-Seq	TRANSCRIPTOMIC	Oligo-dT	PAIRED	515	not available TruSeq Stranded mRNA Sample Preparation kit V2 Guide (with a reduction of the fragmentation time from 8 to 2 minutes)	ILLUMINA	Illumina HiSeq 3000
# 1 (15 fields)
# 42 (60 fields)  *** column no 11 and animal and tissues in column no 2

# $meta
# path	species	labExpId	tissue	animal_old	animal_short	readno
# /work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq.bos_taurus.cd4.cattle1.R1.fastq.gz	bos_taurus	cd4.cattle1	cd4	FR4934530986	cattle1	1
# 83 (7 fields)
# for hic we had this
# path	species	tissue	animal_old	animal_short	animal	readno	Species
# /work/project/fragencode/data/ebi_submission/hic/sequences/hic.capra_hircus.goat2.R1.fastq.gz	capra_hircus	liver	Goat1	goat2	FR250020007415157	R1	Goat
# 25 (8 fields)

# output
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	TITLE	STUDY_REF	DESIGN_DESCRIPTION	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	LIBRARY_LAYOUTNOMINAL_LENGTH	NOMINAL_SDEV	LIBRARY_CONSTRUCTION_PROTOCOL	PLATFORM	INSTRUMENT_MODEL
# SAMEA1088354	INRA_FRAGENCODE_20170716_2_RNASEQ_BOS_0986_CD4	Cattle transcriptome profiling in cd4 cells by the FAANG pilot project FR-AgENCODE	INRA_FRAGENCODE_20180620_RNASEQ	Paired-end RNA-Seq of polyA+ transcripts from cd4 cells on an Illumina HiSeq 3000	INRA_FRAGENCODE_20170716_2_RNASEQ_BOS_0986_CD4	RNA-Seq	TRANSCRIPTOMIC	Oligo-dT	PAIRED	511	not available	TruSeq Stranded mRNA Sample Preparation kit V2 Guide (with a reduction of the fragmentation time from 8 to 2 minutes)	ILLUMINA	Illumina HiSeq 3000
# 1 (15 fields)
# 25 (60 fields)
# 16 (62 fields)


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
    
    # read the rnaseq ena file in order to associate a sample to a lib name and to get the insert size
    # be careful my experiment alias will differ from the one of SF by the addition of _2 
    while (getline < fileRef >0)
    {
	split($0,a,"\t");
	gsub(/20170716/,"20170716_2",a[2]);
	biosamp[a[2]]=a[1];
	size[a[2]]=a[11];
    }
    # print the header of the ena file we want
    print "SAMPLE_DESCRIPTOR", "EXPERIMENT_alias", "TITLE", "STUDY_REF", "DESIGN_DESCRIPTION",	"LIBRARY_NAME", "LIBRARY_STRATEGY", "LIBRARY_SOURCE", "LIBRARY_SELECTION", "LIBRARY_LAYOUT", "NOMINAL_LENGTH", "NOMINAL_SDEV", "LIBRARY_CONSTRUCTION_PROTOCOL", "PLATFORM", "INSTRUMENT_MODEL";
}

# Read the metadata file of fastq files but just read the /1 rows since exact same info in /2 and we need it once
$7==1{
    n=split($1,a,"/");
    split(a[n],b,".R1.fastq.gz");
    # the initial animal id
    n2=split($5,c,"");
    s="";
    for(i=n2-3; i<=n2; i++)
    {
	s=(s)(c[i]);
    }
    auxlib=corr[$2]"_"s"_"corr[$4];
    # !!! need to add _2 after the date to differ from previous submission by SF !!!
    libid="INRA_FRAGENCODE_20170716_2_RNASEQ_"auxlib;
    material=($4=="liver" ? "the liver tissue" : ($4" cells"));
    print biosamp[libid], libid, corr2[$2]" transcriptome profiling in "material" by the FAANG pilot project FR-AgENCODE", "INRA_FRAGENCODE_20180620_RNASEQ", "Paired-end RNA-Seq of polyA+ transcripts from "material" on an Illumina HiSeq 3000", libid, "RNA-Seq", "TRANSCRIPTOMIC", "Oligo-dT", "PAIRED", size[libid], "not available", "TruSeq Stranded mRNA Sample Preparation kit V2 Guide (with a reduction of the fragmentation time from 8 to 2 minutes)", "ILLUMINA", "Illumina HiSeq 3000";
}
