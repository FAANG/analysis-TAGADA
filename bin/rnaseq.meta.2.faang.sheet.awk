# atac.meta.2.faang.sheet.awk
# From 3 information files make a tsv file corresponding to the submission form experiment_faang sheet for hic data
# Input files are
# - $rna1 = the experiment_ena file used for rnaseq by SF in the past in order to get the biosample id (sample_descriptor)
# - $rna2 = the faang experiment file used for rnaseq by SF in the past in order to get the sampling to preparation interval
# - $meta = the metadata file for the biorep
# An improvement would be to have a single such script for hic and atacseq but for this we would need to
# - pass the field no for R1 and R2 or get it from metadata header (9 for atac, 7 for hic)
# - pass the type of attribute we have for R1 and R2 (1 and 2 for atac, R1 and R2 for hic)
# - pass the way the files are ending (_R1.fastq.gz for atac, .R1.fastq.gz for hic)
# - pass the field no for animal init id or get it from metadata  header (6 for both here but only by chance)
# - pass the field no for species or get it from metadata  header (2 for both but again by chance)
# - pass the field no for tissues or get it from metadata  header (4 for atac and 3 for hic)
# !!! changed on June 25th 2018 not to have the same experiment alias as in SF's first submission and for this adding _2 !!!

# We need to put this information in the output tsv file
########################################################
# - SAMPLE_DESCRIPTOR (SAMEA1088320) = OK I have from run sheet and can make the same way
# - EXPERIMENT_alias (INRA_FRAGENCODE_20170716_2_RNASEQ_SUS_1346_CD4) = OK I have from run sheet and can make the same way
#   !!! be careful has to differ from SF's first one by adding _2 after the date !!!
# - assay type  = transcription profiling by high throughput sequencing
# - sample storage = 'frozen, -70 freezer' 
# - sample storage processing = 'cryopreservation, other'
# - sampling to preparation interval = get from /work/project/fragencode/data/ebi_submission/rnaseq/by.SF/experimentfaang_oldtab.tab
# - Units = weeks
# - experimental protocol = ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf
# - extraction protocol	= ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_RNA_extraction_20160504.pdf
# - library preparation location : Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France
# - library preparation location longitude  = 1.499654
# - Units (decimal degrees) = decimal degrees
# - library preparation location latitude = 43.527571
# - Units (decimal degrees) = decimal degrees
# - library preparation date = 2016-07-08T00:00:00 for liver and 2016-08-12T00:00:00 for cd4/8
# - Units (YYYY-MM-DD): YYYY-MM-DD		
# - sequencing location = Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France
# - sequencing location longitude = 1.499654
# - Units = decimal degrees
# - sequencing location latitude = 43.527571
# - Units = decimal degrees
# - sequencing date = 2016-07-08T00:00:00 for liver and 2016-08-12T00:00:00 for cd4/cd8
# - Units (YYYY-MM-DD): YYYY-MM-DD

# example
# cd /work/project/fragencode/data/ebi_submission/hic
# rna1=/work/project/fragencode/data/ebi_submission/rnaseq/by.SF/experiment_ena.tsv
# rna2=/work/project/fragencode/data/ebi_submission/rnaseq/by.SF/experimentfaang_oldtab.tab
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/hic.meta.2.faang.sheet.awk
# meta=/work/project/fragencode/data/ebi_submission/hic/sequences/hic_fastq_biorep_metadata.tsv
# awk -v fileRef1=$rna1 -v fileRef2=$rna2 -f $pgm $meta > hic.experiment_faang.tsv

# $rna1
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	TITLE	STUDY_REF	DESIGN_DESCRIPTION	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	LIBRARY_LAYOUTNOMINAL_LENGTH	NOMINAL_SDEV	LIBRARY_CONSTRUCTION_PROTOCOL	PLATFORM	INSTRUMENT_MODEL
# SAMEA1088320	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	Pig transcriptome profiling in cd4 cells by the FAANG pilot project FR-AgENCODE	INRA_FRAGENCODE_20180228_RNASEQ	Paired-end RNA-Seq of polyA+ transcripts from cd4 cells on an Illumina HiSeq 3000	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	RNA-Seq	TRANSCRIPTOMIC	Oligo-dT	PAIRED	515	not available TruSeq Stranded mRNA Sample Preparation kit V2 Guide (with a reduction of the fragmentation time from 8 to 2 minutes)	ILLUMINA	Illumina HiSeq 3000
# 1 (15 fields)
# 42 (60 fields)

# $rna2
# SAMEA1088320	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	transcription profiling by high throughput sequencing	frozen, -70 freezer	cryopreservation, other	13	weeks	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_RNA_extraction_20160504.pdf	Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France	1.499654	decimal degrees	43.527571	decimal degrees	2016-08-12	YYYY-MM-DD	Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France	1.499654	decimal degrees	43.527571	decimal degrees	2016-08-12T00:00:00	YYYY-MM-DD			
# SAMEA1088325	INRA_FRAGENCODE_20170716_RNASEQ_SUS_0886_CD4	transcription profiling by high throughput sequencing	frozen, -70 freezer	cryopreservation, other	21	weeks	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_RNA_extraction_20160504.pdf	Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France	1.499654	decimal degrees	43.527571	decimal degrees	2016-08-12	YYYY-MM-DD	Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France	1.499654	decimal degrees	43.527571	decimal degrees	2016-08-12T00:00:00	YYYY-MM-DD			
# 42 (47 fields)

# $meta
# path	species	labExpId	tissue	animal_old	animal_short	readno
# /work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq.bos_taurus.cd4.cattle1.R1.fastq.gz	bos_taurus	cd4.cattle1	cd4	FR4934530986	cattle1	1
# 83 (7 fields)
# for hic we had this
# path	species	tissue	animal_old	animal_short	animal	readno	Species
# /work/project/fragencode/data/ebi_submission/hic/sequences/hic.capra_hircus.goat2.R1.fastq.gz	capra_hircus	liver	Goat1	goat2	FR250020007415157	R1	Goat
# 25 (8 fields)

# output
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	assay type	sample storage	sample storage processing	sampling to preparation interval	Units	experimental protocol	extraction protocol	library preparation location	library preparation location longitude	Units	library preparation location latitude	Units	library preparation date	Units	sequencing location	sequencing location longitude	Units	sequencing location latitude	Units	sequencing date	Units
# SAMEA1088354	INRA_FRAGENCODE_20170716_2_RNASEQ_BOS_0986_CD4	transcription profiling by high throughput sequencing	frozen, -70 freezer	cryopreservation, other	33	weeks	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_RNA_extraction_20160504.pdf	Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France	1.499654	decimal degrees	43.527571	decimal degrees	2016-08-12T00:00:00	YYYY-MM-DD	Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France	1.499654	decimal degrees	43.527571	decimal degrees	2016-08-12T00:00:00	YYYY-MM-DD
# 41 (47 fields)
# 1 (48 fields)


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
    
    # read the old rnaseq ena file in order to associate a sample to its  biosample id
    # !!! need to add _2 in the new experiment alias !!!
    while (getline < fileRef1 >0)
    {
	split($0,a,"\t");
	gsub(/20170716/,"20170716_2",a[2]);
	biosamp[a[2]]=a[1];
    }
    # read the old rnaseq faang experiment file in order to associate a sample to its sampling to preparation time interval
    # !!! need to add _2 in the new experiment alias !!!
    while (getline < fileRef2 >0)
    {
	split($0,a,"\t");
	gsub(/20170716/,"20170716_2",a[2]);
	preptime[a[2]]=a[6];
    }

    # print the header of the faang experiment file we want
    print "SAMPLE_DESCRIPTOR", "EXPERIMENT_alias", "assay type", "sample storage", "sample storage processing", "sampling to preparation interval", "Units", "experimental protocol", "extraction protocol", "library preparation location", "library preparation location longitude", "Units", "library preparation location latitude", "Units", "library preparation date", "Units", "sequencing location", "sequencing location longitude", "Units", "sequencing location latitude", "Units", "sequencing date", "Units";
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
    print biosamp[libid], libid, "transcription profiling by high throughput sequencing", "frozen, -70 freezer", "cryopreservation, other", preptime[libid], "weeks", "ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf", "ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_RNA_extraction_20160504.pdf", "Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France", "1.499654", "decimal degrees", "43.527571", "decimal degrees", ($4=="liver" ? "2016-07-08T00:00:00" : "2016-08-12T00:00:00"), "YYYY-MM-DD", "Plateforme Genomique de Toulouse, INRA Castanet-Tolosan, France", "1.499654", "decimal degrees", "43.527571", "decimal degrees", ($4=="liver" ? "2016-07-08T00:00:00" : "2016-08-12T00:00:00"), "YYYY-MM-DD";
}
