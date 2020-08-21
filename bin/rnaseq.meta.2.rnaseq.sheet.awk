# rnaseq.meta.2.rnaseq.sheet.awk
# From 3 information files make a tsv file corresponding to the submission form rnaseq sheet for rnaseq data
# Input files are
# - $rna1 = the experiment_ena file used for an old rnaseq submission by SF for us to get the biosample id (sample_descriptor)
# - $rna2 = the rna-seq file used for and old rnaseq submission by SF for us to get the RIN number
# - $meta = the metadata file for the biorep
# !!! changed on June 25th 2018 not to have the same experiment alias as in SF's first submission and for this adding _2 !!!

# Here we need this information in the output
#############################################
# - SAMPLE_DESCRIPTOR (SAMEA1088320) = OK I have from run sheet and can make the same way
# - EXPERIMENT_alias (INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4) = OK I have from run sheet and can make the same way
#   !!! be careful has to differ from SF's first one by adding _2 after the date !!!
# - experiment target = polyA RNA
# - rna preparation 3' adapter ligation protocol = ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf
# - rna preparation 5' adapter ligation protocol = ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf
# - library generation pcr product isolation protocol = ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf
# - preparation reverse transcription protocol = ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf
# - library generation protocol = ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf
# - read strand = mate 2 sense
# - rna purity - 260:280 ratio = empty
# - rna purity - 260:230 ratio = empty
# - rna integrity number = from /work/project/fragencode/data/ebi_submission/rnaseq/by.SF/RNA-Seq.tsv  (12th field)

# example
# cd /work/project/fragencode/data/ebi_submission/rnaseq
# rna1=/work/project/fragencode/data/ebi_submission/rnaseq/by.SF/experiment_ena.tsv
# rna2=/work/project/fragencode/data/ebi_submission/rnaseq/by.SF/RNA-Seq.tsv
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/rnaseq.meta.2.rnaseq.sheet.awk
# meta=/work/project/fragencode/data/metadata/atacseq/atacseq_combinedreadfile_metadata.tsv
# awk -v fileRef=$rna -f $pgm $meta > rnaseq.rnaseq.tsv

# $rna1
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	TITLE	STUDY_REF	DESIGN_DESCRIPTION	LIBRARY_NAME	LIBRARY_STRATEGY	LIBRARY_SOURCE	LIBRARY_SELECTION	LIBRARY_LAYOUTNOMINAL_LENGTH	NOMINAL_SDEV	LIBRARY_CONSTRUCTION_PROTOCOL	PLATFORM	INSTRUMENT_MODEL
# SAMEA1088320	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	Pig transcriptome profiling in cd4 cells by the FAANG pilot project FR-AgENCODE	INRA_FRAGENCODE_20180228_RNASEQ	Paired-end RNA-Seq of polyA+ transcripts from cd4 cells on an Illumina HiSeq 3000	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	RNA-Seq	TRANSCRIPTOMIC	Oligo-dT	PAIRED	515	not available TruSeq Stranded mRNA Sample Preparation kit V2 Guide (with a reduction of the fragmentation time from 8 to 2 minutes)	ILLUMINA	Illumina HiSeq 3000
# 1 (15 fields)
# 42 (60 fields)

# $rna2
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	experiment target	rna preparation 3' adapter ligation protocol	rna preparation 5' adapter ligation protocol	library generation pcr product isolation protocol	preparation reverse transcription protocol	library generation protocol	read strand	rna purity - 260:280 ratio	rna purity - 260:230 ratio	rna integrity number
# SAMEA1088320	INRA_FRAGENCODE_20170716_RNASEQ_SUS_1346_CD4	polyA RNA	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	mate 2 sense	8.9
# 42 (13 fields)
# 1 (44 fields)  *** last column = 12th one = the rin number

# $meta
# path	species	tissue	animal_old	animal_short	animal	readno	Species
# /work/project/fragencode/data/ebi_submission/rnaseq/sequences/rnaseq.capra_hircus.goat2.R1.fastq.gz	capra_hircus	liver	Goat1	goat2	FR250020007415157	R1	Goat
# 25 (8 fields)

# output
# SAMPLE_DESCRIPTOR	EXPERIMENT_alias	experiment target	rna preparation 3' adapter ligation protocol	rna preparation 5' adapter ligation protocol	library generation pcr product isolation protocol	preparation reverse transcription protocol	library generation protocol	read strand	rna purity - 260:280 ratio	rna purity - 260:230 ratio	rna integrity number
# SAMEA1088354	INRA_FRAGENCODE_20170716_2_RNASEQ_BOS_0986_CD4	polyA RNA	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf	mate 2 sense	9.2
# 41 (13 fields)
# 1 (44 fields)

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
    # read the old rnaseq file in order to associate a sample to its RIN number
    # !!! need to add _2 in the new experiment alias !!!
    while (getline < fileRef2 >0)
    {
	split($0,a,"\t");
	gsub(/20170716/,"20170716_2",a[2]);
	rin[a[2]]=a[12];
    }

    # print the header of the faang experiment file we want
    print "SAMPLE_DESCRIPTOR", "EXPERIMENT_alias", "experiment target", "rna preparation 3' adapter ligation protocol", "rna preparation 5' adapter ligation protocol", "library generation pcr product isolation protocol", "preparation reverse transcription protocol", "library generation protocol", "read strand", "rna purity - 260:280 ratio", "rna purity - 260:230 ratio", "rna integrity number";
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
    proto="ftp://ftp.faang.ebi.ac.uk/ftp/protocols/assays/INRA_SOP_mRNA-seq_long_protocol_20160908.pdf";
    print biosamp[libid], libid, "polyA RNA", proto, proto, proto, proto, proto, "mate 2 sense", "", "", rin[libid];
}
