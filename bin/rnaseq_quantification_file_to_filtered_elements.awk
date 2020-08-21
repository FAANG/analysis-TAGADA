# ~sdjebali/Awk/rnaseq_quantification_file_to_filtered_elements.awk
# = Essential action for element (exon, splice junction, transcript, gene) filtering underlying the summary stat shell script called
#   ~sdjebali/bin/code_for_features_detected_cat2.sh
# Its complication comes from the fact that the script needs to be general enough to deal with
# - any quantification file either gff or bed
# - for gff files files with the element located anywhere in the attribute list of a gff file
# - filtering being done either by values higher than a threshold (max) as for rpkm or lower than a threshold (min) as for iIDR
# - ...

# example
#########
# chr3    GencodeV10      transcript      33537737        33759848        0.000000        -       .       gene_id "ENSG00000163539.11"; transcript_id "ENST00000399362.4"; RPKM1 "0.000000"; RPKM2 "0.000000"; iIDR "NA";
# cd ~sdjebali/ENCODE_AWG/Analyses/SummaryStatistics/Promocells/Annot_Feat_Detected
# zcat /users/rg/projects/encode/scaling_up/whole_genome/encode_DCC_mirror/wgEncodeCshlLongRnaSeq/releaseLatest/wgEncodeCshlLongRnaSeqHwpCellTotalTranscriptGencV10.gtf.gz | awk -v name=transcript -v expid=LID47254-LID47255 -v thresh=0.1 -v ninFile=iIDR -v filetype=gff -v keyelt=transcript_id -v minormax=max -f ~sdjebali/Awk/rnaseq_quantification_file_to_filtered_elements.awk | sort | uniq > Individual/transcript/transcript_detected_in_LID47254-LID47255.tmp.txt 
# ENST00000000233.5
# 33487 (1 fields)

# Note: this script can be made faster if assuming that for a gff file both the element id and the expression value on which to make the filtering
# will always be in the same column in the whole file (here for each row I am looking for those elements everywhere in the attribute field)

{
    if(filetype=="bed")
    {
	elt=$1"_"$2"_"$3"_"$6; 
	if(((minormax=="max")&&($ninFile<=thresh))||((minormax=="min")&&($ninFile>=thresh)))
	{
	    print elt;
	}
    } 
    else
    {
	if($3==name)
	{
	    k=9; 
	    while(k<=(NF-1))
	    {
		if(keyelt!=".")
		{
		    if($k==keyelt)
		    {
			gsub(/;/,"",$(k+1)); 
			gsub(/"/,"",$(k+1)); 
			elt=$(k+1);
		    }
		}
		else
		{
		    elt=$1"_"$4"_"$5"_"$7;
		}
		
		if(($k==ninFile)&&($(k+1)!~/NA/))  # typically for iIDR or for RPKM values in atribute field of a gff file
		{
		    gsub(/;/,"",$(k+1)); 
		    gsub(/"/,"",$(k+1)); 
		    if(((minormax=="max")&&($(k+1)<=thresh))||((minormax=="min")&&($(k+1)>=thresh)))
		    {
			print elt;
		    }
		}
		k+=2;
	    }
	}
    }
}