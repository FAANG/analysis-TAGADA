# ~sdjebali/Awk/intersectBed_pestrandedbam_with_eltgff_to_gff_with_s_as_mappings.awk
# takes as input the result of an intersectBed with -abam $mappings and -b element file in gff (eg genes)
# as well as the element gff file as fileRef file and the mate strand of the data, and produces the same 
# gff file but with addition of the number of mapping that overlap on the sense and the number of mappings
# that overlap the elements on the antisense strand
# This script is used in ~sdjebali/bin/add_sense_antisense_read_counts_to_segments_frombam.sh

# usage
#######
# GFF2GFF=~sdjebali/Awk/gff2gff.awk
# intersectBed -abam $mappings -b $outdir/$b12.ext$toextend.$b22.gff -bed -f 1 -wo | awk -v fileRef=$outdir/$b12.ext$toextend.$b22.gff -v mate_strand=$mate_strand -f intersectBed_pestrandedbam_with_eltgff_to_gff_with_s_as_mappings.awk | awk -f $GFF2GFF > $outdir/$b12.ext$toextend.withreadcount.of.$b22.gff


BEGIN{
    while(getline < fileRef >0)
    {
	line[$1"_"$4"_"$5"_"$7"_"$10]=$0;
    }
} 

{
    coord=$13"_"$16"_"$17"_"$19"_"$22;
    if((mate_strand=="MATE2_SENSE")||(mate_strand=="MATE1_SENSE")||(mate_strand=="MATE_STRAND_CSHL"))  # only possible values to produce anything
    {
	split($4,a,"/");  # only used for MATE2_SENSE and MATE1_SENSE
	if($6==$19)       # same strand between the mapping and the element
	{
	    if(((mate_strand=="MATE2_SENSE")&&(a[2]==2))||((mate_strand=="MATE1_SENSE")&&(a[2]==1))||(mate_strand=="MATE_STRAND_CSHL"))
	    {
		nbsreads[coord]++;
	    }
	    else
	    {
		if(((mate_strand=="MATE2_SENSE")&&(a[2]==1))||((mate_strand=="MATE1_SENSE")&&(a[2]==2))||(mate_strand=="MATE_STRAND_CSHL"))
		{
		    nbasreads[coord]++;
		}
	    }	
	}
	else             # diff strand between the mapping and the element, in the case only MATE2_SENSE and MATE1_SENSE are considered, not MATE_STRAND_CSHL
	{
	    if(((mate_strand=="MATE2_SENSE")&&(a[2]==2))||((mate_strand=="MATE1_SENSE")&&(a[2]==1)))
	    {
		nbasreads[coord]++;
	    }
	    else
	    {
		if(((mate_strand=="MATE2_SENSE")&&(a[2]==1))||((mate_strand=="MATE1_SENSE")&&(a[2]==2)))
		{
		    nbsreads[coord]++;
		}
	    }
	}
    }
}

END{
    for(e in line)
    {
	print line[e], "nbsreads", "\""(nbsreads[e]!="" ? nbsreads[e] : 0)"\"\;", "nbasreads", "\""(nbasreads[e]!="" ? nbasreads[e] : 0)"\"\;";
    }
} 