# ~sdjebali/Awk/make_introns.awk

# starts from an ordered set of exons from transcripts in gtf (order = tr, beg, end) (this file can have additional items)
# and make the corresponding intron file in gtf too.
# Note: introns are made within a given transcript
# Note: assumes the exons are on the same chromosome and strand, which is the case for the annotation but not for the chimeras

# example
# sort -k12,12 -k4,4n -k5,5n schizosaccharomyces_japonicus_yfs275_1_exons.gtf | awk -v fldgn=10 -v fldtr=12 -f ~/Awk/make_introns.awk > schizosaccharomyces_japonicus_yfs275_1_introns.gtf

# schizosaccharomyces_japonicus_yfs275_1_exons.gtf
# supercont1.1  SJ1_FINAL_CALLGENES_1   exon    531     721     .       -       .       gene_id "SJAG_00001"; transcript_id "SJAT_00001";

# schizosaccharomyces_japonicus_yfs275_1_introns.gtf
# supercont1.1  SJ1_FINAL_CALLGENES_1   intron  722     916     .       -       .       gene_id "SJAG_00001"; transcript_id "SJAT_00001";

BEGIN{
    if(fldgn=="")
	fldgn=10;
    if(fldtr=="")
	fldtr=12;
}

$3=="exon"{
    seen[$fldtr]++;
    if(seen[$fldtr]==1)  # if first time we see the current transcript
    {
	prev_end=$5;
	for(i=13; i<=(NF-2); i+=2)
	{
	    if($i!~/exon/)
		info[$fldtr]=(info[$fldtr])($i)(" ")($(i+1))(" ");
	}
	info[$fldtr]=(info[$fldtr])($i)(" ")($(i+1));
    }
    else   # we print the current intron = the one which ends in the current exon and we remember the current exon end for the next intron
    {
	printf("%s\t%s\tintron\t%i\t%i\t.\t%s\t.\tgene_id %s transcript_id %s %s\n", $1, $2, prev_end+1, $4-1, $7, $fldgn, $fldtr, info[$fldtr]);
	prev_end=$5;
    }
}

