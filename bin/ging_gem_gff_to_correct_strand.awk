# ~/Awk/ging_gem_gff_to_correct_strand.awk
# takes as input a gff file mapped by gem out of pe stranded ging rnaseq
# and put the correct strand to the mappings (read /2 was ok but read /1 was
# on the opposite strand of the transcript so need to reverse its strand)


# Assumes an input file like
#############################
# /users/rg/projects/NGS/Projects/ENCODE/hg19main/GingerasCarrie/003WC/work/genome.gtf
# chr12	genome/LID16629_FC61U2UAAXX_1_1.gem.map.gz	paired_read	49330590	49330665	.	+	.	read_id "TUPAC_0006:1:100:10000:10161#0/1"; mismatches "0"; qualities "0/0"; matches "1:0:0";
# where the read id is in $10 

{
    split($10,a,"\""); 
    split(a[2],b,"/"); 
    if(b[2]==1)
    {
	$7=($7=="+" ? "-" : "+");
    } 
    print $0
}