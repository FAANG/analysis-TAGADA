# ~/Awk/extract_coord_from_maf.awk
# - takes as input a .maf file of multiple alignments of one chromosome of a given genome 
#   against many other genomes, typically the UCSC hg19/GRCh37 46way vertebrate genome multiz 
#   alignments (maf format) which can be found here:
#   /seq/genomes/H.sapiens/golden_path_200902/multiz46way/maf/
#   (one file for each chromosome)
# - takes as parameters
#   * a genomic segment of this chromosome for this species (as in bed file = beg is 0 based and end is 1 based)
#   * the name of the assembly of this species
#   * the mode: whole if we want the whole alignment including the wanted segment
#               precise (default) if we only want the part of the alignment corresponding to the input segment 
# - outputs the multiple alignment corresponding to this segment, or the one including this segment in case mode=whole 
# THIS SCRIPT TAKES INTO ACCOUNT THE DELETIONS IN THE REF SPECIES (INDICATED AS -)
# THAT ARE BEFORE THE WANTED BEG AND WITHIN THE WANTED SEQ TO RETRIEVE THE SEQUENCES
# NOTES: 
# - for the moment there is no check whether the chr of the coordinates wanted is the same as the chr of the input file
#   they are thus assumed to be the same
# - the strand is not taken into account since all the alignments are on the + strand of the ref species
# NOTE that generally we do not give as input a complete maf file but only one that totally includes the
# wanted segment

# usage
# awk -v mode=whole|precise(default) -v assembly=ass -v coord=chr_beg_end_strand -f ~/Awk/extract_coord_from_maf.awk file.maf > file.almost.maf

# example
# awk -v mode=whole -v assembly=hg19 -v coord=chr21_14437564_14437573_+ -f ~/Awk/extract_coord_from_maf.awk /seq/genomes/H.sapiens/golden_path_200902/multiz46way/maf/chr21.maf > chr21_14437564_14437573_+.maf

# input
# /seq/genomes/H.sapiens/golden_path_200902/multiz46way/maf/chr21.maf
# ##maf version=1 scoring=autoMZ.v1
# a score=124903.000000
# s hg19.chr21     9411193 153 +  48129895 Gatcttcctccaaagaaattgtagttttcttctggcttagaggtagatcatcttggtccaatc--agac-----tgaaatgccttgaggctagatttcagtctttgtggcagctggtgaatttctagtttgccttttcagctagggattagctttttagg
# s ponAbe2.chrUn 66592936 153 -  72422247 gatctccctccaaagaaattgtagttttcttctggcttagaggtagatcctcttggtccaatc--agac-----tgaaatgccttgaggctagatttcagtttttgtggcagctggtgcatttctagtttgtgttttcagctacggattagctttttagg
# q ponAbe2.chrUn                          999999999999999999999999999999999999999999999999999999999999999--9999-----99999999999999999999999999999999999999999999999999999999999999999999999999999999999999
# i ponAbe2.chrUn N 0 C 0
# s panTro2.chrUn  5526850 152 +  58616431 GATCTCCCTCCAAAGAAATTGTAGTTTTCTTCTCGCTTAGAGGTAGATTCTCTTGGTCCAATC--AGGC-----TGAAATGCCTTGAGGCTAGATTTCAGTTTTTGTGGCAGCTGCTGCATTTCTAGTTTGCCTTT-CAGCTAGGGATTAGCTTTTTAGG
# q panTro2.chrUn                          999999999999999999999999999999999999999999999999999999999999999--9999-----99999999999999999999999999999999999999999999999999999999999999-99999999999999999999999
# i panTro2.chrUn N 0 C 0

# output
# chr21_14437564_14437573_+.maf
# s hg19.chr21                   14437514 59 +  48129895 AATCTTGGATGCAGTTCTTTCTTGTGAAAGAGCAAGGGAACTTAAAAAATATCCCTGTG
# s panTro2.chr18                 1512950 59 +  77261746 AATCTTGGATGCAGTTCTTTCTTGTGAAAGAGCAAGGGAACTTAAAAAATATCCCTGTG
# s ponAbe2.chr2b                 3258874 58 + 135000294 AATCTCAGATGCAG-TCTTTCTTGTGAAAGAGCAAGGGAAGTTAAAAAATATCCCTGTG
# s papHam1.Contig315758                0 59 -      1585 AATCTTGGATGCAGTTCTTCCTTGTGAAAAAGCAAGGGAACTTAAAAAATATCCCTGTG
# s vicPac1.scaffold_4577           69595 40 +     74583 AATTCAGGATCCAGTTGATTCATTC---AGATCAATAGAACTT----------------
# s equCab2.chr29                10923106 59 +  33672925 AATCCAGAATGCAGTTCTTTCATATGAAAGATTAATAGAACTTAAAAAAAATCACTGTG


BEGIN{
    if(mode=="")   # by default we want mode="precise" toextract maf corresponding to given segment
    {
	mode="precise";
    }
    split(coord,coordarr,"_");  # chr_beg_end_strand as in bed file
    searched_chr=coordarr[1];
    searched_beg=coordarr[2]; 
    searched_end=coordarr[3]-1;  # to be truly 0-based
}

{
    if (($1=="s")&&($2==((assembly)(".")(searched_chr)))&&((fincluded(searched_beg,searched_end,$3,$3+$4-1)==1)))
    {
	if(mode=="precise") # then we want only the MA corresponding to the segment
	{
	    split($7,arr,"");

	    # Compute localbeg taking into account the dashes in the ref seq
	    seennt=0;  # non dash - signs that are seen in ref sequence before the begining of the subsequence we want
	    k=1;
	    while((arr[k]!="")&&(seennt<(searched_beg-$3+1)))
	    {
		if(arr[k]!="-")
		{
		    seennt++;
		}
		k++;
	    }
	    localbeg=k-1;  # it is the begining of the sequence we are looking for in the alignment (= with the -)
	    
	    # Compute wantedlength taking into account the dashes in the ref seq
	    seennt=0;  # non dash - signs that are seen in ref sequence between the beg and the end of the seq we want
	    l=k-1;
	    while((arr[l]!="")&&(seennt<(searched_end-searched_beg+1)))
	    {
		if(arr[l]!="-")
		{
		    seennt++;
		}
		l++;
	    }
	    localend=l-1; # it is the end of the sequence we are looking for in the alignment (= with the -) 

	    toprintaux=substr($7,localbeg,(localend-localbeg+1));
	    # the length of an alignment in a species is the number of nt different from dash in this species, 
	    # so here = searched_end-searched_beg+1
	    # the begining of an alignment in a species is the begining when we remove the dashes in this species,
	    # so here = searched_beg
	    toprint=(($1)("\t")($2)("\t")(searched_beg)("\t")(searched_end-searched_beg+1)("\t")($5)("\t")($6)("\t")(toprintaux));
	}
	else   # then we want all the MA including the segment
	{
	    toprint=$0;
	}
	print toprint; 
	printflag=1;
    }
    else
    {
	if($1=="a")
	{
	    printflag=0;
	}
	if(printflag==1&&$1=="s")
	{
	    if(mode=="precise")
	    {
		split($7,arr,"");
		
		# Compute begtoprint taking into account the dashes in the seq of the current species
		seennt=0;  # non dash - signs that are seen in sequence of the current species before the beg of seq we want
		for(i=1; i<localbeg; i++)
		{
		    if(arr[i]!="-")
		    {
			seennt++;
		    }
		}
		begtoprint=$3+seennt;  # this is the start in the sequence of the current species of the alignment provided
		
		# Compute lengthtoprint taking into account the dashes in the seq of the current species
		seennt=0; # non dash - signs that are seen in sequence of the current species between the beg and the end of the seq we want
		toprintaux="";
		for(i=localbeg; i<=localend; i++)
		{
		    if(arr[i]!="-")
		    {
			seennt++;
		    }
		    toprintaux=(toprintaux)(arr[i]); 
		}
		toprint=(($1)("\t")($2)("\t")(begtoprint)("\t")(seennt)("\t")($5)("\t")($6)("\t")(toprintaux));
	    }
	    else
	    {
		toprint=$0;
	    }
	    print toprint;
	}
    }
} 


function fincluded(beg1, end1, beg2, end2)
{
    if((beg2<=beg1)&&(end1<=end2))
	return 1;
    else
	return 0;
}