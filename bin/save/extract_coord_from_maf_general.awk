# ~/Awk/extract_coord_from_maf_general.awk
# same as ~/Awk/extract_coord_from_maf.awk except that it does not require the segment to retrieve
# be totally included in one of the multiple alignments of the maf file: it just needs the wanted 
# segment to overlap the current multiple alignment of the maf file.
# A MAIN DIFFERENCE WITH ~/Awk/extract_coord_from_maf.awk IS THAT HERE THE WANTED SEGMENTS
# IS ALSO INDEXED (at the end by _no) by the coordn of the MA region that we are currently
# extracting it from (and that overlaps it), in order to index the output file in order not
# to erase one by the other in case the wanted segment is overlapping several MA regions.
# But actually this information is not used in this script but just to index the output file
# of this script.
# NOTE that generally we do not give as input a complete maf file but only one that overlaps the
# wanted segment

# Notes:
# - for the moment there is no check whether the chr of the coordinates wanted is the same as the chr of the input file
#   they are thus assumed to be the same
# - the strand is not taken into account since all the alignments are on the + strand of the ref species

# usage
# awk -v mode=whole|precise(default) -v assembly=ass -v coord=chr_beg_end_strand_nomareg -f ~/Awk/extract_coord_from_maf_general.awk file.maf > file.almost.maf

# example
# awk -v mode=whole -v assembly=hg19 -v coord=chr21_14437564_14437573_+_1 -f ~/Awk/extract_coord_from_maf_general.awk /seq/genomes/H.sapiens/golden_path_200902/multiz46way/maf/chr21.maf > chr21_14437564_14437573_+_1.maf

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
    if (($1=="s")&&($2==((assembly)(".")(searched_chr)))&&((foverlap(searched_beg,searched_end,$3,$3+$4-1)==1)))
    {
	if(mode=="precise")  # then we want only the MA or part of MA corresponding to the wanted segment
	{
	    # here the actual searched_beg-searched_end segment in the intersection between
	    # the wanted segment searched_beg-searched_end and the current MA region
	    actual_searched_beg=max(searched_beg,$3);
	    actual_searched_end=min(searched_end,$3+$4-1);
	    split($7,arr,"");

	    # Compute localbeg taking into account the dashes in the ref seq
	    seennt=0;  # non dash - signs that are seen in ref sequence before the begining of the subsequence we want
	    k=1;
	    while((arr[k]!="")&&(seennt<(actual_searched_beg-$3+1)))
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
	    while((arr[l]!="")&&(seennt<(actual_searched_end-actual_searched_beg+1)))
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
	    # so here = actual_searched_end-actual_searched_beg+1
	    # the begining of an alignment in a species is the begining when we remove the dashes in this species,
	    # so here = searched_beg
	    toprint=(($1)("\t")($2)("\t")(actual_searched_beg)("\t")(actual_searched_end-actual_searched_beg+1)("\t")($5)("\t")($6)("\t")(toprintaux));
	}
	else   # then we want all the current MA overlapping the wanted segment
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


function foverlap(beg1,end1,beg2,end2)
{
  return ((end1>=beg2)&&(beg1<=end2)) ? 1 : 0;
}

function min(x,y)
{
  return x <= y ? x : y;
}

function max(x,y)
{
  return x >= y ? x : y;
}
