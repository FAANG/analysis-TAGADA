# make_multiple_maf_out_of_one_maf.awk
# takes as input a maf file of a complete chr as in /seq/genomes/H.sapiens/golden_path_200902/multiz46way/maf
# and as fileRef a bed file corresponding to coord of ALIGNMENTS already present in this file with PRECISE BEG AND STRAND
# and provides as many maf files as there are unique lines in the bed file, corresponding to the blocks
# of corresponding alignments in the huge input maf file. Those file will be indexed by the bed coord of those
# alignments
# Be careful: writes in the MaReg directory (which should thus exist, and means multiple alignment region)
# I have removed -v searched_strand=+ and &&($5==searched_str) down since MA alignments are always
# on the + strand of the ref species

# awk -v refspecies=hg19.chr1 -v fileRef=gencode.v3d.annotation.GRCh37.noncanintrons_24mer_don_acc_separated_nr_DS_coord_chr1_plus_intronsmore24_noncan_MA46species_includingthem.bed -f ~/Awk/make_multiple_maf_out_of_one_maf.awk /seq/genomes/H.sapiens/golden_path_200902/multiz46way/maf/chr1.maf 

# gencode.v3d.annotation.GRCh37.noncanintrons_24mer_don_acc_separated_nr_DS_coord_chr1_plus_intronsmore24_noncan_MA46species_includingthem.bed
# chr1    101705454       101705522       .       .       +

# /seq/genomes/H.sapiens/golden_path_200902/multiz46way/maf/chr1.maf 
# ##maf version=1 scoring=autoMZ.v1
# a score=34237.000000
# s hg19.chr1     10917 479 + 249250621 gagaggcgcaccgcgccggcgcaggcgcagagacacatgctagcgcgtccaggggtggaggcgtggcgcaggcgcagagacgcaagcctacgggcgggggttgggggggcgtgtgttgcaggagcaaagtcgcacggcgccgggctggggcggggggag
# ggtggcgccgtgcacgcgcagaaactcacgtcacggtggcgcggcgcagagacgggtagaacctcagtaatccgaaaagccgggatcgaccgccccttgcttgcagccgggcactacaggacccgcttgctcacggtgctgtgccagggcgccccctgctggcgactagggcaactgcagggctctcttgcttagag
# tggtggccagcgccccctgctggcgccggggcactgcagggccctcttgcttactgtatagtggtggcacgccgcctgctggcagctagggacattgcagggtcctcttgctcaaggtgtagt
# s panTro2.chr15 13606 455 - 100063422 gagaggcgcaccgcgccggcgcag------agacacatactagcgcgtcctgggg-ggaggtgcggcgctgtcgccgagac-cacgcctatgggcgggggttgcgggg-cgcgtggtgcaggagcaaagtcgcacggagcctgtctgggg-----gcag
# gtgggctccgtgcaggcgcagaaacgcacgtcgcggcggcgcggcgcagagacgggtggaacctcagtaatcagaaaagccgggctggaccgccccctgcttgcagccgggcactacaggacccgcttgctcacggtgctctgccagtgcgccccctgctggcaactagggcaactgcagggctctcttgcttagag
# tggtggccagcgggggctgctggctccggggcactgcagggccctcttgcttactgtatagtggtggcacgccgcctgctggcagctagggacattgcagggtcctcttgctc----------
# q panTro2.chr15                       999999999999999999999999------9999999999999999999999999-9999999999999999999999999-99999999999999999999999999-99999999999999999999999999999999999999999-----9999
# 99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999
# 99999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999999----------
# i panTro2.chr15 N 0 C 0

BEGIN{
    while (getline < fileRef >0)
    {
	reg[$2]=$1"_"$2"_"$3"_"$6;
    }
}

{
    if(($1=="s")&&($2==refspecies)&&(reg[$3]!=""))
    {
	# print "I am in the first if (s of the ref species)";
	# print;
	currregion=reg[$3];
	print $0 > "MaReg/"(currregion)"_precise.maf"; 
	printflag=1;
    }
    else
    {
	if($1=="a")
	{
	    # print "I am in the first if of the else (alignment a)";
	    # print;
	    printflag=0;
	}
	if(printflag==1&&$1=="s")
	{
	    # print "I am in the second if of the else (s of other species)";
	    # print;
	    print $0 > "MaReg/"(currregion)"_precise.maf";
	}
    }
}
