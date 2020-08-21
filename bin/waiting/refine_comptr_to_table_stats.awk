# ~/Awk/refine_comptr_to_table_stats.awk
# takes as input the output of ~sdjebali/bin/refine_comptr_output.sh which means a tsv file with for each transcript a specific and a broad class
# and outputs a tsv file good for summary stat table with all numbers and % of total in the right order, meaning
# - all tr
# - annot tr and then split into exact, inclusion, other_spliced, other_monoex
# - extension tr and then split into annot_compatible, other_spliced, other_monoex
# - intergenic tr and then split into spliced and monoex
# - antisense tr and then split into spliced and monoex

# example
# cd ~sdjebali/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Stringtie/WithAnnot
# awk -v lid=stringtie -f ~/Awk/refine_comptr_to_table_stats.awk K562_polya+_wcell_biorep1_star2_stringtie1.0.4_complete_comp_refinedclass_nbex.tsv > K562_polya+_wcell_biorep1_star2_stringtie1.0.4_complete_comp_refinedclass_nbex_sumstat.tsv

# input
# trid  comptrclass     annottrlist     refinedtrclass  nbex
# tSTRG.80746.1   Monoexonic      .       intergenic      n1
# 163350 (5 fields)

# output
# set     total   annot_nb        annot_% exact_nb        exact_% inclusion_nb    inclusion_%     other_spliced_nb        other_spliced_% other_monoex_nb other_monoex_%  extension_nb       extension_%     annot_compatible_nb     annot_compatible_%      other_spliced_nb        other_spliced_% other_monoex_nb other_monoex_%  intergenic_nb   intergenic_%    spliced_nb spliced_%       monoex_nb       monoex_%        antisense_nb    antisense_%     spliced_nb      spliced_%       monoex_nb       monoex_%        unstr_nb        unstr_pcent
# stringtie       163350  93868   57.4643 58810   36.0024 1247    0.763391        1408    0.861953        32403   19.8365 10440   6.39118 8672    5.30885 1767    1.08173 1       0.0
# 00612182        35080   21.4754 1261    0.771962        33819   20.7034 23961   14.6685 3355    2.05387 20606   12.6146 0       0
# 2 (34 fields)

# then to print in 2 tables:
############################
# awk '{s=""; for(i=1; i<=20; i++){s=(s)($i)("\t")}print s}' refine_comptr_output_fortable.tsv
# set     total   annot_nb        annot_% exact_nb        exact_% inclusion_nb    inclusion_%     other_spliced_nb        other_spliced_% other_monoex_nb other_monoex_%  extension_nb       extension_%     annot_compatible_nb     annot_compatible_%      other_spliced_nb        other_spliced_% other_monoex_nb other_monoex_%
# lncrnagal4      5269    1772    33.6307 208     3.94762 414     7.85728 1124    21.3323 26      0.493452        1208    22.9266 426     8.08503 779     14.7846 3       0.0569368
# awk '{s=$1"\t"$2"\t"; for(i=21; i<=NF; i++){s=(s)($i)("\t")}print s}' refine_comptr_output_fortable.tsv
# set     total   intergenic_nb   intergenic_%    spliced_nb      spliced_%       monoex_nb       monoex_%        antisense_nb    antisense_%     spliced_nb      spliced_%       monoex_nb  monoex_%        unstr_nb        unstr_pcent
# lncrnagal4      5268    1762    33.4472 1704    32.3462 58      1.10099 526     9.98481 477     9.05467 49      0.930144        0       0

# check
# from the log file of refine_comptr_output.sh
##############################################
#   93868 annot
#   35080 intergenic
#   23961 antisense
#   10440 extension   
# and from the output of refine_comptr_output.sh
################################################
# cd ~/ENCODE_AWG/Analyses/Tr_build_and_Quantif/Cuff_vs_Stringtie/Stringtie
# awk 'NR>=2{print $2, $4}' K562_polya+_wcell_biorep1_star2_stringtie1.0.4_complete_comp_refinedclass_nbex.tsv | sort | uniq -c | sort -k1,1nr
#   58810 Exact annot
#   33819 Monoexonic intergenic
#   32403 Monoexonic annot
#   20606 Monoexonic antisense
#    8672 Extension extension
#    3355 Intergenic_or_antisense antisense
#    1767 Overlap extension
#    1408 Overlap annot
#    1261 Intergenic_or_antisense intergenic
#    1247 Inclusion annot
#       1 Monoexonic extension

BEGIN{
    OFS="\t";
}

NR>=2{
    split($0,a,"\t"); 
    ntot++; 
    if(a[4]=="unstranded")
    {
	nunstr++
    }
    else
    {
	if(a[4]=="annot")
	{
	    nannot[1]++; 
	    if(a[2]=="Exact")
	    {
		nannot[2]++;
	    }
	    else
	    {
		if(a[2]=="Inclusion")
		{
		    nannot[3]++;
		}
		else
		{
		    if(a[2]=="Overlap")
		    {
			nannot[4]++;
		    }
		    else
		    {
			if(a[2]=="Monoexonic")
			{
			    nannot[5]++;
			}
		    }
		}
	    }
	}
	else
	{
	    if(a[4]=="extension")
	    {
		nextens[1]++
		if(a[2]=="Extension")
		{
		    nextens[2]++
		}
		else
		{
		    if(a[2]=="Overlap")
		    {
			nextens[3]++
		    }
		    else
		    {
			if(a[2]=="Monoexonic")
			{
			    nextens[4]++;
			}
		    }
		}
	    }
	    else
	    {
		if(a[4]=="intergenic")
		{
		    ninter[1]++;
		    if(a[2]=="Intergenic_or_antisense")
		    {
			ninter[2]++;
		    }
		    else
		    {
			if(a[2]=="Monoexonic")
			{
			    ninter[3]++;
			}
		    }
		}
		else
		{
		    if(a[4]=="antisense")
		    {
			nas[1]++;
			if(a[2]=="Intergenic_or_antisense")
			{
			    nas[2]++;
			}
			else
			{
			    if(a[2]=="Monoexonic")
			    {
				nas[3]++
			    }
			}
		    }
		}
	    }
	}
    }
}

END{
    print "set", "total", "annot_nb", "annot_%", "exact_nb", "exact_%", "inclusion_nb", "inclusion_%", "other_spliced_nb", "other_spliced_%", "other_monoex_nb", "other_monoex_%", "extension_nb", "extension_%", "annot_compatible_nb", "annot_compatible_%", "other_spliced_nb", "other_spliced_%", "other_monoex_nb", "other_monoex_%", "intergenic_nb", "intergenic_%", "spliced_nb", "spliced_%", "monoex_nb", "monoex_%", "antisense_nb", "antisense_%", "spliced_nb", "spliced_%", "monoex_nb", "monoex_%", "unstr_nb", "unstr_pcent";
    print lid, ntot, nannot[1], nannot[1]/ntot*100, nannot[2], nannot[2]/ntot*100, nannot[3], nannot[3]/ntot*100, nannot[4], nannot[4]/ntot*100, nannot[5], nannot[5]/ntot*100, nextens[1], nextens[1]/ntot*100, nextens[2], nextens[2]/ntot*100, nextens[3], nextens[3]/ntot*100, nextens[4], nextens[4]/ntot*100, ninter[1], ninter[1]/ntot*100, ninter[2], ninter[2]/ntot*100, ninter[3], ninter[3]/ntot*100, nas[1], nas[1]/ntot*100, nas[2], nas[2]/ntot*100, nas[3], nas[3]/ntot*100, (nunstr!="" ? nunstr : 0), nunstr/ntot*100;
} 