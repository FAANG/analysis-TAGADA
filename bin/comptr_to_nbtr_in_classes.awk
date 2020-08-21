# comptr_to_nbtr_in_classes.awk
# takes 2 files of tr, one with the result of refine comptr script and another one as matrix of expr (tpm) of those
# and computes for a given nb of experiments (i) the number and % of tr of each class (4 broad and 11 refined)
# that are seen with tpm at least 0.1 in at least this nb of experiments (i)

# example
# comptr=/work/project/fragencode/workspace/sdjebali/fragencode/rnaseq/analysis/transcript_models/sus_scrofa/new/all/sus_scrofa_cuff_all_complete_comp_refinedclass_nbex.tsv
# expr=/work/project/fragencode/results/rnaseq/all_quantifications/sus_scrofa/newModel_transcript_tpm_matrix.tsv
# cd /work/project/fragencode/results/rnaseq/multi/analysis/transcript_models
# awk -v sp=sus_scrofa -v i=2 -v fileRef=$comptr -f ~/save/Awk/comptr_to_nbtr_in_classes.awk $expr > tmp

# inputs
# $comptr is like this
# trid	comptrclass	annottrlist	refinedtrclass	nbex
# tTCONS_00000001	Overlap	ENSSSCT00000024788,	extension	n12
# 111185 (5 fields)
# $expr is like this
# 	pig1-cd8	pig1-liver	pig2-cd4	pig2-cd8	pig2-liver	pig3-cd4	pig3-cd8	pig3-liver	pig4-cd4	pig4-cd8	pig4-liver
# TCONS_00000001	0.10	0.82	0.10	0.51	0.75	0.06	0.92	0.85	0.06	0.40	0.69
# 1 (11 fields)
# 111331 (12 fields)

# tmp is like this
# species	nbexp	tottr	totexprtr_nb	totexprtr_pcent	annot_nb	annot_pcent	extension_nb	extension_pcent	antisense_nb	antisense_pcent	intergenic_nb	intergenic_pcent	annot.Exact_nb	annot.Exact_npcent	annot.Inclusion_nb	annot.Inclusion_npcent	annot.Overlap_nb	annot.Overlap_npcent	annot.Monoexonic_nb	annot.Monoexonic_npcent	extension.Extension_nb	extension.Extension_npcent	extension.Overlap_nb	extension.Overlap_npcent	extension.Monoexonic_nb	extension.Monoexonic_npcent	antisense.Intergenic_or_antisense_nb	antisense.Intergenic_or_antisense_npcent	antisense.Monoexonic_nb	antisense.Monoexonic_npcent	intergenic.Intergenic_or_antisense_nb	intergenic.Intergenic_or_antisense_npcent	intergenic.Monoexonic_nb	intergenic.Monoexonic_npcent
# sus_scrofa	2	111037	83662	75.3461	25684	23.131	31764	28.6067	1692	1.52382	24522	22.0845	15222	13.7089	3016	2.71621	6697	6.03132	742	0.668246	8883	8.00004	22678	20.4238	203	0.182822	1094	0.985257	598	0.538559	7728	6.95984	16794	15.1247
# 2 (35 fields)

BEGIN{
    brclass[1]="annot";
    brclass[2]="extension";
    brclass[3]="antisense";
    brclass[4]="intergenic";
    fineclass[1]="annot.Exact";
    fineclass[2]="annot.Inclusion";
    fineclass[3]="annot.Overlap";
    fineclass[4]="annot.Monoexonic";
    fineclass[5]="extension.Extension";
    fineclass[6]="extension.Overlap";
    fineclass[7]="extension.Monoexonic";
    fineclass[8]="antisense.Intergenic_or_antisense";
    fineclass[9]="antisense.Monoexonic";
    fineclass[10]="intergenic.Intergenic_or_antisense";
    fineclass[11]="intergenic.Monoexonic";
    while (getline < fileRef >0)
    {
	OFS="\t";
	if($2!="Unstranded")
	{
	    tr=substr($1,2);
	    broad[tr]=$4;
	    fine[tr]=$4"."$2;
	}
    }
    s=("species")("\t")("nbexp")("\t")("tottr")("\t")("totexprtr_nb")("\t")("totexprtr_pcent")("\t");
    for(k=1; k<=4; k++)
    {
	s=(s)(brclass[k]"_nb")("\t")(brclass[k]"_pcent")("\t");
    }
    for(k=1; k<=10; k++)
    {
	s=(s)(fineclass[k]"_nb")("\t")(fineclass[k]"_pcent")("\t");
    }
    print (s)(fineclass[k]"_nb")("\t")(fineclass[k]"_pcent");
}

NR>=2{
    # if the broad class of the tr exists 
    if(broad[$1]!="")
    {
	ntot++
	n=0;
	for(k=2; k<=NF; k++)
	{
	    if($k>=0.1)
	    {
		n++;
	    }
	}
	if(n>=i)
	{
	    ntotexpr++;
	    nbb[broad[$1]]++;
	    nbf[fine[$1]]++;
	}
    }
}

END{
    s=(sp)("\t")(i)("\t")(ntot)("\t")(ntotexpr)("\t")(ntotexpr/ntot*100)("\t");
    for(k=1; k<=4; k++)
    {
	s=(s)(na(nbb[brclass[k]]))("\t")(na(nbb[brclass[k]])/ntotexpr*100)("\t");
    }
    for(k=1; k<=10; k++)
    {
	s=(s)(na(nbf[fineclass[k]]))("\t")(na(nbf[fineclass[k]])/ntotexpr*100)("\t");
    }
    print (s)(na(nbf[fineclass[k]]))("\t")(na(nbf[fineclass[k]])/ntotexpr*100);
}

function na(x)
{
    return (x!="" ? x : "NA");
}
