# ~sdjebali/Awk/starlog2stats.awk

# usage
# awk -v LID=$lid -f ~sdjebali/Awk/starlog2stats.awk Log.final.out > stats.tsv

# input file is like this one
# /home/sdjebali/fragencode/results/rnaseq/bos_taurus/CD4/MappedData/Bos-taurus-1952-CD4_AATAGT_L005/Log.final.out
     #                           Started job on |	Aug 26 18:43:29
       #                       Started mapping on |	Aug 26 19:06:02
       #                              Finished on |	Aug 26 19:36:47
       # Mapping speed, Million of reads per hour |	42.94

       #                    Number of input reads |	22008729
       #                Average input read length |	299
       #                              UNIQUE READS:
       #             Uniquely mapped reads number |	18183955
       #                  Uniquely mapped reads % |	82.62%
       #                    Average mapped length |	296.05
       #                 Number of splices: Total |	20190548
       #      Number of splices: Annotated (sjdb) |	19198565
       #                 Number of splices: GT/AG |	20093421
       #                 Number of splices: GC/AG |	92879
       #                 Number of splices: AT/AC |	4248
       #         Number of splices: Non-canonical |	0
       #                Mismatch rate per base, % |	0.45%
       #                   Deletion rate per base |	0.02%
       #                  Deletion average length |	2.29
       #                  Insertion rate per base |	0.02%
       #                 Insertion average length |	1.67
       #                       MULTI-MAPPING READS:
       #  Number of reads mapped to multiple loci |	1778300
       #       % of reads mapped to multiple loci |	8.08%
       #  Number of reads mapped to too many loci |	3174
       #       % of reads mapped to too many loci |	0.01%
       #                            UNMAPPED READS:
       # % of reads unmapped: too many mismatches |	0.00%
       #           % of reads unmapped: too short |	9.28%
       #               % of reads unmapped: other |	0.01%
       #                            CHIMERIC READS:
       #                 Number of chimeric reads |	0
#                      % of chimeric reads |	0.00%

# output file is like this
# sample	total_reads	mapped_nb	mapped_pcent	unique_nb	unique_pcent	multi_nb	multi_pcent	uniquely_pcent_oftot	new_splice_sites_pcent
# Bos-taurus-1952-CD4_AATAGT_L005	27248393	26329678	96.6284	25798984	97.9844	530694	2.01557	94.68	7.52036
# note that the second unique number is the pcent of total

BEGIN{OFS="\t";
    print "sample", "total_reads", "mapped_nb", "mapped_pcent", "unique_nb", "unique_pcent", "multi_nb", "multi_pcent", "uniquely_pcent_oftot", "new_splice_sites_pcent";
}


$1=="Number"&&$3=="input"{
    input=$NF;
}

$1=="Uniquely"&&$4=="%"{
    split($NF,a,"%");
    uniquely=a[1];
}

$3=="splices:"&&$4=="Total"{
    totalsplice=$NF;
}

$3=="splices:"&&$4=="Annotated"{
    annotsplice=$NF;
}

$4=="unmapped:"&&$6=="short"{
    short=$NF;
}

$1=="Uniquely"&&$4=="number"{
    uniq=$NF;
}

$1=="Number"&&$6=="multiple"{
    multi=$NF;
}

END{
    OFS="\t";
    print LID, input, (uniq+multi), (uniq+multi)/input*100, uniq, uniq/(uniq+multi)*100, multi, multi/(uniq+multi)*100, uniquely, (totalsplice-annotsplice)/annotsplice*100;
}
