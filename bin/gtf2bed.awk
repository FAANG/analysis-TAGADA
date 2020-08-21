#~/Awk/gtf2bed.awk
#indicates in fld the number of the field that we want to put as the 4th field in the bed file

#in /projects/encode/scaling_up/q2/mRNA_SplESTs
#cat encmrna_mar2008.gtf encsplest_mar2008.gtf | sort -k1,1 -k4,4n -k5,5n | awk -v fld=12 -f ~/Awk/gtf2bed.awk

#cat encmrna_mar2008.gtf encsplest_mar2008.gtf 
#chr7    hg17_all_mrna   exon    115443237       115445468       0.000000        -       .       gene_id "BC045580"; transcript_id "BC045580"; 

#output
#chr7    115443236       115445468       "BC045580";     0.000000        -


# changed in feb 6th 2017 in order to remove the double quotes and ; when $fld is placed in $4

BEGIN{OFS="\t";}


{
    split($fld,a,"\"");
    print $1, $4-1, $5, a[2], ($6!="." ? $6 : 0), $7;  
}

