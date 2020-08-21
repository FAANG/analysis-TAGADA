# gather2.awk
# gathers different pieces of information about atac-seq peaks from 4 livestock species projected to human and merged then
# in particular to a file that is a tsv file with header that has the following
# species	peak.coord	peak.size	peak.align	peak.avg.normaccess	peak.da	peak.class
# will add at the right place (column position)
# - the unique projection status
# - the % of peak size that is aligned when the peak is aligned to human
# - the % of similarity of the aligned peak when the peak is aligned to human
# - the overlap with human dnase peak when the peak is aligned to human
# - the overlap with projected enhancers from 3 tissues when the peak is aligned to human
# - the number of species with an ortholog peak (even if the peak is ambiguous)
# it also adds the number of the class in the existing class (1_tss instead of tss) to be able to plot in ok order with ggplot2
# so the header of the output file will be like this:
#####################################################
# species	peak.coord	peak.size	peak.align	pcent.align	pcent.sim	unique	nborth	over.dhs	over.enh	peak.avg.normaccess	peak.da	peak.class

# in case the main input file does not exist any more then regenerate it like this
# meta=/work/project/fragencode/data/metadata/genome_ourname_ucscname.tsv 
# time awk '{print $1}' $meta | while read sp
# do
#     acc=/work/project/fragencode/results/atacseq/$sp/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.bed
#     da=/work/project/fragencode/results/atacseq/$sp/tissue.peaks.merged/diffcounts.nominsum/data/mergedpeaks.peaknb.allexp.readnb.normcounts.diff.all.bed
#     class=/work/project/fragencode/results/atacseq/$sp/tissue.peaks.merged/mergedpeaks_simpleclass.tsv
#     awk -v fileRef1=$acc -v fileRef2=$da -v fileRef3=$dir/$sp/maf-convert.besthit.psl -v sp=$sp 'BEGIN{OFS="\t"; while (getline < fileRef1 >0){for(i=7; i<=NF; i++){n=NF; acc[$1":"$2"-"$3]+=$i} acc[$1":"$2"-"$3]=acc[$1":"$2"-"$3]/(n-7+1)}  while (getline < fileRef2 >0){da[$1":"$2"-"$3]=1} while (getline < fileRef3 >0){ok[$10]=1}} NR>=2{coord=$1":"$2"-"$3; print sp, coord, $3-$2, (ok[coord]==1 ? "aligned" : "not.aligned"), acc[coord], (da[coord]==1 ? "da" : "non.da"), $4}' $class
# done | awk 'BEGIN{OFS="\t"; print "species", "peak.coord", "peak.size", "peak.align", "peak.avg.normaccess", "peak.da", "peak.class"} {print}' > sp.peak.coord.size.align.avgnormaccess.da.class.tsv
# or better generate this information in the script below

# Example
#########
# cd /work/project/fragencode/results/atacseq/multi/tissue.peaks.merged/project.peaks.to.human
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/gather2.awk
# time awk -v fileRef1=4species.besthit.merged.tsv -v fileRef2=4species.besthit.merged.max1eachspecies.tsv -v fileRef3=4species.besthit.merged.inter.humanpeaks.tsv -v fileRef4=4species.besthit.merged.inter.enh3tissues.tsv -f $pgm sp.peak.coord.size.align.avgnormaccess.da.class.tsv > sp.peak.coord.size.align.palign.psim.over.dhs.enh.avgnormaccess.da.class.tsv
# real	0m40.699s  *** see output below

# inputs
########
# chr1	10264	10433	gallus_gallus:25:2610497-2611185,gallus_gallus:AADN04003056.1:4219-4518	2
# 215620 (5 fields)
# chr1	180083	180276	gallus_gallus:AADN04007855.1:5808-8532:186.54:3.63436:76.7677:intergenic,sus_scrofa:6:63179987-63181712:117.771:8.98551:81.2903:intron,	2
# 212021 (5 fields)
# chr1	10264	10433	gallus_gallus:25:2610497-2611185,gallus_gallus:AADN04003056.1:4219-4518	2	.	-1	-1	.	.	0
# 235532 (11 fields)
# chr1	10264	10433	gallus_gallus:25:2610497-2611185,gallus_gallus:AADN04003056.1:4219-4518	2	.	-1	-1	.	.	0
# 235532 (11 fields)
# species	peak.coord	peak.size	peak.align	peak.avg.normaccess	peak.da	peak.class
# sus_scrofa	1:14949-15262	313	aligned	16.3673	non.da	exon
# 449018 (7 fields)

# output
########
# species	peak.coord	peak.size	peak.align	pcent.align	pcent.sim	unique	nborth	over.dhs	over.enh	peak.avg.normaccess	peak.da	peak.class
# sus_scrofa	1:14949-15262	313	aligned	36.1022	74.3363	1	0	1	1	16.3673	non.da	4_exon
# 449018 (13 fields)

BEGIN{
    OFS="\t";
    sp[1]="bos_taurus";
    sp[2]="capra_hircus";
    sp[3]="gallus_gallus";
    sp[4]="sus_scrofa";
    no["tss"]=1;
    no["tts"]=2;
    no["intron"]=3;
    no["exon"]=4;
    no["intergenic"]=5;
    
    # The first file is to know the ambiguity status of the aligned peaks of each species and the nb of other species with an ortholog peak (even if the peak is ambiguous)
    # chr1	10264	10433	gallus_gallus:25:2610497-2611185,gallus_gallus:AADN04003056.1:4219-4518	2
    while (getline < fileRef1 >0) 
    {
	uniq=1;
	nbsp=0;
	for(i=1; i<=4; i++)
	{
	    nb[sp[i]]=0;
	}
	split($4,a,",");
	k=1;
	while(a[k]!="")
	{
	    split(a[k],b,":");
	    nb[b[1]]++;
	    if(nb[b[1]]==1)
	    {
		nbsp++;
	    }
	    k++;
	}
	i=1;
	while(uniq==1&&i<=4)
	{
	    if(nb[sp[i]]>1)
	    {
		uniq=0;
	    }
	    i++;
	}
	k=1;
	while(a[k]!="")
	{
	    u[a[k]]=uniq;
	    nbothersp[a[k]]=nbsp-1;
	    k++;
	}
    }

    # The second file is to know the % of the initial peak size that is aligned and the % of similarity for the aligned peaks
    # chr1	180083	180276	gallus_gallus:AADN04007855.1:5808-8532:186.54:3.63436:76.7677:intergenic,sus_scrofa:6:63179987-63181712:117.771:8.98551:81.2903:intron,	2
    # where the order of the subfields separated by : is avgnormacces then % size aligned and then % similarity
    while (getline < fileRef2 >0)
    {
	split($4,a,",");
	k=1;
	while(a[k]!="")
	{
	    split(a[k],b,":");
	    palign[b[1]":"b[2]":"b[3]]=b[5];
	    psim[b[1]":"b[2]":"b[3]]=b[6];
	    k++;
	}
    }

    # The third file is to know the overlap with human dnase peak when the peak is aligned to human
    # chr1	10264	10433	gallus_gallus:25:2610497-2611185,gallus_gallus:AADN04003056.1:4219-4518	2	.	-1	-1	.	.	0
    while (getline < fileRef3 >0)
    {
	split($4,a,",");
	k=1;
	while(a[k]!="")
	{
	    if($NF!=0)
		overdhs[a[k]]=1;
	    else
		overdhs[a[k]]=0;
	    k++;
	}
    }
    
    # The fourth file is to know the overlap with projected enhancers from 3 tissues when the peak is aligned to human
    # chr1	10264	10433	gallus_gallus:25:2610497-2611185,gallus_gallus:AADN04003056.1:4219-4518	2	.	-1	-1	.	.	0
    while (getline < fileRef4 >0)
    {
	split($4,a,",");
	k=1;
	while(a[k]!="")
	{
	    if($NF!=0)
		overenh[a[k]]=1;
	    else
		overenh[a[k]]=0;
	    k++;
	}
    }
}

# species	peak.coord	peak.size	peak.align	peak.avg.normaccess	peak.da	peak.class
NR==1{
    print $1, $2, $3, $4, "pcent.align", "pcent.sim", "unique", "nborth", "over.dhs", "over.enh", $5, $6, $7;
}

NR>=2{
    $7=no[$7]"_"$7;
    p=$1":"$2;
    print $1, $2, $3, $4, na(palign[p]), na(psim[p]), na(u[p]), na(nbothersp[p]), na(overdhs[p]), na(overenh[p]), $5, $6, $7;
}

function na(x){
    return (x=="" ? "NA" : x);
}
