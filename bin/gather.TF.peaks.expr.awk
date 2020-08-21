# gather.TF.peaks.expr.awk

# example
# cd /work/project/fragencode/workspace/sdjebali/fragencode/multi/atacseq.rnaseq/tf.in.dapeaks.and.deexpr
# orth=/work/project/fragencode/data/species/multi/homo_sapiens-sus_scrofa/ens91_1to1orth_to_sus_scrofa_hsid_spgnid.tsv
# oktf=/work/project/fragencode/data/species/homo_sapiens/gene_lists/TF/Homo_sapiens_TF_jasparok.txt
# dapeaks=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/diffcounts.cdvsliver/data/mergedpeaks.peaknb.allexp.readnb.normcounts.diff.cd.liver.bed
# degenes=/work/project/fragencode/results/rnaseq/sus_scrofa/diffcounts.cdvsliver/ReferenceModel/data/refgenes.counts.min2tpm0.1.normcounts.diff.cd.liver.bed
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/gather.TF.peaks.expr.awk
# peaks=/work/project/fragencode/results/atacseq/sus_scrofa/tissue.peaks.merged/tfbs/mergedpeaks.interalltf.tsv
# time awk -v fileRef1=$orth -v fileRef2=$oktf -v fileRef3=$dapeaks -v fileRef4=$degenes -f $pgm $peaks > TFname.destatus.nb.cat.dapeakwithTF.nb.pcent.tcell.liver.density.in.dapeaks.tcell.liver.tsv
# real	0m49.247s

# inputs
# ENSG00000275994 ENSSSCG00000018373
# 16391 (2 fields)
# Homo_sapiens	GSX2	ENSG00000180613	Homeobox	ENSP00000319118;ENSP00000483522;	170825
# 360 (6 fields)
# AEMK02000693.1	8381	9237	.	3	+	57.74	121.51	77.25	95.75	62.05	114.16	64.83	94.62	158.97	385.57	267.37	1.9472750971068e-278	4.47191866578495	4.43018532978147
# 9145 (20 fields)
# 10	10115926	10127714	ENSSSCG00000027882	4	+	0	196.02	0	0	231.76	0	0	231.6	0	0	266.72	0	17.1122072489652	7.85279730048079
# 9438 (20 fields)
# RUNX1	TFAP2A	...	ISL2	Hes1	
# 1:14949:15262	0	0	...	0	0
# 1 (479 fields)
# 149333 (480 fields)  *** coord are bed like

# output
# TFname	DEstatus.nb	DEstatus.cat	DE.logfc	nb.tcell.over.peaks.over.it	pcent.tcell.over.peaks.over.it	nb.liver.over.peaks.over.it	pcent.liver.over.peaks.over.itdensity.in.tcell.over.peaks	density.in.liver.over.peaks
# RUNX1	1	tcell.over	-3.3078354180801	112	3.04513	119	2.1767	0.0511147	0.0300983
# 333 (10 fields) 

BEGIN{
    OFS="\t"; 
    namecat[1]="tcell.over"; 
    namecat[2]="liver.over"; 
    namecat[3]="non.diff";  
    while (getline < fileRef1 >0)
    {
	hsid[$2]=$1;
    } 
    while (getline < fileRef2 >0)
    {
	ok[$2]=1; 
	gnname[$3]=$2;
    } 
    while (getline < fileRef3 >0)
    {
	if($(NF-1)<0)
	{
	    cat[$1":"$2":"$3]=1;
	}
	else
	{
	    if($(NF-1)>0)
	    {
		cat[$1":"$2":"$3]=2;
	    }
	}
    } 
    while (getline < fileRef4 >0)
    {
	logfc[gnname[hsid[$4]]]=$(NF-1);
	if($(NF-1)<0)
	{
	    cat[gnname[hsid[$4]]]=1;
	}
	else
	{
	    if($(NF-1)>0)
	    {
		cat[gnname[hsid[$4]]]=2;
	    }
	}
    }
}

NR==1{
    for(i=1; i<=NF; i++)
    {
	if(ok[$i]==1)
	{
	    j++; 
	    idx[j]=i+1; 
	    nameTF[j]=$i; 
	    if(cat[$i]=="")
	    {
		cat[$i]=3;
		logfc[$i]="NA";
	    }
	}
    }
} 

NR>=2{
    split($1,a,":");
    if(cat[$1]=="")
    {
	cat[$1]=3;
    } 
    ntot[cat[$1]]++;
    cumlg[cat[$1]]+=(a[3]-a[2]);

    for(k=1; k<=j; k++)
    {
	nbtf[nameTF[k],namecat[cat[$1]]]+=$(idx[k]);
	if($(idx[k])>0)
	{
	    nb[nameTF[k],namecat[cat[$1]]]++;
	}
    }
} 

END{
    print "TFname", "DEstatus.nb", "DEstatus.cat", "DE.logfc", "nb.tcell.over.peaks.over.it", "pcent.tcell.over.peaks.over.it", "nb.liver.over.peaks.over.it", "pcent.liver.over.peaks.over.it", "density.in.tcell.over.peaks", "density.in.liver.over.peaks"; 
    for(k=1; k<=j; k++)
    {
	print nameTF[k], cat[nameTF[k]], namecat[cat[nameTF[k]]], logfc[nameTF[k]], nb[nameTF[k],"tcell.over"], nb[nameTF[k],"tcell.over"]/ntot[1]*100, nb[nameTF[k],"liver.over"], nb[nameTF[k],"liver.over"]/ntot[2]*100, nbtf[nameTF[k],"tcell.over"]/cumlg[1]*1000, nbtf[nameTF[k],"liver.over"]/cumlg[2]*1000;
    }
}
