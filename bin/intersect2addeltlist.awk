# intersect2addeltlist.awk
# from the output of an intersectBed on a vcf file (in which the combination of the first 5 fields is supposed unique)
# and a bed file of unique elements, add to the vcf file the nb and list of unique elements overlapping
# stranded variable allows to consider the strand or not

# example
# cd /work/project/fragencode/workspace/sdjebali/irsd/data/parkinson/2017/HRC_Imputations/SANGER
# module load bioinfo/bcftools-1.9
# module load bioinfo/bedtools2-2.29.0
# pgm=~/fragencode/tools/multi/Scripts/Awk/intersect2addeltlist.awk
# prom=/work/project/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/jung.ren.2019/HCmerge/pall/HCmerge.pall.score2.prom.uniq.bed.gz
# intersectBed -a allchr.basic.vcf.gz -b $prom -wao | awk -v stranded=0 -f $pgm > allchr.basic.uniqpromlist.tsv 

# time intersectBed -a allchr.basic.vcf.gz -b $prom -wao > allchr.basic.inter.uniqpromHCmerge.tsv
# real	1m9.813s

# $prom
# chr1	943049	965801
# 13319 (3 fields)

# input = intersectBed result
# chr1	13380	rs571093408	C	G	.	PASS	RefPanelAF=7.69941e-05;AN=6046;AC=1;INFO=1	.	-1	-1	0
# chr1	1183858	rs113908945	G	T	.	PASS	RefPanelAF=0.094364;AN=6046;AC=476;INFO=0.930813	chr1	1183838	1209657	1

# allchr.basic.uniqpromlist.tsv output
# same as input but with 2 columns, one with nb of elements in the list and the other one with the list of elements (defined by their chr:beg:end)

BEGIN{
    OFS="\t";
    if(stranded=="")
    {
	stranded=0;
    }
    else
    {
	stranded=1;
    }
}

{
    id=$1":"$2":"$3":"$4":"$5;
    info[id]=$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8;
    if($NF!=0)
    {
	id2=(stranded==0 ? $9":"$10"-"$11 : $9":"$10"-"$11"-"$14);
	nb[id,id2]++;
	if(nb[id,id2]==1)
	{
	    promlist[id]=(promlist[id])(id2)(",");
	    promnb[id]++;
	}
    }
}

END{
    for(s in info)
    {
	print info[s], (promnb[s]!="" ? promnb[s] : 0), (promlist[s]!="" ? promlist[s] : "NA");
    }
}
