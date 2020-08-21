# reformat.truvari.for.tsv.awk

# example
# dir=/work/project/seqoccin/svdetection/nanopore/version.05092019.buildpop/seqoccinsv/compare.lr.sr/bos_taurus/svdetectionevaluation
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/reformat.truvari.for.tsv.awk
# no=1
# rt=short
# pip=cnvpipelines
# st=DEL
# outdir=$dir/allindiv/$rt/$pip/$st
# cd $outdir
# awk -v no=$no -v rt=$rt -v pip=$pip -v st=$st -f $pgm allindiv.$rt.$pip.$st.truvari.out | awk 'BEGIN{OFS="\t"; print "Pipeline", "readtype", "pipeline", "SVtype", "Recall", "Precision", "F1score", "TP", "FP", "FN"} {print}'

# input
# {
#     "TP-base": 2464,
#     "TP-call": 2464,
#     "FP": 29,
#     "FN": 144,
#     "precision": 0.9883674288006418,
#     "recall": 0.9447852760736196,
#     "f1": 0.9660850813565968,    
#     "base cnt": 2608,
#     "call cnt": 2493,
#     "base size filtered": 0,
#     "call size filtered": 3,
#     "base gt filtered": 0,
#     "call gt filtered": 0,
#     "TP-call_TP-gt": 0,
#     "TP-call_FP-gt": 2464,
#     "TP-base_TP-gt": 0,
#     "TP-base_FP-gt": 2464,
#     "gt_precision": 0,
#     "gt_recall": 0,
#     "gt_f1": "NaN"
# }

# output
# 1_short.cnvpipelines	short	cnvpipelines	DEL	0.9447852760736196	0.9883674288006418	0.9660850813565968	2464	29	144


BEGIN{
    OFS="\t";
}

{
    split($2,a,",");
    if($1=="\"TP-base\":")
    {
	tp=a[1];
    }
    else
    {
	if($1=="\"FP\":")
	{
	    fp=a[1];
	}
	else
	{
	    if($1=="\"FN\":")
	    {
		fn=a[1];
	    }
	    else
	    {
		if($1=="\"precision\":")
		{
		    prec=a[1];
		}
		else
		{
		    if($1=="\"recall\":")
		    {
			sn=a[1];
		    }
		    else
		    {
			if($1=="\"f1\":")
			{
			    f1=a[1];
			}
		    }
		}
	    }
	}
    }
}

END{
    print no"_"rt"."pip, rt, pip, st, nn(sn), nn(prec), nn(f1), nn(tp), nn(fp), nn(fn)} function nn(x){return (x=="" ? "NA" : x);
}
