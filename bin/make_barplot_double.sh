#!/bin/bash

# make_barplot_double.sh 
# makes a barplot with two bars, taking a file where column number $2 is the x axis names,
# column number $3 and $4 are the heights of the first and second bar respectively. 
# The argument $5, $6 and $7 are the labels for the x and y axes respectively and are used 
# to name the ps file where we are drawing too. $8 is an optional title

# In ~/ENCODE_AWG/Analyses/Expression/Gencode_Quantification 
# make_barplot_double.sh MaxPlanck_transcripts_rpkm_with_CRG_rpkm_withoutzeros_noapprox_corr_trbiotype_withpropelements.txt 1 3 4 biotype PropElements Correlation biotype_vs_PropElements-Correlation
# make_barplot_double.sh bin_ig_region_cshl_caltech_Ap_cshl_Am_genes_over_long_intergenic_ok_nb.txt 1 2 3 Size_of_Intergenic_Space Number_of_intergenic_regions longPolyA_Detected_Cufflinks longNonPolyA_Detected_Cufflinks 


# postscript(file="MaxPlanck_transcripts_rpkm_with_CRG_rpkm_withoutzeros_noapprox_corr_trbiotype_withpropelements.ps")
# x=read.table("MaxPlanck_transcripts_rpkm_with_CRG_rpkm_withoutzeros_noapprox_corr_trbiotype_withpropelements.txt")[,1]
# y=read.table("MaxPlanck_transcripts_rpkm_with_CRG_rpkm_withoutzeros_noapprox_corr_trbiotype_withpropelements.txt")[,3]
# z=read.table("MaxPlanck_transcripts_rpkm_with_CRG_rpkm_withoutzeros_noapprox_corr_trbiotype_withpropelements.txt")[,4]
# barplot(matrix(c(y,z),2,length(x),byrow=T),beside=T,xpd=F,col=c("blue","orange"),xlab="Biotype",names.arg=x, las=2, cex.axis=0.5)
# legend("center", c("NumberOfElements","Correlation"),fill=c("blue","darkorange"),bty="n")
# dev.off()


echo 'postscript(file="'$5'_vs_'$6'.eps"); x=read.table("'$1'")[,'$2']; y=read.table("'$1'")[,'$3']; z=read.table("'$1'")[,'$4']; barplot(matrix(c(y,z),2,length(x),byrow=T),beside=T,xpd=F,col=c("blue","darkorange"), xlab="'$5'", ylab="'$6'", names.arg=x, cex.axis=2, cex.names=2, main="'$9'", cex.main=2); legend("center", c("'$7'","'$8'"),fill=c("blue","darkorange"),bty="n")' | R --vanilla
# 

# echo 'postscript(file="'$5'_vs_'$6'-'$7'.ps"); x=read.table("'$1'")[,'$2']; y=read.table("'$1'")[,'$3']; z=read.table("'$1'")[,'$4']; barplot(matrix(c(y,z),2,length(x),byrow=T),beside=T,xpd=F,col=c("blue","darkorange"),xlab="'$5'", names.arg=x, las=2, cex.axis=0.5, main="'$8'"); legend("center", c("'$6'","'$7'"),fill=c("blue","darkorange"),bty="n")' | R --vanilla
