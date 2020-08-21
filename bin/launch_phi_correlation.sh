#!/bin/bash

# usage:
# launch_phi_correlation.sh na nb nc nd
# where na, nb, nc, nd is a contingency table entered by column

# compute the phi correlation for a contingency table matrix(c(na,nb,nc,nd),ncol=2) where the 
# matrix is entered per column so that na is the nfeat1feat2, nb is nfeat1barfeat2
# nc is nfeat1feat2bar, nd is nfeat1barfeat2bar

# example:
# launch_phi_correlation.sh 1000 0 0 700
# for a chimeric junction that is present in 1000 samples of a tissue t and not present in 0 sample of this tissue t
# and present in 0 sample of the other tissues and absent in 700 samples of the other tissues

if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ]
then
echo "usage: launch_phi_correlation.sh na nb nc nd" 
echo "where na, nb, nc, nd is a contingency table entered by column"
else
echo 'library(psych); phi(matrix(c('$1','$2','$3','$4'),ncol=2))' | R --vanilla --slave
fi