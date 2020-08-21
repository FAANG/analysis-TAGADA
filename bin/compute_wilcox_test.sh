#!/bin/bash

# usage:
# compute_wilcox_test.sh file field1 field2 
# computes wilcox.test(file[,1],file[,2],alternative="less")

# example:
# In ~/ENCODEExt/CTCF_vs_RF/
# compute_wilcox_test.sh allassign/chr21_GM06990_prop_ctype_3others_crossing_diff.txt 1 2
# W = 13, p-value = 0.08246

echo 'x=read.table("'$1'")[,'$2']; y=read.table("'$1'")[,'$3']; wilcox.test(x,y,alternative="less")' | R --vanilla
