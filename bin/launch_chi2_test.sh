#!/bin/bash
# note: when one has a contigency table he should use chi2 test for independence unless
# the table of expected values contains numbers smaller than 5; in the latter case the test
# to be used is Fisher exact test. Fisher exact test computes the probability of observing
# the result exactly as the one observed, while chi2 test computes the result as extreme as
# the one observed (or more extreme). I.e., chi2 test return cumulative probability
# !!! BE careful: this is a 2 sided test meaning that the pvalue can be very low in case
# the observed value is much higher OR much lower than the expected one

# usage:
# launch_chi2_test.sh nA nB nA&&B nuniv

# compute the chi2 test using the following formula
# fisher.test(matrix(c(n(DNA&&RNA),n(DNA&&not(RNA)),n(not(DNA)&&RNA),n(not(DNA)&&not(RNA))),ncol=2,byrow=T))

# example:
# launch_chi2_test.sh 757 309 30 15736
# in the context of mouse human ortholog genes with mean and variance entropy 0 in all experiments

echo 'chisq.test(matrix(c('$3','$1'-'$3','$2'-'$3',('$4'-'$1')-('$2'-'$3')),ncol=2,byrow=T))' | R --vanilla
