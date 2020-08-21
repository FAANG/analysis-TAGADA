#!/bin/bash
# note: when one has a contigency table he should use chi2 test for independence unless
# the table of expected values contains numbers smaller than 5; in the latter case the test
# to be used is Fisher exact test. Fisher exact test computes the probability of observing
# the result exactly as the one observed, while chi2 test computes the result as extreme as
# the one observed (or more extreme). I.e., chi2 test return cumulative probability
# !!! here I use a one sided test where I want the observed value to be GREATER than the expected one
# !!! this is not possible with chi2 test

# usage:
# launch_fisher_test.sh nA nB nA&&B nuniv

# compute the fisher test using the following formula
# fisher.test(matrix(c(n(DNA&&RNA),n(DNA&&not(RNA)),n(not(DNA)&&RNA),n(not(DNA)&&not(RNA))),ncol=2,byrow=T))

# example:
# launch_fisher_test.sh 75613 24299 695 279719388
# in the context of correlation between gn gn connections found by DNA PET and RNA PET (K562)

if [ ! -n "$1" ] || [ ! -n "$2" ] || [ ! -n "$3" ] || [ ! -n "$4" ]
then
echo "usage: launch_fisher_test.sh nA nB nA&&B nuniv"
else
echo 'fisher.test(matrix(c('$3','$1'-'$3','$2'-'$3',('$4'-'$1')-('$2'-'$3')),ncol=2,byrow=T),alternative="greater")' | R --vanilla --slave
fi