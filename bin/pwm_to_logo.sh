#!/bin/bash

# pwm_to_logo.sh

# example
# cd /work/project/fragencode/workspace/sdjebali/tfbs/jaspar_data/all/Split/sus_scrofa/meme/1
# time pwm_to_logo.sh file_1th_onlymatrix.tsv file_1th_onlymatrix.pdf
# real	0m0.650s
# has produced the file file_1th_onlymatrix.pdf with the logo inside


if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo "Usage: pwm_to_logo.sh pwm.tsv pwm.out" >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a position weight matrix (PWM) of 4 rows and n columns, where n is the length of the motif (usually for TF) with no header" >&2
    echo "- the name of a pdf output file (with extension) where to put the logo of the PWM" >&2
    echo "" >&2
    echo "produces as output:" >&2
    echo "" >&2
    echo "- the file specified as the last argument with the logo of the PWM provided as input" >&2
   exit 1
fi

echo '
library(seqLogo)
pdf("'$2'")
m <- read.table("'$1'")
p <- makePWM(m)
seqLogo(p)
dev.off()
' | R --vanilla
