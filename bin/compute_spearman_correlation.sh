#!/bin/bash

# compute_spearman_correlation.sh file fld1 fld2


if [ -n "$4" ]
then
awk 'NR>=2' $1 > /tmp/coucou.$$.txt 
echo "x=read.table(\"/tmp/coucou."$$".txt\")[,"$2"]; y=read.table(\"/tmp/coucou."$$".txt\")[,"$3"]; cor.test(x,y,method=\"spearman\")" | R --vanilla --slave
rm /tmp/coucou.$$.txt
else
echo 'x=read.table("'$1'")[,'$2']; y=read.table("'$1'")[,'$3']; cor.test(x,y,method="spearman")' | R --vanilla --slave
fi
