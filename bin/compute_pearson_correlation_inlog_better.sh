#!/bin/bash

# compute_pearson_correlation_inlog_better.sh file fld1 fld2

# same as compute_pearson_correlation_inlog.sh except that it does not take out the min in x nor in y
# changed log to log10 on march 25th 2013

# on July 17th 2013 added a 4th param to say if there is a header
# on May 19th 2017 added unicity of temporary file written, tested on a small example

if [ -n "$4" ]
then
awk 'NR>=2' $1 > /tmp/coucou.$$.txt 
echo "x=read.table(\"/tmp/coucou."$$".txt\")[,"$2"]; y=read.table(\"/tmp/coucou."$$".txt\")[,"$3"]; x1=log10(x+10^(-3)); y1=log10(y+10^(-3)); cor.test(x1,y1)" | R --vanilla --slave
rm /tmp/coucou.$$.txt
else
echo 'x=read.table("'$1'")[,'$2']; y=read.table("'$1'")[,'$3']; x1=log10(x+10^(-3)); y1=log10(y+10^(-3)); cor.test(x1,y1)' | R --vanilla --slave
fi
