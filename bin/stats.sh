#!/bin/bash

awk -v no=$2 '{if(no==""&&$1!="."){print $1} else {if($no!="."){print $no}}}' $1 > /tmp/coucoufromsarah.$$.txt
echo "l<-scan(\"/tmp/coucoufromsarah."$$".txt\");summary(l); which.max(l);" | R --vanilla  > /tmp/coucoufromsarahout.$$.txt 2> /tmp/coucoufromsaraherr.$$.txt
cat /tmp/coucoufromsaraherr.$$.txt | awk '{print "#", $0}'
awk '$1=="Min."{print "#", $0; l=NR+1}NR==l{print "#", $0;}NR!=1&&NR==(l+1){print "#", "argmax", $NF}' /tmp/coucoufromsarahout.$$.txt
rm /tmp/coucoufromsarah.$$.txt /tmp/coucoufromsarahout.$$.txt /tmp/coucoufromsaraherr.$$.txt
