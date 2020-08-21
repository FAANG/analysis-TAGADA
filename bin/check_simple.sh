#!/bin/bash
head -n 2 $1 | awk '{print "#", $0}'
awk '{print NF}' $1 | sort -n | uniq -c | awk '{print "#", $1, "("$2" fields)"}'
