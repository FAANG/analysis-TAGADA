#!/usr/bin/env bash

file=$1
select=$2

if [[ $# != 2 ]]; then
echo "USAGE: $0 file.gtf <tags>"
echo "<tags> should be valid gtf keys in the 9th field, comma-separated"
exit 1
fi;



awk -v select=$select 'BEGIN{OFS=FS="\t"; 
split(select, a, ","); 
for(i=1;i<=length(a);i++) {
	if(a[i]~/^[0-9]*$/) {
		fields[a[i]]
	} else {
		tags[a[i]]
	}
}
}
{
for(f=1;f<=NF;f++) {
	if(f in fields) {
		printf $f"\t"
	}
}
split($9,a,"; "); 
for(i=1;i<=length(a);i++) {split(a[i],b," ");gsub(/[";]/,"",b[2]); if(b[1] in tags) {printf b[2]"\t"}}
printf "\n"}' $file | sed 's/\t$//g'
