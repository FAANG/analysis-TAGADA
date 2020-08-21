# read a list of values (one per line at the same column)
# print the mean and the standard deviation

BEGIN{
# if no column is specified (with -v c=N), use the first one
 s=0;s2=0;
 c=((c>0)?c:1);
}
$c==""{print "ERROR: empty column "c" line "NR"!";exit 1}
{
 s += $c;
 s2+= $c * $c;
}
END{
  print s/NR , sqrt ((NR*s2 - s*s) / ((NR*(NR-1))))
}
