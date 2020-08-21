#bed2gff2.awk
# same as bed2gff.awk except that it takes as a parameter the program having generated the features and the type of feature 
# and also enables to print out the score as the 10th field of the gff file in case you specify -v sc=something on the command line

#awk -v sc=1 -v pgm=Affy -v feat=RF -v fldtoadd="expid: 5RACE001 pool: 1 trname: RP1-315G1.5-002 gnname: BIRC4 sourceid: 221485RACE001EX prseq: GTTCTTACCAGACCTCCTCAAGTGAATG" -f ~/Awk/bed2gff2.awk RF/5RACE429knowngenes_5RACE001.bed 

#RF/5RACE429knowngenes_5RACE001.bed 
#chrX    122750680       122750706       POOL1_221485RACE001EX   0       +

#output
#chrX    Affy    RF      122750681       122750680       .       +       .       expid: 5RACE001 pool: 1 trname: RP1-315G1.5-002 gnname: BIRC4 sourceid: 221485RACE001EX prseq: GTTCTTACCAGACCTCCTCAAGTGAATG


($1!="browser")&&($1!="track"){
  if($6=="")
    $6=".";
 
  if(pgm=="")
      pgm=".";
  if(feat=="")
      feat=".";

  if(sc!="")
      printf("%s\t%s\t%s\t%i\t%i\t.\t%s\t.\tname: %s sc: %i %s\n",$1,pgm,feat,$2+1,$3,$6,$4,$5,fldtoadd) 
  else
      printf("%s\t%s\t%s\t%i\t%i\t.\t%s\t.\tname: %s %s\n",$1,pgm,feat,$2+1,$3,$6,$4,fldtoadd)
}
