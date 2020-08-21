# ~/Awk/remove_redund_better.awk
# idem as remove_redund.awk except that it removes the redundancy of all fields specified in fldlist
# Also the keys specified in fldlist will be used for output
# it also adds at the end the number of individual objects of the non redundant object of each line
# with the option -v simple=1 (or any string in fact) it will not report the redundancy factor

#awk -v fldlist=pool:10,tiss:12,asspr:14,assprgn:16,assprgnstr:18 -f ~/Awk/remove_redund_better.awk RF99.7_76_3_full_sort_overuspp_ok_overnegctrol_ok_sorted_closcompexppr_nr.gff  

#RF99.7_76_3_full_sort_overuspp_ok_overnegctrol_ok_sorted_closcompexppr_nr.gff 
#chr21   Cluster nrRF    9884622 9884720 .       .       .       pool_list: 18,26,42,52, tiss_list: 9,9,9,9, asspr_list: chr21_36755144_36755863_m_primer_3race_45,chr21_17903161_17903196_m_primer_3race_0,chr21_32896453_32896548_m_primer_3race_33,chr21_37224421_37224603_m_primer_3race_30, assprgn_list: CLDN14,BTG3,C21orf59,HLCS, assprgnstr_list: -,-,-,-, 

#output
#chr21   Cluster nrRF    9884622 9884720 .       .       .       pool_list: 18,26,42,52, tiss_list: 9,9,9,9, asspr_list: chr21_36755144_36755863_m_primer_3race_45,chr21_17903161_17903196_m_primer_3race_0,chr21_32896453_32896548_m_primer_3race_33,chr21_37224421_37224603_m_primer_3race_30, assprgn_list: CLDN14,BTG3,C21orf59,HLCS, assprgnstr_list: -,-,-,-,  pool_nr: 18:1,26:1,42:1,52:1, tiss_nr: 9:4, asspr_nr: chr21_36755144_36755863_m_primer_3race_45:1,chr21_17903161_17903196_m_primer_3race_0:1,chr21_32896453_32896548_m_primer_3race_33:1,chr21_37224421_37224603_m_primer_3race_30:1, assprgn_nr: CLDN14:1,BTG3:1,C21orf59:1,HLCS:1, assprgnstr_nr: -:4,  totno: 4



BEGIN{split(fldlist,a,",");   #improvement of this is if you do not specify impchar then we print everything that is in the last field (betweem spaces)
 n=1; 
 while(a[n]!="")
   {
     split(a[n],b,":");
     str[n]=b[1];   #name of the field
     fld[n]=b[2];   #position of the field
     n++;
   } 
}

{
  s="";   #s is the string where the string representing the list without redundancy and with the number of occurences will be stored
  for(k=1; k<=n-1; k++)   #for each important feature
    {
      c[k]=$(fld[k]);
      s=(s)(str[k]"_nr: ");
      split(c[k],t,",");   #fld is the number of the input file field where the list with redundacy is

      #initialise auxiliary array NB
      l=1;
      p=0;
      while(t[l]!="")
	{
	  p++;
	  NB[t[l]]=0;
	  l++;
	}
      
      #fill auxiliary array NB
      l=1;  #l is to go over the redundant list t
      i=1;  #i is to go over the non redundant list nrt
      while(t[l]!="")
	{
	  NB[t[l]]++;
	  if(NB[t[l]]==1)
	    {
	      nrt[i]=t[l];
	      i++;
	    }
      l++;
    }
      for(m=1; m<=(i-1); m++)
	{
	    if(simple=="")
		s=(s)(nrt[m])":"(NB[nrt[m]])",";
	    else
		s=(s)(nrt[m])",";	
	}
      s=(s)" ";
    }
  print $0" "(s)" totno: "p;
}

