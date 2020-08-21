# ~sdjebali/Awk/make_gene_pairs_from_annot_better.awk
# same as ~sdjebali/Awk/make_gene_pairs_from_annot.awk except that:
# - it does not hard-code the chr names so it can be applied to any species
# - it produces a 4th file with the gene pairs that lie on different chromosomes (which is a huge file that takes a lot of time to be generated)
# - the 4 files generated only contain gene names (and with their dots), not any other information, just to save space, and are tab separated
# NOTE: should not be run twice in the same dir because produces files with fixed names
# NOTE: the script ~/Awk/real_gene_pairs_5p3p_into_4cat.awk makes 4 files out of real gene pairs which are the 3 above categories and the diff chr one
# Recommendation: gzip the 4 files as soon as they are generated since the 4th one is taking up a lot of space

# example
#########
# cd /no_backup/rg/sdjebali/Chimeras/Benchmark/Data/Make_Chimeras_from_Annot
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/Long/gencode.v19.annotation.long.gtf
# time awk -v fileRef=$annot -f ~sdjebali/Awk/make_gene_pairs_from_annot_better.awk $annot
# real	64m48.329s  ** 1M each of the same chr files and !!! 39G !!! for the diff chr one, and this takes a lot of time as well, need to gzip at the end
# 4.7G in total when zipped, instead of 39G so a lot of space by gzipping
# some checks
#############
# wc -l *
# 1132146444 genepairs_diffchr.tsv
#   29429405 genepairs_nonoverlap_samechr_diffstr.tsv
#   29440465 genepairs_nonoverlap_samechr_samestr_kogxorder.tsv
#   29440465 genepairs_nonoverlap_samechr_samestr_okgxorder.tsv
# 1220456779 total
# check the nb on diff chr is correct (in previous tests I already checked the 3 other nb were correct approx)
# 48807*48806/2 = 1 191 037 221  *** to get all the gene pairs in total but here we only want the ones on diff chr
# to get the ones on the same chr
# annot=/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version19/Long/gencode.v19.annotation.long.gtf
# awk '$3=="gene"{nb[$1]++}END{for(i in nb){s+=((nb[i]*(nb[i]-1))/2)} print s}' $annot
# 58 890 777
# so the ones on diff chr are
# 1191037221-58890777=1132146444  *** ok number


BEGIN{
    OFS="\t";
    while (getline < fileRef >0) 
    {
	nbstr=2;
	str[1]="+";
	str[2]="-";
	if($3=="gene")
	{
	    seenchr[$1]++;
	    if(seenchr[$1]==1)
	    {
		nbchr++;
		chr[nbchr]=$1;
	    }
	}
    }
}

$3=="gene"{
    split($10,a,"\"");
    nb[$1,$7]++; 
    gene[$1,$7,nb[$1,$7]]=a[2]; 
    gbeg[$1,$7,nb[$1,$7]]=$4; 
    gend[$1,$7,nb[$1,$7]]=$5;
} 

END{
    for(k=1; k<=nbchr; k++)
    {
	# for same chr, same strand gene pairs
	for(i=1; i<=nbstr; i++)
	{
	    for(l=1; l<=(nb[chr[k],str[i]]-1); l++)
	    {
		for(m=(l+1); m<=nb[chr[k],str[i]]; m++)
		{
		    if(foverlap(gbeg[chr[k],str[i],l],gend[chr[k],str[i],l], gbeg[chr[k],str[i],m], gend[chr[k],str[i],m])==0)
		    {
			mid1=(gbeg[chr[k],str[i],l]+gend[chr[k],str[i],l])/2; 
			mid2=(gbeg[chr[k],str[i],m]+gend[chr[k],str[i],m])/2; 
			if(((i==1)&&(mid1<mid2))||((i==2)&&(mid1>mid2)))
			{
			    print gene[chr[k],str[i],l], gene[chr[k],str[i],m] > "genepairs_nonoverlap_samechr_samestr_okgxorder.tsv";
			    print gene[chr[k],str[i],m], gene[chr[k],str[i],l] > "genepairs_nonoverlap_samechr_samestr_kogxorder.tsv";
			}
			else
			{
			    if(((i==1)&&(mid1>mid2))||((i==2)&&(mid1<mid2)))
			    {
				print gene[chr[k],str[i],m], gene[chr[k],str[i],l] > "genepairs_nonoverlap_samechr_samestr_okgxorder.tsv";
				print gene[chr[k],str[i],l], gene[chr[k],str[i],m] > "genepairs_nonoverlap_samechr_samestr_kogxorder.tsv";
			    }
			}
		    }
		}
	    }
	}
	# for same chr, different strand gene pairs
	i=1;
	j=2;
	for(l=1; l<=nb[chr[k],str[i]]; l++)
	{
	    for(m=1; m<=nb[chr[k],str[j]]; m++)
	    {
		if(foverlap(gbeg[chr[k],str[i],l],gend[chr[k],str[i],l], gbeg[chr[k],str[j],m], gend[chr[k],str[j],m])==0)
		{
		    print gene[chr[k],str[i],l], gene[chr[k],str[j],m] > "genepairs_nonoverlap_samechr_diffstr.tsv";
		}	
	    }
	}
	# for diff chr gene pairs (longest step and bigger file being generated) (the only condition here is that the chr are different)
	for(k2=1; k2<=nbchr; k2++)
	{
	    if(k2!=k)
	    {
		for(i=1; i<=nbstr; i++)
		{
		    for(i2=1; i2<=nbstr; i2++)
		    {
			for(l=1; l<=nb[chr[k],str[i]]; l++)
			{
			    for(l2=1; l2<=nb[chr[k2],str[i2]]; l2++)
			    {
				if(gene[chr[k],str[i],l] < gene[chr[k2],str[i2],l2])
				{
				    print gene[chr[k],str[i],l], gene[chr[k2],str[i2],l2] > "genepairs_diffchr.tsv";
				}
			    }
			}
		    }
		}
	    }
	}
    }
}


function foverlap(beg1,end1,beg2,end2)
{
  return ((end1>=beg2)&&(beg1<=end2)) ? 1 : 0;
}
