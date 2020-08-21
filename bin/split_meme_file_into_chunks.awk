# ~sdjebali/save/bin/split_meme_file_into_chunks.awk
# splits a meme files into nb chunks (default 10) so as to have nb meme files and be able to run fimo
# nb times and speed up the calculation. It also enables to set the nucleotide composition (default = uniform)

# example
# cd /work/project/fragencode/workspace/sdjebali/tfbs/test/my_data/all
# time awk -v nb=5 -v compo="A 0.29132 C 0.20868 G 0.20868 T 0.29132" -f ~sdjebali/save/bin/split_meme_file_into_chunks.awk JASPAR_CORE_2016_vertebrates.meme
# real	0m0.105s
# will split the meme file into 5 chuncks
# with a chunck starting like this
# MEME version 4

# ALPHABET= ACGT

# strands: + -

# Background letter frequencies (from uniform background):
# A 0.25000 C 0.25000 G 0.25000 T 0.25000 

# MOTIF MA0002.2 RUNX1

# letter-probability matrix: alength= 4 w= 11 nsites= 2000 E= 0
#   0.143500	  0.248000	  0.348000	  0.260500	
#   0.117000	  0.242500	  0.233500	  0.407000	
#   0.061500	  0.536000	  0.074500	  0.328000	
#   0.028500	  0.000000	  0.003500	  0.968000	
#   0.000000	  0.037500	  0.936000	  0.026500	
#   0.043500	  0.063500	  0.035000	  0.858000	
#   0.000000	  0.000000	  0.993500	  0.006500	
#   0.008500	  0.021000	  0.924000	  0.046500	
#   0.005000	  0.200000	  0.125500	  0.669500	
#   0.065500	  0.231500	  0.040500	  0.662500	
#   0.250000	  0.079000	  0.144500	  0.526500	


BEGIN{OFS="\t";
    if(nb=="")
    {
	nb=10;
    }
    if(compo=="")
    {
	compo="A 0.25000 C 0.25000 G 0.25000 T 0.25000";
    }
}

NR<=7{
    hline[NR]=$0;
}

NR>9&&$1=="MOTIF"{
    n++;
    nbl[n]++;
    line[n,nbl[n]]=$0;
}

NR>9&&$1!="MOTIF"{
    nbl[n]++;
    line[n,nbl[n]]=$0;
}

END{
    c=int(n/nb);
    for(i=1; i<=nb; i++)
    {
	for(j=1; j<=7; j++)
	{
	    print hline[j] > "file_"i"th.meme";
	}
	print compo > "file_"i"th.meme";
	print "" > "file_"i"th.meme";
	
	for(k=(i-1)*c+1; k<=i*c; k++)
	{
	    for(j=1; j<=nbl[k]; j++)
	    {
		print line[k,j] > "file_"i"th.meme";
	    }
	}
    }
    for(k=(i-1)*c+1; k<=n; k++)
    {
	for(j=1; j<=nbl[k]; j++)
	{
	    print line[k,j] > "file_"(i-1)"th.meme";
	}	
    }
}
