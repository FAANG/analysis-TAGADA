# ibed_and_cutoff_to_bedpe.awk
# This is to make a complete bedpe file with filtered (cutoff5) connections from promoter capture hic predictions from song et al, nature genetics, 2019 (neurons)
# it takes as input two txt files
# - the main input file is the ibed file that is a kind of bedpe file and with connections with bait/prom first and other end then (see below)
# - the fileRef file is the cutoff file with 3 columns (comma separated list for 1 bp coord of the 1st segment in gx order, comma separated list for 1bp coord of the 2nd segment in gx order and the score)
# here we use the score and the coord that we have 


# Example
# cd ~/fragencode/workspace/sdjebali/irsd/egprediction/from3D/predictions/capturehic/song.shen.2019/GSM3598048_motor
# pgm=~/fragencode/tools/multi/Scripts/Awk/ibed_and_cutoff_to_bedpe.awk
# time zcat GSM3598048_motor.ibed.txt.gz | awk -v fileRef=GSM3598048_motor.cutoff.5.washU.txt -f $pgm > GSM3598048_motor.cutoff.5.washU.bedpe


# GSM3598048_motor.ibed.txt.gz
# bait_chr	bait_start	bait_end	bait_name	otherEnd_chr	otherEnd_start	otherEnd_end	otherEnd_name	N_reads	score
# chr1	903641	927394	baitmap1	chr1	573940	579491	.	1	1.01

# GSM3598048_motor.cutoff.5.washU.txt
# chr1,1041676,1041677	chr1,1048649,1048650	6.1

{
    # reading the ibed file
    nb[$NF]++;
    bait[$NF,nb[$NF]]=$1":"$2":"$3;
    other[$NF,nb[$NF]]=$5":"$6":"$7;
    name[$NF,nb[$NF]]=bait[$NF,nb[$NF]]":"other[$NF,nb[$NF]];
}

END{
    OFS="\t";
    # reading the cutoff file and for each element making a row in the final bedpe file
    while (getline < fileRef >0)
    {
	split($1,a,",");
	split($2,b,",");
	found=0;
	i=1;
	while(found==0&&i<=nb[$3])
	{
	    split(bait[$3,i],c,":");
	    split(other[$3,i],d,":");
	    avgb=(c[2]+c[3])/2;
	    avgo=(d[2]+d[3])/2;
	    if((a[1]==c[1]&&b[1]==d[1]&&abs(a[2]-avgb)<=3&&abs(b[2]-avgo)<=3)||(a[1]==d[1]&&b[1]==c[1]&&abs(a[2]-avgo)<=3&&abs(b[2]-avgb)<=3))
	    {
		found=1;
	    }
	    i++;
	}
	j=i-1;
	print c[1], c[2], c[3], d[1], d[2], d[3], name[$3,j], $3, ".", ".";
    }
}


function abs(x){
    return (x>=0 ? x : -1*x);
}
