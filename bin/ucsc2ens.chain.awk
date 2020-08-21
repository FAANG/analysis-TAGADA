# ucsc2ens.chain.awk
# this is to convert a ucsc chain file into ensembl compatible coordinates, exactly as we do with gff or bed files in ucsc2ens.awk
# except that here the only rows to change are the chain rows with 13 columns that have the 3rd and the 9th columns that contain
# chr name and therefore that should be changed into ensembl name by removing chr when it is there and also if chr is followed by M only
# then also adding T after M (MT instead of chrM)

# !!! for the 13 column rows, the separator is a space, for the 3 column rows, the separator is a tab !!!
# chain 66141295 chrUn_NW_018085100v1 2641821 + 0 1419001 chr10 79102373 - 77246284 78618338 37
# 292	8	10
# ...
# chain 18668689774 chr1 274330532 + 142824 274327149 chr1 315321322 + 224389 308419368 1

# chain 66141295 Un_NW_018085100v1 2641821 + 0 1419001 10 79102373 - 77246284 78618338 37
# 292	8	10
# ...
# chain 18668689774 1 274330532 + 142824 274327149 1 315321322 + 224389 308419368 1

# Example
#########
# cd ~/fragencode/data/species/sus_scrofa/Sscrofa11.1
# pgm=/work/project/fragencode/tools/multi/Scripts/Awk/ucsc2ens.chain.awk
# time zcat ucsc/susScr11ToSusScr3.over.chain.gz | awk -f $pgm | gzip > susScr11ToSusScr3.over.chain.ens.gz


$1!="chain"{
    OFS="\t";
    print $0;
}

$1=="chain"{
    OFS=" ";
    if($3~/^chr/)
    {
	if($3=="chrM")
	{
	    $3="MT";
	}
	else
	{
	    split($3,a,"chr");
	    $3=a[2];
	}
    }
    if($8~/^chr/)
    {
	if($8=="chrM")
	{
	    $8="MT";
	}
	else
	{
	    split($8,a,"chr");
	    $8=a[2];
	}
    }
    print $0;
}
