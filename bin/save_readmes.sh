
# usage
# save_readme.sh sourcedir savedir

# example
# time ~/fragencode/tools/multi/Scripts/Bash/save_readmes.sh ~/fragencode/workspace/sdjebali ~/save/READMEs_frag
# real	0m35.087s
# time ~/fragencode/tools/multi/Scripts/Bash/save_readmes.sh ~/dynagen/sdjebali ~/save/READMEs_dyna
# real	0m3.145s
# time ~/fragencode/tools/multi/Scripts/Bash/save_readmes.sh ~/fragencode/workspace/unipluri ~/save/READMEs_uni
# real	0m7.073s
# time ~/fragencode/tools/multi/Scripts/Bash/save_readmes.sh ~/crct/chimeras ~/save/READMEs_crct
# real	0m0.094s

# Check the inputs do exist
###########################
if [ ! -n "$1" ] || [ ! -n "$2" ]
then
    echo "" >&2
    echo Usage: save_readme.sh sourcedir targetdir >&2
    echo "" >&2
    echo "takes as input:" >&2
    echo "- a source directory under which all the README.sh files will be saved" >&2
    echo "- a target directory under which all the README.sh found in sourcedir will be copied (saved)" >&2
    echo "and copies all the README.sh files that are in the tree structure under sourcedir in a comparable" >&2
    echo "tree structure under targetdir" >&2
    exit 1
fi

# Variable assignment
#####################
source=$1
target=$2

# create the dirs where a file README.sh is on $source
find $source -name README.sh | while read f 
do 
file=${f#$source"/"}
echo $file | awk '{s=""; n=split($1,a,"/"); for(k=1; k<=(n-1); k++){s=(s)(a[k])("/"); print s} }' 
done | while read d
do
mkdir -p $target/$d
done

# copy the readmes in the proper places
find $source -name README.sh | while read f
do 
file=${f#$source"/"}
cd $target/${file%/README.sh}
cp $source/$file .
done

# put the correct rights back
chmod 755 -R $target
