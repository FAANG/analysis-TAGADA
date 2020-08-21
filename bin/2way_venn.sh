#!/bin/bash

# 2way_venn.sh accepts 4 inputs
# - file for 1st set of the venn (string)
# - file for 2nd set of the venn (string)
# - title for the venn (string)
# - file for venn plot (string)

# has been tested on genologin with R 3.3.3 and fragencode configuration
# module use /work/project/fragencode/.local/privatemodules/
# module load fragencode-2017.10.06

# pgm=/work/project/fragencode/tools/multi/Scripts/Bash/2way_venn.sh
# set1=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.callers/connected.components/allrun.minimap.offspring.DEL.2svcallers.cc.of.scsniffles.txt
# set2=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.callers/connected.components/allrun.minimap.offspring.DEL.2svcallers.cc.of.scsvim.txt
# output=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.callers/venn.diagrams/allrun.minimap.offspring.DEL.2svcallers.venn.png
# time $pgm $set1 $set2 "trio1-offspring-allrun-minimap-DEL" $output
# real	0m1.202s


if [ ! -n "$5" ] || [ ! -n "$6" ]
then
    set1=sniffles
    set2=svim
else
    set1=$5
    set2=$6
fi

echo '
     library("futile.logger")
     library("VennDiagram")
     png("'$4'")
     run1=unlist(read.delim("'$1'", sep="\t", h=FALSE))
     run2=unlist(read.delim("'$2'", sep="\t", h=FALSE))
     venn.plot <- venn.diagram(
     x=list( 
     A=run1,
     B=run2
     ),
	main="'$3'",   
	filename = NULL,
	fill=c("#FBB4AE", "#B3CDE3"), 
	alpha=0.50,
	category = c("'$set1'","'$set2'"), #legend
	fontfamily="arial",
	main.cex=2,
	main.fontfamily="arial",
	sub.cex=1.3,
	sub.fontfamily="arial",
	cat.fontfamily=c("arial","arial") 
	)
	grid.draw(venn.plot)
	dev.off()
' | R --vanilla
