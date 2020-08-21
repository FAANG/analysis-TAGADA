#!/bin/bash

# 3way_venn.sh accepts 5 inputs
# - file for 1st set of the venn (string)
# - file for 2nd set of the venn (string)
# - file for 3rd set of the venn (string)
# - title for the venn (string)
# - file for venn plot (string)

# has been tested on genologin with R 3.3.3 and fragencode configuration
# module use /work/project/fragencode/.local/privatemodules/
# module load fragencode-2017.10.06

# pgm=/work/project/fragencode/tools/multi/Scripts/Bash/3way_venn.sh
# set1=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/connected.components/sniffles.minimap.DEL.allrun.cc.of.indoffspring.txt 
# set2=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/connected.components/sniffles.minimap.DEL.allrun.cc.of.indfather.txt 
# set3=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.inds/connected.components/sniffles.minimap.DEL.allrun.cc.of.indmother.txt 
# output=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/venn.diagrams/sniffles.minimap.DEL.allrun.cc.3indiv.venn.png
# time $pgm $set1 $set2 $set3 "trio1-minimap+sniffles-DEL-allrun" $output

# TODO: extract the legend as well

echo '
     library("futile.logger")
     library("VennDiagram")
     png("'$5'")
     run1=unlist(read.delim("'$1'", sep="\t", h=FALSE))
     run2=unlist(read.delim("'$2'", sep="\t", h=FALSE))
     run3=unlist(read.delim("'$3'", sep="\t", h=FALSE))
     venn.plot <- venn.diagram(
     x=list( 
     A=run1,
     B=run2,
     C=run3
     ),
	main="'$4'",   
	filename = NULL,
	fill=c("#FBB4AE", "#B3CDE3", "#CCEBC5"), 
	alpha=0.50,
	category = c("offspring","father","mother"), #legend
	fontfamily="arial",
	main.cex=2,
	main.fontfamily="arial",
	sub.cex=1.3,
	sub.fontfamily="arial",
	cat.fontfamily=c("arial","arial","arial") 
	)
	grid.draw(venn.plot)
	dev.off()
' | R --vanilla
