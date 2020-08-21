#!/bin/bash

# 4way_venn.sh accepts 6 inputs
# - file for 1st set of the venn (string)
# - file for 2nd set of the venn (string)
# - file for 3rd set of the venn (string)
# - file for 4th set of the venn (string)
# - title for the venn (string)
# - file for venn plot (string)

# has been tested on genologin with R 3.3.3 and fragencode configuration
# module use /work/project/fragencode/.local/privatemodules/
# module load fragencode-2017.10.06

# pgm=/work/project/fragencode/tools/multi/Scripts/Bash/4way_venn.sh
# set1=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/connected.components/sniffles.minimap.father.INV.cc.of.run1.txt 
# set2=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/connected.components/sniffles.minimap.father.INV.cc.of.run2.txt 
# set3=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/connected.components/sniffles.minimap.father.INV.cc.of.run3.txt 
# setall=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/connected.components/sniffles.minimap.father.INV.cc.of.runall.txt
# output=/work/project/seqoccin/svdetection/nanopore/bos_taurus/compare.vcfs/vcfs.diff.runs/venn.diagrams/sniffles.minimap.father.INV.cc.3runs.and.merged.venn.png
# time $pgm $set1 $set2 $set3 $setall "trio1-father-minimap+sniffles-INV" $output

# TODO: extract the legend as well

echo '
     library("futile.logger")
     library("VennDiagram")
     png("'$6'")
     run1=unlist(read.delim("'$1'", sep="\t", h=FALSE))
     run2=unlist(read.delim("'$2'", sep="\t", h=FALSE))
     run3=unlist(read.delim("'$3'", sep="\t", h=FALSE))
     runall=unlist(read.delim("'$4'", sep="\t", h=FALSE))
     venn.plot <- venn.diagram(
     x=list( 
     A=run1,
     B=run2,
     C=run3,
     D=runall
     ),
	main="'$5'",   
	filename = NULL,
	fill=c("#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4"), 
	alpha=0.50,
	category = c("run1","run2","run3","runall"), #legend
	fontfamily="arial",
	main.cex=2,
	main.fontfamily="arial",
	sub.cex=1.3,
	sub.fontfamily="arial",
	cat.fontfamily=c("arial","arial","arial","arial") 
	)
	grid.draw(venn.plot)
	dev.off()
' | R --vanilla
