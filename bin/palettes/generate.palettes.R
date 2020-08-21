#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("gplots"))

cbbPalette1 <- rgb(matrix(c(0,0,0,0,73,73,0,146,146,255,109,182,255,182,119,73,0,146,0,109,219,182,109,
255,109,182,255,182,219,255,146,0,0,146,73,0,219,209,0,36,255,36,255,255,109), ncol=3, byrow=T), max=255)

cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#000000", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
palettes = c("rainbow", "terrain.colors", "heat.colors", "topo.colors", "redgreen", "cm.colors")

for (i in 2:15) {
	for (palette in palettes) {
		f = sprintf("%s.%s.txt", palette, i)
		write.table(get(palette)(i), f, quote=F, row.names=F, col.names=F)
	}
}

for (i in 1:nrow(brewer.pal.info)) {
	palette = rownames(brewer.pal.info)[i]
	j = brewer.pal.info[i,1]
	f = sprintf("%s.%s.txt", palette, j)
	write.table(brewer.pal(j, palette), f, quote=F, row.names=F, col.names=F)
}

f = sprintf("cbbPalette.%s.txt", length(cbbPalette))
write.table(cbbPalette, f, quote=F, row.names=F, col.names=F)

f = sprintf("cbbPalette1.%s.txt", length(cbbPalette1))
write.table(cbbPalette1, f, quote=F, row.names=F, col.names=F)

q(save="no")
