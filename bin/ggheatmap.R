#!/usr/bin/env Rscript

# This script was initially written by Alessandra Breschi at the CRG, Barcelona, Spain
# Its initial version is at https://github.com/abreschi/Rscripts

# DEBUG OPTIONS
opt = list()
opt$input_matrix = "~/Documents/blueprint/pilot/Flux/Long/bp.human.long.gene.RPKM.idr_01.thr_0.names_False.tsv"
opt$col_metadata = "~/Documents/blueprint/pilot/bp_rna_dashboard_mmaps.crg.tsv"
opt$col_labels = NULL
opt$row_labels = NULL
opt$colSide_by = "cell"
opt$merge_col_mdata_on = "labExpId"
opt$row_metadata = "/users/rg/projects/encode/scaling_up/whole_genome/Gencode/version15/Long/gen15.gene.super.gene_type.with_header.tsv"
opt$merge_row_mdata_on = "gene"
opt$col_dendro = TRUE
opt$row_dendro = TRUE
opt$dist = "euclidean"
opt$hclust = "complete"


#options(stringsAsFactors=FALSE)


##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(

make_option(c("-i", "--input_matrix"), default="stdin", 
	help="the matrix you want to analyze. \"stdin\" to read from standard input [default=%default]"),

make_option(c("-l", "--log"), action="store_true", default=FALSE, 
	help="apply the log10. NAs are treated as 0s and a pseudocount is added if specified [default=%default]"),

make_option(c("-p", "--pseudocount"), type="double", default=1e-04,
	help="specify a pseudocount for the log [default=%default]"),

make_option(c("--col_metadata"), 
	help="one tsv file with metadata on matrix columns. Can be left empty."),

make_option(c("--merge_col_mdata_on"), default="labExpId",
	help="which field of the metadata corresponds to the column headers? [default=%default]"), 

make_option(c("--row_metadata"), 
	help="one tsv file with metadata on matrix rows. Can be left empty."),

make_option(c("--merge_row_mdata_on"), default="labExpId",
	help="which field of the metadata corresponds to the row names? [default=%default]"), 

make_option(c("--col_labels"), 
	help="Specify the field for the col labels. \"none\" for no col labels. If empty the column headers are used."),

make_option(c("--row_labels"), 
	help="Specify the field for the col labels. \"none\" for no row labels. If empty the row names are used."),

make_option(c("--colSide_by"), 
	help="Specify the field(s), you want the column sides coloured by. If empty no color side is added."),

make_option(c("--colSide_palette"), #default="/users/rg/abreschi/R/palettes/cbbPalette.8.txt",
	help="Palette for colSide colors. Leave empty for color_hue. If multiple palettes are needed, write them comma-separated"),

make_option(c("--rowSide_by"), 
	help="Specify the field(s), you want the row sides coloured by. If empty no color side is added."),

make_option(c("--rowSide_palette"), #default="/users/rg/abreschi/R/palettes/cbbPalette.8.txt",
	help="Palette for rowSide colors"),

make_option(c("--col_dendro"), action="store_true", default=FALSE, 
	help="Print the column dendrogram [default=%default]"),

make_option(c("--row_dendro"), action="store_true", default=FALSE, 
	help="Print the row dendrogram [default=%default]"),

make_option(c("-d", "--dist"), default="euclidean",
	help="distance measure between columns. Choose among <p> (pearson), <s> (spearman), 
    <i> (if a similarity is passed) or <n> (none, if a dissimilarity is passed)
		all methods supported by the function dist(). [default=%default]"),

make_option(c("-c", "--hclust"), default="complete",
	help="hierarchical clustering method. Choose among the method of the function hclust(). [default=%default]"),

make_option(c("-B", "--base_size"), default=16,
	help="The font base size as defined in ggplot2. [default=%default]"),

make_option(c("-W", "--width"), type="integer",
	help="Choose the heatmap width in inches. Default is proportional to the number of columns"),

make_option(c("-H", "--height"), type="integer",
	help="Choose the heatmap height in inches. Default is proportional to the number of rows"),

make_option(c("--matrix_palette"), default="/work/project/fragencode/tools/multi/palettes/terrain.colors.3.txt",
	help="Palette for the heatmap color grandientn [default=%default]"),

make_option(c("--matrix_legend_title"), default="value",
	help="Title for matrix color scale [default=%default]"),

make_option(c("--matrix_fill_limits"), 
	help="Specify limits for the fill scale, e.g. \"\\-1,1\". Escape character for negative numbers [default=%default]"),

make_option(c("-X", "--matrix_text"), default=FALSE, action="store_true",
	help="Add the text to each cell of the matrix [default=%default]"),

#make_option(c("--matrix_legend_breaks"),
#	help="Comma-separated breaks for the color scale"),
#
make_option(c("-o", "--output"), default="ggheatmap.out.pdf",
	help="Output file name, with the extension. [default=%default]"),

make_option(c("-v", "--verbose"), default=FALSE, action="store_true",
	help="Verbose output [default=%default]")

)


parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


#------------#
# LIBRARIES  #
#------------#

cat("Loading libraries... ")
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(RColorBrewer))
cat("DONE\n\n")


# ==========================================
# Function for extracting legend from ggplot
# ==========================================

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}


# ======================
# Plotting variables
# ======================

base_size = opt$base_size
theme_set(theme_grey(base_size))
theme_update(axis.ticks=element_blank())
#theme_update(axis.ticks.margin = unit(0.01, "inch"))
#theme_update(axis.text = element_text(margin=unit(0.01, "inch")))
theme_update(axis.ticks.length = unit(0.01, "inch"))



# ===== #
# BEGIN #
# ===== #


# read table
if (opt$verbose) {cat("Reading input matrix...")}
if (opt$input_matrix == "stdin") {input=file("stdin")} else {input=opt$input_matrix}
m = read.table(input, h=T, sep="\t", comment.char="")
if (opt$verbose) {cat("DONE\n\n")}


# Make valid row and column names
rownames(m) <- make.names(rownames(m))
colnames(m) <- make.names(colnames(m))

# --- read PALETTE files ----

#if (!is.null(opt$rowSide_palette)) {
#	rowSide_palette = as.character(read.table(opt$rowSide_palette, h=F, sep="\t", comment.char="%")$V1)
#	if (opt$verbose) {cat("RowSide Palette:", rowSide_palette, "\n")}
#}

if (!is.null(opt$rowSide_palette)) {
	rowSide_palette_files = strsplit(opt$rowSide_palette, ",")[[1]]
	rowSide_palette = sapply(rowSide_palette_files, function(x)  as.character(read.table(x, h=F, sep="\t", comment.char="%")$V1), simplify=FALSE)
	#if (opt$verbose) {cat("ColSide Palette:", colSide_palette[1], "\n")}
}

if (!is.null(opt$colSide_palette)) {
	colSide_palette_files = strsplit(opt$colSide_palette, ",")[[1]]
	colSide_palette = sapply(colSide_palette_files, function(x)  as.character(read.table(x, h=F, sep="\t", comment.char="%")$V1), simplify=FALSE)
	#if (opt$verbose) {cat("ColSide Palette:", colSide_palette[1], "\n")}
}

matrix_palette = as.character(read.table(opt$matrix_palette, h=F, comment.char="%")$V1)
if (opt$verbose) {cat("Matrix Palette:", matrix_palette, "\n")}

#m = m[1:100,]

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
if (opt$verbose) {
	cat("Columns removed\n")
	sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
}
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# apply the log10 if needed
if (opt$log) {m <- log10(replace(m, is.na(m), 0) + opt$pseudocount)}

# melt the data frame
df = melt(as.matrix(m))

if (opt$verbose) {
	cat("Melt dataframe:\n")
	print(head(df))
}


# --------------- Metadata processing -------------

# read metadata and correct the merging fields
if (!is.null(opt$col_metadata)) {
	col_mdata = read.table(opt$col_metadata, h=T, sep="\t", quote="\"", comment.char="")
	col_mdata[opt$merge_col_mdata_on] <- make.names(col_mdata[,opt$merge_col_mdata_on])
}

if (!is.null(opt$row_metadata)) {
	row_mdata = read.table(opt$row_metadata, h=T, sep="\t", quote="\"", comment.char="")
	row_mdata[opt$merge_row_mdata_on] <- make.names(row_mdata[,opt$merge_row_mdata_on])

}

# --------------
# Palette files
# --------------

# Columns
if (!is.null(opt$colSide_by)) {
	colSide_by = strsplit(opt$colSide_by, ",")[[1]]
	# Check that there are enough color palettes for the column
	if (!is.null(opt$colSide_palette)) {
		if (length(colSide_palette) == 1) {
			colSide_palette = rep(colSide_palette, length(colSide_by))}
		if (length(colSide_palette) >1 && length(colSide_palette) != length(colSide_by))	{
			cat("ERROR: Inconsistent number of column factors and palettes\n");	q(save='no')}
	}
} else {
	colSide_by = NULL
}

# Rows
if (!is.null(opt$rowSide_by)) {
	rowSide_by = strsplit(opt$rowSide_by, ",")[[1]]
	# Check that there are enough color palettes for the column
	if (!is.null(opt$rowSide_palette)) {
		if (length(rowSide_palette) == 1) {
			rowSide_palette = rep(rowSide_palette, length(rowSide_by))}
		if (length(rowSide_palette) >1 && length(rowSide_palette) != length(rowSide_by))	{
			cat("ERROR: Inconsistent number of column factors and palettes\n");	q(save='no')}
	}
} else {
	rowSide_by = NULL
}

#if (!is.null(opt$rowSide_by)) {rowSide_by = strsplit(opt$rowSide_by, ",")[[1]]} else {rowSide_by = NULL}

# read which fields are needed from the metadata

if (!is.null(opt$col_labels) && opt$col_labels != "none") {col_label_fields = strsplit(opt$col_labels,",")[[1]]} else {col_label_fields=NULL}
if (!is.null(opt$row_labels) && opt$row_labels != "none") {row_label_fields = strsplit(opt$row_labels,",")[[1]]} else {row_label_fields=NULL}



col_mdata_header = unique(c(opt$merge_col_mdata_on, colSide_by, col_label_fields))
row_mdata_header = unique(c(opt$merge_row_mdata_on, rowSide_by, row_label_fields))
row_mdata_header = c(opt$merge_row_mdata_on, setdiff(row_mdata_header, intersect(col_mdata_header, row_mdata_header)))


same_mdata = FALSE
# Handle the fact that column metadata and row metadata may be the same file and they may be merged on the field
if (!is.null(opt$row_metadata) & !is.null(opt$col_metadata)) {
	if (opt$row_metadata == opt$col_metadata & opt$merge_row_mdata_on == opt$merge_col_mdata_on) {
		same_mdata = TRUE
	}
}

# Both metadata available
if (same_mdata) {
	# Handle the fact that column metadata and row metadata may be the same file and they may be merged on the field
	mdata_header = union(c(opt$merge_col_mdata_on, colSide_by, col_label_fields),
	c(opt$merge_row_mdata_on, rowSide_by, row_label_fields))
#	col_mdata[opt$merge_col_mdata_on] <- gsub(",", ".", col_mdata[,opt$merge_col_mdata_on])
	df = merge(df, unique(col_mdata[mdata_header]), by.x="Var2", by.y=opt$merge_col_mdata_on)
}else {

	# merge metadata and data (NB: The column Var2 stays)
	if (!is.null(opt$col_metadata)) {
#		col_mdata[opt$merge_col_mdata_on] <- gsub(",", ".", col_mdata[,opt$merge_col_mdata_on])
		df = merge(df, unique(col_mdata[col_mdata_header]), by.x="Var2", by.y=opt$merge_col_mdata_on)
	}
	
	if (!is.null(opt$row_metadata)) {
		row_mdata[opt$merge_row_mdata_on] <- gsub("[-,+]", ".", row_mdata[,opt$merge_row_mdata_on])
		df = merge(df, unique(row_mdata[row_mdata_header]), by.x="Var1", by.y=opt$merge_row_mdata_on)
		print(dim(row_mdata))
		if (length(rownames(m)) != length(intersect(rownames(m), row_mdata[,opt$merge_row_mdata_on]))) {
			cat("ERROR: Not all row names are matched in the row metadata!\n ABORTED\n")
			q(save="no")
		}
	}
}


if (opt$verbose) {cat("merged metadata\n")}
if (opt$verbose) {print(head(df))}

# ---------------- Dendrogram ----------------------

# COLUMNS

if (opt$col_dendro) {
	if (opt$verbose) {cat("column dendrogram... ")}
	if (opt$dist == "p" || opt$dist =="s") {
		colDist = as.dist(1-cor(m, method=opt$dist, use="p"))
	} else if (opt$dist == "i") {
	  colDist = as.dist(max(m) - m)
	} else if (opt$dist == "n") {
	  colDist = as.dist(m)
	} else {
		colDist = dist(t(m), method=opt$dist)
	}
	colHC = hclust(colDist, method=opt$hclust)
	colHC_data = dendro_data(as.dendrogram(colHC))
	col_ggdendro = ggplot(segment(colHC_data))
	col_ggdendro = col_ggdendro + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
	#col_ggdendro = col_ggdendro + coord_flip()
	col_ggdendro = col_ggdendro + scale_x_continuous(expand=c(0.0, 0.5), labels=NULL) 
	col_ggdendro = col_ggdendro + scale_y_continuous(expand=c(0.0, 0.0), labels=NULL)
	col_ggdendro = col_ggdendro + theme(plot.margin=unit(c(0.10, 0.00, 0.00, 0.01), "inch"))
	col_ggdendro = col_ggdendro + theme_dendro()
	col_ggdendro = col_ggdendro + labs(x=NULL, y=NULL)
	if (opt$verbose) {cat("DONE\n")}
}


# ROWS

if (opt$row_dendro) {
	if (opt$verbose) {cat("row dendrogram... ")}
	if (opt$dist == "p" || opt$dist =="s") {
		rowDist = as.dist(1-cor(t(m), method=opt$dist, use="p"))
	} else if (opt$dist == "i") {
	  rowDist = as.dist(max(m) - m)
	} else if (opt$dist == "n") {
	  rowDist = as.dist(m)
	} else {
		rowDist = dist(m, method=opt$dist)
	}
	rowHC = hclust(rowDist, method=opt$hclust)
	rowHC_data = dendro_data(as.dendrogram(rowHC))
	row_ggdendro = ggplot(segment(rowHC_data))
	row_ggdendro = row_ggdendro + geom_segment(aes(x=x, y=y, xend=xend, yend=yend))
	row_ggdendro = row_ggdendro + coord_flip()
	row_ggdendro = row_ggdendro + scale_x_continuous(expand=c(0.0, 0.5), labels=NULL) 
	row_ggdendro = row_ggdendro + scale_y_continuous(expand=c(0.0, 0.0), labels=NULL)
	row_ggdendro = row_ggdendro + theme(plot.margin=unit(c(0.00, 0.10, 0.00, 0.00), "inch"))
	row_ggdendro = row_ggdendro + theme_dendro()
	row_ggdendro = row_ggdendro + labs(x=NULL, y=NULL)
	if (opt$verbose) {cat("DONE\n")}
}



# ------------------- Row and column labels ---------

# Define the row and column labels
if (opt$verbose) {cat("column labels... ")}
if (is.null(opt$col_labels)) {
	col_labels = colnames(m)
} else {
	if (opt$col_labels == "none") {
		col_labels = NULL
	} else {
		col_label_fields[which(col_label_fields == opt$merge_col_mdata_on)] <- "Var2"
		col_labels = apply(df[col_label_fields], 1, paste, collapse=";")
		col_labels <- col_labels[with(df, match(colnames(m), Var2))]
	}
}
if (opt$verbose) {cat("DONE\n")}



if (opt$verbose) {cat("row labels... ")}
if (is.null(opt$row_labels)) {
	row_labels = rownames(m)
} else {
	if (opt$row_labels == "none") {
		row_labels = NULL
	} else {
		# Handle the fact that column metadata and row metadata may be the same file and they may be merged on the field
		if (same_mdata) {
			row_label_fields[which(row_label_fields == opt$merge_row_mdata_on)] <- "Var2"
			row_labels = apply(df[row_label_fields], 1, paste, collapse=";")
			row_labels <- row_labels[with(df, match(rownames(m), Var2))]
		} else {
			row_label_fields[which(row_label_fields == opt$merge_row_mdata_on)] <- "Var1"
			row_labels = apply(df[row_label_fields], 1, paste, collapse=";")
			row_labels <- row_labels[with(df, match(rownames(m), Var1))]
		}
	}
}
if (opt$verbose) {cat("DONE\n")}


if (opt$col_dendro) {
	col_limits = colnames(m)[colHC$order]
	col_labels = col_labels[colHC$order]
} else {
	col_limits = colnames(m)
	col_labels = col_labels
}

if (opt$row_dendro) {
	row_limits = rownames(m)[rowHC$order]
	row_labels = row_labels[rowHC$order]
} else {
	row_limits = rev(rownames(m))
	row_labels = rev(row_labels)
}



# This works in the X11 device
#row_labels_inches = max(strwidth(row_labels, units="in"))
#col_labels_inches = max(strwidth(col_labels, units="in"))

# This works in the pdf device
row_labels_inches = max(strwidth(row_labels, units="in", cex=base_size*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))
col_labels_inches = max(strwidth(col_labels, units="in", cex=base_size*(as.numeric(theme_get()$axis.text$size))*par()$cex/par()$ps))




# ------------------- Matrix plot -------------------

#if (!is.null(opt$matrix_legend_breaks)) {
#	breaks = as.numeric(strsplit(opt$matrix_legend_breaks, ",")[[1]])
#}

matrix_fill_limits = NULL	
if (!is.null(opt$matrix_fill_limits)) {
	opt$matrix_fill_limits = gsub("\\", "", opt$matrix_fill_limits, fixed=TRUE)
	matrix_fill_limits = as.numeric(strsplit(opt$matrix_fill_limits, ",")[[1]])
} 

p1 = ggplot(df, aes(x=Var2, y=Var1))
p1 = p1 + geom_tile(aes(fill=value))
p1 = p1 + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))
p1 = p1 + scale_x_discrete(expand=c(0,0), limits=col_limits, labels=col_labels)
p1 = p1 + scale_y_discrete(expand=c(0,0), limits=row_limits, labels=row_labels)
p1 = p1 + scale_fill_gradientn(colours=matrix_palette, limits=matrix_fill_limits)
#p1 = p1 + scale_fill_gradientn(colours=matrix_palette)
p1 = p1 + theme(plot.margin=unit(c(0.00, 0.00, 0.01, 0.01),"inch"))
p1 = p1 + labs(x=NULL, y=NULL)
p1 = p1 + guides(fill=guide_colourbar(
	title.position="top", 
	direction="horizontal", 
	title.hjust=0, 
	title=opt$matrix_legend_title
#	draw.ulim=FALSE,
#	draw.llim=FALSE
))
#if (!is.null(opt$matrix_legend_breaks)) {
#	p1 = p1 + guides(fill=guide_colourbar(breaks=breaks))
#}
if (opt$matrix_text) {
	p1 = p1 + geom_text( aes( x=Var2, y=Var1, label=value ))
}
p1_legend = g_legend(p1)
p1 = p1 + theme(legend.position = "none")



# -------------------- Column Side Colors ------------

ColSides = list(); ColSide_legends = list()
if (!is.null(opt$colSide_by)) {
	i=1;
	for (colSide in colSide_by) {
		colSide_data = unique(df[c("Var2", colSide)])
		ColSide = ggplot(colSide_data, aes(x=Var2, y="a"))
#		ColSide = ColSide + geom_tile(aes_string(fill=colSide), color="black")
		ColSide = ColSide + geom_tile(aes_string(fill=colSide))
		ColSide = ColSide + scale_x_discrete(limits = col_limits, labels=NULL, expand=c(0,0))
		ColSide = ColSide + scale_y_discrete(labels=NULL, expand=c(0,0))
		if (!is.null(opt$colSide_palette)) {
			ColSide = ColSide + scale_fill_manual(values=colSide_palette[[i]])
		} else {
			ColSide = ColSide + scale_fill_hue()
		}
		ColSide = ColSide + labs(x=NULL, y=NULL)
		ColSide = ColSide + theme(
			plot.margin=unit(c(0.00, 0.00, 0.00, 0.01),"inch"),
			legend.text=element_text(size=0.9*base_size),
			legend.key.size=unit(0.9*base_size, "points")
		)
		ColSide_legends[[i]] = g_legend(ColSide)
		ColSide = ColSide + theme(legend.position="none")
		ColSides[[i]] = ColSide; i=i+1;
	}
}



# ----------------------- Row Side Colors ------------

RowSides = list(); RowSide_legends = list()
if (!is.null(opt$rowSide_by)) {
	i=1;
	for (rowSide in rowSide_by) {
		# Handle the fact that column metadata and row metadata may be the same file and they may be merged on the field
		if (same_mdata) {
			rowSide_data = unique(df[c("Var2", rowSide)])
			RowSide = ggplot(rowSide_data, aes(x=Var2, y="a"))
		} else {
			rowSide_data = unique(df[c("Var1", rowSide)])
			RowSide = ggplot(rowSide_data, aes(x=Var1, y="a"))
		}
		RowSide = RowSide + geom_tile(aes_string(fill=rowSide))
		RowSide = RowSide + scale_x_discrete(limits = row_limits, labels=NULL, expand=c(0,0))
		RowSide = RowSide + scale_y_discrete(labels=NULL, expand=c(0,0))
		if (!is.null(opt$rowSide_palette)) {
#			RowSide = RowSide + scale_fill_manual(values=rowSide_palette)
			RowSide = RowSide + scale_fill_manual(values=rowSide_palette[[i]])
		} else {
			RowSide = RowSide + scale_fill_hue()
		}
		RowSide = RowSide + labs(x=NULL, y=NULL)
		RowSide = RowSide + coord_flip()
		RowSide = RowSide + theme(
			plot.margin=unit(c(0.00, 0.00, 0.00, 0.01),"inch"),
			legend.text=element_text(size=0.9*base_size),
			legend.key.size=unit(0.9*base_size, "points")
		)
		RowSide_legends[[i]] = g_legend(RowSide)
		RowSide = RowSide + theme(legend.position="none")
		RowSides[[i]] = RowSide; i=i+1;
	}
} 





# ============================
# Compose with viewports
# ============================

# >>>>> Matrix viewport <<<<<<<<<<<<<<<<<

matrix_vp_y = max(0, as.numeric(strwidth(rowSide_by, "in")) - col_labels_inches)
matrix_vp_x = max(0, as.numeric(strwidth(colSide_by, "in")) - row_labels_inches)
if (is.null(opt$height)) {matrix_vp_h = base_size/72.27 * nrow(m) + col_labels_inches} else {matrix_vp_h = opt$height}
if (is.null(opt$width)) {matrix_vp_w = base_size/72.27 * ncol(m) + row_labels_inches} else {matrix_vp_w = opt$width}
matrix_vp = viewport(
	y = matrix_vp_y,
	x = matrix_vp_x, 
	h = matrix_vp_h, 
	w = matrix_vp_w, 
	default.units="inch", 
	just=c("left","bottom")
)

# >>>>> Matrix scale viewport <<<<<<<<<<<

#matrix_scale_h = 0.5
#matrix_scale_w = 1
matrix_scale_h = sum(sapply(p1_legend$heights, convertUnit, "in"))
matrix_scale_w = sum(sapply(p1_legend$widths, convertUnit, "in"))
matrix_scale_vp = viewport(
	y = matrix_vp_y + matrix_vp_h + 0.01,
	x = matrix_vp_x + matrix_vp_w + 0.05,
	h = matrix_scale_h,
	w = matrix_scale_w,
	default.units = "inch",
	just = c("left", "bottom")
)


# >>>>>> Column side viewport <<<<<<<<<<<

if (!is.null(opt$colSide_by)){
	ColSide_vps = list(); ColSide_label_vps = list()
	for (i in 1:length(ColSides)) {
		ColSide_vps[[i]] = viewport(
			y = matrix_vp_y + matrix_vp_h + 0.25*(i-1), 
			x = matrix_vp_x + row_labels_inches,
			h = 0.25,
			w = matrix_vp_w - row_labels_inches,
			default.units = "inch",
			just = c("left", "bottom") 
		 )
		ColSide_label_vps[[i]] = viewport(
			y = matrix_vp_y + matrix_vp_h + 0.25*(i-1),
			x = matrix_vp_x + row_labels_inches,
			h = 0.25,
			w = max(row_labels_inches, as.numeric(strwidth(colSide_by, "inch") * (base_size/12) )),
			default.units = "inch",
			just = c("right", "bottom")
		)
	}
}


# >>>>> Row side viewport <<<<<<<<<<<<<<<<

if (!is.null(opt$rowSide_by)){
	RowSide_vps = list(); RowSide_label_vps = list()
	for (i in 1:length(RowSides)) {
		RowSide_vps[[i]] = viewport(
			y = matrix_vp_y + col_labels_inches, 
			x = matrix_vp_x + matrix_vp_w + 0.25*(i-1),
			h = matrix_vp_h - col_labels_inches,
			w = 0.25,
			default.units = "inch",
			just = c("left", "bottom") 
		 )
		RowSide_label_vps[[i]] = viewport(
			y = matrix_vp_y + col_labels_inches,
			x = matrix_vp_x + matrix_vp_w + 0.25*(i-1),
			h = max(col_labels_inches, as.numeric(strwidth(rowSide_by, "inch") * (base_size/12) )),
			w = 0.25,
			default.units = "inch",
			just = c("left", "top")
		)
	}
}


# >>> Column dendrogram viewport <<<<<<<<<

col_dendro_h = 0.1
if (opt$col_dendro) {
	col_dendro_h = 0.5
	colDendro_vp = viewport(
	y = matrix_vp_y + matrix_vp_h + 0.25*length(colSide_by),
	x = matrix_vp_x + row_labels_inches,
	h = col_dendro_h,
	w = matrix_vp_w - row_labels_inches,
	default.units = "inch",
	just = c("left", "bottom")
	)
}


# >>>> Row dendrogram viewport <<<<<<<<<<<<<<<<<<<<<<<<<<<<<

row_dendro_w = 0.1
if (opt$row_dendro) {
	row_dendro_w = 2.0
	rowDendro_vp = viewport(
	y = matrix_vp_y + col_labels_inches,
	x = matrix_vp_x + matrix_vp_w + 0.25*length(rowSide_by), 
	h = matrix_vp_h - col_labels_inches,
	w = row_dendro_w,
	default.units = "inch",
	just = c("left", "bottom")
	)
}


total_h = matrix_vp_y + matrix_vp_h + 0.25*length(colSide_by) + max(col_dendro_h, matrix_scale_h)
total_w = matrix_vp_x + matrix_vp_w + 0.25*length(rowSide_by) + max(row_dendro_w, matrix_scale_w)


# >>>>> Column and row side legends viewport <<<<<<<<<<<<<<<<<

max_legend_width_inch = 0
if (!is.null(opt$colSide_by) || !is.null(opt$rowSide_by)) {
	legend_padding = 0.10
	max_legend_width_inch = max(sapply(c(ColSide_legends, RowSide_legends), function(x) sum(sapply(x$widths, convertUnit, "in"))))
#	legend_height_inch = matrix_vp_h/length(c(colSide_by, rowSide_by))
#	legend_height_inch = total_h/length(c(colSide_by, rowSide_by))
	legend_y = matrix_vp_y + matrix_vp_h
	side_legend_vps = list()
	#for (i in 1:length(c(colSide_by, rowSide_by))) {
	all_side_legends = c(ColSide_legends, RowSide_legends)
	for (i in 1:length(all_side_legends)) {
		legend_height_inch = sum(sapply(all_side_legends[[i]]$heights, convertUnit, "in"))
		legend_width_inch = sum(sapply(all_side_legends[[i]]$widths, convertUnit, "in"))
		side_legend_vp = viewport(
#			y = matrix_vp_y + matrix_vp_h + 0.25*length(colSide_by) + col_dendro_h - legend_height_inch*(i-1),
			x = matrix_vp_x + matrix_vp_w + 0.25*length(rowSide_by) + row_dendro_w,
			y = legend_y - legend_padding,
			h = legend_height_inch,
			w = legend_width_inch,
			default.units = "inch",
			just = c("left", "top")
		)
		legend_y = legend_y - legend_height_inch
		side_legend_vps[[i]] = side_legend_vp
	}
} 

total_w = total_w + max_legend_width_inch


# =======================================
# PRINT PLOT
# =======================================

pdf(opt$output, h = total_h, w=total_w)

#X11(h=total_h, w=total_w)

# Print matrix
print(p1, matrix_vp, newpage=FALSE)

# Print matrix legend
pushViewport(matrix_scale_vp); grid.draw(p1_legend); upViewport()

# Print column side colors
if (!is.null(opt$colSide_by)) {
	for (i in 1:length(ColSide_vps)) {
		print(ColSides[[i]], vp=ColSide_vps[[i]], newpage=FALSE)
		grid.text(colSide_by[i], x=unit(1, "npc"), just="right", vp=ColSide_label_vps[[i]], gp=gpar(face="bold"))
	}
}

# Print row side colors
if (!is.null(opt$rowSide_by)) {
	for (i in 1:length(RowSide_vps)) {
		print(RowSides[[i]], vp=RowSide_vps[[i]], newpage=FALSE)
		grid.text(rowSide_by[i], y=unit(1,"npc"), just="right", rot=90, vp=RowSide_label_vps[[i]], gp=gpar(face="bold"))
	}
}

# Print column dendrogram
if (opt$col_dendro) {print(col_ggdendro, vp=colDendro_vp, newpage=FALSE)}

# Print row dendrogram
if (opt$row_dendro) {print(row_ggdendro, vp=rowDendro_vp, newpage=FALSE)}

# Print column and row side color scales
if (!is.null(opt$colSide_by) || !is.null(opt$rowSide_by)) {
	for (i in 1:length(c(colSide_by, rowSide_by))) {
		all_side_legends = c(ColSide_legends, RowSide_legends)
		pushViewport(side_legend_vps[[i]]); grid.draw(all_side_legends[[i]]); upViewport()
	}
}


dev.off()

file.remove("Rplots.pdf")
q(save="no")
