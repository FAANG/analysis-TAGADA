#!/usr/bin/env Rscript

# This script was initially written by Alessandra Breschi at the CRG, Barcelona, Spain
# Its initial version is at https://github.com/abreschi/Rscripts

options(stringsAsFactors=F)

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))
# suppressPackageStartupMessages(library("funr"))
# exec_dir <- dirname(sys.frame(1)$ofile)
# print(exec_dir)
# args <- commandArgs()
# args <- grep("--file", "rpkm_distribution.R", value=TRUE)
# print(gsub("--file=", "", args))


option_list <- list(

make_option(c("-i", "--input_matrix"), 
	help="the matrix you want to analyze. Can be stdin"),

make_option(c("-l", "--log"), action="store_true", default=FALSE, 
	help="apply the log10 [default=%default]"),

make_option(c("-p", "--pseudocount"), type="double", default = 1e-04,
	help="specify a pseudocount for the log [default=%default]"),

make_option(c("-m", "--metadata"), 
	help="tsv file with metadata on matrix experiment. If empty, column names are used as labels."),

make_option(c("-r", "--representation"), default = "boxplot",
	help="choose the representation <boxplot>, <density>, <histogram>, <boxviol> [default=%default]"),

make_option(c("-o", "--output"), default="out",
	help="additional tags for otuput [default=%default]"),

make_option(c("-d", "--outdir"), default="./",
	help="specificy a directory for the output [default=%default]"),

make_option(c("-f", "--fill_by"), help="choose the color you want to fill by. Leave empty for no fill", type='character'),
make_option(c("-a", "--alpha_by"), help="choose the fator for the transparency in boxplot. Leave empty for no transparency", type="character"),
make_option(c("-v", "--x_title"), help="give a title to the x axis [default=%default]", default="rpkm"), 
make_option(c("-w", "--wrap"), action="store_true", help="say if you want the density plot in facet_wrap [default=%default]", default=FALSE), 
make_option(c("-T", "--title"), default="", 
	help="give a title to the plot [default=%default]"), 

make_option(c("-P", "--palette"), default="palettes/cbbPalette.8.txt", 
	help="palette file name [default=%default]"),

make_option(c("-H", "--height"), default=10, 
	help="height of the plot in inches [default=%default]"),

make_option(c("-W", "--width"), default=10, 
	help="width of the plot in inches [default=%default]"),

make_option(c("-G", "--no_guide"), default=FALSE, action="store_true", 
	help="use this to remove the legend from the plot [default=%default]"),

make_option(c("--verbose"), default=FALSE, action="store_true", 
	help="use this for verbose mode [default=%default]"),

make_option(c("-t", "--tags"), default="labExpId", 
	help="comma-separated field names you want to display in the labels. Leave default for using column names [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


##------------
## LIBRARIES
##------------ 


if (opt$verbose) {cat("Loading libraries... ")}
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plyr))
if (opt$verbose) {cat("DONE\n\n")}


#~~~~~~~~~~~~~~
# BEGIN
#~~~~~~~~~~~~~~

# read the matrix from the command line
if (opt$input_matrix == "stdin") {
	m = read.table(file("stdin"), h=T)
} else {
	m = read.table(opt$input_matrix, h=T)
}

ylab = opt$x_title

# Read palette
palette = read.table(opt$palette, h=F, comment.char="%")$V1
if (opt$verbose) {cat("Palette: ", palette, "\n")}


# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
if (opt$verbose) {sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)}
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# substitute the matrix with its log if required by the user
if (opt$log) {
	m = log10(replace(m, is.na(m), 0) + opt$pseudocount);
}

# prepare data.frame for ggplot
tags = strsplit(opt$tags, ",")[[1]]
df = melt(as.matrix(m), varnames = c("element","labExpId"), value.name="rpkm")

if(opt$verbose) {
	print(head(df))
}

# read the metadata from the metadata file if present
if (!is.null(opt$metadata)) {

	mdata = read.table(opt$metadata, h=T, sep='\t')
	mdata$labExpId <- sapply(mdata$labExpId, function(x) gsub(",", ".", x))
	
	# prepare data.frame for ggplot
	df = merge(unique(mdata[c("labExpId", tags, opt$fill_by, opt$alpha_by)]), df, by="labExpId")
	
	# change color palette in case length is too much
	if (length(unique(df[,opt$fill_by])) > length(palette)) {palette <- rainbow(length(unique(df[,opt$fill_by])))}

	# change the fill_by column to factor
	if (!is.null(opt$fill_by)) {
		df[,opt$fill_by] = as.factor(df[,opt$fill_by])
	}
}


df$labels = apply(df[tags], 1, paste, collapse="_")

if (opt$verbose) {
	print(head(df))
}


###############
# OUTPUT 
###############

output = sprintf("%s/%s.log_%s.psd_%s.%s", opt$outdir, opt$representation, 
ifelse(opt$log, "T", "F"), ifelse(opt$log,opt$pseudocount,"NA"),opt$output)

# plotting...

theme_set(theme_bw(base_size=14))

width = 7; height = 7;

medians = aggregate(rpkm~labels, df, function(x) median(x[is.finite(x)], na.rm=T))
lev = medians[order(medians$rpkm),]$labels
df$labels = factor(df$labels, levels=lev)

if(opt$verbose) {print(medians)}

width = opt$width
height = opt$height
breaks = floor(min(df[is.finite(df[,"rpkm"]),"rpkm"])):ceiling(max(df[is.finite(df[,"rpkm"]),"rpkm"]))


if (opt$representation == "boxplot") {

	gp = ggplot(df, aes_string(x="labels", y="rpkm"))
	gp = gp + geom_boxplot(outlier.size=0.75, size=0.5, aes_string(fill=opt$fill_by, alpha=opt$alpha_by))
	gp = gp + labs(y=ylab, x='', title=opt$title)
	gp = gp + coord_flip()
	gp = gp + scale_fill_manual(values = palette)
	gp = gp + scale_alpha_manual(values=c(.5,1))
	gp = gp + scale_y_continuous(breaks = breaks)
}

if (opt$representation == "boxviol") {
gp = ggplot(df, aes_string(x="labels", y="rpkm"))
gp = gp + geom_violin(size=.4,aes_string(fill=opt$fill_by))
gp = gp + geom_boxplot(alpha=0.2, size=.2)
gp = gp + labs(y=ylab, x='', title=opt$title)
gp = gp + coord_flip()
gp = gp + scale_fill_manual(values = palette)
gp = gp + stat_summary(fun.y="mean", colour="red", geom="point", size = 3)
gp = gp + scale_y_continuous(breaks = breaks)
}

if (opt$representation == "density") {
gp = ggplot(df, aes_string(x="rpkm"))
if (opt$wrap) {
gp = gp + geom_density(aes_string(color=if(is.na(opt$fill_by)){NULL}else{opt$fill_by}))
gp = gp + facet_wrap(~labels)} else {
gp = gp + geom_density( aes_string(group="labels", color=if(is.na(opt$fill_by)){NULL}else{opt$fill_by}) )}
gp = gp + labs(x=ylab, title=opt$title)
gp = gp + scale_color_manual(values = palette)
gp = gp + scale_x_continuous(breaks = breaks)
}

if (opt$representation == "histogram") {
df$labels = apply(df[tags], 1, paste, collapse="\n") # only in this case collapse with new line
gp = ggplot(df, aes_string(x="rpkm"))
gp = gp + geom_histogram(aes_string(fill=if(is.na(opt$fill_by)){NULL}else{opt$fill_by}, y='..count..'), right=TRUE, 
boundary=min(df[is.finite(df[,"rpkm"]),"rpkm"]))
gp = gp + facet_wrap(~labels)
gp = gp + labs(x=ylab, y='', title=opt$title)
gp = gp + scale_fill_manual(values = palette)
gp = gp + scale_x_continuous(breaks = breaks)
}

gp = gp + theme(axis.text=element_text(size=40/log10(length(df$labels))))


if (opt$no_guide) {gp = gp + theme(legend.position="none")}


pdf(sprintf("%s.pdf", output), w=width, h=height); gp; dev.off()

q(save='no')
