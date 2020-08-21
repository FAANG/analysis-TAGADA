#!/usr/bin/env Rscript

# This script was initially written by Alessandra Breschi at the CRG, Barcelona, Spain
# Its initial version is at https://github.com/abreschi/Rscripts

##------------
## LIBRARIES
##------------ 
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(plyr))


options(stringsAsFactors=F)

##################
# OPTION PARSING
##################


option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin", help="the matrix you want to analyze [default=%default]"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-o", "--output"), help="additional tags for otuput", default="out"),
make_option(c("-c", "--color_by"), help="choose the color you want to color by. Leave empty for no color", type='character'),
make_option(c("-y", "--linetype_by"), help="choose the factor you want the linetype by. Leave empty for no linetype", type="character"),
make_option(c("-f", "--file_sel"), help="list of elements of which computing the proportion at each point"),
make_option(c("--out_file"), help="store the coordinates in a file [default=%default]"),
make_option(c("-P", "--palette"), help="file with the colors"),
make_option(c("-t", "--tags"), help="choose the factor by which grouping the lines [default=%default]", default="labExpId")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
#print(opt)



##--------------------##
## CLUSTERING SAMPLES ##
##--------------------##
output = sprintf("TPM_fraction.%s", opt$output)

# 1. read the matrix from the command line
if (opt$input_matrix == "stdin") {inF = file("stdin")} else {inF = opt$input_matrix}
m = read.table(inF, h=T, sep="\t")


# Read color palette if present
if (!is.null(opt$palette)) {palette = read.table(opt$palette, h=F, comment.char="%", sep="\t")$V1}

# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}


# 2. sort all the values for each column
sortm <- apply(m, 2, sort, na.last=T, d=T)
#sortm <- sortm[1:min(which(sortm==0))-1,]

# 3. calculate cumulative sum
cumm <- as.data.frame(apply(sortm, 2, function(x) cumsum(x/sum(x,na.rm=T))))

# 4. read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t')
mdata[,"labExpId"] <- sapply(mdata[,"labExpId"], function(x) gsub(",", ".", x))
if (!is.null(opt$color_by)) {opt$color_by <- strsplit(opt$color_by, ",")[[1]]}

# prepare data.frame for ggplot
df = melt(cumm, variable.name = "labExpId", value.name="TPM_fraction")
fields = unique(c("labExpId", strsplit(opt$tags, ",")[[1]], opt$color_by, opt$linetype_by))
df = merge(unique(mdata[fields]), df, by="labExpId")
df$labels = apply(df[strsplit(opt$tags, ",")[[1]]], 1, paste, collapse="_")
# add a column with the x index
df  = ddply(df, .(labels), transform, x=seq_along(labels), y=sort(TPM_fraction, na.last=T, d=F))

# ===== Add an extra dataframe with the percentage of genes after the union =====

if (!is.null(opt$file_sel)) {

	prop_df = data.frame()
	sel = read.table(opt$file_sel, h=F)[,1]
	# Order the data
	orderm = apply(-m, 2, rank, na.last=T, ties.method='first')
	# NB: orderm has lost the rownames 
	thresholds = c(1:10 %o% 10^(0:4))
	thresholds = thresholds[which(thresholds < nrow(m))]

	props = sapply(thresholds, function(thr) {
		u = rownames(m)[which(rowSums(orderm <= thr) > 0)];
		length(intersect(u, sel))/length(u)
		}
	)
	prop_df = data.frame(x=thresholds, y=props)
}





###############
# OUTPUT 
###############

# coordinate file

if (!is.null(opt$out_file)) {
	write.table(df, file=opt$out_file, row.names=FALSE, quote=FALSE, sep="\t")
}

# plotting...
base_size=16
legend_nrow = 18
#legend_nrow = 4
theme_set(theme_bw(base_size=base_size))
legend_text_inch = theme_get()$legend.text$size * base_size / 72.72
add_w = legend_text_inch * max(nchar(df$labels)) * ceiling(length(levels(as.factor(df$labels)))/legend_nrow)

geom_params = list()
geom_params$size = opt$size
geom_params$alpha = opt$alpha


mapping = list()
mapping <- modifyList(mapping, aes(x=x, y=y, group=labels))

if (!is.null(opt$color_by)) {
	gp_color_by = interaction(df[opt$color_by])
    mapping = modifyList(mapping, aes(color=gp_color_by))
}

if (!is.null(opt$linetype_by)) {
	gp_linetype_by = interaction(df[opt$linetype_by])
	mapping = modifyList(mapping, aes(linetype=gp_linetype_by))
}
	

class(mapping) <- "uneval"

lineLayer <- layer(
        geom = "line",
#       geom_params = geom_params,
        params = geom_params,
        mapping = mapping,
        stat = "identity",
        position = "identity"
)



# GGPLOT

gp = ggplot(df) + lineLayer

#if (!is.null(opt$color_by)) {gp_color_by=interaction(df[opt$color_by])} else {gp_color_by=NULL}
#if (!is.null(opt$linetype_by)) {gp_linetype_by=interaction(df[opt$linetype_by])} else {gp_linetype_by=NULL}
#gp = gp + geom_line(aes(color=gp_color_by, linetype=gp_linetype_by, group=labels))

#gp = gp + scale_color_hue(name=paste(opt$color_by, collapse="."))
if (!is.null(opt$color_by)) {
	if (!is.null(opt$palette)) {
		gp = gp + scale_color_manual(values=palette)
	} else {
		gp = gp + scale_color_hue()
	}
}

gp = gp + labs(y="Fraction of gene TPM", x='Number of genes')
gp = gp + scale_linetype_manual(values=c(2,1))
gp = gp + guides(col = guide_legend(nrow = legend_nrow, title=opt$color_by))
gp = gp + scale_x_log10(expand=c(0,0))
gp = gp + scale_y_continuous(expand=c(0.01,0), limits=c(0,1))
gp = gp + annotation_logticks(sides="b")

if (!is.null(opt$file_sel)) {
gp = gp + geom_point(data=prop_df, aes(x,y), shape=18, size=2)
gp = gp + geom_point(data=prop_df, aes(x,y), shape=18, size=1.7, color='yellow')
}

ggsave(sprintf("%s.pdf",output), h=5, w=6+add_w, title=output)

q(save='no')
