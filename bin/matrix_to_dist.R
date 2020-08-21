#!/usr/bin/env Rscript

# This script was initially written by Alessandra Breschi at the CRG, Barcelona, Spain
# Its initial version is at https://github.com/abreschi/Rscripts


options(stringsAsFactors=F)

# DEBUGGING OPTIONS

opt = list()
opt$input_matrix = "~abreschi/Documents/db/human/gencode18/Flux/encode.promo.human.gene.RPKM.idr_01.thr_0.names_False.tsv" 
opt$verbose = TRUE
opt$log10 = FALSE
opt$pseudocount = 1e-03
opt$cor = "pearson"

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin",
	help="the matrix you want to analyze. \"stdin\" for reading from standard input [default=%default]"),
make_option(c("-v", "--verbose"),action="store_true", default=FALSE, 
	help="if you want detailed output"),
make_option(c("-l", "--log10"), action="store_true", default=FALSE, 
	help="apply the log10"),
make_option(c("-p", "--pseudocount"), type="double", default=0,
	help="specify a pseudocount for the log [default=%default]. NAs are replaced by 0s"),
make_option(c("-c", "--cor"),
	help="choose the correlation method"),
make_option(c("-d", "--dist"), 
	help="choose the distance method"),
make_option(c("-k", "--keep_na"), action="store_true", default=FALSE,
	help="use this if you want to keep the NAs. By default they are converted to 0. [default=%default]"),
make_option(c("-o", "--output"), default="stdout",
	help="a name for the output. \"stdout\" to print on standard output [default=%default]")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options

# Stop if both the correlation method and the distance are different from NULL
if (!is.null(opt$cor) & !is.null(opt$dist)) {
	cat("ERROR: cannot compute both the correlation and the distance\n")
	q(save="no")
}

if (opt$verbose) {print(opt)}

##############
# BEGIN
##############

# read input table
if (opt$input_matrix == "stdin") {input=file("stdin")} else {input=opt$input_matrix}
m = read.table(input, h=T)


# remove potential gene id columns
char_cols <- which(sapply(m, class) == 'character')
if (opt$verbose) {sprintf("WARNING: column %s is character, so it is removed from the analysis", char_cols)}
if (length(char_cols) == 0) {genes = rownames(m)}
if (length(char_cols) != 0) {genes = m[,char_cols]; m = m[,-(char_cols)]}

# Replace NAs with 0s if not asked otherwise
if (!(opt$keep_na)) {m = replace(m, is.na(m), 0)}

# apply the log if required
if (opt$log10) {m = log10(m + opt$pseudocount)}


# compute the correlation
if (!is.null(opt$cor)) {
	df = cor(m, use='p', method=opt$cor)
}

# compute the distance
if (!is.null(opt$dist)) {
	if (opt$dist == "p") {opt$dist <- "pearson"}
	if (opt$dist == "s") {opt$dist <- "spearman"}
	if (opt$dist != "spearman" & opt$dist != "pearson") {
		df = as.matrix(dist(t(m), method=opt$dist))
	} else {
		df = 1 - abs(cor(m, use="p", method=opt$dist))
	}
}

dec = 3

df = round(df, dec)

# print the results
if (opt$output == "stdout") {
	write.table(df, file="", quote=FALSE, sep="\t")
} else {
	write.table(df, file=opt$output, quote=FALSE, sep="\t")
}

q(save="no")
