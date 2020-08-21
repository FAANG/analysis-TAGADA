#!/usr/bin/env Rscript

# This script was initially written by Alessandra Breschi at the CRG, Barcelona, Spain
# Its initial version is at https://github.com/abreschi/Rscripts

# Note for self: used in the senescence project and with R 3.6 on genologin installing in /work2/project/fragencode/.local/R/3.6

# -- Variables --

options(stringsAsFactors=F)
pseudocount = 1e-04
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-u", "--universe"), help="a list of gene identifiers (entrez gene ids), NO header"),
make_option(c("-G", "--genes"), default="stdin",
	help="a list of gene identifiers for the foreground (entrez gene ids), WITHOUT header [default=%default]"),
make_option(c("-c", "--categ"), help="choose the GO category < BP | MF | CC > [default=%default]", default="BP"),
make_option(c("-o", "--output"), help="additional tags for otuput [default=%default]", default="out"),
make_option(c("-d", "--output_dir"), default="./", help="directory for the output [default=%default]"),
make_option(c("-v", "--verbose"), default=FALSE, action="store_true", help="verbose")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}



##------------
## LIBRARIES
##------------

if (opt$verbose) {cat("Loading libraries... ")}


suppressPackageStartupMessages(library("GO.db"))
suppressPackageStartupMessages(library("hgu95av2.db"))
suppressPackageStartupMessages(library("GOstats"))
suppressPackageStartupMessages(library("plyr"))

if (opt$verbose) {cat("DONE\n\n")}

############################
# BEGIN
############################

U = read.table(opt$universe, h=F, col.names='hs')

if (opt$genes == "stdin") {
	G = read.table(file("stdin"), h=F, col.names='hs')
} else {
	G = read.table(opt$genes, h=F, col.names='hs')
}


# I want to create a list of parameters to perform GO enrichment on different gene sets

# take the entrez gene ids for all the orthologous genes which will be my universe (the same for all the sets)
universe = U$hs


if (opt$verbose) {sprintf("%s background genes; %s with a corresponding entrez id", nrow(U), length(unique(universe)))}
# how many genes am I able to map?

createParams = function(x) {
		ann = "hgu95av2.db"
		geneset = x
		sprintf("%s foreground genes; %s with a corresponding entrez id", length(x), length(unique(geneset)))
		pv = 1-(1-0.05)**(1/length(x))
		params = new("GOHyperGParams",
			geneIds = geneset,
			universeGeneIds = universe,
			annotation = ann,
			ontology = opt$categ,
			pvalueCutoff = pv,
			conditional = TRUE,
			testDirection='over')
	return(params)}

res = hyperGTest(createParams(G$hs))

# Reformat the output table
df = summary(res)
df$Pvalue = format(df$Pvalue, digits=1)
df$OddsRatio <- round(df$OddsRatio, 2)
df$ExpCount <- round(df$ExpCount, 2)

# Print output
output = sprintf("%s/%s.%s", opt$output_dir, opt$output, opt$categ)
write.table(df, file=sprintf("%s.tsv", output), quote=F, sep="\t", row.names=F)
htmlReport(res, file=sprintf("%s.html", output))

q(save='no')
