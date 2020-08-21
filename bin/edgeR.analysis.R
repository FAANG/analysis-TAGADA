#!/usr/bin/env Rscript 

# This script was initially written by Alessandra Breschi at the CRG, Barcelona, Spain
# Its initial version is at https://github.com/abreschi/Rscripts

options(stringsAsFactors=F)


#########################
# Default OPT parameters
########################

opt = list()
opt$n_samples = 2
opt$cpm = 1
opt$replace_NAs = TRUE

##################
# OPTION PARSING
##################

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option(c("-i", "--input_matrix"), default="stdin", help="the matrix with READ COUNTS you want to analyze [default=%default]"),
make_option(c("-r", "--replace_NAs"), action="store_true", default=FALSE, help="replace NAs with 0 [default=%default]"),
make_option(c("-m", "--metadata"), help="tsv file with metadata on matrix experiment"),
make_option(c("-f", "--fields"), help="choose the fields you want to use in the differential expression, comma-separated"),
make_option(c("-F", "--formula"), default="~.", help="a formula for the differential expression [default=%default]"),
make_option(c("-T", "--tissue_sp"), default=FALSE, action="store_true", help="Build a linear model with grand mean intercept"),
make_option(c("-C", "--cpm"), help="threshold for cpm [default=%default]", type="integer", default=1),
make_option(c("-s", "--n_samples"), help="minimum number of samples with cpm > \"cpm\" [default=%default]", type="integer", default=2),
make_option(c("-c", "--coefficient"), type="character", help="the coefficient(s) for the comparison"),
make_option(c("-n", "--contrast"), help="Numbers comma-separated giving the vector of contrasts. Use instead of coefficients"),
make_option(c("-o", "--output"), help="additional suffix for output [default=%default]", default='out'),
make_option(c("-d", "--output_dir"), help="choose the output directory [default=%default]", default="./"),
make_option(c("-g", "--genes"), help='a file with a list of genes to filter', type='character'),
make_option(c("-v", "--verbose"), action="store_true", default=FALSE, help="verbose output")
)

parser <- OptionParser(usage = "%prog [options] file", option_list=option_list)
arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options
if (opt$verbose) {print(opt)}


##------------
## LIBRARIES
##------------ 


if(opt$verbose) {cat('Libraries loading... ')}

suppressPackageStartupMessages(library('reshape2'))
suppressPackageStartupMessages(library('ggplot2'))
#suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library('edgeR'))

if (opt$verbose) {cat('DONE\n\n')}

# ==========================================
# Function for loading Rdata
# ==========================================

load_obj <- function(f)
{
    env <- new.env()
    nm <- load(f, env)[1]
    env[[nm]]
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function that estimates the average of gene counts
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mglmAverage = function(object, lev,  prior.count = 0.125, dispersion = 'auto') {
	# Evaluate dispersion argument
    dispersion <- switch(dispersion, common = object$common.dispersion,
        trended = object$trended.dispersion, tagwise = object$tagwise.dispersion,
        auto = getDispersion(object))
    # Get the library size 
    lib.size <- object$samples$lib.size * object$samples$norm.factors
    # Divide each library size by the mean library size across ALL conditions
    # they will be pseudo counts for all genes in each sample
    prior.count <- prior.count * lib.size/mean(lib.size)
    # Take the log of the library size (I don't understand why they add prior.count since it 
    # is negligible compared to the library size)
    offset.aug <- log(lib.size + 2 * prior.count)

    # Indeces of the specified level
    j <- object$samples$group == lev
    # Number of elements in level
    n <- sum(j)
    # Number of elements for which there are counts
    ntags <- nrow(object$counts)
    # Subset the matrix of counts
    y <- object$counts[, j, drop=FALSE]

    mglmAvg = mglmOneGroup(y + matrix(prior.count[j], ntags, n, byrow = TRUE), 
        offset = offset.aug[j], dispersion = dispersion)
	names(mglmAvg) <- rownames(object$counts)
    return(mglmAvg)
}



##-------------------------##
## Differential Expression ##
##-------------------------##

# specify the design to the program
fields = strsplit(opt$fields, ",")[[1]]

# Report error and quit if contrasts and coefficients are both provided
if (!is.null(opt$contrast) & !is.null(opt$coefficient)) {
	cat("ERROR: Provide contrasts OR coefficient!\nABORTED\n")
	q(save='no')
} 
if (is.null(opt$contrast) & is.null(opt$coefficient) & length(fields)>1) {
	cat("ERROR: Provide contrasts or coefficient!\nABORTED\n")
	q(save='no')
}

coeff = NULL
if (!is.null(opt$coefficient)) {
	coeff = eval(parse(text=opt$coefficient))
}

contr = NULL
if (!is.null(opt$contrast)) {
	contr <- eval(parse(text=opt$contrast))
} 

# read the matrix from the command line
#if (opt$input_matrix == "stdin") {inF = file("stdin")} else {inF=opt$input_matrix}
#m = read.table(inF, h=T, sep="\t")
if (opt$input_matrix == "stdin") {
    m = read.table(file("stdin"), h=T)
} else {
    m <- try(load_obj(opt$input_matrix), silent=T)
    if (class(m) == "try-error") {m <- read.table(opt$input_matrix)}
}

if (opt$replace_NAs) {m<-replace(m, is.na(m), 0)}
genes = rownames(m)
m = (apply(m, 2, as.integer))
rownames(m) <- genes

if (opt$verbose) {print(head(m))}


# Keep only selected genes if needed
if (!is.null(opt$genes)) {
fgenes = as.character(read.table(opt$genes, h=F)$V1)
m = m[intersect(rownames(m),fgenes),]
}


# read the metadata from the metadata file
mdata = read.table(opt$metadata, h=T, sep='\t', quote=NULL)
mdata[,"labExpId"] <- gsub("[,-]", ".", mdata[,"labExpId"])
mdata = unique(mdata[,c("labExpId", fields)])
mdata = mdata[match(colnames(m), mdata[,"labExpId"]), c("labExpId", fields)]


# ONLY ONE CONDITION
#if (length(fields) == 1) {
if (length(unique(unlist(mdata[fields]))) == 2) {
condition = factor(mdata[match(colnames(m), mdata[,"labExpId"]), fields])
if (opt$verbose) {print(condition)}

# create count object for edgeR
M = DGEList(counts=na.omit(m), group = condition)
if (opt$verbose) {cat(nrow(M$counts), "/", nrow(m), "did not contain NAs\n\n")}

# This also estimates the library size as the sum of all reads in each sample
# To change the library size
# M$samples$lib.size <- colSums(M$counts)

# TO DO: Add option to provide the library size 


# ~~~~ Normalisation ~~~~~

# Divide by the library size
cpm.m <- cpm(M)
# Filter by cpm and number of samples where the gene is expressed
M <- M[rowSums(cpm.m > opt$cpm) >= opt$n_samples,]
if(opt$verbose) {cat(nrow(M$counts), "passed the required filters\n\n")}

# Normalize by TMM
M <- calcNormFactors(M, method="TMM")

# variance estimation
M <- estimateCommonDisp(M, verbose=T)
M <- estimateTagwiseDisp(M)

# calling differential expression
et <- exactTest(M)
res = topTags(et, n=nrow(et))$table

# Add the estimated averages for each level (the function is a subset of exactTest)
for(lev in levels(condition)) {
#	avg = mglmAverage(M, lev)
#	res[lev] = round(mglmAverage(M, lev)[rownames(res)], digits=2)
	res = cbind(avg=round(mglmAverage(M, lev)[rownames(res)], digits=2), res)
	colnames(res)[which(colnames(res) == "avg")] <- lev
}

# MULTIPLE CONDITIONS

}else{
#design_df = mdata[mdata[,"labExpId"] %in% colnames(m),]
#design_df = droplevels(design_df)
#design_df = design_df[order(design_df$labExpId),]
#design = model.matrix(as.formula(opt$formula), design_df[fields])
design = model.matrix(as.formula(opt$formula), mdata[fields])
#print(colnames(design))
rownames(design) <- mdata[, "labExpId"]

if (opt$tissue_sp) {
	mdata[,fields] <- factor(mdata[,fields])
	contrasts(mdata[,fields]) <- contr.sum(levels(mdata[,fields]))
	design = model.matrix(as.formula(opt$formula), mdata[fields])
	rownames(design) <- mdata[, "labExpId"]
}

print(colnames(design))

cat('\n')


# create count object for edgeR
M = DGEList(counts=na.omit(m))

# normalisation
cpm.m <- cpm(M)
M <- M[rowSums(cpm.m > opt$cpm) >= opt$n_samples,]
#M$samples$lib.size <- colSums(M$counts)
M <- calcNormFactors(M)

# variance estimation
M <- estimateGLMCommonDisp(M, design, verbose=T)
M <- estimateGLMTrendedDisp(M, design)
M <- estimateGLMTagwiseDisp(M, design)


# calling differential expression
if (opt$tissue_sp) {
	fit <- glmQLFit(M, design)
	lrt <- glmQLFTest(fit, coef=coeff, contrast=contr)
} else {
	fit <- glmFit(M, design)
	lrt <- glmLRT(fit, coef=coeff, contrast=contr)
}
res = topTags(lrt, n=nrow(lrt))$table
}

# format the otuput and write to a file
if(length(coeff)<2) {res$logFC <- round(res$logFC, 2)}
res$logCPM <- round(res$logCPM, 2)
res$PValue <- format(res$PValue, digits=2)
res$FDR <- format(res$FDR, digits=2)

output = sprintf("%s/edgeR.cpm%s.n%s.%s", opt$output_dir, opt$cpm, opt$n_samples, opt$output)
write.table(res, file=sprintf("%s.tsv", output), quote=F, sep='\t',row.names=T)
#write.table(cpm.m, file=sprintf("cpm.%s.tsv",opt$output), quote=F, sep='\t',row.names=T)

#------#
# PLOT #
#------#

pdf(sprintf('%s.pdf', output))

if (length(coeff)<2) {
	# MA plot
	gp = ggplot(res, aes(x=logCPM, y=logFC))
	gp = gp + geom_point(shape=1, aes(color=cut(as.double(FDR), breaks=c(0, 0.01, 0.05, 1))))
	gp = gp + scale_color_brewer(name="FDR", palette="Set1")
	gp = gp + geom_hline(yintercept=c(2,-2))
	gp
	
	# Volcano plot
	gp = ggplot(res, aes(x=logFC, y=-log10(as.double(FDR))))
	gp = gp + geom_point(shape=1, aes(color=cut(as.double(FDR), breaks=c(0, 0.01, 0.05, 1))))
	gp = gp + scale_color_brewer(name="FDR", palette="Set1")
	gp = gp + scale_y_log10()
	gp = gp + geom_vline(xintercept=c(2,-2))
	gp
	
}

# plot dispersion estimates
plotBCV(M)
dev.off()










q(save='no')


