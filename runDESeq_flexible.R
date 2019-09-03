# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 7 May 2019
# Title: runDESeq_flexible.R
# Goal: To run DESeq on the quant files
# Usage: Rscript runDESeq.R {analysisID} {/path/to/quant} {tx2gene} {metadata} {dds design}
# ------------------------------------------------------------------
# PREAMBLE
# This script is based off my original deSeq script, but it is more flexible in
# that you can decide to only include certain samples. You can also choose
# which design formula you use (e.g. you could just use ~habitat if you don't
# have any pairs)
# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------
print(paste0("@@@ commencing runDESeq.R analysis at ", date()))
# load arguments and verify correct number
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  cat("ERROR: Incorrect number of arguments!\n")
  cat("usage: Rscript myscript.R arg1 arg2 arg3 arg4 arg5\n")
  cat("example: Rscript runDESeq.R {analysisID} {/path/to/quant} {tx2gene} {metadata} {dds design}\n")
  quit()
}

# process arguments
analysisID <- args[1]
tx2gene <- read.csv(args[3])
metadata <- read.csv(args[4], sep = "\t", row.names = 1)
# if the pair number is numeric, DESeq won't like it, so check and edit if required
if (any(is.numeric(metadata$pair))) {metadata$pair <- paste("pr", metadata$pair, sep="_")}
sampleID <- rownames(metadata) #needed to ensure the txImport files are done in the right order, must be done before the "sub" step
rownames(metadata) <- sub("-", ".", rownames(metadata)) #R doesn't like dashes
quant <- Sys.glob(file.path(args[2],(paste0(sampleID,"*quant.sf")))) # Sys.glob() required for wildcard expansion
design <- args[5]
print("@@@ arguments loaded")

# load packages (and get citations)
bioc_packages <- c("tximport", "sva", "DESeq2", "pheatmap", "apeglm", "genefilter")
r_packages <- c("RColorBrewer", "ggplot2", "pheatmap")
## function to load R packages
baseRpkgTest <- function(x) {
  if (!suppressMessages(require(x,character.only = TRUE, quietly = T))) {
    install.packages(x,dep=TRUE, repos = "https://pbil.univ-lyon1.fr/CRAN/")
    if(!require(x,character.only = TRUE, quietly = T)) stop (paste0(x, "package not found"))
  }
}
## function to load bioconductor packages
biocondpkgTest <- function(x) {
  if (!require(x,character.only = TRUE)) {
    source("http://www.bioconductor.org/biocLite.R")
    biocLite(x)
    if(!require(x,character.only = TRUE)) stop (paste0(x, "bioconductor package not found"))
  }
}
## load packages
for (b_pkg in bioc_packages) {
  biocondpkgTest(b_pkg)
  sink("citations.txt", append=TRUE)
  print(citation(b_pkg))
  sink()
}
for (r_pkg in r_packages) {
  baseRpkgTest(r_pkg)
  sink("citations.txt", append=TRUE)
  print(citation(r_pkg))
  sink()
}
print("@@@ packages loaded - step 1 of 9")

# ------------------------------------------------------------------
# RUN TXIMPORT
# ------------------------------------------------------------------
txi <- tximport(quant, type = "salmon", tx2gene = tx2gene)
colnames(txi$counts) <- rownames(metadata)
print("@@@ tximport complete - step 2 of 9")

# ------------------------------------------------------------------
# TXIMPORT TO DESEQ
# ------------------------------------------------------------------
# get the data into deseq
dds <- DESeqDataSetFromTximport(txi, metadata, design = formula(paste("~",design)))

# prefiltering (optional)
# Benefits: reduces memory use, analysis time (according to DESeq vignette) - not required with my smaller dataset (independent filtering is performed automatically)
#dds <- dds[rowSums(counts(dds)> 1) >= 3, ] # keep genes with >=10 counts and >=3 species' if you prefer to keep genes with total rowsum of >1, use: dds <- dds[ rowSums(counts(dds)) > 1, ]

# save a list of normalised and unnormalised counts for interest
# save list of counts and kept rownames, this will be useful for making a "background set" of expressed genes for enrichment
write.table(counts(dds), file=paste0(analysisID, "_raw-counts.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
dds_sizefact <- estimateSizeFactors(dds)
write.table(counts(dds_sizefact, normalized=TRUE), file=paste0(analysisID, "_norm-counts.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)

print("@@@ data sent to DESeq - step 3 of 9")

# ------------------------------------------------------------------
# PRE-DESEQ DATA VISUALISATION AND EXPLORATION
# ------------------------------------------------------------------
# PREAMBLE
# Data viz is performed on variance-transformed data because the level of
# variance across an RNASeq dataset is not naturally stable across the whole
# range (we see genes with higher counts = greater variance). Transformation
# standardises variance across the range of mean counts for all genes but
# doesn't artificially amplify the variance of *small* count genes (as would
# happen if you just did a log transformation and add a pseudocount to all
# genes). If you have >30 samples, consider using vst instead of rlog
# transformation (p14).
#
# Here we set blind to FALSE because "blind dispersion estimation is not the
# appropriate choice if one expects that many or the majority of genes (rows)
# will have large differences in counts which are explainable by the
# experimental design" (DESeq2 vignette) - as we may expect to be the case with
# cross-species DGE.
# ------------------------------------------------------------------
###############################
## Count data transformation ##
###############################
vsd <- vst(dds, blind=FALSE)

###################################
## Heatmaps - top 100 gene counts##
###################################
select <- order(rowMeans(counts(dds_sizefact,normalized=TRUE)),
                decreasing=TRUE)[1:100]
selectmore <- order(rowMeans(counts(dds_sizefact,normalized=TRUE)),
                decreasing=TRUE)[1:1000]

df <- as.data.frame(colData(dds)[,c("pair","habitat")])
pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df,
         cutree_cols = 3, cutree_rows = 2,
         main = "100 genes with highest row means",
         filename=paste0(analysisID, "_pheatmap_top100genes.pdf"))

# or a prettier one, for species with sex data
if ("A" %in% levels(metadata$sex) || "U" %in% levels(metadata$sex)) {
  df2 <- as.data.frame(colData(dds)[,c("pair","habitat","sex")])
  ann_colors = list(
    habitat = c(epi = "black", hypo = "grey96"),
    sex = c(F = "deeppink1", M = "dodgerblue1", U = "gray43"))
  ###pair = c(pr_10 = "lightblue4", pr_11 = "red", pr_13 = "lightcoral", pr_15 = "yellow3", pr_16 = "palegreen3", pr_17 = "blue", pr_2 = "skyblue2", pr_3 = "thistle3", pr_5 = "turquoise", pr_9 = "orange2"))
  pheatmap(assay(vsd)[select,], cluster_rows=TRUE, show_rownames=FALSE,
           cluster_cols=TRUE, annotation_col=df2, annotation_colors = ann_colors,
           cutree_rows = 2,
           main = "100 genes with highest row means",
           fontsize = 6,
           filename=paste0(analysisID, "_pheatmap_top100genes_v2.pdf"))
}

#########################################
## Heatmaps - deviation between samples##
#########################################
# Heatmap of relative VST-transformed values across samples. Pair and habitat
# information are shown with colored bars at the top of the heatmap. Blocks of
# genes that covary across samples.

topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("pair","habitat")])
pheatmap(mat, annotation_col = anno,
         cluster_cols = TRUE, cluster_rows = TRUE,
         main = "50 most variable genes - coloured by deviation from gene average",
         fontsize = 6,
         filename=paste0(analysisID, "_pheatmap_50variablegenesv1.pdf"))

#############################
## Sample distance heatmap ##
#############################
## First make a sample distance matrix
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( vsd$pair, vsd$habitat, sep = " - " )
colnames(sampleDistMatrix) <- paste( vsd$pair, vsd$habitat, sep = " - " )
## Now make heatmaps
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main=paste0(analysisID, " sample distance"),
         filename=paste0(analysisID, "_pheatmap_sampDist.pdf"),
         col=colors)

# or a prettier one
colors2 <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(255)
pheatmap(sampleDistMatrix,
         cellwidth = 5, cellheight = 5, fontsize = 6,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         main=paste0(analysisID, " sample distance"),
         filename=paste0(analysisID, "_pheatmap_sampDist_v2.pdf"),
         col=colors2)

# or look at the top 300 genes
  topVarGenes_300 <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 300)
  mat_300  <- assay(vsd)[ topVarGenes_300, ]
  mat_300  <- mat_300 - rowMeans(mat_300)
  pheatmap(mat_300, annotation_col = anno,
           main = "300 most variable genes - coloured by deviation from gene average",
           filename=paste0(analysisID, "_pheatmap_300variablegenes_v1.pdf"))

  # now make the same 300 and 100 graphs but reorder the matrices first so the hypo and epi are grouped
  ## get desired order
nameorder <- rownames(metadata[order(metadata$habitat),])
  mat_reorder  <- assay(vsd)[ topVarGenes, nameorder]

  ## top 100
  # or a prettier one
  pheatmap(mat_reorder, annotation_col = df2, annotation_colors = ann_colors,
           cluster_cols = FALSE, cluster_rows = TRUE,
           main = "50 most variable genes - coloured by deviation from gene average (grouped by habitat)",
           fontsize = 6,
           filename=paste0(analysisID, "_pheatmap_50variablegenesv3.pdf"))
  ## top 1000
  topVarGenes_300_reorder <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 300)
  mat_300_reorder  <- assay(vsd)[ topVarGenes_300_reorder, nameorder]
  mat_300_reorder  <- mat_300_reorder - rowMeans(mat_300_reorder)
  pheatmap(mat_300_reorder, annotation_col = anno,
           main = "300 most variable genes - coloured by deviation from gene average (ordered by habitat)",
           filename=paste0(analysisID, "_pheatmap_300variablegenes_v2.pdf"))

##################################
## Principal Component Analysis ##
##################################
# First define functions (adapting DESeq2's plotPCA but giving more PCs)
# function to perform PCA analysis
fullerPCA <- function(object, intgroup = "condition", ntop = 500, returnData = TRUE) {
  # usage: pcatable <- fullerPCA({vsd}, c("pair", "habitat")
  # NB: commands modified from DESeq2:::plotPCA.DESeqTransform
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup,
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
d <- data.frame(pca$x, group = group, intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
}
# function to get % variance explained by each PC
fullerPCA.pcvar <- function(object, intgroup = "condition", ntop = 500, returnData = TRUE) {
  # usage: pcatable.pcvar <- fullerPCA.pcvar({vsd}, c("pair", "habitat"))
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  return(percentVar)
}
# function to graph PCs of your choosing (PC1-10 possible)
plotMyPCA <- function (PCX, PCY) {
  # usage: someplot <- plotMyPCA("PC1", "PC8")
  require(ggplot2)
  ggplot(pcatable, aes(x = pcatable[[PCX]], y = pcatable[[PCY]], color = pair, shape = habitat)) +
    geom_point(size =3) +
    xlab(paste0("PC", which(colnames(pcatable)==PCX), ": ", round((pcatable.pcvar[which(colnames(pcatable)==PCX)])*100), "% variance")) +
    ylab(paste0("PC", which(colnames(pcatable)==PCY), ": ", round((pcatable.pcvar[which(colnames(pcatable)==PCY)])*100), "% variance")) +
    scale_color_brewer(palette="Paired") +
    ggtitle(paste0("PCA: ", analysisID, " (", paste0("PC", which(colnames(pcatable)==PCX)), " vs ", paste0("PC", which(colnames(pcatable)==PCY)), ")")) +
    theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) +
    #stat_ellipse(type = "t") +
    coord_fixed()
}

# Now run fullerPCA, fullerPCA.pcvar, plot graphs
pcatable <- fullerPCA(vsd, c("pair", "habitat"))
write.table(pcatable, file=paste0(analysisID, "_pcatable.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
pcatable.pcvar <- fullerPCA.pcvar(vsd, c("pair", "habitat"))
write.table(pcatable.pcvar, file=paste0(analysisID, "_pcatable_pcvar.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
myplot <- plotMyPCA("PC1","PC2") # change if you want other PCs
ggsave(paste0(analysisID, "_pca.pdf"))
myplot.ellipse <- myplot + stat_ellipse(type = "t")
ggsave(paste0(analysisID, "_pca_ellipse.pdf"))
# play with colours and labels
myplot.changecolour <- myplot + geom_text(aes(label=rownames(pcatable), color=pair),hjust="inward", vjust="inward", check_overlap = TRUE, size = 3) + geom_point(aes(color=habitat), size = 3)
ggsave(paste0(analysisID, "_pca_v2.pdf"))

# Finally let's see a pairwise comparison of the first five PCs
## Based on Kevin Blighe https://www.biostars.org/p/282685/ (accessed 10.05.19)
nPCs <- ncol(pcatable)-4
if (nPCs > 10) {
  par(cex=1.0, cex.axis=0.8, cex.main=0.8)
  colours <- c('black', 'red')[unclass(pcatable$habitat)]
  pdf(paste0(analysisID, "_PCA1to5.pdf"))
  pairs(pcatable[,1:5], col=colours, main="Principal components analysis bi-plot\nPCs 1-5, red = hypo", pch=16, upper.panel = NULL)
  dev.off()
  pdf(paste0(analysisID, "_PCA6to10.pdf"))
  pairs(pcatable[,6:10], col=colours, main="Principal components analysis bi-plot\nPCs 6-10, red = hypo", pch=16, upper.panel = NULL)
  dev.off()
} else {
  par(cex=1.0, cex.axis=0.8, cex.main=0.8)
  colours <- c('black', 'red')[unclass(pcatable$habitat)]
 # pdf(paste0(analysisID, "_PCA1toN.pdf"))
  pairs(pcatable[,1:nPCs], col=colours, main="Principal components analysis bi-plot\nPCs 1-N, red = hypo", pch=16, upper.panel = NULL)
  dev.off()
}

### let's adapt and make some prettier ones and some that look at sex
if ("A" %in% levels(metadata$sex) || "U" %in% levels(metadata$sex)) {
  plotMyPCA2 <- function (PCX, PCY) {
    # usage: someplot <- plotMyPCA("PC1", "PC8")
    require(ggplot2)
    ggplot(pcatable, aes(x = pcatable[[PCX]], y = pcatable[[PCY]], color = pair, shape = habitat)) +
      geom_point(size =3) +
      xlab(paste0("PC", which(colnames(pcatable)==PCX), ": ", round((pcatable.pcvar[which(colnames(pcatable)==PCX)])*100), "% variance")) +
      ylab(paste0("PC", which(colnames(pcatable)==PCY), ": ", round((pcatable.pcvar[which(colnames(pcatable)==PCY)])*100), "% variance")) +
      #scale_color_manual(values=c("lightblue4", "red", "lightcoral", "yellow3", "palegreen3", "blue","skyblue2","thistle3","turquoise","orange2")) +
      ggtitle(paste0("PCA: ", analysisID, " (", paste0("PC", which(colnames(pcatable)==PCX)), " vs ", paste0("PC", which(colnames(pcatable)==PCY)), ")")) +
      theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) +
      #stat_ellipse(type = "t") +
      coord_fixed()
  }
  myplot2 <- plotMyPCA2("PC1","PC2") # change if you want other PCs
  myplot2
  ggsave(paste0(analysisID, "_pca_v2.pdf"))

  pcatable.sex <- fullerPCA({vsd}, c("sex", "habitat"))
  pcatable.pcvar.sex <- fullerPCA.pcvar({vsd}, c("sex", "habitat"))
}

if ("A" %in% levels(metadata$sex) || "U" %in% levels(metadata$sex)) {
  plotMyPCA3 <- function (PCX, PCY) {
    # usage: someplot <- plotMyPCA("PC1", "PC8")
    require(ggplot2)
    ggplot(pcatable.sex, aes(x = pcatable.sex[[PCX]], y = pcatable.sex[[PCY]], color = sex, shape = habitat)) +
      geom_point(size =3) +
      xlab(paste0("PC", which(colnames(pcatable.sex)==PCX), ": ", round((pcatable.pcvar.sex[which(colnames(pcatable.sex)==PCX)])*100), "% variance")) +
      ylab(paste0("PC", which(colnames(pcatable.sex)==PCY), ": ", round((pcatable.pcvar.sex[which(colnames(pcatable.sex)==PCY)])*100), "% variance")) +
      scale_color_manual(values=c("deeppink1", "dodgerblue1", "gray43")) +
      ggtitle(paste0("PCA: ", analysisID, " (", paste0("PC", which(colnames(pcatable.sex)==PCX)), " vs ", paste0("PC", which(colnames(pcatable.sex)==PCY)), ")")) +
      theme(text = element_text(size=10), axis.text.x = element_text(angle=90, hjust=1)) +
      #stat_ellipse(type = "t") +
      coord_fixed()
  }
  myplot3 <- plotMyPCA3("PC1","PC2") # change if you want other PCs
  myplot3
  ggsave(paste0(analysisID, "_pca_v3.pdf"))
}

### and some more PC1-5 matrices

# colours_hab <- c('black', 'red')[unclass(pcatable$habitat)]
# #colourful_pairs <- c("lightblue4", "red", "lightcoral", "yellow3", "palegreen3", "blue","skyblue2","thistle3","turquoise","orange2")[unclass(pcatable$pair)]
# #boringcolours_pairs <- c("gray", "red", "gray", "gray", "gray", "blue","gray","gray","gray","gray")[unclass(pcatable$pair)]
# pdf(paste0(analysisID, "_PCA1to5_v2.pdf"))
# pairs(pcatable[,1:5], main="Principal components analysis bi-plot\nPCs 1-5,", pch=16, upper.panel = NULL)
# dev.off()
#
# pdf(paste0(analysisID, "_PCA1to5_v3.pdf"))
# pairs(pcatable[,1:5], col=boringcolours_pairs, main="Principal components analysis bi-plot\nPCs 1-5,", pch=16, upper.panel = NULL)
# dev.off()

#par(cex=1.0, cex.axis=0.8, cex.main=0.8)
if ("A" %in% levels(metadata$sex) || "U" %in% levels(metadata$sex)) {
  colours_sex <- c("deeppink1", "dodgerblue1", "gray43")[unclass(pcatable.sex$sex)]
  pdf(paste0(analysisID, "_PCA1to5_sex1-5.pdf"))
  pairs(pcatable.sex[,1:5], col=colours_sex, main="Principal components analysis bi-plot\nPCs 1-5, by sex", pch=16, upper.panel = NULL)
  dev.off()
  pdf(paste0(analysisID, "_PCA1to5_sex6-10.pdf"))
  pairs(pcatable.sex[,6:10], col=colours_sex, main="Principal components analysis bi-plot\nPCs 6-10, by sex", pch=16, upper.panel = NULL)
  dev.off()
}

###############
## MDS Plots ##
###############
mdsData <- cbind(data.frame(cmdscale(sampleDistMatrix)), as.data.frame(colData(vsd)))
mymds <- ggplot(mdsData, aes(X1,X2,color=pair,shape=habitat)) + geom_point(size=3) + ggtitle(paste0("MDS - ", analysisID))
ggsave(paste0(analysisID, "_mds.jpeg"))

# MDS: sex + habitat
if ("A" %in% levels(metadata$sex) || "U" %in% levels(metadata$sex)) {
  mymds_S.H1 <- ggplot(mdsData, aes(X1,X2,color=habitat,shape=sex)) + geom_point(size=3) + ggtitle(paste0("MDS - ", analysisID))
  ggsave(paste0(analysisID, "_mds_sexHabitat1.jpeg"))
}

if ("A" %in% levels(metadata$sex) || "U" %in% levels(metadata$sex)) {
  mymds_S.H2 <- ggplot(mdsData, aes(X1,X2,color=sex,shape=habitat)) + geom_point(size=3) + ggtitle(paste0("MDS - ", analysisID)) + scale_color_manual(values=c("deeppink1", "dodgerblue1", "gray43"))
  ggsave(paste0(analysisID, "_mds_sexHabitat2.jpeg"))
}

# MDS : sex + pair
if ("A" %in% levels(metadata$sex) || "U" %in% levels(metadata$sex)) {
  mymds_S.P1 <- ggplot(mdsData, aes(X1,X2,color=pair,shape=sex)) + geom_point(size=3) + ggtitle(paste0("MDS - ", analysisID))
  ggsave(paste0(analysisID, "_mds_sexPair1.jpeg"))
}

if ("A" %in% levels(metadata$sex) || "U" %in% levels(metadata$sex)) {
  mymds_S.P2 <- ggplot(mdsData, aes(X1,X2,color=pair,shape=sex)) +
    geom_point(size=3, alpha = 0.9) +
    ggtitle(paste0("MDS - ", analysisID))
  #scale_color_manual(values=c("lightblue4", "red", "lightcoral", "yellow3", "palegreen3", "blue","skyblue2","thistle3","turquoise","orange2"))
  ggsave(paste0(analysisID, "_mds_sexPair2.jpeg"))
}
print("@@@ pre-DESeq data viz complete - step 4 of 9")

# ------------------------------------------------------------------
# Run DESeq (not considering batch effects)
# ------------------------------------------------------------------
ddsDE <- DESeq(dds)
res <- results(ddsDE, contrast=c("habitat","hypo","epi"))
res.05 <- results(ddsDE, contrast=c("habitat","hypo","epi"), alpha = 0.05)
res.FC1 <- results(ddsDE, contrast=c("habitat","hypo","epi"), lfcThreshold=1)
res.05.FC1 <- results(ddsDE, contrast=c("habitat","hypo","epi"), alpha = 0.05 , lfcThreshold=1)
res.05.FC.5 <- results(ddsDE, contrast=c("habitat","hypo","epi"), alpha = 0.05 , lfcThreshold=0.5)

# get a summary
sink(paste0(analysisID, "_resultsSumm_noSVA.txt"), append=TRUE, split=TRUE)
print("default parameters")
summary(res)
print("alpha = 0.05")
summary(res.05)
print("LFC = 1")
summary(res.FC1)
print("alpha = 0.05 and LFC = 1")
summary(res.05.FC1)
print("alpha = 0.05 and LFC = 0.5")
summary(res.05.FC.5)
sink()
print("@@@ standard DESeq analysis complete - step 5a of 9")

# ------------------------------------------------------------------
# Export data (not considering batch effects)
# ------------------------------------------------------------------
# write the non-batch-effect results
res.tab <- as.data.frame(subset(res, padj < 0.1))
write.table(res.tab, file=paste0(analysisID, "_results-default.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
res.tab.05 <- as.data.frame(subset(res, padj < 0.05))
write.table(res.tab.05, file=paste0(analysisID, "_results-alpha.05.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
res.tab.FC1 <- as.data.frame(subset(res.FC1, padj < 0.1))
write.table(res.tab.FC1, file=paste0(analysisID, "_results-LFC1.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
res.tab.05.FC1 <- as.data.frame(subset(res.05.FC1, padj < 0.05))
write.table(res.tab.05.FC1, file=paste0(analysisID, "_results-LFC1-alpha05.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
res.tab.05.FC.5 <- as.data.frame(subset(res.05.FC.5, padj < 0.05))
write.table(res.tab.05.FC.5, file=paste0(analysisID, "_results-LFC0.5-alpha05.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)

print("@@@ standard DESeq analysis results written to file - step 5b of 9")

# ------------------------------------------------------------------
# POST-DESEQ DATA VISUALISATION AND EXPLORATION (not considering batch effects)
# ------------------------------------------------------------------
# first we run lfcShrink to shrink log2 FC  using apeglm
resShr <- lfcShrink(ddsDE, coef="habitat_hypo_vs_epi", type="apeglm")

##############
## MA Plots ##
##############
# The plot visualizes the differences between measurements taken in two
# samples, by transforming the data onto M (log ratio) and A (mean average)
# scales, then plotting these values.
#
# The log2 fold change for a particular comparison is plotted on the y-axis and
# the average of the counts normalized by size factor is shown on the x-axis.
# Each gene is represented with a dot. Genes with an adjusted p value below a
# threshold (here 0.1, the default) are shown in red.
jpeg(paste0(analysisID, "_MA.jpg"))
plotMA(resShr, ylim = c(-5, 5), main = paste0(analysisID, ": MA"))
dev.off()

###########################
## Histogram of p-values ##
###########################
# This plot is best formed by excluding genes with very small counts, which
# otherwise generate spikes in the histogram.
jpeg(paste0(analysisID, "_pvalhisto.jpg"))
hist(resShr$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
dev.off()

print("@@@ post-DESeq data viz complete - step 6 of 9")

# # ------------------------------------------------------------------
# # Incorporate batch effects
# # ------------------------------------------------------------------
# ddsDEcounts <- counts(ddsDE, normalized=TRUE)
# ddsDEgetrows <- rowMeans(ddsDEcounts) > 1
# ddsDEcounts <- ddsDEcounts[ddsDEgetrows,]
# modelH1 <- model.matrix(design = formula(paste("~",myformula)), colData(ddsDE))
# modelH0 <- model.matrix(~ 1, colData(ddsDE))
#
# # calculate number of latent variables
# n.sv = num.sv(ddsDEcounts, modelH1, method="be")
# if (n.sv == 0) {
#   cat("WARNING: 0 variables predicted!\n")
#   n.sv <- 5
#   cat("WARNING: n.sv manually set to 5 variables\n")
#   #quit()
# }
# print(paste0("@@@ analysing ", n.sv, "  variables"))
# # perform SVA estimation
# correctionSVA <- svaseq(ddsDEcounts, modelH1, modelH0, n.sv=n.sv)
# print("@@@ batch effects estimated - step 7a of 9")
#
# ########################
# ## Plot Batch Effects ##
# ########################
# # define a function to plot SVA variables on known metadata
# plotSVABatch <- function (condition) {
#   # usage1: jpeg(paste0(analysisID, "_SVAgroupings_morphospecies.jpg"))
#   # usage2: plotBatch("condition")
#   # usage3: dev.off()
#   par(mfrow=c(5,1),mar=c(3,5,3,1))
#   stripchart(correctionSVA$sv[,1] ~ dds[[condition]],vertical=TRUE,main=paste0("SV1 - ", analysisID, " by ", condition))
#   abline(h=0)
#   stripchart(correctionSVA$sv[,2] ~ dds[[condition]],vertical=TRUE,main=paste0("SV2 - ", analysisID, " by ", condition))
#   abline(h=0)
#   stripchart(correctionSVA$sv[,3] ~ dds[[condition]],vertical=TRUE,main=paste0("SV3 - ", analysisID, " by ", condition))
#   abline(h=0)
#   stripchart(correctionSVA$sv[,4] ~ dds[[condition]],vertical=TRUE,main=paste0("SV4 - ", analysisID, " by ", condition))
#   abline(h=0)
#   stripchart(correctionSVA$sv[,5] ~ dds[[condition]],vertical=TRUE,main=paste0("SV5 - ", analysisID, " by ", condition))
#   abline(h=0)
# }
#
# # plot all metadata variables
# for (condition in colnames(colData(dds))) {
#   jpeg(paste0(analysisID, "_SVAgroupings_", condition, ".jpg"))
#   plotSVABatch(condition)
#   dev.off()
# }
# print("@@@ batch effects plotted - step 7b of 9")
#

# # ------------------------------------------------------------------
# # Run DESeq (incorporating batch effects)
# # ------------------------------------------------------------------
#
# ddsSVA <- ddsDE
# ddsSVA_partial <- ddsDE
# ddsSVA$SV1 <- correctionSVA$sv[,1]
# ddsSVA$SV2 <- correctionSVA$sv[,2]
# ddsSVA$SV3 <- correctionSVA$sv[,3]
# ddsSVA$SV4 <- correctionSVA$sv[,4]
# ddsSVA$SV5 <- correctionSVA$sv[,5]
#
# ddsSVA_partial$SV1 <- correctionSVA$sv[,1]
# ddsSVA_partial$SV2 <- correctionSVA$sv[,2]
#
# # we will make ddsSVA (incorporating all 5 predicted sources of batch effects) and ddsSVA_partial (incorporating the first two)
# ## first ddsSVA
# design(ddsSVA) <- ~ pair + SV1 + SV2 + SV3 + SV4 + SV5 + habitat
# ddsSVA <- DESeq(ddsSVA)
# resSVA <- results(ddsSVA, contrast=c("habitat","hypo","epi"))
# resSVA.05 <- results(ddsSVA, contrast=c("habitat","hypo","epi"), alpha = 0.05)
# resSVA.FC1 <- results(ddsSVA, contrast=c("habitat","hypo","epi"), lfcThreshold=1)
# resSVA.05.FC1 <- results(ddsSVA, contrast=c("habitat","hypo","epi"), alpha = 0.05 , lfcThreshold=1)
#
# ## get a summary
# sink(paste0(analysisID, "_resultsSumm_SVAfull.txt"), append=TRUE, split=TRUE)
# print("default parameters")
# summary(resSVA)
# print("alpha = 0.05")
# summary(resSVA.05)
# print("LFC = 1")
# summary(resSVA.FC1)
# print("alpha = 0.05 and LFC = 1")
# summary(resSVA.05.FC1)
# sink()
#
# ## and now ddsSVA_partial
# design(ddsSVA_partial) <- ~ pair + SV1 + SV2 + habitat
# ddsSVA_partial <- DESeq(ddsSVA_partial)
# resSVA_partial <- results(ddsSVA_partial, contrast=c("habitat","hypo","epi"))
# resSVA_partial.05 <- results(ddsSVA_partial, contrast=c("habitat","hypo","epi"), alpha = 0.05)
# resSVA_partial.FC1 <- results(ddsSVA_partial, contrast=c("habitat","hypo","epi"), lfcThreshold=1)
# resSVA_partial.05.FC1 <- results(ddsSVA_partial, contrast=c("habitat","hypo","epi"), alpha = 0.05 , lfcThreshold=1)
#
# ## get a summary
# sink(paste0(analysisID, "_resultsSumm_SVApartial.txt"), append=TRUE, split=TRUE)
# print("default parameters")
# summary(resSVA_partial)
# print("alpha = 0.05")
# summary(resSVA_partial.05)
# print("LFC = 1")
# summary(resSVA_partial.FC1)
# print("alpha = 0.05 and LFC = 1")
# summary(resSVA_partial.05.FC1)
# sink()
#
# print("@@@ DESeq2 run incorporating batch effects - results run - step 8 of 9")

# # ------------------------------------------------------------------
# # Export data (incorporating batch effects)
# # ------------------------------------------------------------------
#
# # full model effects
# resSVA.tab <- as.data.frame(subset(resSVA, padj < 0.1))
# write.table(resSVA.tab, file=paste0(analysisID, "_SVA_results-default.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
# resSVA.tab.05 <- as.data.frame(subset(resSVA.05, padj < 0.05))
# write.table(resSVA.tab.05, file=paste0(analysisID, "_SVA_results-alpha.05.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
# resSVA.tab.FC1 <- as.data.frame(subset(resSVA.FC1, padj < 0.1))
# write.table(resSVA.tab.FC1, file=paste0(analysisID, "_SVA_results-LFC1.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
# resSVA.tab.05.FC1 <- as.data.frame(subset(resSVA.05.FC1, padj < 0.05))
# write.table(resSVA.tab.05.FC1, file=paste0(analysisID, "_SVA_results-LFC1-alpha.05.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
#
# #partial (SV1+2) effcts
# resSVA_partial.tab <- as.data.frame(subset(resSVA_partial, padj < 0.1))
# write.table(resSVA_partial.tab, file=paste0(analysisID, "_SVApartial_results-default.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
# resSVA_partial.tab.05 <- as.data.frame(subset(resSVA_partial.05, padj < 0.05))
# write.table(resSVA_partial.tab.05, file=paste0(analysisID, "_SVApartial_results-alpha.05.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
# resSVA_partial.tab.FC1 <- as.data.frame(subset(resSVA_partial.FC1, padj < 0.1))
# write.table(resSVA_partial.tab.FC1, file=paste0(analysisID, "_SVApartial_results-LFC1.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
# resSVA_partial.tab.05.FC1 <- as.data.frame(subset(resSVA_partial.05.FC1, padj < 0.05))
# write.table(resSVA_partial.tab.05.FC1, file=paste0(analysisID, "_SVApartial_results-LFC1-alpha.05.txt"), sep="\t", quote = FALSE, col.names = NA, row.names = TRUE)
#
# print("@@@ exported SVA DESeq2 results - step 8 of 9")

# # ------------------------------------------------------------------
# # POST-DESEQ DATA VISUALISATION AND EXPLORATION (incorporating batch effects)
# # ------------------------------------------------------------------
# # first we run lfcShrink to shrink log2 FC  using apeglm
# resShrSVA <- lfcShrink(ddsSVA, coef="habitat_hypo_vs_epi", type="apeglm")
# resShrSVA_partial <- lfcShrink(ddsSVA_partial, coef="habitat_hypo_vs_epi", type="apeglm")
#
# ##############
# ## MA Plots ##
# ##############
# # SVA full
# jpeg(paste0(analysisID, "_SVA_MA.jpg"))
# plotMA(resShrSVA, ylim = c(-5, 5), main = paste0(analysisID, ": MA (SVA full)"))
# dev.off()
# # SVA partial
# jpeg(paste0(analysisID, "_SVApartial_MA.jpg"))
# plotMA(resShrSVA_partial, ylim = c(-5, 5), main = paste0(analysisID, ": MA (SVA partial)"))
# dev.off()
#
# ###########################
# ## Histogram of p-values ##
# ###########################
# jpeg(paste0(analysisID, "_SVA_pvalhisto.jpg"))
# hist(resShrSVA$pvalue[resShrSVA$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white", main = paste0(analysisID, ": pvalues (SVA full)"))
# dev.off()
# jpeg(paste0(analysisID, "_SVApartial_pvalhisto.jpg"))
# hist(resShrSVA_partial$pvalue[resShrSVA_partial$baseMean > 1], breaks = 0:20/20, col = "grey50", border = "white", main = paste0(analysisID, ": pvalues (SVA full)"))
# dev.off()
#
# print("@@@ post-DESeq2 (batch effect version) data viz complete - step 9 of 9")

# ------------------------------------------------------------------
# WRAPPING UP
# ------------------------------------------------------------------

# save useful R objects to file
save(dds, file = paste0(analysisID, "_dds.Rdata"))

# Write session info to file
writeLines(capture.output(sessionInfo()), paste0(analysisID, "_sessionInfo.txt"))

# Acknowledgements
## NB we already have the loaded package citations from the start
sink("citations.txt", append=TRUE)
print("#base R package")
citation()

print(paste0("@@@ runDESeq.R analysis complete at ", date()))
