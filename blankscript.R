# R SCRIPT
# ------------------------------------------------------------------
# Author: 
# Date: 
# Title: 
# Goal: 
# Usage: Rscript enrichment.R {analysisID} {genesOfInterest} {backgroundGenes} {annotations} {descriptions}
# ------------------------------------------------------------------
# USAGE
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# INITIALISATION
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 5) {
  cat("ERROR: 5 arguments expected\n")
  cat("example: Rscript pfam_enrichment.R {analysisID} {genesOfInterest} {backgroundGenes} {annotations} {descriptions}\n")
  quit()
}

analysisID <- args[1]
myInterestGenes <- args[2]
myBackgroundGenes <- args[3]
myAnnotations <- args[4]
myDescriptions <- args[5]
# load packages
bioc_packages <- c("DESeq2", "pheatmap", "apeglm", "genefilter")
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
}
for (r_pkg in r_packages) {
  baseRpkgTest(r_pkg)
}

# checks if your outdir ends in / and adds one if not
if (endsWith(outdir, "/") == FALSE) {
  outdir <- paste0(outdir, "/", sep="")
}





sessionInfo()
