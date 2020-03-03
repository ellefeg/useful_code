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

bioc_packages <- c("multtest")
## function to load bioconductor packages
biocondpkgTest <- function(x) {
  if (!suppressMessages(require(x,character.only = TRUE, quietly = T))) {
    source("http://www.bioconductor.org/biocLite.R")
    biocLite(x)
    if(!require(x,character.only = TRUE, quietly = T)) stop (paste0(x, "bioconductor package not found"))
  }
}
## load packages
for (b_pkg in bioc_packages) {
  biocondpkgTest(b_pkg)
  sink("citations.txt", append=TRUE)
  print(citation(b_pkg))
  sink()
}
