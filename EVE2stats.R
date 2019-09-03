# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 18 July 2019
# Title: EVE2stats.R
# Goal:
# Usage: Rscript EVE2stats.R workdirPath lrtPath geneIDPath outputName
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
	cat("usage: Rscript EVE2stats.R workdirPath lrtPath geneIDPath outputName\n")
	quit()
}

workdirPath <- args[1]   #where the output will be saved
lrtPath <- args[2]       #unmodified LRT output from EVE, single line, space delimited
geneIDPath <- args[3]    #gene IDs, same order as LRT file, one per row
graphName <- "LRT values of in/significant genes"     #desired name for LRT significance graph
outputName <- args[4]    #desired name for output files


# ------------------------------------------------------------------
# PREAMBLE
# ------------------------------------------------------------------

#title: "Custom R script to process results of EVE comparative transcriptomics software"
#author: Laura Grice, after "Brauer, C.J., Unmack, P.J., and Beheregaray, L.B. (2017). Comparative ecological transcriptomics and the contribution of gene expression to the evolutionary potential of a threatened fish. Mol Ecol 26, 6841???6856."
#date: "2017 and 30.01.19"
#output:
#  html_notebook:
#    theme: cerulean
#    toc: yes


# Introduction

## Aim
#This script was originally included as Appendix S2 from : "Brauer, C.J., Unmack, P.J., and Beheregaray, L.B. (2017). Comparative ecological transcriptomics and the contribution of gene expression to the evolutionary potential of a threatened fish. Mol Ecol 26, 6841???6856"

#The authors first run EVE on expression and phylogenetic data:

#**"...A recent extension of such ANOVA-based methods is the Expression Variance and Evolution Model (EVE) (Rohlfs & Nielsen, 2015). EVE models gene expression as a quantitative trait across a phylogeny; considering the ratio (b) of among-lineage expression divergence to within-lineage expression diversity ... The expectation is that b should be consistent among the majority of genes that have undergone similar evolutionary and demographic processes, but higher for genes with more variance within than among lineages, and lower for genes with more variance among compared to within lineages. Here, the EVE model was used to parameterize the ratio b across the five populations of N. australis, to identify transcripts potentially under divergent selection for expression level (low b) or transcripts with high expression plasticity (high b). The model utilizes gene expression data and a phylogeny as input."**

#Note that EVE is currently installed on Asellus (`~/scripts/software/EVE/EVEmodel`) and can be run according to the README file in the same directory.

#Brauer et al. then take the LRT test output and process it using a slightly different version of the commands contained in this R notebook:

#**"For each transcript (i), maximum-likelihood values were calculated and a likelihood ratio test (LRT) was performed to assess the null hypothesis that b-i is equal to b for all transcripts. Under the null model, the LRT statistic follows a v21 distribution (Rohlfs & Nielsen, 2015) and a custom R script (Appendix S2, Supporting information) was used to identify candidate transcripts where the LRT statistic deviated from this distribution at a FDR of 10%."**

## Input
#This R notebook takes the EVE LRT test results, calculates FDR p-values, pulls out #candidate genes for statistically significant divergence (between species) or diversity #(within species), and plots a graph of the divergence results.

#To run this script you will need to enter the following information when you run the command

#*lrtPath* - The path to your EVE LRT output. Do not modify this file - i.e. it should be a space-delimited list on a single line. See the EVE README file for how to create this file.
#*geneIDPath* - make a file listing the ID of every gene in the LRT file, in the same order, with one gene per row. Give the path to this file.
#*graphName* - This title will be shown on a graph of the in/significant LRT divergence results
#*outputName* - This prefix will be given to all output files
#*workdirPath* - Your output files will be saved here

## Output
#2x divergence files:
#  1. *_EVE_diverge_Pvalues.txt - a table of genes of interest, LRT value, divergence p-value and FDR value
#  2. *_EVE_diverge_FDRgenes.txt - the same table as above, containing only genes with FDR < threshold.
#2x diversity files:
#  3. *_EVE_diverse_Pvalues.txt - as above, but for diversity p-value/FDR
#  4. *_EVE_diverse_FDRgenes.txt - diversity values for genes below FDR threshold
#1x graph:
#  5. <output>_EVE_diverge_result.pdf - a PDF showing the LRT values of all genes, coloured red if they show significant evidence of divergence between species

#If no genes are significant, the *FDRgenes.txt file/s will be empty.

#The function EVEresults takes teh LRT values and generates p-values and FDR statistics. The first steps of the function do the statistical heavy-lifting. P-values are calculated with a chi-squared test with 1 degree of freedom, because "the alternate hypothesis has one additional degree of freedom as compared to the null hypothesis, [so] the asymtotic distribution for the LR test statistic under the null hypothesis is ????2 with one degree of frededom (LRT??i?????shared ~ ????12)")"

#define EVEresults function
EVEresults <- function(genes, LRT, df, title, rate, output){
  # LG: get p-values and FDR
  P_diverge <- 1 - pchisq(LRT, df = df)
  P_diverse <- pchisq(LRT, df = df)
  FDR_diverge <- p.adjust(P_diverge, "fdr")
  FDR_diverse <- p.adjust(P_diverse, "fdr")

  # LG: Write input data, p-value and FDR to table
  EVE_diverge <- as.data.frame(cbind(genes, LRT, P_diverge, FDR_diverge))
  EVE_diverse <- as.data.frame(cbind(genes, LRT, P_diverse, FDR_diverse))
  colnames(EVE_diverge) <- c("Gene", "LRT", "P", "FDR")
  colnames(EVE_diverse) <- c("Gene", "LRT", "P", "FDR")

  # LG: Filter to keep only genes where FDR < rate (e.g. 0.1)
  EVE_diverge.sub <- subset(EVE_diverge, EVE_diverge$FDR < rate)
  EVE_diverge.sub <- as.data.frame(EVE_diverge.sub)
  colnames(EVE_diverge.sub)<- c("Gene", "LRT", "P", "FDR")
  EVE_diverse.sub <- subset(EVE_diverse, EVE_diverse$FDR < rate)
  EVE_diverse.sub <- as.data.frame(EVE_diverse.sub)
  colnames(EVE_diverse.sub)<- c("Gene", "LRT", "P", "FDR")

  # LG: Print some stats on the screen
  FDR_diverge.len <- length(EVE_diverge.sub[,1])
  EVE_diverge.len <- length(EVE_diverge[,1])
  FDR_diverse.len <- length(EVE_diverse.sub[,1])
  EVE_diverse.len <- length(EVE_diverse[,1])
  print(paste0(FDR_diverge.len, " genes have higher variance among than within lineages at ", rate, "FDR."))
  print(paste0(FDR_diverse.len, " genes have higher variance within than among lineages at ", rate, "FDR."))
  diverge_result <- list(EVE_diverge, EVE_diverge.sub)
  diverse_result <- list(EVE_diverse, EVE_diverse.sub)
  wd <- getwd()
  print(paste0("Writing results to ", wd))

  # LG: Write results to table
  # EDIT LG: Previously an if-else statement where results written to .xslx files unless 5000+ genes, in which case csv files written. I'd prefer to have tab-delim files for everything.
  write.table(diverge_result[1], file=paste0(output,"_EVE_diverge_Pvalues.txt"), sep="\t")
  write.table(diverge_result[2], file=paste0(output,"_EVE_diverge_FDRgenes.txt"), sep="\t")
  write.table(diverse_result[1], file=paste0(output,"_EVE_diverse_result_Pvalues.txt"), sep="\t")
  write.table(diverse_result[2], file=paste0(output,"_EVE_diverse_result_FDRgenes.txt"), sep="\t")

  # LG: To make a graph, work out the axis scale in multiples of 500, based on number of genes
    # LG: (10 genes = scale of 500, 501 genes = scale of 1000)
    # LG: ceiling(index.lim/500)*500 = [(# genes / 500) --> round to nearest number] * 500
    # LG: Axis on multiples of 500 is annoying for low gene #, I changed to 50
    # LG: If you change this number back to 500, change the "axis" line below so that 2x "by=25" are now "by=1000"
  index.lim <- length(EVE_diverge[,1])
  index.lim <- ceiling(index.lim/50)*50

  # LG: Plot LRT values, red for significant, black for not significant
    # LG: Note that xact = "n" suppresses the axis labels for the X-axis - we set this in "axis" next
    # LG: axis sets the X-axis, this originally said "by=1000"
  fig <- plot.default(EVE_diverge$LRT, col=ifelse(EVE_diverge$FDR < rate, "red", "black"), main = title, xlim = c(0, index.lim), xaxt = "n")
  axis(side=1,at=pretty(seq(0, index.lim, by=25)),labels=pretty(seq(0, index.lim, by=25)))
  pdf(file = paste0(output,"_EVE_diverge_result.pdf"), height = 8.27, width = 11.69)
  # LG: Same thing, to view in the R window
    plot.default(EVE_diverge$LRT, col=ifelse(EVE_diverge$FDR < rate, "red", "black"), main = title, xlim = c(0, index.lim), xaxt = "n")
  axis(side=1,at=pretty(seq(0, index.lim, by=25)),labels=pretty(seq(0, index.lim, by=25)))
  dev.off()
  print("Finished!")
  return(diverge_result)
  return(diverse_result)
  return(fig)
}


#The following commands assume that your LRT file is exactly as produced by EVE - that is, #space-delimited on a single line. If you have modified the file so it is one value per #row, you will need to change the lrtValues assignment as follows:
#`lrtValues<-read.table("betaTestLRTs_<suffix>.res")`
#`lrtValues<-as.matrix(lrttab)`

#Your geneID file, however, should be one gene per row, in the same order as the LRT file.

# load data files}
# NB: this bit added by LG
lrtValues<-t(read.table(lrtPath))
geneID<-(read.table(geneIDPath))
dfValue <- 1 # degrees of freedom, as specified in Brauer et al. script
FDRrate <- 0.1 # FDR threshold


#run EVEresults function}
MyEVEresults <- EVEresults(genes=geneID, LRT=lrtValues, df=dfValue, title=graphName, rate=FDRrate, output=outputName)


#session info
sessionInfo()
