# R SCRIPT
# ------------------------------------------------------------------
# Author: Laura Grice
# Date: 17 April 2019
# Title: txSalmon2geneSalmon.R
# Goal: To collapse transcriptome-level Salmon output to gene family-level values, ready for txImport + DESeq2
# ------------------------------------------------------------------
# USAGE
# Usage: Rscript txSalmon2geneSalmon.R {QuantPath} {metadata} {fnodes} {outdir path}
# NB: You should run ALL samples at once (they are all required to calculate EffectiveLengths)
## QuantPath: /path/to/directory (that holds quant.sf)
## Metadata: The easiest thing is to use the metadata file you will use for DESeq2, but in its simplest form it can be a single-row file listing $(basename --suffix=.quant.sf {my.quant.sf})
## fnodes: tab-delim file, no header, col1 = family ID, col2 = txID (all samples combined)
## outdir path: /path/to/output/dir
# ------------------------------------------------------------------
# ACKNOWLEDGEMENTS
# This script adapts the "transcriptome to gene" summarisation methods of txImport
# Paper: Soneson, Love, Robinson (2015). F1000Research. http://dx.doi.org/10.12688/f1000research.7563.1
# R scripts on Github: tximport.R, summarizeToGene.R and helper.R on https://github.com/mikelove/tximport/tree/master/R
# ------------------------------------------------------------------
# PITFALLS AND TO DO LIST
# 1. if you end the quantpath with a "/" it won't run properly
# 2. if you end the outdir path with a "/" it will save the output to one dir higher than you want
# 3. if you're going to run tximport and deseq anyway, it'd be more efficient to require a tx2gene file than an fnodes file
# 4. the merge fnodes + EffLenTx step is very slow
# 5. remove the random print status updates
# ------------------------------------------------------------------

## ----load arguments------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  cat("usage: Rscript txSalmon2geneSalmon.R arg1 arg2 arg3 arg4\n")
  quit()
}

QuantPath <- args[1]
metadataPath <- args[2]
fnodesPath <- args[3]
outDir <- args[4]

print("arguments loaded - step 1 of 9")

## ----load user input-----------------------------------------------------
outDir <- if(endsWith(outDir, "/") == FALSE){
  paste(outDir, "/", sep="")
} #make sure the outdir ends in a / and fix if required
  metadata <- read.table(metadataPath, sep = "\t", row.names=1, header = TRUE)
  sampleID <- rownames(metadata)
  quant <- Sys.glob(file.path(QuantPath,(paste0(sampleID,"*quant.sf")))) # Sys.glob() required for wildcard expansion
  fnodes <- read.table(fnodesPath, header = FALSE, row.names = 2)
  colnames(fnodes) <- c("Family")

print("user input loaded - step 2 of 9")

## ----load some functions-------------------------------------------------
#### ----function to add calculated values into an empty matrix------------
# usage: CoolMatrix <- addToMatrix(CoolMatrix, TPM, sampleID[i])
addToMatrix <- function(someMatrix, someInput, someDataname) {
  colnames(someInput) <- someDataname
  someMatrix <- merge(someMatrix, someInput, by = 0, all.x = TRUE)
  rownames(someMatrix) = someMatrix$Row.names
  someMatrix$Row.names <- NULL
  return(someMatrix)
}
#### ----function to add values into final sample-specific "quant.sf" files--
# Usage: myOutput <- createFakeQuantSF("ColName in e.g. TPMMat")
# Usage: myOutput <- createFakeQuantSF(sampleID[i])
createFakeQuantSF <- function(sampleOfInterest) {
  # make empty matrix for data
  quantSF_template <- matrix(nrow = length(usedFams), ncol = 5)
  colnames(quantSF_template) <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads")
  quantSF_template[,1] <- usedFams
  # slot in the data
  quantSF_template[,"Length"] <- LengthMat[,sampleOfInterest]
  quantSF_template[,"EffectiveLength"] <- EffLenMat[,sampleOfInterest]
  quantSF_template[,"TPM"] <- TPMMat[,sampleOfInterest]
  quantSF_template[,"NumReads"] <- NumReadsMat[,sampleOfInterest]
  quantSF_template <- as.data.frame(quantSF_template)
  return(quantSF_template)
}

print("functions loaded - step 3 of 9")

## ----get list of gene families in quant files----------------------------
# Here we find the full list of gene families present across our datasets
# This is important if we are using 0:N orthologues (i.e. families present in SpeciesX but not SpeciesY)
library(data.table)
usedIDs <- list(c())
for (i in seq_along(quant)) {
  names <- fread(quant[i], select = c(1:1), header = TRUE)
  usedIDs <- rbind(usedIDs, names)
} # loop gets "Names" for each quant.sf file, adds to a list
usedIDs <- usedIDs[!duplicated(usedIDs),]
usedIDs <- as.data.frame(usedIDs$Name)
colnames(usedIDs) <- c("transcriptID")
usedIDs <- merge(fnodes, usedIDs, by.x = 0, by.y = "transcriptID", all.x = FALSE, all.y = TRUE)
usedFams <- levels(droplevels(usedIDs$Family)) #droplevels needed or you get weird "all NA" samples

print("gene families done - step 4 of 9")

## ----make empty matrices-------------------------------------------------
# Make empty matrices ready to be filled with the results
## We need to make matrices with NO columns in case we are using 0:N orthologues
## This method allows each quant.sf table to have a different length.
## Unfortunately, this tactic (using merge) is slower than that used in txImport.R (lines 381-383)
matTemplate <- matrix(ncol = 0, nrow = length(usedFams))
rownames(matTemplate) <- usedFams

# Each of these matrices is named based on the quant.sf nomenclature, not the txImport nomenclature (where "length" = "effectiveLength")
LengthMat <- matTemplate
EffLenMat <- matTemplate
TPMMat <- matTemplate
NumReadsMat <- matTemplate

print("empty matrices made - step 5 of 9")

## ----analyse the data----------------------------------------------------
for (i in seq_along(quant)) {
  # read in the data and add family IDs
  tab <- read.table(quant[i], header = TRUE, row.names = 1)
  tab <- merge(tab, fnodes, by = 0, all = FALSE)
  rownames(tab)=tab$Row.names
  tab$Row.names <- NULL
  # Now calculate those values which can be made independently for each library
  TPM <- rowsum(tab$TPM, tab$Family)
  NumReads <- rowsum(tab$NumReads, tab$Family)
  Length <- rowsum(tab$Length, tab$Family)
  weightedLength <- rowsum(tab$TPM * tab$EffectiveLength, tab$Family)
  simpleEffLen <- weightedLength / TPM
  # Now we slot the data into the matrices, using the function `addToMatrix` above:
  ## NB: For TPMMat, NumReadsMat and LengthMat, we could theoretically slot directly into the blank quant.sf templates
  ## But we need some of these values to calculate the pooled EffectiveLength values
  ## So we follow the way it works in txImport
  TPMMat <- addToMatrix(TPMMat, TPM, sampleID[i])
  NumReadsMat <- addToMatrix(NumReadsMat, NumReads, sampleID[i])
  LengthMat <- addToMatrix(LengthMat, Length, sampleID[i])
  EffLenMat <- addToMatrix(EffLenMat, simpleEffLen, sampleID[i])
}

print("data analysed - step 6 of 9")

## ----calculate genewise average (for families with no TPM data)----------
# Pull out Salmon's EffectiveLength values for all samples
EfLenTx <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Name", "EffectiveLength"))
print("step 7a of 7d - make empty matrix")

for (i in seq_along(quant)) {
  eflens <- fread(quant[i], select = c(1:1, 3:3), header = TRUE)
  EfLenTx <- rbind(EfLenTx, eflens)
} # loop gets "Names" and "EffLen" for each quant.sf file, adds to a list
print("step 7b of 7d - rbind")
EfLenTxGn <- merge(fnodes, EfLenTx, by.x = 0, by.y = "Name", all.x = FALSE, all.y = TRUE)
print("step 7c of 7d - merge")
# Now calculate family-wise means (across all samples)
EfLen_GenewiseMeans <- aggregate(EfLenTxGn[,3], list(EfLenTxGn$Family), mean)
rownames(EfLen_GenewiseMeans) = EfLen_GenewiseMeans$Group.1
EfLen_GenewiseMeans$Group.1 <- NULL
colnames(EfLen_GenewiseMeans) <- c("AveMean")
print("step 7d of 7d - genewise means")

print("genewise ave for effective lengths calculated - step 7 of 9")

## ----calculate alternative EffLen values---------------------------------
# now we look at each row with NaNs
# if ALL values for a given gene are NaN, we slot in the EfLen_GenewiseMeans value we just calculated
# if only SOME values are NaN, we can calculate their EffLens from the OTHER samples
nanRows <- which(apply(EffLenMat, 1, function(row) any(is.nan(row)))) #doesn't find NA only NaN
for (i in nanRows) {
  if (all(is.na(EffLenMat[i,]))) {
    EffLenMat[i,] <- EfLen_GenewiseMeans[i,]
  } else {
    fix <- which(is.na(EffLenMat[i,])) #which rows we want to fix
    EffLenMat[i,fix] <- as.numeric(exp(rowMeans(log(EffLenMat[i,]), na.rm = TRUE)))
   }
}

print("EffLens fixed - step 8 of 9")

## ----slot data into fake quantSF files-----------------------------------
# Make a dir to save the data
mainDir <- "."
outDir_path <- file.path(mainDir, outDir)
if (!dir.exists(outDir_path)){
	print("creating")
	dir.create(outDir_path)
} else {
    print("outDir_path already exists!")
}

for (i in seq_along(sampleID)){
  MyFakeQuantSF <- createFakeQuantSF(sampleID[i])
  write.table(MyFakeQuantSF, file = paste0(outDir, "/", sampleID[i], "_FAM.quant.sf"), sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE )
}

print("data slotted in - FINISHED 9 of 9")

## ----session info--------------------------------------------------------
sessionInfo()
