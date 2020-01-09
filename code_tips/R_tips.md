# Working with real files

Getting files into R
```r
# tell R where a file is
myfilepath <- "/path/to/file"
data <- read.delim(myfilepath, header = T, comment.char = "", row.names = 1, stringsAsFactors = F, header = TRUE)
# convert missing data (NAs) to 0
data <- data[apply(data, 1, function(x) sum(is.na(x)))  == 0,]
```

Getting files out of R:
```r
write.table(data, file = "/path/to/file.tab", sep = "\t", quote = F) 
```

# Look at an R object

```r
# what sort of data is it?
typeof(data)
# how many rows?
ncol(data)
# do we have 7 rows?
ncol(data) == 7
# look at the top 3 rows (plus colnames)
data[1:3,]
# look at the first 3 columns (plus rownames)
data[,1:3]
# look at a specific column
data$firstcolname
# look at a specific value in a specific column
data$firstcolname[2]
# compare this value to an expectation
data$firstcolname[2] > 1
# how many rows are less than 0.05?
sum(df$column < 0.05)
```

# Manipulate a dataset
```r
# pull out a smaller part of a table
datasmall <- data[1:3,1:3]
# convert to a matrix
datasmallM <- as.matrix(datasmall)
# calculate row-wise mean (requires matrix)
apply(datasmallM, 1, function(x) mean(x, na.rm=T))
## and convert back to a matrix
as.matrix(apply(datasmallM, 1, function(x) mean(x, na.rm=T)))
# rename column data$stupidname to data$goodname
colnames(data)[colnames(data)=="stupidname"] <- "goodname"
# rename all columns
colnames(avecounts_AAD3) <- c("name", "self", "nonself")
# remove a column
data$badColumn <- NULL
# calculate a log
log(datasmall)
# merge two dataframes
merge(data1, data2, by = "sharedColID")
## merge by rownames
merge(data1, data2, by=0)
## NB this will "break" the row.names - which will be moved to their own column called Row.names, to fix:
rownames(de)=de$Row.names
```

A useful blog post about how you can (and can't) modify R objects by running functions on them: https://stackoverflow.com/questions/15497947/call-by-reference-in-r-using-function-to-modify-an-object

Filter a dataframe to only keep rows that equal a given value in a given column
```{r}
table[table$col==val,]
# or if col is a variable
table[table[[col]]==val,]
```

## Adding new data to existing data frames

The `mutate` function (`library(dplyr)`) can calculate new column data from existing columns. See [this useful tutorial](https://rstudio-pubs-static.s3.amazonaws.com/116317_e6922e81e72e4e3f83995485ce686c14.html#/5).

```{r}
# example: make a new column called "FavSpecies", where favourite species are annotated with T (and others with F)
## This might be useful for differential colouring in graphs
### test the if/else command:
### ifelse(data$species == "FavSpeciesA" | data$speices == "FavSpeciesB", T, F)
data <- mutate(data, FavSpecies = ifelse(data$species == "FavSpeciesA" | data$speices == "FavSpeciesB", T, F))

# example: make new columns data$LogLength and data$LogHeight from data$Length and data$Height
data2 <- mutate(data, LogLength = log(Length))
## don't forget to use data2 as input here!
data2 <- mutate(data2, LogHeight = log(Height))
```

Sum all values in Column X that are united by a shared value in Column Y (e.g. sum TPMs of all isoforms of GeneFamilyY)
```{r}
rowsum(data$TPM, FAMID)
```

## Loops

The most basic forloop. 
```{r}
for (i in 1:10) {
paste(i)
}
```

Loop through samples, perform an analysis, and save the results to a vector. Then paste together with some other data:
```{r}
for i i
```

Concatenate two tables together vertically. Example: Two tables, each with the columns "Year" and "MaxTemp".

```{r}
#data1 and data2 must have the same columns
rbind(data1, data2)
```

Concatenate two tables together horizontally. Example: Two tables, each with the column "Year" (and the same years in each), and one with the column "MaxTemp" and the other with the column "MinTemp"
```{r}
cbind(data1, data2)
```

Looping over objects:

myobjects = list("1" = avecounts_AAD3, "2" = rawcounts_AAD3, "3" = rawcounts_AAD3rep1, "4" = rawcounts_AAD3rep2, "5" = rawcounts_AAD3rep3, "6" = rawcounts_AAD3rep4)

for (i in 1:length(myobjects)) {
  l <- lm(formula = log10nonself ~ log10self, data = myobjects[[i]])
  print(l)
}

# useful packages
```r
library(ggplot2)
library(RColorBrewer) # colour palettes (`paired` is good for lots of samples)
library(wesanderson) # colour palettes
library(dplyr) #for the mutate function
```

# useful resources

[the ggplot flipbook](https://evamaerey.github.io/ggplot_flipbook/ggplot_flipbook_xaringan.html#1) - shows step-by-step construction and modification of ggplot graphs

Dummy R data, see all pre-installed datasets by typing `data()`. Load a given dataset by typing `data(DatasetName)`. WorldPhones is a nice small one.

# graphing

[R Graph Gallery](https://www.r-graph-gallery.com/), a website with templates/tutorials for different types of R graph

[A useful pagge about plotitng multiple grarpphs . on the sae axis](https://www.r-bloggers.com/ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/)

A very basic scatterplot in ggplots:
```{r}
ggplot(my_data, aes(x = var1, y = var2, color = var3)) +
  geom_point() +
  ggtitle("My Title") +
  xlab("the x label") +
  ylab("the y label") +
  labs(x = "the x label", y = "the y label", col = "legend title")
```

We can add any graphical element to our graph with a new line of command (which must be ended with a `+` except for the last line in the ggplot command). Here are some useful graphical elements.
* A trendline, without confidence intervals: `geom_smooth(method = "lm", se = F)
* use specific colours for each category: `scale_color_manual(values = c("colour1", "colour2", "colour3")`
    + Note that the first colour in the list will be given to the first variable in `Levels` for `data$VarToColour`, and so on.
* set x axis: `xlim(c(0,150)`
    + Specify only the lower bound: `xlim(c(0,NA)`
* Alternative title method: `labs(title = "My Title", subtitle = "This is a graph")`
* White background on graph (not grey): `theme_bw()`


Convert a ggplot graph into an interactive graph

```{r}
library(ggplot2)
library(plotly)
# where anothersample is something we want in our hover text that is NOT xvar/yvar/colorvar/fillvar
# this ggplot command can be as complex as you like
p <- ggplot(data, aes(x = xvar, y = yvar, color = colorvar, shape = fillvar, label = anothersample))
ggplotly(p, tooltip = "anothersample")
```

# Working with R notebooks

Convert an R notebook to a script

```{r}
library(knitr)
purl("/path/to/notebook.Rmd", documentation = 1, output = "/path/to/output.R")
# documentation = 0 (no text included), 1 (headers converted to comments), 2 (everything converted to comments)
```

# Some PCA tools you can run
```{r}
## plots from https://www.biostars.org/p/282685/
# first make your PCA data: (where object = vsd or rlog object)
rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)

# Plot proportion of variance explained by each PC
barplot((pcatable.pcvar*100), cex.names=1, xlab=paste("Principal component (PC), 1-", length(pcatable$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))

# pairs plots
par(cex=1.0, cex.axis=0.8, cex.main=0.8)
pairs(pcatable$x[,1:5], col="black", main="Principal components analysis bi-plot\nPCs 1-5", pch=16)
##pairs(pcatable$x[,6:10], col="black", main="Principal components analysis bi-plot\nPCs 6-10", pch=16)

# bi-plots
par(mar=c(4,4,4,4), mfrow=c(1,3), cex=1.0, cex.main=0.8, cex.axis=0.8)
##Plots scatter plot for PC 1 and 2
plot(pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(percentVar[1], 2), "%"), ylab=paste("PC2, ", round(percentVar[2], 2), "%"))
points(pca$x, col="black", pch=16, cex=1)
##Plots scatter plot for PC 1 and 3
plot(pca$x[,1], pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(percentVar[1], 2), "%"), ylab=paste("PC3, ", round(percentVar[3], 2), "%"))
points(pca$x[,1], pca$x[,3], col="black", pch=16, cex=1)
##Plots scatter plot for PC 2 and 3
plot(pca$x[,2], pca$x[,3], type="n", main="Principal components analysis bi-plot", xlab=paste("PC2, ", round(percentVar[2], 2), "%"), ylab=paste("PC3, ", round(percentVar[3], 2), "%"))
points(pca$x[,2], pca$x[,3], col="black", pch=16, cex=1)

# tri-plot
require(scatterplot3d)
par(mar=c(4,4,4,4), cex=1.0, cex.main=0.8, cex.axis=0.8)
scatterplot3d(pca$x[,1:3], angle=-40, main="", color="black", pch=17, xlab=paste("PC1, ", round(percentVar[1], 2), "%"), ylab=paste("PC2, ", round(percentVar[2], 2), "%"), zlab=paste("PC3, ", round(percentVar[3], 2), "%"), grid=FALSE, box=FALSE)
source('http://www.sthda.com/sthda/RDoc/functions/addgrids3d.r')
addgrids3d(pca$x[,1:3], grid = c("xy", "xz", "yz"))

## NB: this is an R package called PCAtools! load this and try
```

# Replicating historical plots with ggplots
https://www.statswithmatt.com/post/ggplot2-meets-w-e-b-du-bois/
