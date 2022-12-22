# How to plot multi-panel figures in ggplot

```{r}
# load the modules
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(scales)
```

```{r}
# load the inputs and outputs
df <- readRDS("/path/to/SeuratObject.RDS")
Idents(df) <- "myclusters"
DefaultAssay(df) <- "Spatial"

plot_param <- read.delim("/path/to/plotting_parameters_EachTissueSeparate.txt", row.names = 2)

outdir <- "/path/to/outdir/"

sig_all <- readRDS("/path/to/genesignatures.RDS")
```

# calculate gene signature scores

```{r}
# for "all signatures"
for (i in 1:length(sig_all)) {
  df <- AddModuleScore(df, features = sig_all[i], assay = "Spatial", name = names(sig_all)[i])
}
```

# plot in spatialfeatureplot

## first subset df into separate objects for each sample

```{r}
# get cell names
B1_cells <- names(which(df$tissues == "B1"))
B2_cells <- names(which(df$tissues == "B2"))
D1_cells <- names(which(df$tissues == "D1"))
WTA_cells <- names(which(df$tissues == "A_WT"))
WTB_cells <- names(which(df$tissues == "B_WT"))
WTC_cells <- names(which(df$tissues == "C_WT"))
WTD_cells <- names(which(df$tissues == "D_WT"))

# subset
B1 <- subset(df, cells = B1_cells)
B2 <- subset(df, cells = B2_cells)
D1 <- subset(df, cells = D1_cells)
WTA <- subset(df, cells = WTA_cells)
WTB <- subset(df, cells = WTB_cells)
WTC <- subset(df, cells = WTC_cells)
WTD <- subset(df, cells = WTD_cells)

# remove empty image slots from each new file
B1@images$ArrayA <- NULL
B1@images$ArrayD <- NULL
B1@images$WT_A <- NULL
B1@images$WT_B <- NULL
B1@images$WT_C <- NULL
B1@images$WT_D <- NULL

B2@images$ArrayA <- NULL
B2@images$ArrayD <- NULL
B2@images$WT_A <- NULL
B2@images$WT_B <- NULL
B2@images$WT_C <- NULL
B2@images$WT_D <- NULL

D1@images$ArrayA <- NULL
D1@images$ArrayB <- NULL
D1@images$WT_A <- NULL
D1@images$WT_B <- NULL
D1@images$WT_C <- NULL
D1@images$WT_D <- NULL

WTA@images$ArrayA <- NULL
WTA@images$ArrayB <- NULL
WTA@images$ArrayD <- NULL
WTA@images$WT_B <- NULL
WTA@images$WT_C <- NULL
WTA@images$WT_D <- NULL

WTB@images$ArrayA <- NULL
WTB@images$ArrayB <- NULL
WTB@images$ArrayD <- NULL
WTB@images$WT_A <- NULL
WTB@images$WT_C <- NULL
WTB@images$WT_D <- NULL

WTC@images$ArrayA <- NULL
WTC@images$ArrayB <- NULL
WTC@images$ArrayD <- NULL
WTC@images$WT_A <- NULL
WTC@images$WT_B <- NULL
WTC@images$WT_D <- NULL

WTD@images$ArrayA <- NULL
WTD@images$ArrayB <- NULL
WTD@images$ArrayD <- NULL
WTD@images$WT_A <- NULL
WTD@images$WT_B <- NULL
WTD@images$WT_C <- NULL
```

```{r}
# generate a full-page figure for each gene in "sig_all"
for (i in paste0(names(sig_all), "1")) {
  # determine the maximum and minimum values across all datasets so that all plots share a common scale
  # the best way to do this depends on your data
    temp <- FeaturePlot(df, features = i)  
    minval <- min(temp$data[,i])
    maxval <- max(temp$data[,i])
    minrange <- round(min(abs(minval), abs(maxval)), 2)
    maxrange <- round(max(abs(minval), abs(maxval)), 2)
    midrange <- minrange/2
  
  # now plot a graph for each sample separately, and save each plot to an object
  # B1
  p_B1 <- SpatialFeaturePlot(B1, features = i, pt.size.factor = plot_param["B1","spot_size"]) +
      theme(aspect.ratio = plot_param["B1","ratio"]) +
    scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))
  # B2
  p_B2 <- SpatialFeaturePlot(B2, features = i, pt.size.factor = plot_param["B2","spot_size"]) +
      theme(aspect.ratio = plot_param["B2","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))

    # D1
  p_D1 <- SpatialFeaturePlot(D1, features = i, pt.size.factor = plot_param["D1","spot_size"]) +
      theme(aspect.ratio = plot_param["D1","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))

  # WTA
  p_WTA <- SpatialFeaturePlot(WTA, features = i, pt.size.factor = plot_param["WTA","spot_size"]) +
      theme(aspect.ratio = plot_param["WTA","ratio"]) +
     scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))

  # WTB
  p_WTB <- SpatialFeaturePlot(WTB, features = i, pt.size.factor = plot_param["WTB","spot_size"]) +
      theme(aspect.ratio = plot_param["WTB","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))

  # WTC
  p_WTC <- SpatialFeaturePlot(WTC, features = i, pt.size.factor = plot_param["WTC","spot_size"]) +
      theme(aspect.ratio = plot_param["WTC","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))

  # WTD
  p_WTD <- SpatialFeaturePlot(WTD, features = i, pt.size.factor = plot_param["WTD","spot_size"]) +
      theme(aspect.ratio = plot_param["WTD","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))
  
  # in this case, the legend is the same for every figure, so I want to just include it once
  # first I generate a plot again (it doesn't matter which dataset I use, because they all have the same legend) and focus on making the legend pretty
  # I'm not actually going to include this in my final figure
  p_WTD_legend <- SpatialFeaturePlot(WTD, features = i, pt.size.factor = plot_param["WTD","spot_size"]) +
      theme(aspect.ratio = plot_param["WTD","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(angle = -90, size = 6, hjust = 0))
  
  # now I extract the legend from the plot I just made; this one WILL go in the final figure
  p_legend <- as_ggplot(get_legend(p_WTD_legend))
  
  # plot all graphs together
  # NoLegend() is a wrapper function from the Seurat package, but it works on all ggplot objects
  # if you don't have this package you can generate your plots above but remove the legend manually (google how to do it)
  p_all <- ggarrange(p_B1 + NoLegend(),
            p_D1 + NoLegend(),
            p_B2 + NoLegend(),
            p_WTA + NoLegend(),
            p_WTB + NoLegend(),
            p_WTC + NoLegend(),
            p_WTD + NoLegend(),
            p_legend,
            labels = c("d1", "d3", "d7", "WT_A", "WT_B", "WT_C", "WT_D"),
            font.label = list(size = 8),
            ncol = 3, nrow = 3)
  
  # add title to the graph, this is optional
  # I've done gsub... to automatically generate the title but you could just enter free text here
    title <- text_grob(gsub("1$", "", i) %>% gsub("_", " ", .), size = 8, face = "bold")
    mygenes <- sig_all[gsub("1$", "", i)][[1]][1:10]
    subtitle <- text_grob(paste0(paste(mygenes, collapse = ", "), "....."), size = 8, face = "italic")
    annotate_figure(p_all, top = title, bottom = subtitle)
  
  ggsave(filename = paste0(outdir, "modulescore_", i, ".pdf"), width = 5, height = 5)
  ggsave(filename = paste0(outdir, "modulescore_", i, ".jpg"), width = 5, height = 5)

  }
```

# same thing but lesion zones

## first extract the tissues

```{r}
plotting_param_lesion <- read.delim("/Volumes/SciDrive/Archive/SCI/MND_WT/1_IntegrateWithSCI/GetPlottingCoordinates/outdir_lesionOnly/plotting_parameters_LesionZones.txt", row.names = 2)

# get cells
B1_cells_lesion <- names(which(df$tissues == "B1" & df$lesion_zones5 == "lesion"))
B2_cells_lesion <- names(which(df$tissues == "B2" & df$lesion_zones5 == "lesion"))
D1_cells_lesion <- names(which(df$tissues == "D1" & df$lesion_zones5 == "lesion"))
temp_all_cells <- c(B1_cells_lesion, B2_cells_lesion, D1_cells_lesion, WTA_cells, WTB_cells, WTC_cells, WTD_cells)

# subset
B1l <- subset(df, cells = B1_cells_lesion)
B2l <- subset(df, cells = B2_cells_lesion)
D1l <- subset(df, cells = D1_cells_lesion)
df_l <- subset(df, cells = temp_all_cells)

# remove images
B1l@images$ArrayA <- NULL
B1l@images$ArrayD <- NULL
B1l@images$WT_A <- NULL
B1l@images$WT_B <- NULL
B1l@images$WT_C <- NULL
B1l@images$WT_D <- NULL

B2l@images$ArrayA <- NULL
B2l@images$ArrayD <- NULL
B2l@images$WT_A <- NULL
B2l@images$WT_B <- NULL
B2l@images$WT_C <- NULL
B2l@images$WT_D <- NULL

D1l@images$ArrayA <- NULL
D1l@images$ArrayB <- NULL
D1l@images$WT_A <- NULL
D1l@images$WT_B <- NULL
D1l@images$WT_C <- NULL
D1l@images$WT_D <- NULL
```

```{r message=FALSE, warning=FALSE}
dir.create(file.path(outdir, "spatial_plots_lesion"), showWarnings = FALSE)

for (i in paste0(names(sig_all), "1")) {
    # calculate the range
    temp <- FeaturePlot(df_l, features = i)  
    minval <- min(temp$data[,i])
    maxval <- max(temp$data[,i])
    minrange <- round(min(abs(minval), abs(maxval)), 2)
    maxrange <- round(max(abs(minval), abs(maxval)), 2)
    midrange <- minrange/2
  
  # B1
  p_B1 <- SpatialFeaturePlot(B1l, features = i, pt.size.factor = plotting_param_lesion["B1","spot_size"]) +
      theme(aspect.ratio = plotting_param_lesion["B1","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))

  # B2
  p_B2 <- SpatialFeaturePlot(B2l, features = i, pt.size.factor = plotting_param_lesion["B2","spot_size"]) +
      theme(aspect.ratio = plotting_param_lesion["B2","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))
  
    # D1
  p_D1 <- SpatialFeaturePlot(D1l, features = i, pt.size.factor = plotting_param_lesion["D1","spot_size"]) +
      theme(aspect.ratio = plotting_param_lesion["D1","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))
  
  # WTA
  p_WTA <- SpatialFeaturePlot(WTA, features = i, pt.size.factor = plotting_param_lesion["WTA","spot_size"]) +
      theme(aspect.ratio = plotting_param_lesion["WTA","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))

  # WTB
  p_WTB <- SpatialFeaturePlot(WTB, features = i, pt.size.factor = plotting_param_lesion["WTB","spot_size"]) +
      theme(aspect.ratio = plotting_param_lesion["WTB","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))

  # WTC
  p_WTC <- SpatialFeaturePlot(WTC, features = i, pt.size.factor = plotting_param_lesion["WTC","spot_size"]) +
      theme(aspect.ratio = plotting_param_lesion["WTC","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))
  
    # WTD
  p_WTD <- SpatialFeaturePlot(WTD, features = i, pt.size.factor = plotting_param_lesion["WTD","spot_size"]) +
      theme(aspect.ratio = plotting_param_lesion["WTD","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange))
  
  # legend
  p_WTD_legend <- SpatialFeaturePlot(WTD, features = i, pt.size.factor = plot_param["WTD","spot_size"]) +
      theme(aspect.ratio = plot_param["WTD","ratio"]) +
      scale_fill_gradientn(colours = c("#2166AC", "#2166AC", "#92C5DE", "#F7F7F7", "#F4A582", "#B2182B", "#B2182B"),
                         breaks = c(-maxrange, -minrange, -midrange, 0, midrange, minrange, maxrange),
                         limits = c(-maxrange, maxrange)) +
    theme(legend.title = element_blank(),
          legend.text = element_text(angle = -90, size = 6, hjust = 0))
  
  p_legend <- as_ggplot(get_legend(p_WTD_legend))
  
  # plot all graphs together
  p_all <- ggarrange(p_B1 + NoLegend(),
            p_D1 + NoLegend(),
            p_B2 + NoLegend(),
            p_WTA + NoLegend(),
            p_WTB + NoLegend(),
            p_WTC + NoLegend(),
            p_WTD + NoLegend(),
            p_legend,
            labels = c("d1", "d3", "d7", "WT_A", "WT_B", "WT_C", "WT_D"),
            font.label = list(size = 8),
            ncol = 3, nrow = 3)
  # add title
    title <- text_grob(gsub("1$", "", i) %>% gsub("_", " ", .), size = 8, face = "bold")
    mygenes <- sig_all[gsub("1$", "", i)][[1]][1:10]
    subtitle <- text_grob(paste0(paste(mygenes, collapse = ", "), "....."), size = 8, face = "italic")
    annotate_figure(p_all, top = title, bottom = subtitle)

  ggsave(filename = paste0(outdir, "spatial_plots_lesion/",  "modulescore_", i, ".pdf"), width = 5, height = 5)
  ggsave(filename = paste0(outdir, "spatial_plots_lesion/",  "modulescore_", i, ".jpg"), width = 5, height = 5)
  
  }
```
