20th January 2022

Seurat issue: https://github.com/satijalab/seurat/issues/5141

I have updated Seurat/some R packages and I am now getting an issue where the SpatialDimPlot and/or SpatialFeaturePlot are forced into a square. This is a problem for subset data, because it distorts the shape of the H+E images.

Here is a reproducible example, using the stxKidney data:

```
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

kidney <- LoadData("stxKidney")
# subset an arbitrary set of points to make a rectangle (my real data is rectangular)
keep <- rownames(kidney@images$image@coordinates)[kidney@images$image@coordinates$imagecol > 5000 & kidney@images$image@coordinates$imagecol < 7000]
kidney_sub <- subset(kidney, cells = keep)
SpatialDimPlot(kidney_sub, crop = TRUE) + ggtitle("crop = true")
SpatialDimPlot(kidney_sub, crop = FALSE) + ggtitle("crop = false")
```
![image](https://user-images.githubusercontent.com/22607689/150287747-4fb85251-af3a-42af-bf76-01739192c8c2.png)

![image](https://user-images.githubusercontent.com/22607689/150287770-eb21a555-5773-440a-ae5d-973809a2272d.png)

However, my desired output (and what Seurat used to do) is to produce something like this (NB: I have cropped this manually as an example):

![image](https://user-images.githubusercontent.com/22607689/150287858-037b5f86-4ba9-4015-933e-589bdd285b8b.png)

This issue arose for me late 2021. I must have updated/downgraded some R packages because at the start of 2022 the images looked as per the desired image above. However, yesterday I had to update some packages and downgrade others, and the issue is back. I am not sure which package is causing the problem. I have installed the most recent CRAN version of Seurat.

I have found a hack/workaround to produce an image like what I want, though it is a bit annoying when you have multiple slices:

```
coord <- GetTissueCoordinates(object = kidney_sub@images$image)
# calculate the aspect ratio of rows to columns
myratio <- (max(coord$imagerow) - min(coord$imagerow)) / (max(coord$imagecol) - min(coord$imagecol))
# force the image into the right aspect ratio
SpatialDimPlot(kidney_sub, crop = TRUE, pt.size.factor = 5) + theme(aspect.ratio = myratio)
```

I'm a bit nervous this isn't perfectly scaling the image, but from comparing my data to the plots I had previously generated with Seurat (when the cropping worked normally) it seems to be OK:
![image](https://user-images.githubusercontent.com/22607689/150287968-5d3b7ba7-49dc-44be-ae03-15e34fbe8be9.png)


From digging into the code, I'm guessing the issue is arising from the aspect.ratio part of SingleSpatialPlot().

```
sessionInfo()
R version 4.0.3 (2020-10-10)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS  11.6.1

Matrix products: default
LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib

locale:
[1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] dplyr_1.0.7                patchwork_1.1.1            ggplot2_3.3.5              stxKidney.SeuratData_0.1.0
[5] stxBrain.SeuratData_0.1.1  ifnb.SeuratData_3.1.0      SeuratData_0.2.1           SeuratObject_4.0.4        
[9] Seurat_4.0.6              

loaded via a namespace (and not attached):
  [1] Rtsne_0.15            colorspace_2.0-2      deldir_1.0-6          ellipsis_0.3.2        ggridges_0.5.3       
  [6] rstudioapi_0.13       spatstat.data_2.1-2   farver_2.1.0          leiden_0.3.9          listenv_0.8.0        
 [11] ggrepel_0.9.1         fansi_1.0.2           codetools_0.2-18      splines_4.0.3         knitr_1.37           
 [16] polyclip_1.10-0       jsonlite_1.7.3        ica_1.0-2             cluster_2.1.2         png_0.1-7            
 [21] uwot_0.1.11           shiny_1.7.1           sctransform_0.3.3     spatstat.sparse_2.1-0 compiler_4.0.3       
 [26] httr_1.4.2            assertthat_0.2.1      Matrix_1.4-0          fastmap_1.1.0         lazyeval_0.2.2       
 [31] cli_3.1.0             later_1.3.0           htmltools_0.5.2       tools_4.0.3           igraph_1.2.11        
 [36] gtable_0.3.0          glue_1.6.0            RANN_2.6.1            reshape2_1.4.4        rappdirs_0.3.3       
 [41] Rcpp_1.0.8            scattermore_0.7       vctrs_0.3.8           nlme_3.1-155          lmtest_0.9-39        
 [46] xfun_0.29             stringr_1.4.0         globals_0.14.0        mime_0.12             miniUI_0.1.1.1       
 [51] lifecycle_1.0.1       irlba_2.3.5           goftest_1.2-3         future_1.23.0         MASS_7.3-55          
 [56] zoo_1.8-9             scales_1.1.1          spatstat.core_2.3-2   promises_1.2.0.1      spatstat.utils_2.3-0 
 [61] parallel_4.0.3        RColorBrewer_1.1-2    yaml_2.2.1            reticulate_1.23       pbapply_1.5-0        
 [66] gridExtra_2.3         rpart_4.1-15          stringi_1.7.6         rlang_0.4.12          pkgconfig_2.0.3      
 [71] matrixStats_0.61.0    evaluate_0.14         lattice_0.20-45       ROCR_1.0-11           purrr_0.3.4          
 [76] tensor_1.5            labeling_0.4.2        htmlwidgets_1.5.4     cowplot_1.1.1         tidyselect_1.1.1     
 [81] parallelly_1.30.0     RcppAnnoy_0.0.19      plyr_1.8.6            magrittr_2.0.1        R6_2.5.1             
 [86] generics_0.1.1        DBI_1.1.2             pillar_1.6.4          withr_2.4.3           mgcv_1.8-38          
 [91] fitdistrplus_1.1-6    survival_3.2-13       abind_1.4-5           tibble_3.1.6          future.apply_1.8.1   
 [96] crayon_1.4.2          KernSmooth_2.23-20    utf8_1.2.2            spatstat.geom_2.3-1   plotly_4.10.0        
[101] rmarkdown_2.11        grid_4.0.3            data.table_1.14.2     digest_0.6.29         xtable_1.8-4         
[106] tidyr_1.1.4           httpuv_1.6.5          munsell_0.5.0         viridisLite_0.4.0    
```
