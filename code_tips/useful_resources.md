note to self so far i've done the bookmarks tabs from chrome for DGE and Bioinformtics

# RNA-seq analysis

## Differential gene expression

### Official vignettes by creators

* [RNA-seq workflow](http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html) R vignette, includes data exploration, DESeq, batch effects etc.
* [DESeq2 vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
* [txImport vignette](https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html)
* [R vignette - salmon to DGE](https://www.bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html) by Mike Love et al.
* [edgeR user guide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) (PDF)

### Other tutorials

* [RNA-seq tutorial with DESeq2](http://folk.uio.no/jonbra/MBV-INF4410_2017/exercises/2017-12-07_R_DESeq2_exercises_with_results.html) an online tutorial
* [Differential expression tutorial in EdgeR](http://www.compbio.dundee.ac.uk/user/pschofield/Projects/teaching_pg/workshops/biocDGE.html) By Pietà Schofield.
* [RNASeq tutorial - preprocessing data](https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/02_Preprocessing_Data.nb.html) Cruk suummer schoolschool
* [DGE with limma-voom](https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html) course material from UC davis.

## Brownian motion and OU modelling

* [Various online course content](http://www.eve.ucdavis.edu/~wainwrightlab/Roi/Site/Teaching.html) BY Sam Price and Roi Holzman, UC Davis. Contains OUCH tutorial.
* [evee tools](https://evee-tools.github.io/) scripts associated with Chen et al. 2019 OU paper

# Single cell RNASeq (scRNASeq)

* Seurat has many [tutorials](https://satijalab.org/seurat/vignettes.html); pf particular interest is the [basic clustering vignette](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html), the [label transfer vignette](https://satijalab.org/seurat/v3.1/integration.html) and the [control vs. stimulated PBMCs vignette](https://satijalab.org/seurat/v3.1/immune_alignment.html). Also potentially of interest are the data downloads linked with the various vignettes, e.g. the [Microwell Mouse Cell Atlas data](https://satijalab.org/seurat/v3.1/mca.html). See also this useful ["cheat sheet" of Seurat commands](https://satijalab.org/seurat/essential_commands.html).

* Ming Tang is a bioinformatician at Harvard. He often discusses scRNASeq on [twitter](https://twitter.com/tangming2005). He also has a [blog](https://divingintogeneticsandgenomics.rbind.io/) and a very useful [scRNAseq tutorial](https://crazyhottommy.github.io/scRNA-seq-workshop-Fall-2019/index.html). This tute runs over the same basic content as the Seurat vignettes, but also includes useful links, explanations and "DIY versions" of different plots. He also has his own long list of scRNASeq resources on [github](https://github.com/crazyhottommy/scRNAseq-analysis-notes).

* The Broad Institute has published a long set of notes for an scRNASeq workshop. Of particular interest is [section 9.4]((https://broadinstitute.github.io/2019_scWorkshop/data-wrangling-scrnaseq.html#beginning-with-seurat-httpsatijalab.orgseurat) which discuss Seurat.

* Another [Seurat tutorial](https://www.fimm.fi/sites/default/files/Seurat-guideline-10x.pdf) from FIMM (Institute for Molecular Medicine Finland); this tutorial is not very detailed.

* Yet [another tutorial](https://nbisweden.github.io/excelerate-scRNAseq/session-integration/Data_Integration.html); this one involves data integration.

* [And another](https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html)

* An extensive [Introduction to single-cell RNA-seq analysis](http://barc.wi.mit.edu/education/hot_topics/scRNAseq_March2019/SingleCellRNAseq.pdf) from MIT (powerpoint PDF) which introduces the theory (with a bit of code) through scRNASeq technology, experimental design, analysis pipeline, pseudotime...

* UMAP algorithm resources: [How Exactly UMAP Works](https://towardsdatascience.com/how-exactly-umap-works-13e3040e1668 ) blog post, [Understanding UMAP](https://pair-code.github.io/understanding-umap/) overview, introduction to the UMAP algorithm at [SciPy 2018](https://www.youtube.com/watch?v=nq6iPZVUxZU) by the original author.

* Dave Tang is a computational biologist in Perth at UWA. He has written about scRNASeq [on his computational biology and genomics blog](https://davetang.org/muse/category/single-cell-2/). However, these posts (especially the ones about [Seurat](https://davetang.org/muse/2017/08/01/getting-started-seurat/) and [Monocle](https://davetang.org/muse/2017/10/01/getting-started-monocle/) are quite old (2017-2018).

# Other transcriptomics

* [Grouper: Clustering contigs](https://github.com/COMBINE-lab/grouper)


# General coding resources

* [Search engine for symbols](http://symbolhound.com/)
* [regex tester 1](https://regex101.com/) and [regex tester 2](https://www.regextester.com/)
* [markkdown generator](https://dillinger.io/)
* [markdown table generator](https://www.tablesgenerator.com/markdown_tables)
* [GUI to help change your bashrc profile](http://ezprompt.net/)
* [Tao Te Programming](https://www.burns-stat.com/documents/books/tao-te-programming/table-of-contents/) Patrick Burns

## Bash

* [Explain Shell](https://explainshell.com/)
* [Shell Check](https://www.shellcheck.net/)
* various bash scripting guides:
    + [Handy bash one-liners](https://github.com/onceupon/Bash-Oneliner/blob/master/README.md)
    + [how to bash](http://tldp.org/HOWTO/Bash-Prog-Intro-HOWTO.html#toc8)
    + [bash templates](https://www.networkworld.com/article/2694348/unix--scripting-with-templates.html)
    + [writing robust bash scrippts](https://www.davidpashley.com/articles/writing-robust-shell-scripts/)
    + [advanced bash scripting guide](http://www.tldp.org/LDP/abs/html/index.html)
    + [coding style guide](https://github.com/robbyrussell/oh-my-zsh/wiki/Coding-style-guide)

## R

* [R for Data Science](https://r4ds.had.co.nz/) G Grolemund and H Wickham

* [R inferno book](https://www.burns-stat.com/pages/Tutor/R_inferno.pdf) Patric Burns 2011

* [Impatient R book](https://www.burns-stat.com/documents/tutorials/impatient-r/)

* [Data transformation in R and dplyr](https://craig.rbind.io/post/2019-12-30-asgr-2-1-data-transformation-part-1/)

* [A Scientist's Guide to R](https://craig.rbind.io/post/2019-05-17-asgr-basic-workflow/)

### Writing nice code

* [Principles for structuring R projects](https://chrisvoncsefalvay.com/2018/08/09/structuring-r-projects/)

* [More principles for structuring R projects](https://nicercode.github.io/blog/2013-04-05-projects/)

* [Software Carpentry: R code best practices](https://swcarpentry.github.io/r-novice-inflammation/06-best-practices-R/)

### R data viz

* Data viz in R [lecture](http://datacarpentry.org/semester-biology/materials/ggplot/) and [exercises](http://datacarpentry.org/semester-biology/assignments/r-datavis/) and [a textbook](https://r4ds.had.co.nz/data-visualisation.html)

* [R Graph Gallery](https://www.r-graph-gallery.com/)

* [R colour picker](https://deanattali.com/blog/colourpicker-ggmarginal-gadgets/)

## Bioinformatics

* [Handy bash commands for fasta files](https://www.biostars.org/p/17680/)
* [blastP](https://blast-ncbi-nlm-nih-gov.docelec.univ-lyon1.fr/Blast.cgi?PROGRAM=blastp&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome)

## Statistics

* Online course [Improving your statistical inferences](https://www.coursera.org/learn/statistical-inferences)

* [How PCA works](https://divingintogeneticsandgenomics.rbind.io/post/pca-in-action/) by Ming Tang

* The [StatQuest Youtube page](https://www.youtube.com/channel/UCtYLUTtgS3k1Fg4y5tAhLbw) has a lot of statistical explanation videos, and in particular a playlist of videos about [high throughput sequencing](https://www.youtube.com/playlist?list=PLblh5JKOoLUJo2Q6xK4tZElbIvAACEykp). For instance, see his short videos about [PCA](https://www.youtube.com/watch?v=HMOI_lkzW08) or [tSNE](https://www.youtube.com/watch?v=NEaUSP4YerM).


## Maths

* [Equation solver](https://www.symbolab.com/solver/equation-calculator/) takes an equation, solves it, and shows you the steps

