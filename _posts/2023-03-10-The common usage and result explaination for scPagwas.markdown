---
layout: post
title:  "The common usage and result explaination for scPagwas!"
date:   2023-03-10 19:31:29 +0900
categories: Sample
---
## Preprocessed data
This is a real example for the artical. We use BMMC scRNA-seq data and monocyte count trait to test scPagwas.
The processed gwas data can be download from [here](https://1drv.ms/t/s!As-aKqXDnDUHi6sx7Hqblj2Sgl7P8w?e=cb5Ihf).

BMMC example scRNA-seq data can be download from [here](https://1drv.ms/u/s!As-aKqXDnDUHi9pNciEgQFbe-AHgLw?e=4JiHRw).

## 1. Compute the singlecell and celltype result for monocytecount trait

In this example, we run the scPagwas for usual steps, both running singlecell and celltype process.
```ruby
library(scPagwas)
Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data ="monocytecount_prune_gwas_data.txt",
                     Single_data ="Seu_Hema_data.rds",
                     output.prefix="monocytecount_scPagwas",
                     output.dirs="monocytecount_bmmc",
                     #seruat_return=T,
                     Pathway_list=Genes_by_pathway_kegg,
                     assay="RNA",
                     singlecell=T, 
                     celltype=T,
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)

```

## 2. Result interpretation

There are two types of result, Seruat format return result and files output;

### 2.1 Pagwas : return result for seruat_return=TRUE

```ruby
> print(Pagwas)
An object of class Seurat 
16488 features across 35582 samples within 3 assays 
Active assay: RNA (15860 features, 5000 variable features)
 2 other assays present: scPagwasPaHeritability, scPagwasPaPca
 3 dimensional reductions calculated: pca, tsne, umap

> names(Pagwas@misc)
[1] "Pathway_list"                  "pca_cell_df"                  
[3] "lm_results"                    "bootstrap_results"            
[5] "scPathways_rankPvalue"         "gene_heritability_correlation"
```

Returns a Seruat data with entries(seruat_return=T):

1.  In this Seruat result, two new assays were added:

    -   **scPagwasPaPca**: An assay for S4 type of data; the svd score result for pathways in each cells;

    -   **scPagwasPaHeritability**: An assay for S4 type of data; the gPas score matrix for pathways in each cells;


2.  In meta.data, four columns were added:

    -   **scPagwas.TRS.Score1**: inheritance related score, enrichment socre for genetics top genes;

    -   **scPagwas.gPAS.score**: Inheritance pathway regression effects score for each cells;

    -   **CellScalepValue**: Ranked CellScalepValue for each cells;

    -   **CellScaleqValue**: Ranked CellScaleqValue for each cells, adjust p value.


3.  A new element names misc is added in result, `Pagwas@misc` is a list data including:

    -   **Pathway_list**: filtered pathway gene list;

    -   **pca_cell_df**: a pahtway and celltype data frame for pathway svd(1'st pca) result for each celltype;

    -   **lm_results**: the regression result for each celltype;

    -   **gene_heritability_correlation**: heritability correlation value for each gene;

    -   **scPathways_rankPvalue**: p values for each pathway in each single cell;

    -   **bootstrap_results**: The bootstrap data frame results for celltypes including bootstrap pvalue and confidence interval.

### 2.1 Pagwas : output files result

In **monocytecount_bmmc** result document folder, several result files are including:

-   **scPagwas.run.log** : The running time log for scPagwas;

-   **\*\_parameters.txt** : The parameters log for scPagwas;

-   **\*\_cellytpes_bootstrap_results.csv** : The result of cellytpes for bootstrap p results;

-   **\*\_gene_heritability_correlation.csv** : The result of gene heritability correlation with gPAS score.

-   **\*\_Pathway_singlecell_lm_results.txt** : The regression result for all pahtway and single cell matrix;

-   **\*\_singlecell_Pathways_rankPvalue.csv** : The pathway pvalue for eache singlecell;

-   **\*\_singlecell_scPagwas_score_pvalue.Result.csv** : The data frame result for each cell inlcuding scPagwas.TRS.Score, scPagwas.gPAS.score, pValueHighScale, qValueHighScale;

**scPagwas_cache** is a temporary folder to save the SOAR data, you can remove it when finish the scPagwas.
Sometimes, If scPagwas is failed of running for some reason, you should remove this folder or run this `SOAR::Remove(SOAR::Objects())` to remove the object in this folder.

## 3. Result Visualization

### 3.1 Visualize the scPagwas Score of single cell data in UMAP or TSNE plot.

1) DimPlot of singlecell data.

```ruby
 require("RColorBrewer")
 require("Seurat")
 require("SeuratObject")
 require("ggsci")
 require("dplyr")
 require("ggplot2")
 require("ggpubr")
 #check the objects
color26 <- c("#D9DD6B","#ECEFA4","#D54C4C","#8D2828","#FDD2BF","#E98580","#DF5E5E","#492F10","#334257","#476072","#548CA8",
"#00A19D","#ECD662","#5D8233","#284E78","#3E215D","#835151","#F08FC0","#C6B4CE","#BB8760","#FFDADA","#3C5186",
"#558776","#E99497","#FFBD9B","#0A1D37")

png("monocyte_bmmc_dimplot_umap.png",width = 600, height = 600)
 Seurat::DimPlot(Pagwas,pt.size=1,reduction="tsne",label = T, repel=TRUE)+
 scale_colour_manual(name = "celltypes", values = color26)+
 umap_theme()+ggtitle("Monocyte BMMC")+labs(x="TSNE",y="")+theme(aspect.ratio=1)
dev.off()

```
![alt_text](/public/img/monocyte_bmmc_dimplot_umap.png)

scPagwas.TRS.Score1 and positive(p<0.05) cells showing in dimplot.

```ruby
 scPagwas_Visualization(Single_data=Pagwas,
                        p_thre = 0.05,
                        FigureType = "tsne",
                        width = 7,
                        height = 7,
                        lowColor = "white", 
                        highColor = "red",
                        output.dirs="figure",
                        size = 0.5,
                        do_plot = T)
```

2) scPagwas.gPAS.score dimplot. 

![alt_text](/public/img/scPagwas.gPAS.score_tsne.png)
3) scPagwas.TRS.Score dimplot. 

![alt_text](/public/img/scPagwas.TRS.Score_tsne.png)
4) scPagwas_CellScaleqValue dimplot. 
![alt_text](/public/img/scPagwas_CellScaleqValue0.05_tsne.png)
Positive cells(q value<0.05): red dot; Negative cells: other dot. 

### 3.2 Plot the barplot of the proportion of positive Cell. 

Plot the barplot of the proportion of positive Cells in celltypes:
```ruby
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                         var_ident="celltypes",
                         var_group="positiveCells",
                         vec_group_colors=c("#E8D0B3","#7EB5A6"),
                         do_plot = T)
```
![alt_text](/public/img/plot_bar_positie_nagtive.png)
Plot the barplot of the proportion of celltypes in positive and negative Cells:
```ruby
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                              var_ident="positiveCells",
                              var_group="celltypes",
                              p_thre = 0.01,
                              vec_group_colors=NULL,
                              f_color=colorRampPalette(brewer.pal(n=10, name="RdYlBu")),
                              do_plot = T)
```
![alt_text](/public/img/plot_bar_positie_nagtive2.png)
### 3.3 Visualize the heritability correlated Pathways

Plot the specific genetics pathway for each celltypes
```ruby
 library(tidyverse)
 library("rhdf5")
 library(ggplot2)
 library(grDevices)
 library(stats)
 library(FactoMineR)
 library(scales)
 library(reshape2)
 library(ggdendro)
 library(grImport2)
 library(gridExtra)
 library(grid)
 library(sisal)

 source(system.file("extdata", "plot_scpathway_contri_dot.R", package = "scPagwas"))
png("plot_scpathway_dot.png",width = 1000, height = 600)
plot_scpathway_dot(Pagwas=Pagwas,
                   celltypes=c("01_HSC","02_Early.Eryth","05_CMP.LMPP","11_CD14.Mono.1","12_CD14.Mono.2","13_CD16.Mono","15_CLP.2","19_CD8.N","20_CD4.N1","21_CD4.N2"),
                   topn_path_celltype=5,
                   filter_p=0.05,
                   max_logp=15,
                   display_max_sizes=F,
                   size_var ="logrankPvalue" ,
                   col_var="proportion",
                   shape.scale = 8,
                   cols.use=c("lightgrey", "#E45826"),
                   dend_x_var = "logrankPvalue",
                   dist_method="euclidean",
                   hclust_method="ward.D",
                   do_plot = T,
                   #figurenames = "Pathway_plot.pdf",
                   width = 7,
                   height = 7)
dev.off()
```
![alt_text](/public/img/plot_scpathway_dot.png)
### 3.4 Plot the heritability correlated genes

```ruby
pdf("heritability_cor_scatterplot.pdf")
heritability_cor_scatterplot(gene_heri_cor=Pagwas@misc$gene_heritability_correlation,
                             topn_genes_label=10,
                             color_low="#035397",
                             color_high ="#F32424",
                             color_mid = "white",
                             text_size=2,
                             do_plot=T,
                             max.overlaps =20,
                             width = 7,
                             height = 7)
dev.off()
```
![alt_text](/public/img/heritability_cor_scatterplot.png)
### 3.5 Show expression of the top heritability correlation genes in celltypes
```ruby
library("ggsci")
library("Seurat")
top5genes<-rownames(Pagwas@misc$gene_heritability_correlation)[order(Pagwas@misc$gene_heritability_correlation,decreasing = T)[1:5]]
pdf("plot_vln_Corgenes.pdf",width = 6, height =7)
plot_vln_Corgenes(seurat_obj=Pagwas,
             assay="RNA", slot="data",
             var_group="celltypes",
             vec_features=top5genes,
             vec_group_colors= color26,
             do_plot = T
             )
dev.off()
```
![alt_text](/public/img/plot_vln_Corgenes.png)
### 3.6 celltypes bootstrap_results reuslt 

Barplot for celltypes 
```ruby
Bootstrap_P_Barplot(p_results=Pagwas@misc$bootstrap_results$bp_value[-1],
                    p_names=rownames(Pagwas@misc$bootstrap_results)[-1],
                    figurenames = "Bootstrap_P_Barplot.pdf",
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = "monocytecount_bmmc")
```
![alt_text](/public/img/Bootstrap_P_Barplot.png)
### Forest plot for estimate value. 
```ruby
pdf("estimate_Plot.pdf")
Bootstrap_estimate_Plot(Pagwas=Pagwas,
                        width = 9,
                        height = 7,
                        do_plot=T)
dev.off()
```
![alt_text](/public/img/estimate_Plot_00.png)

### sessionInfo
```ruby
R version 4.2.2 (2022-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

attached base packages:
[1] stats     graphics  grDevices datasets  utils     methods   base     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.9          compiler_4.2.2      RColorBrewer_1.1-3  progressr_0.11.0   
 [5] tools_4.2.2         digest_0.6.30       evaluate_0.18       Rtsne_0.16         
 [9] lattice_0.20-45     rlang_1.0.6         Matrix_1.5-3        cli_3.3.0          
[13] rstudioapi_0.14     yaml_2.3.6          parallel_4.2.2      xfun_0.34          
[17] fastmap_1.1.0       knitr_1.40          SOAR_0.99-11        bigreadr_0.2.4     
[21] bigassertr_0.1.5    globals_0.16.1      grid_4.2.2          data.table_1.14.4  
[25] listenv_0.8.0       future.apply_1.10.0 RcppAnnoy_0.0.20    parallelly_1.32.1  
[29] RANN_2.6.1          rmarkdown_2.18      ROCR_1.0-11         codetools_0.2-18   
[33] htmltools_0.5.3     MASS_7.3-58.1       future_1.29.0       renv_0.16.0        
[37] KernSmooth_2.23-20 
```

