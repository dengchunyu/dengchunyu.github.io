---
layout: post
title:  "Running scPagwas steps by SubFunctions"
date:   2023-03-10 19:31:29 +0900
categories: RoutineUse
---

## Three routine running situations
There are some different situations for running scPagwas!

### 1.1. Run both singlecell and celltypes functions.  
```ruby
library(scPagwas)
Pagwas<-scPagwas_main(Pagwas =NULL,
                     gwas_data =system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"),
                     Single_data =system.file("extdata", "scRNAexample.rds", package = "scPagwas"),
                     output.prefix="Test",
                     output.dirs="Test",
                     Pathway_list=Genes_by_pathway_kegg,
                     assay="RNA",
                     singlecell=T, 
                     celltype=T,
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
print(Pagwas)
An object of class Seurat 
39249 features across 1000 samples within 4 assays 
Active assay: SCT (18715 features, 2000 variable features)
 3 other assays present: RNA, scPagwasPaHeritability, scPagwasPaPca
 2 dimensional reductions calculated: pca, umap
```
Data in the scPagwas return result.  
```ruby
names(Pagwas@meta.data)
# [1] "orig.ident"          "nCount_RNA"          "nFeature_RNA"        "percent.mt"         
# [5] "S.Score"             "G2M.Score"           "Phase"               "old.ident"          
# [9] "nCount_SCT"          "nFeature_SCT"        "SCT_snn_res.0.5"     "seurat_clusters"    
#[13] "SCT_snn_res.0.4"     "SCT_snn_res.0.8"     "anno"                "scPagwas.TRS.Score1"
#[17] "scPagwas.gPAS.score" "CellpValue"          "CellqValue" 
```
- The colomn "scPagwas.TRS.Score1","scPagwas.gPAS.score" is the trait relavent score and genetic pathway score for each cells.
- The colomn "CellpValue" , "CellqValue" is the statistics p value for trait relavent score for each cells.

```ruby
names(Pagwas@misc)
#[1] "Pathway_list"                  "pca_cell_df"                  
#[3] "lm_results"                    "bootstrap_results"            
#[5] "scPathways_rankPvalue"         "gene_heritability_correlation"
```
- "bootstrap_results" is the cell type results; "scPathways_rankPvalue" is the p value for each pahtway in each celltype; "gene_heritability_correlation" is the correlation result for each gene for genetic pathway score, the top genes were used for calculating the trait relavent score.
- "Pathway_list" ,"pca_cell_df" and "lm_results" are the intermediate results for scPagwas.  

```ruby
names(Pagwas@assays)
#[1] "RNA"                    "SCT"                    "scPagwasPaHeritability"
#[4] "scPagwasPaPca"
```
the intermediate results "scPagwasPaHeritability" and "scPagwasPaPca" are the additional assays for scPagwas. 
- "scPagwasPaHeritability" is the Heritability score for each pahtway in each cell.
- "scPagwasPaPca"  is the pca score for each pahtway in each cell.

### 1.2. Only run celltypes functions.

The executive programs for celltypes and single cell can running separatly, If you only want to know the celtypes, set the `celltype=T` and `singlecell=F`.
It can save a lot of times, for omitting many running processes.

```ruby
library(scPagwas)
Pagwas_celltypes<-scPagwas_main(Pagwas =NULL,
                     gwas_data =system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"),
                     Single_data =system.file("extdata", "scRNAexample.rds", package = "scPagwas"),
                     output.prefix="Test",
                     output.dirs="Test",
                     Pathway_list=Genes_by_pathway_kegg,
                     assay="RNA",
                     singlecell=T, 
                     celltype=T,
                     block_annotation = block_annotation,
                     chrom_ld = chrom_ld)
```

The reuslt of celltypes is list format(not seurat format).  
```ruby
names(Pagwas_celltypes)
# [1] "Celltype_anno"     "data_mat"          "VariableFeatures" 
# [4] "merge_scexpr"      "rawPathway_list"   "Pathway_list"     
# [7] "pca_scCell_mat"    "pca_cell_df"       "snp_gene_df"      
#[10] "lm_results"        "bootstrap_results
```

- bootstrap_results: The bootstrap data frame results for celltypes including bootstrap pvalue and confidence interval.
- pca_scCell_mat : a pahtway and cell data matrix for pathway svd(1'st pca) result for each cell;
- pca_cell_df : a pahtway and celltype data matrix for pathway svd(1'st pca) result for each celltype;

- Other elements are the intermediate data. 

### 1.3. Continue to run singlecell functions  
Because the parameters and input data are the same with celltypes function, we can take the celltypes result as the input data for single cell function. The advantage is there is no need to run the `svd` code block, saving a lot of time. 

```ruby
Pagwas_singlecell<-scPagwas_main(Pagwas =Pagwas_celltypes, 
                                 gwas_data =system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"),
                                 Single_data =system.file("extdata", "scRNAexample.rds", package = "scPagwas"),
                                 output.prefix="Test",
                                 output.dirs="Test",
                                 Pathway_list=Genes_by_pathway_kegg,
                                 assay="RNA",
                                 singlecell=T, 
                                 celltype=T,
                                 block_annotation = block_annotation,
                                 chrom_ld = chrom_ld)
print(Pagwas_singlecell)
An object of class Seurat 
39249 features across 1000 samples within 4 assays 
Active assay: SCT (18715 features, 2000 variable features)
 3 other assays present: RNA, scPagwasPaHeritability, scPagwasPaPca
 2 dimensional reductions calculated: pca, umap
```
The result is seurat format.  

## 2. Running scPagwas step by step

To utilize and understand scPagwas better, we run scPagwas flexibly step by step. In usual, we needn't run these sub-functions one by one.Most of the steps can founded in 'scPagwas_main' functions.  

We use an example provided by scPagwas package.

### 2.1 Single data input

Obtain the single cell expression matrix and filter the cluster for minimum cells.

```ruby
library(scPagwas)
Pagwas <- list()
Single_data <- readRDS(system.file("extdata", "scRNAexample.rds", package = "scPagwas"))
Pagwas <- Single_data_input(
      Pagwas = Pagwas,
      assay = "RNA",
      Single_data = Single_data,
      Pathway_list = Genes_by_pathway_kegg
    )
Single_data <- Single_data[, colnames(Pagwas$data_mat)]
names(Pagwas)
[1] "Celltype_anno"    "data_mat"         "VariableFeatures" "merge_scexpr"
```

The result list including single cell expression matrix and cell types mean expression matrix.

### 2.2 Run pathway pca score 

Obtain the pathway pca score matrix.

```ruby
Pagwas <- Pathway_pcascore_run(
        Pagwas = Pagwas,
        Pathway_list = Genes_by_pathway_kegg
      )
names(Pagwas)
[1] "Celltype_anno"    "data_mat"         "VariableFeatures" "merge_scexpr"    
[5] "rawPathway_list"  "Pathway_list"     "pca_scCell_mat"   "pca_cell_df" 
```

The result list including `pca_scCell_mat` single cell pca score matrix and `pca_cell_df` cell types pca score matrix.

### 2.3 GWAS summary data input

Read the GWAS summary data, remove the MHC and sex chromosome and filtered the maf of SNP.

```ruby
gwas_data <- bigreadr::fread2(system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"))
Pagwas <- GWAS_summary_input(
    Pagwas = Pagwas,
    gwas_data = gwas_data,
    maf_filter = 0.1
  )
names(Pagwas)
[1] "Celltype_anno"    "data_mat"         "VariableFeatures" "merge_scexpr"    
[5] "rawPathway_list"  "Pathway_list"     "pca_scCell_mat"   "pca_cell_df"     
[9] "gwas_data"    
```
### 2.4 Mapping Snps to Genes

We set the `marg` is 10KB, means the position for SNP is less 10KB distance to TSS of gene.

```ruby
Pagwas$snp_gene_df <- SnpToGene(
        gwas_data = Pagwas$gwas_data,
        block_annotation = block_annotation,
        marg = 10000
      )
```
SNP to gene dataframe:  
![alt_text](/public/img/snp_gene_df.png) 

### 2.5 Pathway-SNP annotation  
Mapping SNPs to pahtways and getting block data. 
```ruby
Pagwas <- Pathway_annotation_input(
      Pagwas = Pagwas,
      block_annotation = block_annotation
    )
names(Pagwas)
[1] "Celltype_anno"    "data_mat"         "VariableFeatures" "merge_scexpr"    
[5] "rawPathway_list"  "Pathway_list"     "pca_scCell_mat"   "pca_cell_df"     
[9] "gwas_data"        "snp_gene_df"      "pathway_blocks" 
```
### 2.6 Link the pathway blocks to pca score 
Also run the regression function for single cell.   
```ruby
Pagwas <- Link_pathway_blocks_gwas(
      Pagwas = Pagwas,
      chrom_ld = chrom_ld,
      singlecell = T,
      celltype = T,
      backingpath="./temp")
```
The pathway regression result for single cell.   

### 2.7 Perform regression for celltypes 
Run the regression function for celltypes. 
```ruby
Pagwas$lm_results <- Pagwas_perform_regression(Pathway_ld_gwas_data = Pagwas$Pathway_ld_gwas_data)
Pagwas <- Boot_evaluate(Pagwas, bootstrap_iters = 200, part = 0.5)
Pagwas$Pathway_ld_gwas_data <- NULL
```

Show the regression result for celltypes.
```ruby
names(Pagwas)
Pagwas$lm_results
knitr::kable(Pagwas$bootstrap_results)
```
![alt_text](/public/img/bootstrap_result.png)
### 2.8 Construct the scPagwas score
The gPAS scPagwas score mainly to deal with the single-cell regression result.
```ruby
Pagwas <- scPagwas_perform_score(
      Pagwas = Pagwas,
      remove_outlier = TRUE
    )
names(Pagwas)
 [1] "Celltype_anno"          "data_mat"               "VariableFeatures"      
 [4] "merge_scexpr"           "rawPathway_list"        "Pathway_list"          
 [7] "pca_scCell_mat"         "pca_cell_df"            "snp_gene_df"           
[10] "Pathway_sclm_results"   "lm_results"             "bootstrap_results"     
[13] "Pathway_single_results" "scPathways_rankPvalue"  "scPagwas.gPAS.score" 
```
### 2.9 Get the gene heritability correlation 
Run heritability correlation for all genes.
```ruby
    Pagwas$gene_heritability_correlation <- scGet_gene_heritability_correlation(scPagwas.gPAS.score=Pagwas$scPagwas.gPAS.score,
                                                                                data_mat=Pagwas$data_mat)

```
### 2.10 Calculate the TRS score for top genes. 

Calculate the TRS score for top genes by `AddModuleScore` and running the p value for each cell by `rankPvalue`. 

```ruby
scPagwas_topgenes <- names(Pagwas$gene_heritability_correlation[order(Pagwas$gene_heritability_correlation, decreasing = T), ])[1:1000]

Single_data <- Seurat::AddModuleScore(Single_data, assay = "RNA", list(scPagwas_topgenes), name = c("scPagwas.TRS.Score"))
#run pvalue
CellScalepValue <- scGene_scaleP(
      Single_mat = t(data.matrix(
        Seurat::GetAssayData(Single_data, assay = 'RNA')[scPagwas_topgenes, ]
      ))
    )
    
```
All these sub-functions for scPagwas can running dependently, but need to run orderly. The `scPagwas_main` function is a wrapper functions for these sub-functions.


## 3. Result Visualization

### 3.1 Visualize the scPagwas Score of single cell data in UMAP or TSNE plot.

1) DimPlot of singlecell data.
```ruby
 suppressMessages(require("RColorBrewer"))
 suppressMessages(require("Seurat"))
 suppressMessages(require("SeuratObject"))
 suppressMessages(require("ggsci"))
 suppressMessages(require("dplyr"))
 suppressMessages(require("ggplot2"))
 suppressMessages(require("ggpubr"))
 Seurat::DimPlot(Pagwas,pt.size=1,reduction="umap",label = T, repel=TRUE)+
 umap_theme()+ggtitle("Test")+labs(x="UMAP",y="")+theme(aspect.ratio=1)
```
![alt_text](/public/img/steps_Figure1.png)  

scPagwas.TRS.Score1 and positive(p<0.05) cells showing in dimplot.  

```ruby
 scPagwas_Visualization(Single_data=Pagwas_singlecell,
                        p_thre = 0.05,
                        FigureType = "umap",
                        width = 7,
                        height = 7,
                        lowColor = "white", 
                        highColor = "red",
                        output.dirs="figure",
                        size = 0.5,
                        do_plot = T)
```
![alt_text](/public/img/steps_figure2.png)
![alt_text](/public/img/steps_figure3.png)
![alt_text](/public/img/steps_figure4.png)

### 3.2 Plot the barplot of the proportion of positive Cell. 

Plot the barplot of the proportion of positive Cells in celltypes:
```ruby
plot_bar_positie_nagtive(seurat_obj=Pagwas_singlecell,
                         var_ident="anno",
                         var_group="positiveCells",
                         vec_group_colors=c("#E8D0B3","#7EB5A6"),
                         do_plot = T)
```
![alt_text](/public/img/steps_figure5.png)  
Plot the barplot of the proportion of celltypes in positive and negative Cells:

```ruby
plot_bar_positie_nagtive(seurat_obj=Pagwas_singlecell,
                              var_ident="positiveCells",
                              var_group="anno",
                              p_thre = 0.01,
                              vec_group_colors=NULL,
                              f_color=colorRampPalette(brewer.pal(n=10, name="RdYlBu")),
                              do_plot = T)
```
![alt_text](/public/img/steps_figure6.png)  

### 3.3 Plot the heritability correlated genes

```ruby
heritability_cor_scatterplot(gene_heri_cor=Pagwas_singlecell@misc$gene_heritability_correlation,
                             topn_genes_label=10,
                             color_low="#035397",
                             color_high ="#F32424",
                             color_mid = "white",
                             text_size=2,
                             do_plot=T,
                             max.overlaps =20,
                             width = 7,
                             height = 7)
```
![alt_text](/public/img/steps_figure7.png)  

### 3.4 Show expression of the top heritability correlation genes in celltypes

```ruby
suppressMessages(library("ggsci"))
suppressMessages(library("Seurat"))
top5genes<-rownames(Pagwas_singlecell@misc$gene_heritability_correlation)[order(Pagwas_singlecell@misc$gene_heritability_correlation,decreasing = T)[c(4,30,36,55)]]

plot_vln_Corgenes(seurat_obj=Pagwas_singlecell,
             assay="RNA", slot="data",
             var_group="anno",
             vec_features=top5genes,
             do_plot = T
             )
```
![alt_text](/public/img/steps_figure9.png)  

### 4.6 celltypes bootstrap_results reuslt 
Barplot for celltypes 

```ruby
Bootstrap_P_Barplot(p_results=Pagwas_singlecell@misc$bootstrap_results$bp_value[-1],
                    p_names=rownames(Pagwas_singlecell@misc$bootstrap_results)[-1],
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = "Test")
```
![alt_text](/public/img/steps_figure8.png)  

Forest plot for estimate value. 
```ruby
Bootstrap_estimate_Plot(Pagwas=Pagwas_singlecell,
                        width = 9,
                        height = 7,
                        do_plot=T)

```
![alt_text](/public/img/steps_figure10.png)
