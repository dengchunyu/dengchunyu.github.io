---
layout: post
title:  "Running scPagwas steps by SubFunctions"
date:   2023-03-10 19:31:29 +0900
categories: Sample
---

There are some different situations for running scPagwas!
## 1. Only run celltypes functions.

The executive programs for celltypes and single cell are independent, If you only want to know the celtypes, set the `celltype=T` and `singlecell=F`.
The advantages is save a lot of times, for celltype only can omit many running processes.

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

## 2. Continue to run singlecell functions  
Because the parameters and input data are the same with celltypes function, we can take the celltypes result as the input data for single cell function. The advantage is there is no need to 
run the `svd` code block, save a lot of time. 

```ruby
system.time(
Pagwas_singlecell<-scPagwas_main(Pagwas =NULL, 
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
)

```
The result is seurat format.

## 3. Running scPagwas step by step

To utilize and understand scPagwas better, we run scPagwas flexibly step by step. In usual, we needn't run these sub-functions one by one.

We use an example provided by scPagwas package.

### 3.1 Single data input

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
```

The result list including single cell expression matrix and cell types mean expression matrix.

### 3.2 Run pathway pca score 

Obtain the pathway pca score matrix.

```ruby
Pagwas <- Pathway_pcascore_run(
        Pagwas = Pagwas,
        Pathway_list = Genes_by_pathway_kegg
      )
names(Pagwas)
```


The result list including `pca_scCell_mat` single cell pca score matrix and `pca_cell_df` cell types pca score matrix.

### 3.3 GWAS summary data input

Read the GWAS summary data, remove the MHC and sex chromosome and filtered the maf of SNP.

```ruby
gwas_data <- bigreadr::fread2(system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"))
Pagwas <- GWAS_summary_input(
    Pagwas = Pagwas,
    gwas_data = gwas_data,
    maf_filter = 0.1
  )
names(Pagwas)
```
