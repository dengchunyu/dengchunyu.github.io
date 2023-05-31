---
layout: page
title: About
permalink: /about/
---

<img src="/public/img/logo.png" width="10%" style="display: block; margin: auto;" />

**scPagwas** employing the polygenic regression model to uncover trait-relevant cell subpopulations by incorporating pathway activity transformed scRNA-seq data with genome-wide association studies (GWAS) data.

<img src="/public/img/Figure 1_v3_00.png" width="60%" style="display: block; margin: auto;" />

The methodology and benchmarking performance are described in: 

> Polygenic regression uncovers trait-relevant cellular contexts through pathway activation transformation of single-cell RNA sequencing data.(medRxiv.2023) 

Code for reproducing the analysis from the paper is available [here](https://github.com/dengchunyu/scPagwas_reproduce). 
You can install the released version of scPagwas from [github](https://github.com/sulab-wmu/scPagwas) with: 

```ruby
#install some dependence packages
install.packages("Seurat")
install.packages("ggpubr")
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("GenomicRanges")
BiocManager::install("IRanges")

devtools::install_github("sulab-wmu/scPagwas")
```
## Usage 
quick-start example: 
```ruby
library(scPagwas)
 #1.start to run the wrapper functions for example.
 Pagwas_data<-scPagwas_main(Pagwas = NULL,
                     gwas_data =system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"), # The GWAS Summary statistics files 
                     Single_data =system.file("extdata", "scRNAexample.rds", package = "scPagwas"),# scRNA-seq data in seruat format with "RNA" assays and normalized.
                     output.prefix="test", # the prefix name for output files
                     output.dirs="scPagwastest_output",# the directory file's name for output
                     block_annotation = block_annotation,# gene position in chromosome is provided by package.
                     assay="RNA", # the assays for scRNA-seq data to use.
                     Pathway_list=Genes_by_pathway_kegg,# pathway list is provided by package, including gene symbols.
                     n.cores=1,
                     iters_singlecell = 100,
                     chrom_ld = chrom_ld,# The LD data is provided by package.
                     singlecell=T, # Whether to run the singlecell process.
                     celltype=T# Whether to run the celltype process.
)
```
