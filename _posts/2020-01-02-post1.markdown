---
layout: post
title: Data_input_preproccess_for_scPagwas
date: 2020-01-02 19:20:23 +0900
category: sample
---
# Data input
There were three input resources for scPagwas: scRNA-seq dataset(seruat format), a GWAS summary dataset(txt file) on a given phenotype and an extensive panel of pathways or functional gene sets(gene symbol list).

## 1.Single cell data Input

### 1.1.Downloading scRNA-seq dataset

The example scRNA-seq data for T cells of melanoma was downloaded from the GEO [GSE115978](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115978)
Two data were download from this website:
Single cell annotations data: [GSE115978_cell.annotations.csv.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_cell.annotations.csv.gz)

Single cell count data:[GSE115978_counts.csv.gz](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115978/suppl/GSE115978_counts.csv.gz)
### 1.2.Progressing scRNA-seq dataset

-   1). The scRNA-seq dataset should have Idents for cell types or clusters.
-   2). The scRNA-seq dataset should be normalized.

```ruby
##read the data
library(Seurat)
library(SeuratObject)
counts <- read.csv("GSE115978_counts.csv.gz",row.names=1)
Anno<- read.csv("GSE115978_cell.annotations.csv.gz")
##create the SeuratObject
Single_data<-Seurat::CreateSeuratObject(
  counts,
  assay = "RNA",
  meta.data=Anno
)

Idents(Single_data)<-Single_data$BioClassification

Single_data <- ScaleData(Single_data)
Single_data <- NormalizeData(Single_data, normalization.method = "LogNormalize", scale.factor = 10000)
```

## 2.GWAS summary data Input

GWAS Summary statistics are download from figshare:
[pgc.bip.full.2012-04.txt.gz](https://doi.org/10.6084/m9.figshare.14671995)

### 2.1.Read and progrocess the example GWAS Summary statistics file

In R environment, select the specific columns and output the result.
```ruby
#read the raw data
GWAS_raw <-read_table("finngen_R7_C3_GBM_EXALLC.gz")
#select the specific columns
GWAS_raw<-GWAS_raw[,c(1,2,3,4,5,6,7,9,10,11,12,13)]
colnames(GWAS_raw)<-c("chrom","pos","REF","ALT","rsid","nearest_genes","pval","beta","se","maf")
x<-tolower(x)
x<-str_replace(x,".gz","")
#Output the data
write.table(GWAS_raw,file="finngen_R7_C3_GBM_EXALLC.txt",row.names=F,quote=F)
```

There were more than 10M snp in the GWAS summary datasets. To facilitate the calculation of scPagwas, we prune the snp by plink.
Here, we 

`awk  '{print $3 }' ${i}.txt  > /share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/IEU_opengwas_project/tempfile/${i}_SNP_list.txt`

The GWAS Summary statistics file need to be processed into a "txt" file including six coloumn and tab-delimited.
```ruby
library(scPagwas)
gwas_data <- bigreadr::fread2("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/finngen_R7_C3_GBM_EXALLC.gz")

gwas_data <- bigreadr::fread2(system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"))
knitr::kable(head(gwas_data))
```
