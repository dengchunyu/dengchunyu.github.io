---
layout: post
title: Data input and preproccess for scPagwas
date: 2023-03-10 19:20:23 +0900
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
![alt_text](/public/img/figure1.1.png)

### 2.2.Take the monocyte count as real example and deal with it step by step:

Preprogress the raw monocytecount GWAS Summary statistics file.
Extract the useful information from raw files.
Download the Raw gwas GWAS Summary statistics file: <https://gwas.mrcieu.ac.uk/files/ieu-b-31/ieu-b-31.vcf.gz>
```ruby
library(readr)
############read the raw vcf file
monocytecount<-"ieu-b-31.vcf.gz"
GWAS_raw <-read_table2(monocytecount,comment = "#")
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a type of file:split the FORMAT and IEU coloumn, which including different information.
value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- unlist(strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T))
  value <- unlist(strsplit(GWAS_raw$IEU[i], split = ":",fixed=T))
  names(value)<-index
  value<-value[c("ES","SE","LP","AF","SS")]
  return(value)
  })

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=5)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 colnames(gwas_data)
#[1] "CHROM" "POS"   "ID"    "REF"   "ALT"   "ES"    "SE"    "LP"    "AF"   
#[10] "SS"
#change the colnames:
colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf","N")
#select the useful information for scPagwas.
gwas_data<-gwas_data[,c("chrom","pos","rsid","beta","se","maf")]
write.table(gwas_data,file="monocytecount_gwas_data.txt",row.names=F,quote=F)
```

## 3.Pathway gene list

There are some processed pathway gene list provided by scPagwas package: all these pathway list are including gene symbols.

-   Genes_by_pathway_kegg
-   genes.by.celltype.pathway
-   genes.by.gpbp.pathway
-   genes.by.reactome.pathway
-   genes.by.regulatory.pathway
-   genes.by.tft.pathway
-   genes.by.hallmark.pathway
-   genes.by.immunologic.pathway
-   genes.by.immunesigdb.pathway

There are some steps for get these pathway gene list: There is no need to repeat running these code.

### KEGG pathway list
```ruby
library(KEGGREST)
pathways.list <- keggList("pathway", "hsa")
# Pull all genes for each pathway
pathway.codes <- sub("path:", "", names(pathways.list))	
genes.by.pathway_kegg <- sapply(pathways.list,
                           function(pwid){
                             pw <- keggGet(pwid)
                             if (is.null(pw[[1]]$GENE)) return(NA)
                             pw2 <- pw[[1]]$GENE[c(FALSE,TRUE)] # may need to modify this to c(FALSE, TRUE) for other organisms
                             pw2 <- unlist(lapply(strsplit(pw2, split = ";", fixed = T), function(x)x[1]))
                             return(pw2)
                           })
```
### Other pathway list

These pathways data are download from [GSEA](https://www.gsea-msigdb.org/gsea/index.jsp) website.
```
x <- readLines("c2.cp.reactome.v7.5.1.symbols.gmt")
res <- strsplit(x, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
genes.by.reactome.pathway <- lapply(res, "[", -c(1:2))
save(genes.by.reactome.pathway,file = "genes.by.reactome.pathway.RData")
```
Note:

1.  Sometimes, there is no need to change the pathway list frequently once you choose a fitable pathway list.
2.  The summed number of genes in pathway list should not too much to cost a long time or memory.Once your pathway list is too big, you can choose some sub-pathway suitable for you analysisi or remove some redundance pathways.

### Reduce pathway gene list

We set the `genes.by.reactome.pathway` for example, there are 1615 gene list and 89476 genes, a mount of resouces and time will be cost when running scPagwas, we provid `reduce_pathway` to reduce the pathway list:

-   `pathway_seed` choose a list of pathway names as seed, which are some pahtways cann't be remove.
-   `remove_proporion` The propotion of duplicated between seed pathway and the others

```ruby
reduce_genes.by.reactome.pathway<-scPagwas::reduce_pathway(pathway_seed=names(genes.by.reactome.pathway)[sample(1:length(genes.by.reactome.pathway),500)],
                                                 pathway_list=genes.by.reactome.pathway,
                                                 remove_proporion=0.6)
length(reduce_genes.by.reactome.pathway)
#[1] 1045
length(unlist(reduce_genes.by.reactome.pathway))
#[1] 34980

##C7: immunologic signature gene sets
x <- readLines("c7.all.v2022.1.Hs.symbols.gmt")
res <- strsplit(x, "\t")
names(res) <- vapply(res, function(y) y[1], character(1))
genes.by.immunologic.pathway <- lapply(res, "[", -c(1:2))

###The length of Pathway will not lager than 1000 and not smaller than 100.
reduce_genes.by.immunologic.pathway<-genes.by.immunologic.pathway[unlist(lapply(genes.by.immunologic.pathway,function(x) length(x)>100))]
reduce_genes.by.immunologic.pathway<-genes.by.immunologic.pathway[unlist(lapply(genes.by.immunologic.pathway,function(x) length(x)<1000))]

set.seed(1234)
for (i in 1:30) {
  reduce_genes.by.immunologic.pathway<-scPagwas::reduce_pathway(
	pathway_seed=names(reduce_genes.by.immunologic.pathway)[sample(1:length(reduce_genes.by.immunologic.pathway),300)],
	pathway_list=reduce_genes.by.immunologic.pathway,
	remove_proporion=0.1)
}
```

## 4.Gene block annotation

Gene block annotation in chromosome is need to prepare for scPagwas.
scPagwas can provide a block annotation data for protein-coding genes.

File downloaded from MAGMA website [MAGMA | CTG (cncr.nl)](https://ctg.cncr.nl/software/magma)

Here is the procedure of obtaining the data:
```ruby
library("rtracklayer")
gtf_df<- rtracklayer::import("gencode.v34.annotation.gtf.gz")
gtf_df <- as.data.frame(gtf)
gtf_df <- gtf_df[,c("seqnames","start","end","type","gene_name")]
gtf_df <- gtf_df[gtf_df$type=="gene",]
block_annotation<-gtf_df[,c(1,2,3,5)]
colnames(block_annotation)<-c("chrom", "start","end","label")
```
## 5.LD data

The 1,000 Genomes Project Phase 3 Panel was applied to calculate the linkage disequilibrium (LD) among SNPs extracted from GWAS summary statistics.
the processed LD data are show here:

We use `vcftools` and `plink` to deal with the 1,000 Genomes genotypes.


`./vcftools --vcf ./1000genomes_all_genotypes.vcf --plink-tped --out ./1000genomes_all_genotypes
./plink --tfile ./1000genomes_all_genotypes --recode --out ./1000genomes_all_genotypes
./plink --map ./1000genomes_all_genotypes.map --ped ./1000genomes_all_genotypes.ped --allow-no-sex --autosome --r2 --ld-window-kb 1000 --ld-window-r2 0.2 --out ./ld_1000genome`

R environment to print out the scPagwas-needed data.
```ruby
covid_ld<-read.delim("./ld_1000genome.ld")
#remove sex chrome
covid_ld<-covid_ld[!(covid_ld$ %in% 23),]
colnames(covid_ld)[7]<-"R"
#print out the result in chrom number
lapply(unique(covid_ld$CHR_A), function(i){
  a<-data.table(covid_ld[covid_ld$CHR_A == i,])
  file_name <- paste0("./LD/",i,".Rds")
  saveRDS(a, file = file_name)
})
#integrate the data
chrom_ld<-lapply(as.character(1:22),function(chrom){
  chrom_ld_file_path <- paste(ld_folder, '/', chrom, '.Rds', sep = '')
 ld_data <- readRDS(chrom_ld_file_path)[(R**2 > r2_threshold), .(SNP_A, SNP_B, R)]
  return(ld_data)
})
```
