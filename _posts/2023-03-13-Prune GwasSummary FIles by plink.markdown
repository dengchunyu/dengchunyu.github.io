---
layout: post
title:  "Prune GwasSummary FIles by plink"
date:   2023-03-13 19:31:29 +0900
categories: DataPrepare
---
```{r setup}
library(scPagwas)
```


# Prune the gwas summary statistics file by plink

The markdown is Under editing!!
## Reprocess the gwas data
```ruby
  GWAS_raw <-read_table(paste("./",x,".gz",sep=""))
  GWAS_raw<-GWAS_raw[,c(1,2,3,4,5,6,9,10)]
  colnames(GWAS_raw)<-c("chrom","pos","REF","ALT","rsid","nearest_genes","beta","se")
  GWAS_raw$maf<-0.1
  
  write.table(GWAS_raw,file=paste("./",tolower(x),".txt",sep=""),row.names=F,quote=F)
```

## Extract the snp data
```ruby
mkdir tempfile
for i in finngen_r7_c3_breast_exallc
do
awk  '{print $5 }' ${i}.txt  > ./tempfile/${i}_SNP_list.txt
done

```


```ruby
for i in $(seq 1 22)  
do 
echo $i
plink 
--bfile ./02_partitioned_LD_score_estimation/1000G_EUR_Phase3_plink/1000G.EUR.QC.$i 
--extract ./tempfile/${j}_SNP_list.txt 
--noweb --make-bed --out ./tempfile/1000G.EUR.QC.${j}_${i}_filtered
done 
```

## Filter LD data
```ruby
for i in $(seq 1 22)  
do 
echo $i
plink 
--bfile ./tempfile/1000G.EUR.QC.${j}_${i}_filtered 
--indep-pairwise 50 5 0.8 
--out  ./tempfile/${j}_${i}_plink_prune_EUR_filtered_LD0.8
done 
```


## Integrate the result of snp
```ruby
cat [${j}]*.prune.in > ${j}_EUR_LD0.8.prune
```

## Merge the snp and gwas data
```ruby
library(readr)
library(dplyr)


gwas<-read_table(paste0("./",i,".txt"))
SNP_prune<- read_table(paste0("./finngen_R7/tempfile/",i,"_EUR_LD0.8.prune"))
SNP_prune<-SNP_prune[!duplicated(unlist(SNP_prune)),]
colnames(SNP_prune)<-"rsid"
#### Left Join using inner_join function 
gwas= gwas %>% inner_join(SNP_prune,by="rsid")
print(nrow(gwas))
write.table(gwas,file= paste0("./",i,"_prune.txt"),,row.names=F,quote=F)
        print(i)

```
