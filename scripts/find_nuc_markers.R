args <- commandArgs(trailingOnly = TRUE)
library(tidyverse)
library(DESeq2)


setwd(args[1])
tile_list = read.table("tiles_info",sep=",")
my_dataset <- read.table(paste0("./dge/",tile_list[1,1],".cytoplasmic"), sep="\t")
colnames(my_dataset)<- c("gene",paste0(tile_list[1,1],".cyto"))
my_data <- read.table(paste0("./dge/",tile_list[1,1],".nuclear"), sep="\t")
colnames(my_data)<- c("gene",paste0(tile_list[1,1],".nuc"))
my_dataset <- merge(my_dataset,my_data,all=T,by="gene")


for (mytile in tile_list[-1,1]) {
	my_data <- read.table(paste0("dge/",mytile,".cytoplasmic"), sep="\t")
	colnames(my_data)<- c("gene",paste0(mytile,".cyto"))
	my_dataset <- merge(my_dataset,my_data,all=T,by="gene")
	my_data <- read.table(paste0("dge/",mytile,".nuclear"), sep="\t")
	colnames(my_data)<- c("gene",paste0(mytile,".nuc"))
	my_dataset <- merge(my_dataset,my_data,all=T,by="gene")
	}


metadata <- as.data.frame(rep(c("cyto","nuc"),length(tile_list[,1])))
rownames(metadata) <- colnames(my_dataset[,-1])
colnames(metadata)<- "condition"

rownames(my_dataset) <- my_dataset[,1]
my_dataset<- my_dataset[,-1]
my_dataset[is.na(my_dataset)]<- 0
dds <- DESeqDataSetFromMatrix(countData = my_dataset,colData = metadata,design = ~ condition)
dds$condition <- relevel(dds$condition, ref = "cyto")
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
write.table(resOrdered,"deseq2_results.txt")
res$padj[is.na(res$padj) ] <- 1
nuclear_markers <- rownames(res[res$padj<0.05 & res$log2FoldChange > 2,])
write.table(nuclear_markers,file="nuclear_markers",row.names=F,quote=F,col.names=F)

