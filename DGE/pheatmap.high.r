#!/usr/bin/env Rscript

library(pheatmap)
library(grid)

pdf("FigureS3.pdf")


data = read.table("zscore.skin.high.tsv",sep="\t",header=TRUE,row.names = 1,,check.names=F)
bk<-c(seq(min(data),0,length.out = 50),seq(0.01,max(data),length.out = 50))
annotation <- data.frame(Group = factor(1:29 <=26, labels = c("SMG8/9-", "SMG+")))
rownames(annotation) <- colnames(data) 

pheatmap(data,cluster_rows = FALSE,cluster_cols = FALSE,cellwidth = 15,
		 border=FALSE,
		 show_rownames=F,
		 color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
		 breaks = bk,
		 annotation = annotation,annotation_legend = FALSE)
grid.text(levels(annotation$Group), x=c(0.86,0.4),y=c(0.976,0.976), gp=gpar(fontsize=10))

data = read.table("zscore.blood.high.tsv",sep="\t",header=TRUE,row.names = 1,,check.names=F)
bk<-c(seq(min(data),0,length.out = 50),seq(0.01,max(data),length.out = 50))
annotation <- data.frame(Group = factor(1:39 <=34, labels = c("SMG8/9-", "SMG+")))
rownames(annotation) <- colnames(data) 

pheatmap(data,cluster_rows = FALSE,cluster_cols = FALSE,cellwidth = 12,
		 border=FALSE,
		 show_rownames=F,
		 color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),colorRampPalette(colors = c("white","red"))(length(bk)/2)),
		 breaks = bk,
		 annotation = annotation,annotation_legend = FALSE)
grid.text(levels(annotation$Group), x=c(0.87,0.45),y=c(0.976,0.976), gp=gpar(fontsize=10))
