#!/usr/bin/env Rscript
library(ggfortify)
library(ggrepel)
library(ggplot2)
library(grid)
library(gridExtra)
pdf('PCA.pdf',width=13,height=6)

cts = read.table("blood.input.tsv",header=TRUE,sep="\t")

#sample <- c("10DG0934", "11DG0060", "11DG0268",
#"11DG0165",	"11DG0840",	"12DG1794",	"14DG2098",	"15DG2154",	"15DG2530",	"16DG0559",
#"16DG0676",	"16DG0790",	"16DG1353",	"18DG0180",	"18DG0295",	"19DG0151",
#"13DG2283",	"14DG2019",	"16DG0144",	"16DG0518",	"16DG0932",	"17DG0349",
#"18DG0348", "18DG0464F","18DG0603F","19DG0230",
#"19DG0152F"," 19DG1391F", "19DG2599F")

sample = colnames(cts)[1:39]

t1 <- t(data.frame(cts,row.names=1))
t2 <- as.data.frame(t1,row.names=F)
t3 <- as.data.frame(cbind(sample,t2))
cts <- data.frame(t3,row.names=1)
end = ncol(cts)-2
df=as.data.frame(lapply(cts[1:end],as.numeric))
#df = log2(df)
#df <- df[1:13619]
#cts[13620:13620]
pca_res <- prcomp(df, scale. = TRUE)
#plot(pca_res)
#ggbiplot(pca_res)+theme_bw()+


#grid.newpage()
#pushViewport(viewport(layout = grid.layout(1, 2)))
#vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)


p1 = autoplot(pca_res,data = cts,colour = 'Group')+theme_bw()+
geom_text_repel(label = cts$Label,size=5)+
labs(tag = "B")+
scale_color_discrete(name = "Group")+
theme(axis.text = element_text(size=15),axis.title = element_text(size=15))+
theme(legend.text=element_text(size=12))



cts2 = read.table("skin.input.tsv",header=TRUE,sep="\t")

sample2 = colnames(cts2)[1:29]

t4 <- t(data.frame(cts2,row.names=1))
t5 <- as.data.frame(t4,row.names=F)
t6 <- as.data.frame(cbind(sample2,t5))
cts2 <- data.frame(t6,row.names=1)
end2 = ncol(cts2)-2
df2=as.data.frame(lapply(cts2[1:end2],as.numeric))

pca_res2 <- prcomp(df2, scale. = TRUE)
p3 = autoplot(pca_res2, x=1,y=2,data = cts2,colour = 'Group')+theme_bw()+
geom_text_repel(label = cts2$Label,size=5)+
labs(tag = "A")+
scale_color_discrete(name = "Group")+
theme(axis.text = element_text(size=15),axis.title = element_text(size=15))+
theme(legend.position='none')
#p4 = autoplot(pca_res, x=3,y=4,scale = 1,data = cts,colour = 'Group')+theme_bw()+
#geom_text_repel(label = cts$Label,size=2.5)+
#labs(tag = "B")+
#scale_color_discrete(name = "Group")+
#theme(legend.position = 'none')+
#theme(axis.text = element_text(size=10),axis.title = element_text(size=10))


#new('ggmultiplot', plots = list(p1, p2),nrow=2,ncol=2)

grid.arrange(p3, p1, ncol = 2, widths = c(5,6))
#print(p1, vp = vplayout(1, 2))
#print(p3, vp = vplayout(1, 1))


#p3=autoplot(pca_res, x=1,y=3,scale = 1,data = cts,colour = 'Group')+theme_bw()+
#biplot(pca_res,scale=0,loadings=FALSE,colour = 'Group')+theme_bw()+
#pairsplot(pca_res, components = getComponents(pca_res, c(1:10)))+theme_bw()+
#geom_text_repel(label = cts$Label)+
#scale_color_discrete(name = "Group")+
#theme(legend.position = c(0.86,0.9),legend.text=element_text(size=12))+
#theme(axis.text = element_text(size=15),axis.title = element_text(size=15))
