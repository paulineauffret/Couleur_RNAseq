#!/usr/bin/Rscript
#Differential Analysis Flesh RNA-seq 14 nov 2017
#see https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.pdf

#INPUT = HTseq-count output files (txt files, one per sample, raw) + sample sheet with samples metadata (txt file, rows = samples ; col = attributes)
#sessionInfo()

#Loading libraries
library(RSvgDevice)
library(DESeq2)
library(tximport)
library(gplots)
library(ggplot2)
library(heatmap)
library(RColorBrewer)
library(VennDiagram)
library(readODS)
library(assertthat)
library(scales)
library(WGCNA)
library(genefilter)
library(ggrepel)
library(rafalib)

session_info()
###################################################################################################
#
#     Set parameters
#
###################################################################################################
rm(list=ls()) 
setwd()           #set the working directory to the output directory
curdir=           #path to current directory containing sample sheet and count files from step 04
sampleSheet=paste(curdir,"/SampleSheet.txt",sep="")

threshold=10 #min expression value accross samples

#Loading transcriptomes annotations
annot41k <- read.csv("",h=T, sep="\t") ; head(annot41k) ; dim(annot41k)
annotFuc <- read.csv("",h=T, sep="\t") ; head(annotFuc) ; dim(annotFuc)
annotFuc <- annotFuc[,c(1,2)] ; head(annotFuc)
colnames(annotFuc) <- c("transcript_name","fucata_hit")
head(annotFuc) ; dim(annotFuc)

#Source functions file
source("05_running_de_analysis_deseq2_functions.r")
###################################################################################################
#
#     DESeq2 analysis from HTSeq-count data
#
###################################################################################################
#Imports samples infos
samplesInfo <- read.table(sampleSheet, h=T) ; head(samplesInfo) ; dim(samplesInfo)
target_file <- samplesInfo[,c(3,1,5)]
colnames(target_file) <- c("label","files","color") ; head(target_file) ; dim(target_file)

#Creates DEseq object
sampleFiles <- grep("F4_F256_q5_f2.bam_htseq-count.txt",list.files(curdir),value=TRUE)
sampleCondition <- samplesInfo
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles)
sampleTable <- merge(sampleTable,samplesInfo,by.x="sampleName",by.y="Sample_name")
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = curdir, design= ~Color) ; ddsHTSeq


#Trims too low transcripts
ddsHTSeq <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > threshold, ] ; ddsHTSeq

#DEseq object
dds <- DESeq(ddsHTSeq, fitType="parametric")

#Plot dispersion plot
tiff("plotDispersion.tiff")
plotDispEsts(dds)
dev.off()

###################################################################################################
#     PCA 
## PC1 et PC2
rld <- rlog(ddsHTSeq, blind=FALSE)
data <- plotPCA(rld, intgroup=c("Color"), returnData=TRUE)#, ntop=100000)
percentVar <- round(100*attr(data, "percentVar"))
tiff("PCA_PC1_PC2.tiff")
ggplot(data, aes(PC1, PC2, color=Color)) +
  geom_point(size=2.5) +#, fill="black", stroke=0.8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + ggtitle("Principal Component Analysis plot") +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  #scale_shape_manual(values=c(21,24)) +
  scale_color_manual(values=c("white","black","orange")) +
  geom_text(aes(label=colData(ddsHTSeq)$Id),hjust=0.5, vjust=1.8, size=3)
dev.off()

##PC1 et PC3
data <- plotPCA.san(rld, intgroup=c("Color"), returnData=TRUE)#, ntop=100000)
percentVar <- round(100*attr(data, "percentVar"))
tiff("PCA_PC1_PC3.tiff")
ggplot(data, aes(PC1, PC3, color=Color)) +
  geom_point(size=2.5) +#, fill="black", stroke=0.8) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar[3],"% variance")) +
  coord_fixed() + ggtitle("Principal Component Analysis plot") +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  #scale_shape_manual(values=c(21,24)) +
  scale_color_manual(values=c("white","black","orange")) +
  geom_text(aes(label=colData(ddsHTSeq)$Id),hjust=0.5, vjust=1.8, size=3)
dev.off()

###################################################################################################
#     Heatmap
#Selects only most abundant transcripts
nCounts <- counts(dds, normalized=TRUE)
select <- order(rowMeans(nCounts),decreasing=TRUE)[1:75]

#Selects corresponding norm counts
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)[select,]

#Gets the metadata
colnames(log2.norm.counts) <- colData(dds)$Id
df <- as.data.frame(colData(dds)["Color"])
rownames(df) <- colData(dds)$Id

#Heatmap parameters
dist1 <- "euclidean"
clust <- "average"
ann_colors = list(Color = c(Albinos = "white", Black = "black", Orange= "orange"))

## 'Overwrite' default draw_colnames with your own version 
assignInNamespace(x="draw_colnames", value="draw_colnames_45",ns=asNamespace("pheatmap"))

#Draw heatmap
tiff("heatmap.tiff")
pheatmap(log2.norm.counts,
         clustering_distance_cols = dist1, 
         clustering_method = clust, 
         cluster_rows=FALSE,
         cluster_cols=TRUE, 
         annotation_col=df,
         show_rownames=TRUE,
         fontsize_row=5,
         fontsize_col=10,
         fontsize=8,
         annotation_colors = ann_colors, 
         main=paste("Clustered heatmap of 75 most abundant genes\n",dist1," distance with ",clust, " clustering method",sep=""),
         las=1)
dev.off()

###################################################################################################
#     Clustering
#Normalized counts
nt <- normTransform(dds)
log2.norm.counts <- assay(nt)

#Gets the metadata
colnames(log2.norm.counts) <- colData(dds)$Id
df <- as.data.frame(colData(dds)["Color"])
rownames(df) <- colData(dds)$Id

#Heatmap parameters
dist1 <- "euclidean"
clust <- "average"

distance_matrix <- dist(t(log2.norm.counts), method = dist1)
hh <- hclust(distance_matrix, method = clust)
colo <- c(rep("darkorange2",4),rep("grey20", 4),rep("antiquewhite3",4))
tiff("clustering.tiff")
myplclust(hh, labels=rownames(df), lab.col=colo)
dev.off()

###################################################################################################
#     Heatmap sample to sample
distsRL <- dist(t(assay(rld)))
mat<- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds), paste(Id, Color , sep=' : '))
hc <- hclust(distsRL)
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
colo <- c(rep("darkorange2",4),rep("grey20", 4),rep("antiquewhite3",4))
tiff("heatmap_samplebysample.tiff")
heatmap.2(mat, Rowv=as.dendrogram(hc),
          symm=TRUE, trace='none',
          col=rev(hmcol),margin=c(10, 10), colRow=colo, colCol=colo)
dev.off()


###################################################################################################
#     DE analysis
res <- results(dds, pAdjustMethod = "BH", alpha=0.01) 
res

#1-1 comparisons
res_alb_ora <- results(dds,contrast=c("Color","Albinos","Orange"), pAdjustMethod = "BH", alpha=0.01)
res_alb_bla <- results(dds,contrast=c("Color","Albinos","Black"), pAdjustMethod = "BH", alpha=0.01)
res_ora_bla <- results(dds,contrast=c("Color","Orange","Black"), pAdjustMethod = "BH", alpha=0.01)

###################################################################################################
#     MA plots 
tiff("MAplot_res_alb_ora.tiff")
plotMA.DESeqResults(res_alb_ora, main="Albinos vs. Orange", ylim=c(-20,20),alpha=0.01)
dev.off()
tiff("MAplot_res_alb_bla.tiff")
plotMA.DESeqResults(res_alb_bla, main="Albinos vs. Black", ylim=c(-20,20),alpha=0.01)
dev.off()
tiff("MAplot_res_ora_bla.tiff")
plotMA.DESeqResults(res_ora_bla, main="Orange vs. Black", ylim=c(-20,20),alpha=0.01)
dev.off()

#Selects results with p-value <= 0.05 et |log2FC| >= 2
res_alb_ora <-res_alb_ora[which(res_alb_ora$padj<=0.01),] ; res_alb_ora<-res_alb_ora[which(abs(res_alb_ora$log2FoldChange) >= 2),] ; dim(res_alb_ora)
res_alb_bla <-res_alb_bla[which(res_alb_bla$padj<=0.01),] ; res_alb_bla<-res_alb_bla[which(abs(res_alb_bla$log2FoldChange) >= 2),] ; dim(res_alb_bla)
res_ora_bla <-res_ora_bla[which(res_ora_bla$padj<=0.01),] ; res_ora_bla<-res_ora_bla[which(abs(res_ora_bla$log2FoldChange) >= 2),] ; dim(res_ora_bla)

###################################################################################################
#     Venn diagrams

inter_alb_ora_alb_bla <- intersect(rownames(res_alb_ora), rownames(res_alb_bla)) ; length(inter_alb_ora_alb_bla)
tiff("alb_bla_vennD")
draw.pairwise.venn(length(rownames(res_alb_ora)),length(rownames(res_alb_bla)),length(inter_alb_ora_alb_bla),category=c("Albinos vs. Orange","Albinos vs. Black"),
                   fill=c("grey","black"), cat.col=c("grey","black"), cat.cex=0.8, cat.pos=c(3,6), cat.dist=c(0.02,0.02),
                   cat.fontfamily=c("sans","sans"))
dev.off()

tiff("ora_bla_vennD")
inter_alb_ora_ora_bla <- intersect(rownames(res_alb_ora), rownames(res_ora_bla)) ; length(inter_alb_ora_ora_bla)
draw.pairwise.venn(length(rownames(res_alb_ora)),length(rownames(res_ora_bla)),length(inter_alb_ora_ora_bla),category=c("Albinos vs. Orange","Orange vs. Black"),
                   fill=c("grey","orange"), cat.col=c("grey","orange"), cat.cex=0.8, cat.pos=c(3,6), cat.dist=c(0.02,0.02),
                   cat.fontfamily=c("sans","sans"))
dev.off()

tiff("ora_bla.tiff")
inter_alb_bla_ora_bla <- intersect(rownames(res_alb_bla), rownames(res_ora_bla)) ; length(inter_alb_bla_ora_bla)
draw.pairwise.venn(length(rownames(res_alb_bla)),length(rownames(res_ora_bla)),length(inter_alb_bla_ora_bla),category=c("Albinos vs. Black","Orange vs. Black"),
                   fill=c("black","orange"), cat.col=c("black","orange"), cat.cex=0.8, cat.pos=c(3,6), cat.dist=c(0.02,0.02),
                   cat.fontfamily=c("sans","sans"))
dev.off()


temp1 <- intersect(rownames(res_alb_ora), rownames(res_alb_bla)) ; length(temp1)
temp2 <- intersect(rownames(res_alb_ora), rownames(res_ora_bla)) ; length(temp2)
temp3 <- intersect(rownames(res_alb_bla), rownames(res_ora_bla)) ; length(temp3)
temp123 <- intersect(temp1,temp2); length(temp123)
tiff("vennD.tiff")
draw.triple.venn(length(rownames(res_alb_ora)),length(rownames(res_alb_bla)),length(rownames(res_ora_bla)),length(temp1),length(temp3),length(temp2),length(temp123), 
                 category=c(paste("Albinos vs. Orange\nn=",length(rownames(res_alb_ora)),sep=""),paste("Albinos vs. Black\nn=",length(rownames(res_alb_bla)),sep=""),paste("Orange vs. Black\nn=",length(rownames(res_ora_bla)),sep="")),fill=c("grey","black","orange"), 
                 cat.col=c("grey","black","orange"), cat.fontfamily=c("sans","sans","sans"),margin=0.05,cat.dist=c(0.07,0.07,0.05))
dev.off()


###################################################################################################
#     Merge DE genes and annotations
res_alb_ora <- as.data.frame(res_alb_ora)
res_alb_bla <- as.data.frame(res_alb_bla)
res_ora_bla <- as.data.frame(res_ora_bla)

res_alb_ora$tr_name <- rownames(res_alb_ora)
res_alb_bla$tr_name <- rownames(res_alb_bla)
res_ora_bla$tr_name <- rownames(res_ora_bla)

res_alb_ora_annot <- merge(res_alb_ora,annot41k,by.x="tr_name",by.y="Name", all.x=TRUE, all.y=FALSE) ; head(res_alb_ora_annot)
res_alb_ora_annot_fuc <- merge(res_alb_ora_annot,annotFuc,by.x="tr_name",by.y="transcript_name", all.x=TRUE, all.y=FALSE) ; head(res_alb_ora_annot_fuc)

res_alb_bla_annot <- merge(res_alb_bla,annot41k,by.x="tr_name",by.y="Name", all.x=TRUE, all.y=FALSE) ; head(res_alb_bla_annot)
res_alb_bla_annot_fuc <- merge(res_alb_bla_annot,annotFuc,by.x="tr_name",by.y="transcript_name", all.x=TRUE, all.y=FALSE) ; head(res_alb_bla_annot_fuc)

res_ora_bla_annot <- merge(res_ora_bla,annot41k,by.x="tr_name",by.y="Name", all.x=TRUE, all.y=FALSE) ; head(res_ora_bla_annot)
res_ora_bla_annot_fuc <- merge(res_ora_bla_annot,annotFuc,by.x="tr_name",by.y="transcript_name", all.x=TRUE, all.y=FALSE) ; head(res_ora_bla_annot_fuc)

f1 <- res_ora_bla_annot_fuc$Fullname
f2 <- res_alb_bla_annot_fuc$Fullname
f3 <- res_alb_ora_annot_fuc$Fullname

temp1 <- intersect(f3, f2) ; length(temp1)
temp2 <- intersect(f3, f1) ; length(temp2)
temp3 <- intersect(f2, f1) ; length(temp3)
temp123 <- intersect(temp1,temp2); length(temp123)
draw.triple.venn(length(rownames(res_alb_ora)),length(rownames(res_alb_bla)),length(rownames(res_ora_bla)),length(temp1),length(temp3),length(temp2),length(temp123), 
                 category=c(paste("Albinos vs. Orange\nn=",length(rownames(res_alb_ora)),sep=""),paste("Albinos vs. Black\nn=",length(rownames(res_alb_bla)),sep=""),paste("Orange vs. Black\nn=",length(rownames(res_ora_bla)),sep="")),fill=c("grey","black","orange"), 
                 cat.col=c("grey","black","orange"), cat.fontfamily=c("sans","sans","sans"),margin=0.05,cat.dist=c(0.07,0.07,0.05))
dev.off()

###################################################################################################
#     Export annotated DE genes lists
write.table(res_alb_ora_annot,"res_alb_ora_pv0.05_FC2_annot",row.names = FALSE, quote=F, sep="\t")
write.table(res_alb_bla_annot,"res_alb_bla_pv0.05_FC2_annot",row.names = FALSE, quote=F, sep="\t")
write.table(res_ora_bla_annot,"res_ora_bla_pv0.05_FC2_annot",row.names = FALSE, quote=F, sep="\t")

write.table(res_alb_ora_annot_fuc,"res_alb_ora_pv0.05_FC2_annot_fuc",row.names = FALSE, quote=F, sep="\t")
write.table(res_alb_bla_annot_fuc,"res_alb_bla_pv0.05_FC2_annot_fuc",row.names = FALSE, quote=F, sep="\t")
write.table(res_ora_bla_annot_fuc,"res_ora_bla_pv0.05_FC2_annot_fuc",row.names = FALSE, quote=F, sep="\t")

###################################################################################################
#     Common DE genes
res_alb_bla_res_alb_ora <- res_alb_bla_annot_fuc[which(res_alb_bla_annot_fuc$tr_name%in%res_alb_ora_annot_fuc$tr_name),] ; head(res_alb_bla_res_alb_ora) ; dim(res_alb_bla_res_alb_ora)
res_alb_bla_res_ora_bla <- res_alb_bla_annot_fuc[which(res_alb_bla_annot_fuc$tr_name%in%res_ora_bla_annot_fuc$tr_name),] ; head(res_alb_bla_res_ora_bla) ; dim(res_alb_bla_res_ora_bla)
res_alb_ora_res_ora_bla <- res_alb_bla_annot_fuc[which(res_alb_ora_annot_fuc$tr_name%in%res_ora_bla_annot_fuc$tr_name),] ; head(res_alb_ora_res_ora_bla) ; dim(res_alb_ora_res_ora_bla)

res_alb_ora_res_ora_bla_res_alb_bla <- res_alb_bla_res_alb_ora[which(res_alb_bla_res_alb_ora$tr_name%in%res_ora_bla_annot_fuc$tr_name),]  ; head(res_alb_ora_res_ora_bla_res_alb_bla) ; dim(res_alb_ora_res_ora_bla_res_alb_bla)

write.table(res_alb_bla_res_alb_ora,"res_alb_bla_res_alb_ora_pv0.05_FC2_annot_fuc",row.names = FALSE, quote=F, sep="\t")
write.table(res_alb_bla_res_ora_bla,"res_alb_bla_res_ora_bla_pv0.05_FC2_annot_fuc",row.names = FALSE, quote=F, sep="\t")
write.table(res_alb_ora_res_ora_bla,"res_alb_ora_res_ora_bla_pv0.05_FC2_annot_fuc",row.names = FALSE, quote=F, sep="\t")

write.table(res_alb_ora_res_ora_bla_res_alb_bla,"res_alb_ora_res_ora_bla_res_alb_bla_pv0.05_FC2_annot_fuc",row.names = FALSE, quote=F, sep="\t")

