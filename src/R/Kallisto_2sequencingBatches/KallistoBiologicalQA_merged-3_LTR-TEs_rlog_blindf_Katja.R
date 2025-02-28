#' ---
#' title: "Somatic Embryogenesis Project - Biological QA (rlog)"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("~/Git/UPSCb/projects/spruce-somatic-embryogenesis")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="~/Git/UPSCb/projects/spruce-somatic-embryogenesis")
#' ```


#' Load libraries
#'```{r drop tcltk, echo=FALSE}
#' options(gsubfn.engine = "R")
#'```
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(hyperSpec))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(scatterplot3d))
suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(Mfuzz))

#' Source some helper functions
source("~/Git/UPSCb/src/R/plot.multidensity.R")
source("~/Git/UPSCb/src/R/featureSelection.R")

#' Create a palette
pal <- brewer.pal(10,"Paired")
pal2 <- c(brewer.pal(12, "Paired"), "#000000")
hpal <- colorRampPalette(c("blue","white","red"))(100)

brewer.pal(9, "Greys")
#' Register the default plot margin
mar <- par("mar")


#' # Raw data
#' ## Loading
#' Read the sample information
samples <- read.csv("doc/technical_samples-3_genes_and_TEs.csv")

#' Read the rlog normalised counts
rlog.counts_genes_and_TEs <- read.csv("library-size-corrected_blindFALSE-rlog_gene-expression_data-3_genes_and_TEs.csv", row.names = 1)
colnames(rlog.counts_genes_and_TEs) <- sub("^X", "", colnames(rlog.counts_genes_and_TEs))
rlog.counts_genes_and_TEs <- as.matrix(rlog.counts_genes_and_TEs)

#' Match the names
samples <- samples[match(samples$ID, colnames(rlog.counts_genes_and_TEs)), ]
samples$Stages <- as.factor(samples$Stages)

#' # Data normalisation 
#' ## Preparation
#'  For visualization, the data was submitted to a rlog
#' transformation using DESeq2. The dispersion was not estimated independently
#' of the sample tissue and replicate (blind = FALSE).  


#' ## PCA
pc <- prcomp(t(rlog.counts_genes_and_TEs))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Stages)],
              pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples$Stages)],
       legend=levels(samples$Stages))
par(mar=mar)

#' ### 1st and 2nd by Stage
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples$Stages)],
     main="1st and 2nd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples$Stages)],
       legend=levels(samples$Stages))
par(mar=mar)

#' ### 2nd and 3rd by Stage
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(samples$Stages)],
     main="2nd and 3rd PC",
     pch=19)
legend("top",pch=19,
       col=pal[1:nlevels(samples$Stages)],
       legend=levels(samples$Stages))
par(mar=mar)

#' Make a selection of genes
res.sel_genes_and_TEs <- sapply(1:10,function(i){sum(featureSelect(rlog.counts_genes_and_TEs, conditions = samples$Stages, exp = i,nrep = 2))})

plot(res.sel_genes_and_TEs)

sel_genes_and_TEs <- featureSelect(rlog.counts_genes_and_TEs,
                                   conditions = samples$Stages,
                                   exp = 5,
                                   nrep = 2)

message(sprintf("%s genes are selected",sum(sel_genes_and_TEs)))

#' Heatmap with the selection of TEs
heatmap.2(rlog.counts_genes_and_TEs[sel_genes_and_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(samples$Batch, "-", samples$Stages),
          trace = "none"
)

#' Write it out
write.csv(rlog.counts_genes_and_TEs,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/library-size-corrected_blindFALSE-rlog_gene-expression_data-3_genes_and_TEs_R.csv")
  
#' ## PCA of only TEs (subsetted data)
pc <- prcomp(t(rlog.counts_genes_and_TEs[grepl("^MA_", rownames(rlog.counts_genes_and_TEs)) == FALSE, ]))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples$Stages)],
              pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples$Stages)],
       legend=levels(samples$Stages))
par(mar=mar)

#' ### 1st and 2nd by Stages
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples$Stages)],
     main="1st and 2nd PC",
     pch=19)
legend("bottomright",pch=19,
       col=pal[1:nlevels(samples$Stages)],
       legend=levels(samples$Stages))
par(mar=mar)

#' ### 1st and 2nd by Batch
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples$Batch)],
     main="1st and 2nd PC",
     pch=19)
legend("bottomright",pch=19,
       col=pal[1:nlevels(samples$Batch)],
       legend=levels(samples$Batch))
par(mar=mar)

#' ### 2nd and 3rd by Stage
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(samples$Stages)],
     main="2nd and 3rd PC",
     pch=19)
legend("bottomright",pch=19,
       col=pal[1:nlevels(samples$Stages)],
       legend=levels(samples$Stages))
par(mar=mar)

#' ### 2nd and 3rd by Batch
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(samples$Batch)],
     main="2nd and 3rd PC",
     pch=19)
legend("bottomleft",pch=19,
       col=pal[1:nlevels(samples$Batch)],
       legend=levels(samples$Batch))
par(mar=mar)


#' # Heatmap  
#' select all TEs
rlog.counts_TEs <- rlog.counts_genes_and_TEs[grepl("^MA_", rownames(rlog.counts_genes_and_TEs)) == FALSE, ]

heatmap.2(rlog.counts_TEs, 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(samples$Batch, "-", samples$Stages),
          trace = "none"
          )

#' Make a selection of TEs  
#' TEs with expression value of 1
sel_TEs <- featureSelect(rlog.counts_TEs,
                                   conditions = samples$Stages,
                                   exp = 1,
                                   nrep = 2)

res.sel_TEs <- sapply(1:10,function(i){sum(featureSelect(rlog.counts_TEs,conditions = samples$Stages,exp = i,nrep = 2))})

plot(res.sel_TEs)

message(sprintf("%s genes are selected",sum(sel_TEs)))

#' Heatmap with the selection of TEs
heatmap.2(rlog.counts_TEs[sel_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(samples$Batch, "-", samples$Stages),
          trace = "none"
)

#' TEs with expression value of 2
sel_TEs <- featureSelect(rlog.counts_TEs,
                         conditions = samples$Stages,
                         exp = 2,
                         nrep = 2)

message(sprintf("%s genes are selected",sum(sel_TEs)))

heatmap.2(rlog.counts_TEs[sel_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(samples$Batch, "-", samples$Stages),
          trace = "none"
)

#' ### Hierarchical clustering
# sample.hc <- hclust(d = pearson.dist(t(scale(rlog.counts[sel,]))),method = "ward.D")
# plot(sample.hc,labels=paste(samples$Batch,samples$Stages,sep="-"))

#' ### Heatmap
# gene.hc <- hclust(d = pearson.dist(t(scale(t(rlog.counts[sel,])))),method = "ward.D")
# ord <- order(samples$Stages)
# pdf(file="analysis/kallisto/heatmap.pdf",width=24,height=24)
# heatmap.2(
# rlog.counts[sel,ord][gene.hc$order,],
# labRow = NA,trace = "none",
# scale = "row",Rowv=FALSE,
# dendrogram="none",Colv=FALSE,
# labCol = paste(samples$Batch,
# samples$Stages,sep="-")[ord],
# col=hpal)


#' # Differential expression for genes and TEs

#' Load DESeq object
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/DESeqDataSet-3_genes_and_TEs_R_blindf.rda")

#' Differential Expression
dds.genes_and_TEs <- DESeq(dds.counts_genes_and_TEs)


#' Get DE genes and TEs between the developmental stages
alpha=0.01
lfc=0.5
#for all results?
allres_genes_and_TEs <- lapply(2:8,function(i){
  res_genes_and_TEs <- results(dds.genes_and_TEs,c("Stages",as.character(i),as.character(i-1)))
  return(rownames(res_genes_and_TEs[res_genes_and_TEs$padj <= alpha & abs(res_genes_and_TEs$log2FoldChange) >= lfc & ! is.na(res_genes_and_TEs$padj) ,]))
})

#' Get names of DE TEs
allres_TEs <- lapply(allres_genes_and_TEs, function(x) {
  x[grepl("^MA_", x) == FALSE]
})

#' Plot number of DE genes between the stages
nr_DE_TEs <- sapply(allres_TEs, function(x) length(x))
nr_DE_copia <- sapply(allres_TEs, function(x) {
  length(x[grepl("^C", x)])
})
nr_DE_gypsy <- sapply(allres_TEs, function(x) {
  length(x[grepl("^G", x)])
})

barplot(nr_DE_TEs, names.arg = c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8"), col = "darkblue", ylim = c(0,50))
barplot2(rbind(nr_DE_copia, nr_DE_gypsy), 
         beside = TRUE, 
         names.arg = c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8"), 
         ylim = c(0,30), 
         col = c("lightgreen", "lightblue"),
         legend.text = c("Copia LTR-TEs", "Gypsy LTR-TEs"))

#' Make a heatplot with all the DE TEs
names_DE_TEs <- unique(unlist(allres_TEs))
names_DE_copia <- unique(unlist(sapply(allres_TEs, function(x) {
  x[grepl("^C", x)]
})))
names_DE_gypsy <- unique(unlist(sapply(allres_TEs, function(x) {
  x[grepl("^G", x)]
})))

#' Heatmap
heatmap.2(rlog.counts_genes_and_TEs[names_DE_TEs, ], 
          scale = "row",
          labCol = paste0(samples$Batch, "-", samples$Stages),
          trace = "none"
)

# Copia
heatmap.2(rlog.counts_genes_and_TEs[names_DE_copia, ], 
          scale = "row",
          labCol = paste0(samples$Batch, "-", samples$Stages),
          trace = "none"
)

# Gypsy
heatmap.2(rlog.counts_genes_and_TEs[names_DE_gypsy, ], 
          scale = "row",
          labCol = paste0(samples$Batch, "-", samples$Stages),
          trace = "none"
)


#' Plot counts of DE TEs  
# calculate mean of biological replicates
mean_counts <- do.call(
  cbind,
  lapply(split.data.frame(t(rlog.counts_genes_and_TEs),
                          samples$Stages),
         colSums))

# plot gypsy and copia elements
matplot(t(mean_counts[names_DE_copia, ]), type = "l", xlab = "Stages", ylab = "mean normalised counts")
matplot(t(mean_counts[names_DE_gypsy, ]), type = "l", xlab = "Stages", ylab = "mean normalised counts")

#' How many TEs is DE in the first and second big change in the development?
commonTEs_2and6 <- intersect(allres_TEs[[2]], allres_TEs[[6]])
length(commonTEs_2and6)

# plot TEs that are DE expressed in stages 2-3 and 6-7
matplot(t(mean_counts[commonTEs_2and6, ]),
        type = "l",
        xlab = "Stages", 
        ylab = "mean normalised counts",
        col = pal2[1:nlevels(factor(commonTEs_2and6))],
        lwd = 2
        )
plot.new()
legend("center", commonTEs_2and6, col = pal2[1:nlevels(factor(commonTEs_2and6))], lwd = 2, bty = "n")

heatmap.2(rlog.counts_genes_and_TEs[commonTEs_2and6, ], 
          scale = "row",
          labCol = paste0(samples$Batch, "-", samples$Stages),
          trace = "none"
)

#' DE data of TEs  
#' 3 vs. 2
res_TEs_3vs2 <- results(dds.genes_and_TEs, c("Stages", "3", "2"), alpha = alpha)
res_TEs_3vs2 <- res_TEs_3vs2[res_TEs_3vs2$padj <= alpha & abs(res_TEs_3vs2$log2FoldChange) >= lfc & ! is.na(res_TEs_3vs2$padj), ]
res_TEs_3vs2 <- res_TEs_3vs2[grepl("^MA", rownames(res_TEs_3vs2)) == FALSE, ]
res_TEs_3vs2 <- res_TEs_3vs2[order(res_TEs_3vs2$log2FoldChange), ]

#' 7 vs. 6
res_TEs_7vs6 <- results(dds.genes_and_TEs, c("Stages", "7", "6"), alpha = alpha)
res_TEs_7vs6 <- res_TEs_7vs6[res_TEs_7vs6$padj <= alpha & abs(res_TEs_7vs6$log2FoldChange) >= lfc & ! is.na(res_TEs_7vs6$padj), ]
res_TEs_7vs6 <- res_TEs_7vs6[grepl("^MA", rownames(res_TEs_7vs6)) == FALSE, ]
res_TEs_7vs6 <- res_TEs_7vs6[order(res_TEs_7vs6$log2FoldChange), ]


#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
