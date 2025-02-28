#' ---
#' title: "Somatic Embryogenesis Project Kallisto Biological QA_"
#' author: "Nicolas Delhomme & Camilla Canovi"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/src/R")
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
samples <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/samples-3.csv")
samples$Batch <- factor(sprintf("B%d",grepl("P7614",samples$ScilifeID)+1))
samples$Stages <- factor(samples$Stages)

#' ## Kallisto
orig_all <- list.files("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/kallisto_LTR-TEs", 
                    recursive = TRUE, 
                    pattern = "abundance.tsv",
                    full.names = TRUE)

#' name them
names(orig_all) <- sub("*_sortmerna_trimmomatic","",
                   sapply(strsplit(orig_all, "/"), .subset, 9))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
orig <- orig_all[intersect(names(orig_all), samples$ScilifeID)]
samples <- samples[match(names(orig),samples$ScilifeID),]

#' Read the expression at the transcript level
tx <- suppressMessages(tximport(files = orig, 
                                type = "kallisto", 
                                txOut = TRUE))
kg <- round(tx$counts)


#' ## Combining samples - Genes and TEs (kg)
samples$ID <- sub("_L00[1,2]", "",
                  samples$ScilifeID)

counts_genes_and_TEs <- do.call(
  cbind,
  lapply(split.data.frame(t(kg),
                          samples$ID),
         colSums))

csamples <- samples[,-1]
csamples <- csamples[match(colnames(counts_genes_and_TEs),csamples$ID),]
write.csv(csamples,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3_genes_and_TEs_R.csv")


#' # QC
#' Check how many genes are never expressed
sel <- rowSums(counts_genes_and_TEs) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts_genes_and_TEs),digits=1),
        sum(sel),
        nrow(counts_genes_and_TEs))

#' ## Coverage
#' The cumulative gene coverage is as expected
plot(density(log10(rowMeans(counts_genes_and_TEs))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by stages
plot.multidensity(lapply(1:ncol(counts_genes_and_TEs),function(k){log10(counts_genes_and_TEs)[,k]}),
                  col=pal[as.integer(csamples$Stages)],
                  legend.x="topright",
                  legend=levels(csamples$Stages),
                  legend.col=pal[1:nlevels(csamples$Stages)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' Colored by Batch
plot.multidensity(lapply(1:ncol(counts_genes_and_TEs),function(k){log10(counts_genes_and_TEs)[,k]}),
                  col=pal[as.integer(csamples$Batch)],
                  legend.x="topright",
                  legend=levels(csamples$Batch),
                  legend.col=pal[1:nlevels(csamples$Batch)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' ## Raw data export
dir.create(file.path("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis","kallisto_LTR-TEs"),showWarnings=FALSE)
# ??Permission denied?? write.csv(counts_genes_and_TEs,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/raw-unnormalised-gene-expression_data-3_genes_and_TEs_R_blindf.csv")

#' # Data normalisation 
#' ## Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate  

#' Create the dds object
rownames(csamples) <- csamples$ID
dds.counts_genes_and_TEs <- DESeqDataSetFromMatrix(
  countData = counts_genes_and_TEs,
  colData = csamples[,c("SubmittedID","Stages","Description","Batch","ID")],
  design = ~Batch+Stages)

#' Save it
# ??Permission?? save(dds.counts_genes_and_TEs,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/DESeqDataSet-3_genes_and_TEs_R_blindf.rda")

#' Check the size factors (i.e. the sequencing library size effect)
dds.counts_genes_and_TEs <- estimateSizeFactors(dds.counts_genes_and_TEs)
sizes.counts_genes_and_TEs <- sizeFactors(dds.counts_genes_and_TEs)
names(sizes.counts_genes_and_TEs) <- colnames(counts_genes_and_TEs)
pander(sizes.counts_genes_and_TEs)
boxplot(sizes.counts_genes_and_TEs, main="Sequencing libraries size factor")

#' ## Variance stabilising transformation (blind = TRUE)
vsd.counts_genes_and_TEs <- varianceStabilizingTransformation(dds.counts_genes_and_TEs, blind=TRUE)
vst.counts_genes_and_TEs <- assay(vsd.counts_genes_and_TEs)
vst.counts_genes_and_TEs <- vst.counts_genes_and_TEs - min(vst.counts_genes_and_TEs)

#' Validate the VST 
meanSdPlot(vst.counts_genes_and_TEs[rowSums(counts_genes_and_TEs)>0,])

#' ## PCA
pc <- prcomp(t(vst.counts_genes_and_TEs))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(csamples$Stages)],
              pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 1st and 2nd by Stage
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(csamples$Stages)],
     main="1st and 2nd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 2nd and 3rd by Stage
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(csamples$Stages)],
     main="2nd and 3rd PC",
     pch=19)
legend("top",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' Make a selection of genes
res.sel_genes_and_TEs <- sapply(1:10,function(i){sum(featureSelect(vst.counts_genes_and_TEs,conditions = csamples$Stages,exp = i,nrep = 2))})

plot(res.sel_genes_and_TEs)

sel_genes_and_TEs <- featureSelect(vst.counts_genes_and_TEs,
                                   conditions = csamples$Stages,
                                   exp = 5,
                                   nrep = 2)

message(sprintf("%s genes are selected",sum(sel_genes_and_TEs)))

#' Heatmap with the selection of TEs
heatmap.2(vst.counts_genes_and_TEs[sel_genes_and_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
)

#' Write it out
# ??Permission ?? write.csv(vst.counts_genes_and_TEs,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs//library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_genes_and_TEs_R.csv")
  
#' ## PCA of only TEs (subsetted data)
pc <- prcomp(t(vst.counts_genes_and_TEs[grepl("^MA_", rownames(vst.counts_genes_and_TEs)) == FALSE, ]))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(csamples$Stages)],
              pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 1st and 2nd by Stages
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(csamples$Stages)],
     main="1st and 2nd PC",
     pch=19)
legend("bottomright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 1st and 2nd by Batch
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(csamples$Batch)],
     main="1st and 2nd PC",
     pch=19)
legend("bottomright",pch=19,
       col=pal[1:nlevels(csamples$Batch)],
       legend=levels(csamples$Batch))
par(mar=mar)

#' ### 2nd and 3rd by Stage
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(csamples$Stages)],
     main="2nd and 3rd PC",
     pch=19)
legend("bottomright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 2nd and 3rd by Batch
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(csamples$Batch)],
     main="2nd and 3rd PC",
     pch=19)
legend("bottomleft",pch=19,
       col=pal[1:nlevels(csamples$Batch)],
       legend=levels(csamples$Batch))
par(mar=mar)


#' # Heatmap  
#' select all TEs
vst.counts_TEs <- vst.counts_genes_and_TEs[grepl("^MA_", rownames(vst.counts_genes_and_TEs)) == FALSE, ]

heatmap.2(vst.counts_TEs, 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
          )

#' Make a selection of TEs  
#' TEs with expression value of 1
sel_TEs <- featureSelect(vst.counts_TEs,
                                   conditions = csamples$Stages,
                                   exp = 1,
                                   nrep = 2)

res.sel_TEs <- sapply(1:10,function(i){sum(featureSelect(vst.counts_TEs,conditions = csamples$Stages,exp = i,nrep = 2))})

plot(res.sel_TEs)

message(sprintf("%s genes are selected",sum(sel_TEs)))

#' Heatmap with the selection of TEs
heatmap.2(vst.counts_TEs[sel_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
)

#' TEs with expression value of 2
sel_TEs <- featureSelect(vst.counts_TEs,
                         conditions = csamples$Stages,
                         exp = 2,
                         nrep = 2)

message(sprintf("%s genes are selected",sum(sel_TEs)))

heatmap.2(vst.counts_TEs[sel_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
)

#' ## Variance stabilising transformation (blind = FALSE)
vsd.counts_genes_and_TEs <- varianceStabilizingTransformation(dds.counts_genes_and_TEs, blind=FALSE)
vst.counts_genes_and_TEs <- assay(vsd.counts_genes_and_TEs)
vst.counts_genes_and_TEs <- vst.counts_genes_and_TEs - min(vst.counts_genes_and_TEs)

#' Validate the VST 
meanSdPlot(vst.counts_genes_and_TEs[rowSums(counts_genes_and_TEs)>0,])

#' ## PCA
pc <- prcomp(t(vst.counts_genes_and_TEs))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(csamples$Stages)],
              pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 1st and 2nd
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(csamples$Stages)],
     main="1st and 2nd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 2nd and 3rd
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(csamples$Stages)],
     main="2nd and 3rd PC",
     pch=19)
legend("top",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' Make a selection of genes
sel_genes_and_TEs <- featureSelect(vst.counts_genes_and_TEs,
                                   conditions = csamples$Stages,
                                   exp = 5,
                                   nrep = 2)

res.sel_genes_and_TEs <- sapply(1:10,function(i){sum(featureSelect(vst.counts_genes_and_TEs,conditions = csamples$Stages,exp = i,nrep = 2))})

plot(res.sel_genes_and_TEs)

message(sprintf("%s genes are selected",sum(sel_genes_and_TEs)))

#' Heatmap with the selection of genes
heatmap.2(vst.counts_genes_and_TEs[sel_genes_and_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
)



#' Write it out
# ??Permission ?? write.csv(vst.counts_genes_and_TEs,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs//library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_genes_and_TEs_R.csv")

# PCA of only TEs (subsetted data)

#' ## PCA
pc <- prcomp(t(vst.counts_genes_and_TEs[grepl("^MA_", rownames(vst.counts_genes_and_TEs)) == FALSE, ]))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(csamples$Stages)],
              pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 1st and 2nd by Stages
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(csamples$Stages)],
     main="1st and 2nd PC",
     pch=19)
legend("bottomright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 1st and 2nd by Batch
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(csamples$Batch)],
     main="1st and 2nd PC",
     pch=19)
legend("bottomright",pch=19,
       col=pal[1:nlevels(csamples$Batch)],
       legend=levels(csamples$Batch))
par(mar=mar)

#' ### 2nd and 3rd by Stage
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(csamples$Stages)],
     main="2nd and 3rd PC",
     pch=19)
legend("bottomright",pch=19,
       col=pal[1:nlevels(csamples$Stages)],
       legend=levels(csamples$Stages))
par(mar=mar)

#' ### 2nd and 3rd by Batch
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(csamples$Batch)],
     main="2nd and 3rd PC",
     pch=19)
legend("bottomleft",pch=19,
       col=pal[1:nlevels(csamples$Batch)],
       legend=levels(csamples$Batch))
par(mar=mar)

#' # Heatmap  
#' select only TEs
vst.counts_TEs <- vst.counts_genes_and_TEs[grepl("^MA_", rownames(vst.counts_genes_and_TEs)) == FALSE, ]

heatmap.2(vst.counts_TEs, 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
)

#' Make a selection of TEs
sel_TEs <- featureSelect(vst.counts_TEs,
                         conditions = csamples$Stages,
                         exp = 1,
                         nrep = 2)

res.sel_TEs <- sapply(1:10,function(i){sum(featureSelect(vst.counts_TEs,conditions = csamples$Stages,exp = i,nrep = 2))})

plot(res.sel_TEs)

message(sprintf("%s genes are selected",sum(sel_TEs)))

#' Heatmap with the selection of TEs
heatmap.2(vst.counts_TEs[sel_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
)

#' ### Hierarchical clustering
# sample.hc <- hclust(d = pearson.dist(t(scale(vst.counts[sel,]))),method = "ward.D")
# plot(sample.hc,labels=paste(csamples$Batch,csamples$Stages,sep="-"))

#' ### Heatmap
# gene.hc <- hclust(d = pearson.dist(t(scale(t(vst.counts[sel,])))),method = "ward.D")
# ord <- order(csamples$Stages)
# pdf(file="analysis/kallisto/heatmap.pdf",width=24,height=24)
# heatmap.2(
# vst.counts[sel,ord][gene.hc$order,],
# labRow = NA,trace = "none",
# scale = "row",Rowv=FALSE,
# dendrogram="none",Colv=FALSE,
# labCol = paste(csamples$Batch,
# csamples$Stages,sep="-")[ord],
# col=hpal)


#' # Differential expression for genes and TEs

#' Load DESeq object
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/DESeqDataSet-3_genes_and_TEs_R_blindf.rda")

#' Differential Expression
dds_genes_and_TEs <- DESeq(dds.counts_genes_and_TEs)


#' Get DE genes and TEs between the developmental stages
alpha=0.01
lfc=0.5
#for all results?
allres_genes_and_TEs <- lapply(2:8,function(i){
  res_genes_and_TEs <- results(dds_genes_and_TEs,c("Stages",as.character(i),as.character(i-1)))
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
heatmap.2(vst.counts_genes_and_TEs[names_DE_TEs, ], 
          scale = "row",
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
)

# Copia
heatmap.2(vst.counts_genes_and_TEs[names_DE_copia, ], 
          scale = "row",
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
)

# Gypsy
heatmap.2(vst.counts_genes_and_TEs[names_DE_gypsy, ], 
          scale = "row",
          labCol = paste0(csamples$Batch, "-", csamples$Stages),
          trace = "none"
)


#' Plot counts of DE TEs  
# calculate mean of biological replicates
mean_counts <- do.call(
  cbind,
  lapply(split.data.frame(t(vst.counts_genes_and_TEs),
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

#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
