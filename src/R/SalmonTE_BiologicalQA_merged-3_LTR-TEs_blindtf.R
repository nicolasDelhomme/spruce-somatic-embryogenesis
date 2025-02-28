#' ---
#' title: "Somatic Embryogenesis Project SalmonTE Biological QA"
#' author: "Katja Stojkovič"
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

#' ## SalmonTE
orig_counts <- read.csv("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmonte/EXPR.csv", row.names = 1)

#' set sample names as in "samples"
colnames(orig_counts) <- sub("*_sortmerna_trimmomatic","", colnames(orig_counts))
colnames(orig_counts) <- sub("^X", "", colnames(orig_counts))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
counts <- orig_counts[intersect(names(orig_counts), samples$ScilifeID)]
samples <- samples[match(names(counts),samples$ScilifeID),]

#' Round the counts
counts <- round(counts)

#' ## Combining samples - Genes and TEs (kg)
samples$ID <- sub("_L00[1,2]", "",
                  samples$ScilifeID)

counts_biorep <- do.call(
  cbind,
  lapply(split.data.frame(t(counts),
                          samples$ID),
         colSums))

samples_biorep <- samples[,-1]
samples_biorep <- samples_biorep[match(colnames(counts_biorep),samples_biorep$ID),]
#write.csv(samples_biorep,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3_genes_and_TEs_R.csv")

#' # QC
#' Check how many TEs are never expressed
sel <- rowSums(counts_biorep) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts_biorep),digits=1),
        sum(sel),
        nrow(counts_biorep))
#' These are:
which(sel)

#' ## Coverage
#' The cumulative gene coverage
plot(density(log10(rowMeans(counts_biorep))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by stages
plot.multidensity(lapply(1:ncol(counts_biorep),function(k){log10(counts_biorep)[,k]}),
                  col=pal[as.integer(samples_biorep$Stages)],
                  legend.x="topright",
                  legend=levels(samples_biorep$Stages),
                  legend.col=pal[1:nlevels(samples_biorep$Stages)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per TE raw counts (log10)")

#' Colored by Batch
plot.multidensity(lapply(1:ncol(counts_biorep),function(k){log10(counts_biorep)[,k]}),
                  col=pal[c(2,5)][as.integer(samples_biorep$Batch)],
                  legend.x="topright",
                  legend=levels(samples_biorep$Batch),
                  legend.col=pal[c(2,5)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' ## Raw data export
dir.create(file.path("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis","salmonTE"),showWarnings=FALSE)
write.csv(counts_biorep,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/salmonTE/raw-unnormalised-TE-expression_data-3.csv")

#' # Data normalisation 
#' ## Preparation
#' For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate  

#' Create the dds object
rownames(samples_biorep) <- samples_biorep$ID
dds.counts_biorep <- DESeqDataSetFromMatrix(
  countData = counts_biorep,
  colData = samples_biorep[,c("SubmittedID","Stages","Description","Batch","ID")],
  design = ~Batch+Stages)

#' Save it
save(dds.counts_biorep,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/salmonTE/DESeqDataSet-3_TEs.rda")

#' Check the size factors (i.e. the "mapping" library size effect)
dds.counts_biorep <- estimateSizeFactors(dds.counts_biorep)
sizes.counts_biorep <- sizeFactors(dds.counts_biorep)
names(sizes.counts_biorep) <- colnames(counts_biorep)
pander(sizes.counts_biorep)
boxplot(sizes.counts_biorep, main="'mapping' libraries size factor")

#' ## Variance stabilising transformation (blind = TRUE)
vsd.counts_biorep <- varianceStabilizingTransformation(dds.counts_biorep, blind=TRUE)
vst.counts_biorep <- assay(vsd.counts_biorep)
vst.counts_biorep <- vst.counts_biorep - min(vst.counts_biorep)

#' Validate the VST 
meanSdPlot(vst.counts_biorep[rowSums(counts_biorep)>0,])

#' ## Rlog transformation (blind = TRUE)
rld.counts_biorep <- rlog(dds.counts_biorep, blind=TRUE)
rlt.counts_biorep <- assay(rld.counts_biorep)
rlt.counts_biorep <- rlt.counts_biorep - min(rlt.counts_biorep)

#' Validate the rlog transformation
meanSdPlot(rlt.counts_biorep[rowSums(counts_biorep)>0,])

#' Results from the vst and rlog method are similar. I will continue with rlog method, as the method performs a bit better then vst
#' and the transformation is a bit more compact.  
#' 

#' ## PCA
pc <- prcomp(t(rlt.counts_biorep))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples_biorep$Stages)],
              pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Stages)],
       legend=levels(samples_biorep$Stages))
par(mar=mar)

#' ### 1st and 2nd by Batch
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples_biorep$Batch)],
     main="1st and 2nd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Batch)],
       legend=levels(samples_biorep$Batch))
par(mar=mar)

#' ### 2nd and 3rd by Batch
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(samples_biorep$Batch)],
     main="2nd and 3rd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Batch)],
       legend=levels(samples_biorep$Batch))
par(mar=mar)


#' ### 1st and 2nd by Stage
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples_biorep$Stages)],
     main="1st and 2nd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Stages)],
       legend=levels(samples_biorep$Stages))
par(mar=mar)

#' ### 2nd and 3rd by Stage
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(samples_biorep$Stages)],
     main="2nd and 3rd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Stages)],
       legend=levels(samples_biorep$Stages))
par(mar=mar)

#' ## Heatmap  
#' Define new palette:
hpal <- brewer.pal(n = 11, name = "RdBu")
#' Heatmap with all the TEs  
heatmap.2(rlt.counts_biorep, 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(samples_biorep$Batch, "-", samples_biorep$Stages),
          trace = "none", 
          col = hpal,
          ColSideColors = pal[as.integer(samples_biorep$Stages)]
)

#' Make a selection of TEs
res.sel_TEs <- sapply(1:10,function(i){
  sum(featureSelect(rlt.counts_biorep,conditions = samples_biorep$Stages,exp = i,nrep = 2))
  })

plot(res.sel_TEs)

sel_TEs <- featureSelect(rlt.counts_biorep,
                                   conditions = samples_biorep$Stages,
                                   exp =8,
                                   nrep = 2)

message(sprintf("%s TEs are selected",sum(sel_TEs)))

#' Heatmap with the selection of TEs
heatmap.2(rlt.counts_biorep[sel_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(samples_biorep$Batch, "-", samples_biorep$Stages),
          trace = "none", 
          col = hpal,
          ColSideColors = pal[as.integer(samples_biorep$Stages)]
)

#' ## Hierarchical clustering
sample.hc <- hclust(d = pearson.dist(t(scale(rlt.counts_biorep))),method = "ward.D")
plot(sample.hc,labels=paste(samples_biorep$Batch,samples_biorep$Stages,sep="-"))

#' Write the normalised data
write.csv(rlt.counts_biorep,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/salmonTE/library-size-corrected_blindTRUE-rlog_TE-expression_data-3.csv")
  

#' ## Rlog transformation (blind = FALSE)
rld.counts_biorep <- rlog(dds.counts_biorep, blind=FALSE)
rlt.counts_biorep <- assay(rld.counts_biorep)
rlt.counts_biorep <- rlt.counts_biorep - min(rlt.counts_biorep)

#' Validate the rlog
meanSdPlot(rlt.counts_biorep[rowSums(counts_biorep)>0,])

#' ## PCA
pc <- prcomp(t(rlt.counts_biorep))

percent <- round(summary(pc)$importance[2,]*100)

#' ### 3 first dimensions
mar=c(5.1,4.1,4.1,2.1)
scatterplot3d(pc$x[,1],
              pc$x[,2],
              pc$x[,3],
              xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
              ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
              zlab=paste("Comp. 3 (",percent[3],"%)",sep=""),
              color=pal[as.integer(samples_biorep$Stages)],
              pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Stages)],
       legend=levels(samples_biorep$Stages))
par(mar=mar)

#' ### 1st and 2nd by Batch
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples_biorep$Batch)],
     main="1st and 2nd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Batch)],
       legend=levels(samples_biorep$Batch))
par(mar=mar)

#' ### 2nd and 3rd by Batch
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(samples_biorep$Batch)],
     main="2nd and 3rd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Batch)],
       legend=levels(samples_biorep$Batch))
par(mar=mar)


#' ### 1st and 2nd by Stage
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples_biorep$Stages)],
     main="1st and 2nd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Stages)],
       legend=levels(samples_biorep$Stages))
par(mar=mar)

#' ### 2nd and 3rd by Stage
plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(samples_biorep$Stages)],
     main="2nd and 3rd PC",
     pch=19)
legend("topright",pch=19,
       col=pal[1:nlevels(samples_biorep$Stages)],
       legend=levels(samples_biorep$Stages))
par(mar=mar)

#' ## Heatmap  
#' 
#' Heatmap with all the TEs  
heatmap.2(rlt.counts_biorep, 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(samples_biorep$Batch, "-", samples_biorep$Stages),
          trace = "none", 
          col = hpal,
          ColSideColors = pal[as.integer(samples_biorep$Stages)]
)

#' Make a selection of TEs
res.sel_TEs <- sapply(1:10,function(i){
  sum(featureSelect(rlt.counts_biorep,conditions = samples_biorep$Stages,exp = i,nrep = 2))
})

plot(res.sel_TEs)

sel_TEs <- featureSelect(rlt.counts_biorep,
                         conditions = samples_biorep$Stages,
                         exp =8,
                         nrep = 2)

message(sprintf("%s TEs are selected",sum(sel_TEs)))

#' Heatmap with the selection of TEs
heatmap.2(rlt.counts_biorep[sel_TEs, ], 
          scale = "row", 
          labRow = NULL, 
          labCol = paste0(samples_biorep$Batch, "-", samples_biorep$Stages),
          trace = "none", 
          col = hpal,
          ColSideColors = pal[as.integer(samples_biorep$Stages)]
)

#' Write the normalised data
write.csv(rlt.counts_biorep,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/salmonTE/library-size-corrected_blindFALSE-rlog_TE-expression_data-3.csv")

#' ## Hierarchical clustering
sample.hc <- hclust(d = pearson.dist(t(scale(rlt.counts_biorep))),method = "ward.D")
plot(sample.hc,labels=paste(samples_biorep$Batch,samples_biorep$Stages,sep="-"))


#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
