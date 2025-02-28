#' ---
#' title: "Somatic Embryogenesis Project Kallisto Biological QA"
#' author: "Nicolas Delhomme & Iryna Shutava"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")
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
hpal <- colorRampPalette(c("blue","white","red"))(100)

#' Register the default plot margin
mar <- par("mar")

#' # Raw data
#' ## Loading
#' Read the sample information
samples <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/samples.csv")
samples$Batch <- factor(sprintf("B%d",grepl("P7614",samples$ScilifeID)+1))
samples$Stages <- factor(samples$Stages)

#' ## Kallisto
orig <- list.files("kallisto", 
                    recursive = TRUE, 
                    pattern = "abundance.tsv",
                    full.names = TRUE)

#' name them
names(orig) <- sub("*_sortmerna_trimmomatic","",
                   sapply(strsplit(orig, "/"), .subset, 2))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
samples <- samples[match(names(orig),samples$ScilifeID),]

#' Read the expression at the transcript level
tx <- suppressMessages(tximport(files = orig, 
                                type = "kallisto", 
                                txOut = TRUE))
kg <- round(tx$counts)

#' # QC
#' Check how many genes are never expressed
sel <- rowSums(kg) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(kg),digits=1),
        sum(sel),
        nrow(kg))

#' ## Coverage
#' The cumulative gene coverage is as expected
plot(density(log10(rowMeans(kg))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by stages
plot.multidensity(lapply(1:ncol(kg),function(k){log10(kg)[,k]}),
                  col=pal[as.integer(samples$Stages)],
                  legend.x="topright",
                  legend=levels(samples$Stages),
                  legend.col=pal[1:nlevels(samples$Stages)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' Colored by Batch
plot.multidensity(lapply(1:ncol(kg),function(k){log10(kg)[,k]}),
                  col=pal[as.integer(samples$Batch)],
                  legend.x="topright",
                  legend=levels(samples$Batch),
                  legend.col=pal[1:nlevels(samples$Batch)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' ## Raw data export
dir.create(file.path("analysis","kallisto"),showWarnings=FALSE)
write.csv(kg,file="analysis/kallisto/raw-unormalised-gene-expression_data.csv")

#' # Data normalisation 
#' ## Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' ## # Setup
#' Create the dds object
dds.kg <- DESeqDataSetFromMatrix(
  countData = kg,
  colData = data.frame(batch=samples$Batch,stage=samples$Stages),
  design = ~ batch + stage)

#' Save it
save(dds.kg,file="analysis/kallisto/DESeqDataSet.rda")

#' Check the size factors (i.e. the sequencing library size effect)
dds.kg <- estimateSizeFactors(dds.kg)
sizes.kg <- sizeFactors(dds.kg)
names(sizes.kg) <- colnames(kg)
pander(sizes.kg)
boxplot(sizes.kg, main="Sequencing libraries size factor")

#' ## Variance stabilising transformation
#' ### blind
vsd.kg <- varianceStabilizingTransformation(dds.kg, blind=TRUE)
vst.kg <- assay(vsd.kg)
vst.kg <- vst.kg - min(vst.kg)

#' Validate the VST 
meanSdPlot(vst.kg[rowSums(kg)>0,])

#' Write it out
write.csv(vst.kg,file="analysis/kallisto/library-size-corrected_blind-variance-stabilised_gene-expression_data.csv")

#' # QC on the normalised data
#' 
#' ## PCA
pc <- prcomp(t(vst.kg))
  
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

#' ### 1st and 2nd
plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(samples$Stages)],
     main="1st and 2nd PC",
     pch=19)
legend("top",pch=19,
       col=pal[1:nlevels(samples$Stages)],
       legend=levels(samples$Stages))

#' ### 2nd and 3rd
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

#' ## Heatmap
#' ### Feature selection
sel <- featureSelect(vst.kg,conditions = samples$Stages,exp = 5,nrep = 2)
message(sprintf("%s genes are selected for the heatmap",sum(sel)))

#' ### Hierarchical clustering
sample.hc <- hclust(d = pearson.dist(t(scale(vst.kg[sel,]))),method = "ward.D")
plot(sample.hc,labels=paste(samples$Batch,samples$Stages,sep="-"))

#' ### Heatmap
gene.hc <- hclust(d = pearson.dist(t(scale(t(vst.kg[sel,])))),method = "ward.D")
ord <- order(samples$Stages)
pdf(file="analysis/kallisto/heatmap.pdf",width=24,height=24)
heatmap.2(
  vst.kg[sel,ord][gene.hc$order,],
  labRow = NA,trace = "none",
  scale = "row",Rowv=FALSE,
  dendrogram="none",Colv=FALSE,
  labCol = paste(samples$Batch,
                 samples$Stages,sep="-")[ord],
  col=hpal)
dev.off()

#' # Combining samples
#' 
#' ## VST
#' 
#' ## QC
#' 
#' ### PCA
#' 
#' ### Heatmap
#' 
#' # Conclusion
#' 
#' WRITE ME UP
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
