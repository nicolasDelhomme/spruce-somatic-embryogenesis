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
samples <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/samples-3.csv")
samples$Batch <- factor(sprintf("B%d",grepl("P7614",samples$ScilifeID)+1))
samples$Stages <- factor(samples$Stages)
levels(samples$Stages) <- c("S1","S2","S3","S4","S5","S6","S7","S8")

#' ## Kallisto
orig_all <- list.files("kallisto", 
                    recursive = TRUE, 
                    pattern = "abundance.tsv",
                    full.names = TRUE)

#' name them
names(orig_all) <- sub("*_sortmerna_trimmomatic","",
                   sapply(strsplit(orig_all, "/"), .subset, 2))

#' Reorder the sample data.frame according to the way
#' the results were read in and accounting for tech reps
orig <- orig_all[intersect(names(orig_all), samples$ScilifeID)]
samples <- samples[match(names(orig),samples$ScilifeID),]

#' Read the expression at the transcript level
tx <- suppressMessages(tximport(files = orig, 
                                type = "kallisto", 
                                txOut = TRUE))
kg <- round(tx$counts)

#' # Combining samples
samples$ID <- sub("_L00[1,2]", "",
                  samples$ScilifeID)

counts <- do.call(
  cbind,
  lapply(split.data.frame(t(kg),
                          samples$ID),
         colSums))

csamples <- samples[,-1]
csamples <- csamples[match(colnames(counts),csamples$ID),]
write.csv(csamples,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3.csv")

#' # QC
#' Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ## Coverage
#' The cumulative gene coverage is as expected
plot(density(log10(rowMeans(counts))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by stages
plot.multidensity(lapply(1:ncol(counts),function(k){log10(counts)[,k]}),
                  col=pal[as.integer(csamples$Stages)],
                  legend.x="topright",
                  legend=levels(csamples$Stages),
                  legend.col=pal[1:nlevels(csamples$Stages)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' Colored by Batch
plot.multidensity(lapply(1:ncol(counts),function(k){log10(counts)[,k]}),
                  col=pal[as.integer(csamples$Batch)],
                  legend.x="topright",
                  legend=levels(csamples$Batch),
                  legend.col=pal[1:nlevels(csamples$Batch)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' ## Raw data export
dir.create(file.path("analysis","kallisto"),showWarnings=FALSE)
write.csv(counts,file="analysis/kallisto/raw-unormalised-gene-expression_data-3_blindf.csv")

#' # Data normalisation 
#' ## Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' ## # Setup
#' Create the dds object
rownames(csamples) <- csamples$ID
dds.counts <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = csamples[,c("SubmittedID","Stages","Description","Batch","ID")],
  design = ~Batch+Stages)

#' Save it
save(dds.counts,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/DESeqDataSet-3_blindf.rda")

#' Check the size factors (i.e. the sequencing library size effect)
dds.counts <- estimateSizeFactors(dds.counts)
sizes.counts <- sizeFactors(dds.counts)
names(sizes.counts) <- colnames(counts)
pander(sizes.counts)
boxplot(sizes.counts, main="Sequencing libraries size factor")

#' ## Variance stabilising transformation
#' ### blind
vsd.counts <- varianceStabilizingTransformation(dds.counts, blind=FALSE)
vst.counts <- assay(vsd.counts)
vst.counts <- vst.counts - min(vst.counts)

#' Validate the VST 
meanSdPlot(vst.counts[rowSums(counts)>0,])

#' Write it out
write.csv(vst.counts,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3.csv")

#' Make a selection of genes
rangeFeatureSelect(as.matrix(vst.counts), conditions = csamples$Stages, nrep = 2)
rangeFeatureSelect(as.matrix(vst.counts), conditions = csamples$Stages, nrep = 3)

#check with 2replicates
#sel <- featureSelect(as.matrix(vst.counts),conditions = csamples$Stages,exp = 2,nrep = 2)
#message(sprintf("%s genes are selected for the heatmap",sum(sel)))

#selection with 3replicates
sel <- featureSelect(as.matrix(vst.counts),conditions = csamples$Stages,exp = 2,nrep = 3)
message(sprintf("%s genes are selected for the heatmap",sum(sel)))

csamples <- csamples[match(csamples$ID, colnames(vst.counts)), ]


#' Rename genes and samples
rownames(vst.counts[sel,]) <- sub("\\.1$","",rownames(vst.counts[sel,]))

colData(dds.counts)$SeidrID <- paste(colData(dds.counts)$Batch,colData(dds.counts)$Stages,sep="_")

colnames(vst.counts) <- colData(dds.counts)$SeidrID

ord <- order(colData(dds.counts)$Stages,colData(dds.counts)$Batch)

#' Create the anova table
anv.table <- data.frame(Experiment=1,Perturbations=NA,PerturbationLevels=NA,
                        Treatment=NA,DeletedGenes=NA,OverexpressedGenes=NA,
                        Time=colData(dds.counts)$Stages[ord],
                        Repeat=unlist(lapply(lapply(table(colData(dds.counts)$Stages[ord]),":",1),rev)))

write.table(anv.table,quote=FALSE,sep="\t",row.names=FALSE,file="analysis/vala/data/NetworkE_chip_features.tsv")

#' Write out a selection of genes
write.table(vst.counts[sel,ord],file="analysis/vala/data/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection.tsv", sep = "\t", quote = FALSE)
vst.selection3 <- read.table("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/vala/data/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection.tsv",
                    row.names=1)
#' # QC on the normalised data
#' 
#' ## PCA
pc <- prcomp(t(vst.counts))
  
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

#' ## Heatmap
#' ### Feature selection
sel <- featureSelect(vst.counts,conditions = csamples$Stages,exp = 5,nrep = 2)
message(sprintf("%s genes are selected for the heatmap",sum(sel)))

#' ### Hierarchical clustering
sample.hc <- hclust(d = pearson.dist(t(scale(vst.counts[sel,]))),method = "ward.D")
plot(sample.hc,labels=paste(csamples$Batch,csamples$Stages,sep="-"))

#' ### Heatmap
gene.hc <- hclust(d = pearson.dist(t(scale(t(vst.counts[sel,])))),method = "ward.D")
ord <- order(csamples$Stages)
pdf(file="analysis/kallisto/heatmap.pdf",width=24,height=24)
heatmap.2(
  vst.counts[sel,ord][gene.hc$order,],
  labRow = NA,trace = "none",
  scale = "row",Rowv=FALSE,
  dendrogram="none",Colv=FALSE,
  labCol = paste(csamples$Batch,
                 csamples$Stages,sep="-")[ord],
  col=hpal)



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
