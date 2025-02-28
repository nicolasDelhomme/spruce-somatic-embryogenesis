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

#' # Genes and TEs (kg)
#' # Combining samples
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

#' ## # Setup
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

#' ## Variance stabilising transformation
#' ### blind
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


#' Write it out
# ??Permission ?? write.csv(vst.counts_genes_and_TEs,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs//library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_genes_and_TEs_R.csv")

#' Make a selection of genes
sel_genes_and_TEs <- featureSelect(vst.counts_genes_and_TEs,conditions = csamples$Stages,exp = 2,nrep = 2)

res.sel_genes_and_TEs <- sapply(1:10,function(i){sum(featureSelect(vst.counts_genes_and_TEs,conditions = csamples$Stages,exp = i,nrep = 2))})

plot(res.sel_genes_and_TEs)

message(sprintf("%s genes are selected for the heatmap",sum(sel_genes_and_TEs)))

#' Rename genes and samples
rownames(vst.counts_genes_and_TEs[sel_genes_and_TEs,]) <- sub("\\.1$","",rownames(vst.counts_genes_and_TEs[sel_genes_and_TEs,]))

colData(dds.counts_genes_and_TEs)$SeidrID <- paste(colData(dds.counts_genes_and_TEs)$Batch,colData(dds.counts_genes_and_TEs)$Stages,sep="_")

colnames(vst.counts_genes_and_TEs) <- colData(dds.counts_genes_and_TEs)$SeidrID

ord_genes_and_TEs <- order(colData(dds.counts_genes_and_TEs)$Stages,colData(dds.counts_genes_and_TEs)$Batch)

#' Create the anova table
anv.table <- data.frame(Experiment=1,Perturbations=NA,PerturbationLevels=NA,
                        Treatment=NA,DeletedGenes=NA,OverexpressedGenes=NA,
                        Time=colData(dds.counts_genes_and_TEs)$Stages[ord_genes_and_TEs],
                        Repeat=unlist(lapply(lapply(table(colData(dds.counts_genes_and_TEs)$Stages[ord_genes_and_TEs]),":",1),rev)))

# ??Permission?? write.table(anv.table,quote=FALSE,sep="\t",row.names=FALSE,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs//NetworkE_chip_features_genes_and_TEs_R.tsv")

#' Write out a selection of genes
# ??Permission?? write.table(vst.counts_genes_and_TEs[sel_genes_and_TEs,ord_genes_and_TEs],file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs//library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection_genes_and_TEs_R.tsv", sep = "\t", quote = FALSE)



names_LTR <- rownames(kg)[grepl("^MA_", rownames(kg)) == FALSE]
kg_LTR <- kg[names_LTR, ]
#kg_LTR <- kg_LTR[,8:55]
#' # TEs (kg_LTR)



#' # Combining samples
samples$ID <- sub("_L00[1,2]", "",
                  samples$ScilifeID)

#samples_B2 <- samples[samples$Batch == "B2", ]
  
  
counts_TEs <- do.call(
  cbind,
  lapply(split.data.frame(t(kg_LTR),
                          samples$ID),
         colSums))

csamples <- samples[,-1]
csamples <- csamples[match(colnames(counts_TEs),csamples$ID),]
# ??Permission?? write.csv(csamples,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3_TEs_R.csv")



#' # QC
#' Check how many genes are never expressed
sel <- rowSums(counts_TEs) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts_TEs),digits=1),
        sum(sel),
        nrow(counts_TEs))

#' ## Coverage
#' The cumulative gene coverage is as expected
plot(density(log10(rowMeans(counts_TEs))),col=pal[1],
     main="gene mean raw counts distribution",
     xlab="mean raw counts (log10)")

#' The same is done for the individual
#' samples colored by stages
plot.multidensity(lapply(1:ncol(counts_TEs),function(k){log10(counts_TEs)[,k]}),
                  col=pal[as.integer(csamples$Stages)],
                  legend.x="topright",
                  legend=levels(csamples$Stages),
                  legend.col=pal[1:nlevels(csamples$Stages)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' Colored by Batch
plot.multidensity(lapply(1:ncol(counts_TEs),function(k){log10(counts_TEs)[,k]}),
                  col=pal[as.integer(csamples$Batch)],
                  legend.x="topright",
                  legend=levels(csamples$Batch),
                  legend.col=pal[1:nlevels(csamples$Batch)],
                  legend.lwd=2,
                  main="sample raw counts distribution",
                  xlab="per gene raw counts (log10)")

#' ## Raw data export
dir.create(file.path("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs"),showWarnings=FALSE)
# ??Permission?? write.csv(counts_TEs,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/raw-unnormalised-gene-expression_data-3_TEs_R_blindf.csv")

#' # Data normalisation 
#' ## Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' ## # Setup
#' Create the dds object
rownames(csamples) <- csamples$ID
dds.counts_TEs <- DESeqDataSetFromMatrix(
  countData = counts_TEs,
  colData = csamples[,c("SubmittedID","Stages","Description","Batch","ID")],
  design = ~Batch+Stages)

#' Save it
# ??Permission?? save(dds.counts_TEs,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs//DESeqDataSet-3_TEs_R_blindf.rda")

#' Check the size factors (i.e. the sequencing library size effect)
dds.counts_TEs <- estimateSizeFactors(dds.counts_TEs)
sizes.counts_TEs <- sizeFactors(dds.counts_TEs)
names(sizes.counts_TEs) <- colnames(counts_TEs)
pander(sizes.counts_TEs)
boxplot(sizes.counts_TEs, main="Sequencing libraries size factor")

#' ## Variance stabilising transformation
#' ### blind
vsd.counts_TEs <- varianceStabilizingTransformation(dds.counts_TEs, blind=FALSE)
vst.counts_TEs <- assay(vsd.counts_TEs)
vst.counts_TEs <- vst.counts_TEs - min(vst.counts_TEs)

#' Validate the VST 
meanSdPlot(vst.counts_TEs[rowSums(counts_TEs)>0,])

#' ## PCA
 pc <- prcomp(t(vst.counts_TEs))

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
 legend("bottomleft",pch=19,
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
 legend("bottomleft",pch=19,
 col=pal[1:nlevels(csamples$Stages)],
 legend=levels(csamples$Stages))
 par(mar=mar)

#' Write it out
# ??Permission?? write.csv(vst.counts_TEs,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_TEs_R.csv")

#' Make a selection of genes
sel_TEs <- featureSelect(vst.counts_TEs,conditions = csamples$Stages,exp = 2,nrep = 2)

res.sel_TEs <- sapply(2:8,function(i){sum(featureSelect(vst.counts_TEs,conditions = csamples$Stages,exp = i,nrep = 2))})

plot(res.sel_TEs)

message(sprintf("%s genes are selected for the heatmap",sum(sel_TEs)))

#' Rename genes and samples
rownames(vst.counts_TEs[sel_TEs,]) <- sub("\\.1$","",rownames(vst.counts_TEs[sel_TEs,]))

colData(dds.counts_TEs)$SeidrID <- paste(colData(dds.counts_TEs)$Batch,colData(dds.counts_TEs)$Stages,sep="_")

colnames(vst.counts_TEs) <- colData(dds.counts_TEs)$SeidrID

ord_TEs <- order(colData(dds.counts_TEs)$Stages,colData(dds.counts_TEs)$Batch)

#' Create the anova table
anv.table <- data.frame(Experiment=1,Perturbations=NA,PerturbationLevels=NA,
                        Treatment=NA,DeletedGenes=NA,OverexpressedGenes=NA,
                        Time=colData(dds.counts_TEs)$Stages[ord_TEs],
                        Repeat=unlist(lapply(lapply(table(colData(dds.counts_TEs)$Stages[ord_TEs]),":",1),rev)))

# write.table(anv.table,quote=FALSE,sep="\t",row.names=FALSE,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs//NetworkE_chip_features_TEs_R.tsv")

#' Write out a selection of genes
# write.table(vst.counts_TEs[sel_TEs,ord_TEs],file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection_TEs_R.tsv", sep = "\t", quote = FALSE)


#' #### Soft clustering (Mfuzz)
suppressPackageStartupMessages(library(Mfuzz))
suppressPackageStartupMessages(library(org.At.tair.db))

#' Genes and TEs

#' Load DESeq object
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/DESeqDataSet-3_genes_and_TEs_R_blindf.rda")

#' # Differential Expression
dds_genes_and_TEs <- DESeq(dds.counts_genes_and_TEs)

#' Create the eSet
vst_genes_and_TEs.g <- read.table("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection_genes_and_TEs_R.tsv",
                    row.names=1)
alpha=0.01
lfc=0.5
#for all results?
allres_genes_and_TEs <- lapply(2:8,function(i){
  res_genes_and_TEs <- results(dds_genes_and_TEs,c("Stages",as.character(i),as.character(i-1)))
  return(rownames(res_genes_and_TEs[res_genes_and_TEs$padj <= alpha & abs(res_genes_and_TEs$log2FoldChange) >= lfc & ! is.na(res_genes_and_TEs$padj) ,]))
})

sel_genes_and_TEs <- unique(unlist(allres_genes_and_TEs))
eset_genes_and_TEs <- ExpressionSet(as.matrix(vst_genes_and_TEs.g[rownames(vst_genes_and_TEs.g) %in% sel_genes_and_TEs,]))
#' Standardise
eset_genes_and_TEs.s <- standardise(eset_genes_and_TEs)
#' Average the replicates (mean)
samples <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3_genes_and_TEs_R.csv")
eset_genes_and_TEs.m <- ExpressionSet(
  sapply(split.data.frame(t(exprs(eset_genes_and_TEs.s)),
                          samples$Stages),
         colMeans))
#' Estimate the fuzzification
m_genes_and_TEs <- mestimate(eset_genes_and_TEs.m)
#' Find the clusters (8 is based on the previous heatmap)
cl_genes_and_TEs <- mfuzz(eset_genes_and_TEs.m,centers=24,m=m_genes_and_TEs)
# cl$cluster1
#' There are a number of clusters that behave similarly
par(mar=c(0.1,0.1,0.1,0.1))
mfuzz.plot(eset_genes_and_TEs.m,
           cl=cl_genes_and_TEs,
           mfrow=c(6,4),
           time.labels = colnames(eset_genes_and_TEs.m),
           new.window=FALSE)

#' TEs


#' Load DESeq object
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/DESeqDataSet-3_TEs_R_blindf.rda")

#' # Differential Expression
dds_TEs <- DESeq(dds.counts_TEs)

#' Create the eSet
vst_TEs.g <- read.table("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection_TEs_R.tsv",
                     row.names=1)
 alpha=0.01
 lfc=0.5
#for all results?
 allres_TEs <- lapply(2:8,function(i){
   res_TEs <- results(dds_TEs,c("Stages",as.character(i),as.character(i-1)))
   return(rownames(res_TEs[res_TEs$padj <= alpha & abs(res_TEs$log2FoldChange) >= lfc & ! is.na(res_TEs$padj) ,]))
 })

 sel_TEs <- unique(unlist(allres_TEs))
 eset_TEs <- ExpressionSet(as.matrix(vst_TEs.g[rownames(vst_TEs.g) %in% sel_TEs,]))
#' Standardise
 eset_TEs.s <- standardise(eset_TEs)
#' Average the replicates (mean)
 samples_TEs <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3_TEs_R.csv")
 samples_TEs <- samples_TEs[ , -1]
 samples_TEs$Stages <- as.factor(samples_TEs$Stages)
 eset_TEs.m <- ExpressionSet(
   as.matrix(sapply(split.data.frame(t(exprs(eset_TEs.s)),
                           samples_TEs$Stages),
          colMeans),nrow=1))
#' A single value is present
par(mfrow=c(1,1),mar=c(5.1,4.1,4.1,2.1))
plot(exprs(eset_TEs.m),type="l")
 
#' samples colored by stages
sel_TEs_16 <- allres_TEs[[2]]
counts_DE_TEs <- vst.counts_genes_and_TEs[sel_TEs_16, ]
plot.default(counts_DE_TEs)
plot.default(counts_DE_TEs[1:8])
plot.default(counts_DE_TEs[9:16])


#' # Heatmap

heatmap.2(counts_TEs, scale = "row", labRow = NULL, labCol = paste0(csamples$Batch, "-", csamples$Stages))
 
#' ## Heatmap
#' ### Feature selection
# sel <- featureSelect(vst.counts,conditions = csamples$Stages,exp = 5,nrep = 2)
# message(sprintf("%s genes are selected for the heatmap",sum(sel)))

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


#' WRITE ME UP
#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
