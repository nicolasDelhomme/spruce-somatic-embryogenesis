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

#' ## Combining samples
samples$ID <- sub("_L00[1,2]", "",
                  samples$ScilifeID)

counts <- do.call(
  cbind,
  lapply(split.data.frame(t(kg),
                          samples$ID),
         colSums))

csamples <- samples[,-1]
csamples <- csamples[match(colnames(counts),csamples$ID),]
write.csv(csamples,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3_genes_and_TEs.csv")



#' ## QC
#' Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ### Coverage
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

#' ### Raw data export
dir.create(file.path("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis","kallisto_LTR-TEs"),showWarnings=FALSE)
write.csv(counts,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/raw-unormalised-gene-expression_data-3_genes_and_TEs_blindf.csv")

#' ## Data normalisation 
#' ### Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' ### Setup
#' Create the dds object
rownames(csamples) <- csamples$ID
dds.counts_genes_and_TEs <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = csamples[,c("SubmittedID","Stages","Description","Batch","ID")],
  design = ~Batch+Stages)

#' Save it
save(dds.counts_genes_and_TEs,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/DESeqDataSet-3_genes_and_TEs_blindf.rda")

#' Check the size factors (i.e. the sequencing library size effect)
dds.counts_genes_and_TEs <- estimateSizeFactors(dds.counts_genes_and_TEs)
sizes.counts_genes_and_TEs <- sizeFactors(dds.counts_genes_and_TEs)
names(sizes.counts_genes_and_TEs) <- colnames(counts)
pander(sizes.counts_genes_and_TEs)
boxplot(sizes.counts_genes_and_TEs, main="Sequencing libraries size factor")

#' ### Rlog transformation
#' #### blind
rlogd.counts_genes_and_TEs <- rlogTransformation(dds.counts_genes_and_TEs, blind=FALSE)
rlog.counts_genes_and_TEs <- assay(rlogd.counts_genes_and_TEs)
rlog.counts_genes_and_TEs <- rlog.counts_genes_and_TEs - min(rlog.counts_genes_and_TEs)

#' Validate the VST 
meanSdPlot(rlog.counts_genes_and_TEs[rowSums(counts)>0,])

#' Write it out
write.csv(rlog.counts_genes_and_TEs,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/library-size-corrected_blindFALSE-rlog_gene-expression_data-3_genes_and_TEs.csv")



#' # TEs (kg_LTR)  
names_LTR <- rownames(kg)[grepl("^MA_", rownames(kg)) == FALSE]
kg_LTR <- kg[names_LTR, ]

#' ## Combining samples
samples$ID <- sub("_L00[1,2]", "",
                  samples$ScilifeID)

counts <- do.call(
  cbind,
  lapply(split.data.frame(t(kg_LTR),
                          samples$ID),
         colSums))

csamples <- samples[,-1]
csamples <- csamples[match(colnames(counts),csamples$ID),]
write.csv(csamples,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3_LTR-TEs.csv")



#' ## QC
#' Check how many genes are never expressed
sel <- rowSums(counts) == 0
sprintf("%s%% percent (%s) of %s genes are not expressed",
        round(sum(sel) * 100/ nrow(counts),digits=1),
        sum(sel),
        nrow(counts))

#' ### Coverage
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

#' ### Raw data export
dir.create(file.path("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis","kallisto_LTR-TEs"),showWarnings=FALSE)
write.csv(counts,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/raw-unnormalised-gene-expression_data-3_LTR-TEs_blindf.csv")

#' ## Data normalisation 
#' ### Preparation
#'  For visualization, the data is submitted to a variance stabilization
#' transformation using DESeq2. The dispersion is estimated independently
#' of the sample tissue and replicate

#' ### Setup
#' Create the dds object
rownames(csamples) <- csamples$ID
dds.counts_TEs <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = csamples[,c("SubmittedID","Stages","Description","Batch","ID")],
  design = ~Batch+Stages)

#' Save it
save(dds.counts_TEs,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/DESeqDataSet-3_LTR-TEs_blindf.rda")

#' Check the size factors (i.e. the sequencing library size effect)
dds.counts_TEs <- estimateSizeFactors(dds.counts_TEs)
sizes.counts_TEs <- sizeFactors(dds.counts_TEs)
names(sizes.counts_TEs) <- colnames(counts)
pander(sizes.counts_TEs)
boxplot(sizes.counts_TEs, main="Sequencing libraries size factor")

#' ### Rlog transformation
#' #### blind
rlogd.counts_TEs <- rlogTransformation(dds.counts_TEs, blind=FALSE)
rlog.counts_TEs <- assay(rlogd.counts_TEs)
rlog.counts_TEs <- rlog.counts_TEs - min(rlog.counts_TEs)

#' Validate the VST 
meanSdPlot(rlog.counts_TEs[rowSums(counts)>0,])

#' Write it out
write.csv(rlog.counts_TEs,file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/library-size-corrected_blindFALSE-rlog_gene-expression_data-3_LTR-TEs.csv")
