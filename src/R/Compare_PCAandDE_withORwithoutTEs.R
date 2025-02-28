#' title: "Compare analysis of genes with or without TEs"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
#' # Aim  
#' Get PCA plots for publications. 
#' Check if anything changes if analysis is done with or without TEs 
#' - probably not as there are only ~300 TEs and many more genes (~66000). 
#' For the purpose of this project they should be excluded, but it might be good for future 
#' analysis to leave them in, as it would be easy to compare analysis between them.  
#' 
#'  Setup  
#' ## Libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))

#' ## DESeq object
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes+TEs.rda")
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes.rda")

#' ## Functions
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
    if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
    if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
    if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
} 

#' ## Graphics
pal <- brewer.pal(12,"Paired")
mar <- par("mar")

#' # Analysis of dataset with only genes included  
#' ## VST  
vsd_genes <- varianceStabilizingTransformation(dds_genes, blind=TRUE)
vst_genes <- assay(vsd_genes)
vst_genes <- vst_genes - min(vst_genes)

vsda_genes <- varianceStabilizingTransformation(dds_genes, blind=FALSE)
vsta_genes <- assay(vsda_genes)
vsta_genes <- vsta_genes - min(vsta_genes)

#' ## PCA  
pc_g <- prcomp(t(vsta_genes))
percent_g <- round(summary(pc_g)$importance[2,]*100)

barplot(percent_g, ylim = c(0,40), main = "% of variation explained by PCs")
plot(cumsum(percent_g), 
     xlab = "Number of PCs", 
     ylab = "Percent variation", 
     main = "Variation explained by PCs (cummulative)")

plot(pc_g$x[,1],
     pc_g$x[,2],
     xlab=paste("Comp. 1 (",percent_g[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent_g[2],"%)",sep=""),
     col=pal[as.integer(as.factor(colData(dds_genes)$Stages))],
     main="1st and 2nd PC",
     pch=19)
legend("top",pch=19,
       col=pal[1:nlevels(as.factor(colData(dds_genes)$Stages))],
       legend=levels(as.factor(colData(dds_genes)$Stages)))

plot(pc_g$x[,2],
     pc_g$x[,3],
     xlab=paste("Comp. 2 (",percent_g[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent_g[3],"%)",sep=""),
     col=pal[as.integer(as.factor(colData(dds_genes)$Stages))],
     main="2nd and 3rd PC",
     pch=19)
legend("top",pch=19,
       col=pal[1:nlevels(as.factor(colData(dds_genes)$Stages))],
       legend=levels(as.factor(colData(dds_genes)$Stages)))

#' ## DE
dds_genes <- DESeq(dds_genes)

res_list_genes <- lapply(2:8, function(i){
    results(dds_genes, c("Stages", paste0("S", as.character(i)), paste0("S", as.character(i-1))), 
            filter = rowMedians(counts(dds_genes)))
})

names(res_list_genes) <- lapply(2:8, function(i){
    paste0("res_", i, "vs", i-1)
})

#' Filter
res_sig_genes <- lapply(res_list_genes, sigDeg)

#' Count all, up-regulated and down-regulated DE genes

nr_DEG <- sapply(res_sig_genes, function(x){
    all <- nrow(x)
    up <- nrow(sigDeg(x, genes = "up"))
    down <- nrow(sigDeg(x, genes = "down"))
    return(cbind(all, up, down))
})
rownames(nr_DEG) <- c("all", "up", "down")

#' Number of all DEGs

barplot(nr_DEG[1,],
        beside = TRUE, 
        main = "Number of all DEGs",
        xlab = "Stage comparison", 
        ylab = "Number of DEGs", 
        names.arg = sub("res_", "", colnames(nr_DEG)), 
        ylim = c(0,20000))

#' # Analysis of genes and TEs together  
#' ## VST  
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)
vst <- assay(vsd)
vst <- vst - min(vst)

vsda <- varianceStabilizingTransformation(dds, blind=FALSE)
vsta <- assay(vsda)
vsta <- vsta - min(vsta)

#' ## PCA  
pc <- prcomp(t(vsta))
percent <- round(summary(pc)$importance[2,]*100)

barplot(percent, ylim = c(0,40), main = "% of variation explained by PCs")
plot(cumsum(percent), 
     xlab = "Number of PCs", 
     ylab = "Percent variation", 
     main = "Variation explained by PCs (cummulative)")

plot(pc$x[,1],
     pc$x[,2],
     xlab=paste("Comp. 1 (",percent[1],"%)",sep=""),
     ylab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     col=pal[as.integer(as.factor(colData(dds)$Stages))],
     main="1st and 2nd PC",
     pch=19)
legend("top",pch=19,
       col=pal[1:nlevels(as.factor(colData(dds)$Stages))],
       legend=levels(as.factor(colData(dds)$Stages)))

plot(pc$x[,2],
     pc$x[,3],
     xlab=paste("Comp. 2 (",percent[2],"%)",sep=""),
     ylab=paste("Comp. 3 (",percent[3],"%)",sep=""),
     col=pal[as.integer(as.factor(colData(dds)$Stages))],
     main="2nd and 3rd PC",
     pch=19)
legend("top",pch=19,
       col=pal[1:nlevels(as.factor(colData(dds)$Stages))],
       legend=levels(as.factor(colData(dds)$Stages)))

#' ## DE
dds <- DESeq(dds)

res_list <- lapply(2:8, function(i){
    results(dds, c("Stages", paste0("S", as.character(i)), paste0("S", as.character(i-1))), 
            filter = rowMedians(counts(dds)))
})

names(res_list) <- lapply(2:8, function(i){
    paste0("res_", i, "vs", i-1)
})

#' Filter
res_sig <- lapply(res_list, sigDeg)

#' Count all, up-regulated and down-regulated DE genes

nr_DEGT <- sapply(res_sig, function(x){
    all <- nrow(x)
    up <- nrow(sigDeg(x, genes = "up"))
    down <- nrow(sigDeg(x, genes = "down"))
    return(cbind(all, up, down))
})
rownames(nr_DEGT) <- c("all", "up", "down")

#' Number of all DEGs

barplot(nr_DEGT[1,],
        beside = TRUE, 
        main = "Number of all DEGs",
        xlab = "Stage comparison", 
        ylab = "Number of DEGs and TEs", 
        names.arg = sub("res_", "", colnames(nr_DEGT)), 
        ylim = c(0,20000))

#' Compare number of DEGs in both dds objects
res_sig_subGenes <- lapply(res_sig, function(x){
    x[grepl("MA_", rownames(x)), ]
})

nr_subDEGs <- sapply(res_sig_subGenes, function(x){
    all <- nrow(x)
    up <- nrow(sigDeg(x, genes = "up"))
    down <- nrow(sigDeg(x, genes = "down"))
    return(cbind(all, up, down))
})
rownames(nr_subDEGs) <- c("all", "up", "down")

nr_subDEGs
nr_DEG

#' There is subtle change in number of DEGs, whether TEs are included in the analysis or not, 
#' though it is so small, that it is negligible.