#' ---
#' title: "Differential expression of TEs in spruce somatic embryogenesis (SalmonTE)"
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

#' Load libraries
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(gplots))

pal <- brewer.pal(12, "Paired")

mar <- par("mar")

#' # Differential expression

#' Load DESeq object
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/salmonTE/DESeqDataSet-3_TEs.rda")

#' Differential Expression
dds <- DESeq(dds.counts_biorep)

#' Get DE TEs between the developmental stages
alpha=0.01
lfc=0.5

# Extract names
allres <- lapply(2:8,function(i){
  res <- results(dds,c("Stages",as.character(i),as.character(i-1)))
  #return(rownames(res[res$padj <= alpha & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) ,]))
})

sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
    if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
    if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
    if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
}

allres <- lapply(allres, sigDeg)
res_up <- lapply(allres, sigDeg, genes = "up")
res_down <- lapply(allres, sigDeg, genes = "down")

#' Plot number of DE TEs between the stages
nr_DE_TEs <- elementNROWS(allres)

# all DE TEs
barplot(nr_DE_TEs, 
        names.arg = c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8"), 
        col = pal[1], 
        ylim = c(0,25))

# up- and down-regulated TEs
par(mar = c(5.1, 6.1, 5.1, 2.1))
barplot(rbind(elementNROWS(res_up), elementNROWS(res_down)),
        beside = TRUE,
        col = pal[c(3,2)],
        names.arg = c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8"), 
        ylim = c(0,15),
        cex.axis = 1.4, 
        cex.names = 1.4)
mtext(side=1, line=3, "Stages", cex=1.6)
mtext(side=2, line=4, "Number of DE", cex=1.6)
mtext(side=2, line=2.5, "LTR retrotransposons", cex=1.6)
mtext(side=3, line=3, "Number of up- and down-regulated", font=2, cex=1.8)
mtext(side=3, line=1, "LTR retrotransposons", font=2, cex=1.8)

legend("top", bty = "n",
       fill = pal[c(3,2)],
       legend=c("up-regulated", "down-regulated"), cex = 1.4)

par(mar=mar)

# DE Copia and Gypsy TEs
nr_DE_copia <- sapply(allres, function(x) {
    TEname <- rownames(x)
    sum(grepl("^C", TEname), na.rm = TRUE)
})

nr_DE_gypsy <- lapply(allres, function(x) {
    TEname <- rownames(x)
    sum(grepl("^G", TEname), na.rm = TRUE)
})

barplot2(rbind(nr_DE_copia, nr_DE_gypsy), 
         beside = TRUE, 
         names.arg = c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8"), 
         ylim = c(0,12), 
         col = pal[3:4],
         legend.text = c("Copia LTR-TEs", "Gypsy LTR-TEs"))

#' Make a heatplot with all the DE TEs
names_DE_TEs <- unique(unlist(allres))
names_DE_copia <- unique(unlist(sapply(allres, function(x) {
  x[grepl("^C", x)]
})))
names_DE_gypsy <- unique(unlist(sapply(allres, function(x) {
  x[grepl("^G", x)]
})))

#' Read in rlog transformed data, blind = TRUE
rlt.counts_biorep <- read.csv("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/salmonTE/library-size-corrected_blindTRUE-rlog_TE-expression_data-3.csv", row.names = 1)
colnames(rlt.counts_biorep) <- sub("^X", "", colnames(rlt.counts_biorep))

#' Read in sample info
samples_biorep <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3_genes_and_TEs_R.csv")

#' Heatmap
heatmap.2(as.matrix(rlt.counts_biorep[names_DE_TEs, ]), 
          scale = "row",
          labCol = paste0(samples_biorep$Batch, "-", samples_biorep$Stages),
          trace = "none",
          ColSideColors = pal[as.integer(samples_biorep$Stages)]
          )

# Copia
heatmap.2(as.matrix(rlt.counts_biorep[names_DE_copia, ]), 
          scale = "row",
          labCol = paste0(samples_biorep$Batch, "-", samples_biorep$Stages),
          trace = "none",
          ColSideColors = pal[as.integer(samples_biorep$Stages)]
)

# Gypsy
heatmap.2(as.matrix(rlt.counts_biorep[names_DE_gypsy, ]), 
          scale = "row",
          labCol = paste0(samples_biorep$Batch, "-", samples_biorep$Stages),
          trace = "none",
          ColSideColors = pal[as.integer(samples_biorep$Stages)]
)


#' Plot counts of DE TEs  
# calculate mean of the counts in biological replicates (counts per stage for each TE)
mean_counts <- do.call(
  cbind,
  lapply(split.data.frame(t(rlt.counts_biorep),
                          samples_biorep$Stages),
         colSums))

# plot counts of DE gypsy and copia elements
matplot(t(mean_counts[names_DE_copia, ]), type = "l", xlab = "Stages", ylab = "mean normalised counts", main = "Expression profile of DE Copia TEs")
matplot(t(mean_counts[names_DE_gypsy, ]), type = "l", xlab = "Stages", ylab = "mean normalised counts", main = "Expression profile of DE Gypsy TEs")

#' How many TEs is DE between stages 2 and 3 as well as between stages 6 and 7?
commonTEs_2and6 <- intersect(allres[[2]], allres[[6]])
length(commonTEs_2and6)

# plot those TEs that are DE expressed between stages 2-3 and also betweeen 6-7
matplot(t(mean_counts[commonTEs_2and6, ]),
        type = "l",
        xlab = "Stages", 
        ylab = "mean normalised counts",
        col = pal[1:nlevels(factor(commonTEs_2and6))],
        lwd = 2
)
legend("top", commonTEs_2and6, col = pal[1:nlevels(factor(commonTEs_2and6))], lwd = 2, bty = "n")

#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
