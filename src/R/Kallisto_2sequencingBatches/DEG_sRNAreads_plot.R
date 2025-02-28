#' Aim: prepare figure with number of DEGs and number of sRNA reads of different length between each stage of SE
#' 
#' # Get nr of DEGs

#' Set the working dir
setwd("~/Git/UPSCb/projects/spruce-somatic-embryogenesis")

suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(RColorBrewer))

pal <- brewer.pal(8, "Set2")

#' Load DESeq object
load("DESeqDataSet-3_S.rda")

#' Differential Expression
dds <- DESeq(dds.counts)

#' Obtain the results 
res_2vs1_med <- results(dds,c("Stages","S2","S1"), filter = rowMedians(counts(dds)))
res_3vs2_med <- results(dds,c("Stages","S3","S2"), filter = rowMedians(counts(dds)))
res_4vs3_med <- results(dds,c("Stages","S4","S3"), filter = rowMedians(counts(dds)))
res_5vs4_med <- results(dds,c("Stages","S5","S4"), filter = rowMedians(counts(dds)))
res_6vs5_med <- results(dds,c("Stages","S6","S5"), filter = rowMedians(counts(dds)))
res_7vs6_med <- results(dds,c("Stages","S7","S6"), filter = rowMedians(counts(dds)))
res_8vs7_med <- results(dds,c("Stages","S8","S7"), filter = rowMedians(counts(dds)))

#' Filter DEGs by lfc and padj
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
  if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
  if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
  if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
} 

res_list <- list(res_2vs1_med, res_3vs2_med, res_4vs3_med, res_5vs4_med, res_6vs5_med, res_7vs6_med, res_8vs7_med)
names(res_list) <- c("res_2vs1", "res_3vs2", "res_4vs3", "res_5vs4", "res_6vs5", "res_7vs6", "res_8vs7")
res_sig_list <- lapply(res_list, sigDeg)

#' Number of DEGs
nr_DEgenes <- (elementNROWS(res_sig_list))
names(nr_DEgenes) <- sub("res_", "", names(res_sig_list))

#' # Get normalised counts of sRNA reads of different length  
#' Load normalised counts of redundant reads in all the stages
load("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/mapping_stats/LibSizeNormNrRedReads_byStage.Rda")

#' Load counts of non-redundant sequences in different stages
load("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/mapping_stats/NrNonredReads_byStage.Rda")

#' # Plot  
#' 
#' Make a plot with double y axis, one for nr of DEGs and the other one for nr of redundant sRNA reads  
#' 
#' ## with redundant sRNA reads
par(mar = c(5, 5, 4, 5))

pdf("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DEGs_redundantsRNAreads_SE.pdf", width = 10, height = 6)

barplot2(nr_DEgenes, 
         beside = TRUE, 
         border = "grey", 
         xlab = "Stage",
         xaxt = "n",
         ylab = "Number of DE genes", 
         ylim = c(0,20000),
         main = "Number of DEGs and sRNA reads in SE")

par(new = TRUE)

# code for the sRNA plot:
plot(red_SE_norm_byStage$count.18.sum,
type = "l",
lwd = 2,
col = pal[1],
yaxt = "n",
ylim = c(0, 400000),
ylab = "",
xlab = "",
bty = "u")

lines(red_SE_norm_byStage$count.19.sum, col = pal[4], lwd = 2)
lines(red_SE_norm_byStage$count.20.sum, col = pal[7], lwd = 2)
lines(red_SE_norm_byStage$count.21.sum, col = pal[2], lwd = 2)
lines(red_SE_norm_byStage$count.22.sum, col = pal[6], lwd = 2)
lines(red_SE_norm_byStage$count.23.sum, col = pal[5], lwd = 2)
lines(red_SE_norm_byStage$count.24.sum, col = pal[3], lwd = 2)

legend("topleft",
       legend = paste(str_extract(colnames(red_SE_norm_byStage[ , 22:28]), "[0-9]+"), "nt"),
       col = pal[c(1,4,7,2,6,5,3)],
       lwd = 2,
       bty = "n", cex = 0.8, title = "sRNA size")

axis(side = 4)
mtext("Redundant mapped sRNA reads (rpm)", side = 4, line = 3)

dev.off()

#' ## with non-redundant sRNA reads  
#' 
#' Calculate % of non-redundant reads
nonred_SE_norm_byStage <- nonred_SE_byStage/nonred_SE_byStage$count.all.sum*100

par(mar = c(5, 5, 4, 5))

pdf("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DEGs_nonredundantsRNAreads_SE.pdf", width = 10, height = 6)

barplot2(nr_DEgenes, 
         beside = TRUE, 
         border = "grey", 
         xlab = "Stage",
         xaxt = "n",
         ylab = "Number of DE genes", 
         ylim = c(0,20000),
         main = "Number of DEGs and sRNA reads in SE")

par(new = TRUE)

# code for the sRNA plot:
plot(nonred_SE_norm_byStage$count.18.sum,
     type = "l",
     lwd = 2,
     col = pal[1],
     yaxt = "n",
     ylim = c(0, 80),
     ylab = "",
     xlab = "",
     bty = "u")

lines(nonred_SE_norm_byStage$count.19.sum, col = pal[4], lwd = 2)
lines(nonred_SE_norm_byStage$count.20.sum, col = pal[7], lwd = 2)
lines(nonred_SE_norm_byStage$count.21.sum, col = pal[2], lwd = 2)
lines(nonred_SE_norm_byStage$count.22.sum, col = pal[6], lwd = 2)
lines(nonred_SE_norm_byStage$count.23.sum, col = pal[5], lwd = 2)
lines(nonred_SE_norm_byStage$count.24.sum, col = pal[3], lwd = 2)

legend("topleft",
       legend = paste(str_extract(colnames(nonred_SE_norm_byStage[ , 22:28]), "[0-9]+"), "nt"),
       col = pal[c(1,4,7,2,6,5,3)],
       lwd = 2,
       bty = "n", cex = 0.8, title = "sRNA size")

axis(side = 4)
mtext("% of non-redundant mapped sRNA reads", side = 4, line = 3)

dev.off()
