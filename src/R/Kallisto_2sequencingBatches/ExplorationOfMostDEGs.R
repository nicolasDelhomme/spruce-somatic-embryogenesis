#' ---
#' title: "Exploration of the most DE genes in spruce SE"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---

#' # Setup  
#' # Environment
#' Set the working dir
setwd("~/Git/UPSCb/projects/spruce-somatic-embryogenesis")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="~/Git/UPSCb/projects/spruce-somatic-embryogenesis")
#' ```
#' Libs
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(Glimma))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(stringr))

#' Load DESeq object
load("DESeqDataSet-3_S.rda")

#' Setup graphics
pal=brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
mar <- par("mar")

#' # Analysis

#' ## Differential Expression
dds <- DESeq(dds.counts)

#' Extract differentially expressed genes (DEGs) 
resultsNames(dds)

res_2vs1 <- results(dds,c("Stages","S2","S1"))
res_3vs2 <- results(dds,c("Stages","S3","S2"))
res_4vs3 <- results(dds,c("Stages","S4","S3"))
res_5vs4 <- results(dds,c("Stages","S5","S4"))
res_6vs5 <- results(dds,c("Stages","S6","S5"))
res_7vs6 <- results(dds,c("Stages","S7","S6"))
res_8vs7 <- results(dds,c("Stages","S8","S7"))

#' ## Filter DEGs by log2fc and FDR  
#'   
#' Select significant genes (FDR < 0.01, |log2FC| => 0.5)
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
  if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
  if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
  if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
} 

res_list <- list(res_2vs1, res_3vs2, res_4vs3, res_5vs4, res_6vs5, res_7vs6, res_8vs7)
names(res_list) <- c("res_2vs1", "res_3vs2", "res_4vs3", "res_5vs4", "res_6vs5", "res_7vs6", "res_8vs7")
res_sig_list <- lapply(res_list, sigDeg)

#' Number of DEGs between the stages:
elementNROWS(res_sig_list)
barplot(elementNROWS(res_sig_list),
        ylim = c(0,20000))

#' ## Explore baseMean of DEGs  
#'  
#' baseMean = mean of normalised counts for all samples
#' 
#' ### baseMean of genes with the highest abs(log2fc)  
#' Arrange DEGs by highest absolute log2fc value and print their baseMean
res_sig_list_byLog2fc <- lapply(names(res_sig_list), function(x){
  res_sig_list[[x]][order(abs(res_sig_list[[x]][,"log2FoldChange"]), decreasing = TRUE), ]
})
names(res_sig_list_byLog2fc) <- names(res_sig_list)

sapply(res_sig_list_byLog2fc, function(x){
  round(x[1:20, "baseMean"], 1)
})

#' Some highly differentially expressed genes have very low baseMean, the value can be as low as 1.  
#' Plot expression profile of some of these genes:  
# plot profile of genes with baseMean less than 10 among top 5 genes ordered by log2fc in comparison of different stages
lapply(names(res_sig_list_byLog2fc), function(x){
  first5 <- res_sig_list_byLog2fc[[x]][1:5, ]
  genes <- rownames(first5[first5$baseMean < 10, ])
  print(x)
  lapply(genes, function(gene){
    plotCounts(dds, gene, "Stages")
    plotCounts(dds, gene, "Batch")
  })
})

# this way you can obtain count data for a gene: plotCounts(dds, "MA_170445g0010.1", "Stages", returnData = TRUE)

#' ### Distribution of baseMean of DE genes
boxplot(sapply(res_sig_list, function(x){x[,"baseMean"]}),
        log = "y",
        ylab = "baseMean (log)",
        xlab = "comparison of 2 consecutive stages of somatic embryogenesis",
        main = "Distribution of baseMean of DEGs between two consecutive stages of SE",
        notch = TRUE,
        col = pal[3])

#' ### Is there a correlation between baseMean and padj?
lapply(names(res_sig_list), function(x){
  plot(res_sig_list[[x]][order(res_sig_list[[x]]$padj), "padj"],
       res_sig_list[[x]][order(res_sig_list[[x]]$padj), "baseMean"],
       xlab = "padj",
       ylab = "baseMean",
       main = x)
})

#' There are genes with high baseMean, enriched in low padj, but as expected there is no correlation.  
#' NOTE: it would be better to plot it on -log scale (as the example later in the document)
#' 
#' ### What would be an appropriate threshold to filter the genes?  
#' 
#' Number of DEGs at different cutoff values of baseMean  
#' Check wider values of baseMean  
# define baseMean
bM <- seq(0, 300, by = 10)

# extract nr of DEGs at different baseMean cutoffs
nrDEGs_baseMean <- sapply(bM, function(i){
  sapply(res_sig_list, function(x){
    unique(elementNROWS(x[x[, "baseMean"] > i, ]))
  })
})

colnames(nrDEGs_baseMean) <- bM

# plot
matplot(t(nrDEGs_baseMean), 
        type = "l",
        col = pal[1:7],
        lty = 1,
        lwd = 2,
        xaxt = "n",
        xlab = "baseMean cutoff value",
        ylab = "number of DEGs",
        main = "Number of DEGs between different stages of SE at different baseMean cutoff values")
axis(1, at =  c(1:length(bM)), labels = bM)
legend("topright", legend = rownames(nrDEGs_baseMean), col = pal[1:7], lty = 1, lwd = 2)

#' Check narrower range of baseMean  
# define baseMean
bM <- seq(0, 30, by = 2)

# extract nr of DEGs at different baseMean cutoffs
nrDEGs_baseMean <- sapply(bM, function(i){
  sapply(res_sig_list, function(x){
    unique(elementNROWS(x[x[, "baseMean"] > i, ]))
  })
})

colnames(nrDEGs_baseMean) <- bM

# plot
matplot(t(nrDEGs_baseMean), 
        type = "l",
        col = pal[1:7],
        lty = 1,
        lwd = 2,
        xaxt = "n",
        xlab = "baseMean cutoff value",
        ylab = "number of DEGs",
        main = "Number of DEGs between different stages of SE at different baseMean cutoff values")
axis(1, at =  c(1:length(bM)), labels = bM)
legend("bottomright", legend = rownames(nrDEGs_baseMean), col = pal[1:7], lty = 1, lwd = 2)

#' ### Diagnostics for independent filtering  
#' 
#' Examples taken from: Diagnostics for independent filtering: choosing filter statistic and cutoff (W. Huber)  
#'   
#'   
#' Replicate diagnostics plots from above reference for the problematic time points in our experiment (3vs2, 7vs6) 

with(res_3vs2, plot(rank(baseMean)/length(baseMean), -log10(pvalue), pch=16, cex=0.45))

with(res_7vs6, plot(rank(baseMean)/length(baseMean), -log10(pvalue), pch=16, cex=0.45))

library(genefilter)
theta = seq(from=0, to=0.5, by=0.1)
pBH = filtered_p(filter=res_3vs2$baseMean, test=res_3vs2$pvalue, theta=theta, method="BH")
rejection_plot(pBH, at="sample",
               xlim=c(0, 0.5), ylim=c(0, 30000), 
               xlab="FDR cutoff (Benjamini & Hochberg adjusted p-value)", main="3vs2")

pBH = filtered_p(filter=res_7vs6$baseMean, test=res_7vs6$pvalue, theta=theta, method="BH")
rejection_plot(pBH, at="sample",
               xlim=c(0, 0.5), ylim=c(0, 40000), 
               xlab="FDR cutoff (Benjamini & Hochberg adjusted p-value)", main="7vs6")

# the plot shows the number of rejections (i. e. genes detected as differentially expressed) 
# as a function of the FDR threshold (x-axis) and the filtering cutoff θ (line colours, specifies
# quantiles of the distribution of the filter statistic)


theta = seq(from=0, to=0.8, by=0.02)
rejBH = filtered_R(alpha=0.1, filter=res_3vs2$baseMean, test=res_3vs2$pvalue, theta=theta, method="BH")
plot(theta, rejBH, type="l",xlab=expression(theta), ylab="number of rejections", main = "3vs2")

rejBH = filtered_R(alpha=0.1, filter=res_7vs6$baseMean, test=res_7vs6$pvalue, theta=theta, method="BH")
plot(theta, rejBH, type="l",xlab=expression(theta), ylab="number of rejections", main = "7vs6")

# If we select a fixed cutoff for the adjustedp-values, we can also look more closely at the relationship 
# between the fraction of null hypotheses filtered and the total number of discoveries.Thefiltered_Rfunction 
# wraps filtered_pand just returns rejection counts. It requires you to choose a particular p-value cutoff, 
# specified through the argument alpha.


filterChoices = data.frame(`mean` = rowMeans(counts(dds)),
                           `median`  = rowMedians(counts(dds)))

rejChoices = sapply(filterChoices, function(f){
  filtered_R(alpha=0.1, filter=f, test=res_3vs2$pvalue, theta=theta, method="BH")
})
myColours = c("red", "blue")
matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,xlab=expression(theta), ylab="number of rejections")
legend("bottomleft", legend=colnames(filterChoices), fill=myColours)

rejChoices = sapply(filterChoices, function(f){
  filtered_R(alpha=0.1, filter=f, test=res_7vs6$pvalue, theta=theta, method="BH")
})
matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,xlab=expression(theta), ylab="number of rejections")
legend("bottomleft", legend=colnames(filterChoices), fill=myColours)

# The number of rejections at FDR=10% as a function of θ for a number of different choices of the filter statistic.  


#' ### use median instead of mean, while extracting results
res_2vs1_med <- results(dds,c("Stages","S2","S1"), filter = rowMedians(counts(dds)))
res_3vs2_med <- results(dds,c("Stages","S3","S2"), filter = rowMedians(counts(dds)))
res_4vs3_med <- results(dds,c("Stages","S4","S3"), filter = rowMedians(counts(dds)))
res_5vs4_med <- results(dds,c("Stages","S5","S4"), filter = rowMedians(counts(dds)))
res_6vs5_med <- results(dds,c("Stages","S6","S5"), filter = rowMedians(counts(dds)))
res_7vs6_med <- results(dds,c("Stages","S7","S6"), filter = rowMedians(counts(dds)))
res_8vs7_med <- results(dds,c("Stages","S8","S7"), filter = rowMedians(counts(dds)))

res_list_med <- list(res_2vs1_med, res_3vs2_med, res_4vs3_med, res_5vs4_med, res_6vs5_med, res_7vs6_med, res_8vs7_med)
names(res_list_med) <- c("res_2vs1", "res_3vs2", "res_4vs3", "res_5vs4", "res_6vs5", "res_7vs6", "res_8vs7")
res_sig_list_med <- lapply(res_list_med, sigDeg)

#' Arrange DEGs by highest/lowest log2fc and print their baseMean
res_sig_list_med_byLog2fc <- lapply(names(res_sig_list_med), function(x){
  res_sig_list_med[[x]][order(abs(res_sig_list_med[[x]][,"log2FoldChange"]), decreasing = TRUE), ]
})
names(res_sig_list_med_byLog2fc) <- names(res_sig_list_med)

sapply(res_sig_list_med_byLog2fc, function(x){
  round(x[1:20, "baseMean"], 1)
})

#' ### DEGs with the lowest baseMean
sapply(res_sig_list_med_byLog2fc, function(x){
  round(min(x[,"baseMean"]), 1)
})

#' Plot counts of the genes with the lowest baseMean in each comparison
lapply(res_sig_list_med_byLog2fc, function(x){
  print(x[x$baseMean == min(x[,"baseMean"]), ])
  plotCounts(dds, rownames(x[x$baseMean == min(x[,"baseMean"]), ]), "Stages")
})


#' This looks much better now. There are still genes with very low baseMean among DEGs, but the expression profile seems clear
#' as values in replicates are more similar than before.  
#' 
#' ### Number of DEGs
elementNROWS(res_sig_list_med)
barplot(elementNROWS(res_sig_list_med), ylim = c(0,20000))

#' Counts of DEGs are slightly lower, but not much. Altogether there is less genes:  
#' number
sum(elementNROWS(res_sig_list))-sum(elementNROWS(res_sig_list_med))
elementNROWS(res_sig_list)-elementNROWS(res_sig_list_med)
#' percent
(sum(elementNROWS(res_sig_list))-sum(elementNROWS(res_sig_list_med)))/sum(elementNROWS(res_sig_list))*100
(elementNROWS(res_sig_list)-elementNROWS(res_sig_list_med))/elementNROWS(res_sig_list)*100

#' Are DEGs the same (except the ones that are missing) or are they different?
#' 
#' Genes that are common to both sets of results or different between them
excluded_DEGs <- sapply(seq_along(res_sig_list), function(x){
  DEGs_mean <- rownames(res_sig_list[[x]])
  DEGs_median <- rownames(res_sig_list_med[[x]])
  mean_only <- length(setdiff(DEGs_mean, DEGs_median))
  common <- length(intersect(DEGs_mean, DEGs_median))
  median_only <- length(setdiff(DEGs_median, DEGs_mean))
  nr_mean <- length(DEGs_mean)
  nr_median <- length(DEGs_median)
  return(c(nr_mean, nr_median, mean_only, common, median_only))
})

excluded_DEGs <- as.data.frame(t(excluded_DEGs))
rownames(excluded_DEGs) <- names(res_sig_list)
colnames(excluded_DEGs) <- c("nr_DEGs_mean", "nr_DEGs_median", "mean_only", "common", "median_only")
# sum of genes that are changing (in one or the other dataset)
excluded_DEGs$all_changing_DEGs <- apply(excluded_DEGs, 1, function(x){
  sum(x["mean_only"], x["median_only"])
})

excluded_DEGs

#' There are also genes that are found only in the set of results for which median was used as a filtering method.
#' True number of different genes between two different set of results is in the last column of the table above.  

#' Are the "problematic genes" excluded from the list of DEGs, when we use median instead of mean?  
# get names of the genes, which were plotted as the bad examples in the beginnig of the script
problematic_genes <- lapply(names(res_sig_list_byLog2fc), function(x){
  first5 <- res_sig_list_byLog2fc[[x]][1:5, ]
  rownames(first5[first5$baseMean < 10, ])
  })
# is any of those genes still present among the DEGs?
lapply(seq_along(res_sig_list_med), function(x){
  any(rownames(res_sig_list_med[[x]]) %in% problematic_genes[[x]])
})

#' ### Distribution of baseMean of DEGs filtered by median.
boxplot(sapply(res_sig_list_med, function(x){x[,"baseMean"]}),
        log = "y",
        ylab = "baseMean (log)",
        xlab = "comparison of 2 consecutive stages of somatic embryogenesis",
        main = "Distribution of baseMean of DEGs between two consecutive stages of SE",
        notch = TRUE,
        col = pal[3])

#' It looks the same as before.  

#' ### Is there any correlation between baseMean?
lapply(names(res_sig_list_med), function(x){
  plot(res_sig_list_med[[x]][order(res_sig_list_med[[x]]$padj), "padj"],
       res_sig_list_med[[x]][order(res_sig_list_med[[x]]$padj), "baseMean"],
       xlab = "padj",
       ylab = "baseMean",
       main = x)
})


