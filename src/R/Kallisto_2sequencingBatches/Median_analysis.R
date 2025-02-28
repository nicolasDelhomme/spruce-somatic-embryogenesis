#' ---
#' title: "Spruce-somatic-embryogenesis Median Analysis"
#' author: "Camilla Canovi"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup  
#' # Environment
#' Set the working dir
setwd("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/src/R")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/src/R")
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
suppressPackageStartupMessages(library(VennDiagram))
suppressPackageStartupMessages(library(Mfuzz))
suppressPackageStartupMessages(library(org.At.tair.db))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(treemap))
suppressPackageStartupMessages(library(pheatmap))

#' Helper files
suppressMessages(source("~/Git/UPSCb/src/R/densityPlot.R"))
suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))
suppressMessages(source("~/Git/UPSCb/src/R/plot.multidensity.R"))
suppressMessages(source("~/Git/UPSCb/src/R/gopher.R"))

#' Load DESeq object
load("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/DESeqDataSet-3_S.rda")

#' Setup graphics
pal=brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' ```{r drop limma,echo=FALSE}
#' detach("package:limma")
#' ```
 
#' # Process

#' # Differential Expression
dds <- DESeq(dds.counts)

# Dispersion Estimation
#
# The dispersion estimation is adequate
plotDispEsts(dds)

#'  ## Obtain the results 
resultsNames(dds)
     
res_2vs1_med <- results(dds,c("Stages","S2","S1"), filter = rowMedians(counts(dds)))
res_3vs2_med <- results(dds,c("Stages","S3","S2"), filter = rowMedians(counts(dds)))
res_4vs3_med <- results(dds,c("Stages","S4","S3"), filter = rowMedians(counts(dds)))
res_5vs4_med <- results(dds,c("Stages","S5","S4"), filter = rowMedians(counts(dds)))
res_6vs5_med <- results(dds,c("Stages","S6","S5"), filter = rowMedians(counts(dds)))
res_7vs6_med <- results(dds,c("Stages","S7","S6"), filter = rowMedians(counts(dds)))
res_8vs7_med <- results(dds,c("Stages","S8","S7"), filter = rowMedians(counts(dds)))

#<<<<<<< HEAD
#' all the other comparisons
res_3vs1_med <- results(dds,c("Stages","S3","S1"), filter = rowMedians(counts(dds)))
res_4vs1_med <- results(dds,c("Stages","S4","S1"), filter = rowMedians(counts(dds)))
res_5vs1_med <- results(dds,c("Stages","S5","S1"), filter = rowMedians(counts(dds)))
res_6vs1_med <- results(dds,c("Stages","S6","S1"), filter = rowMedians(counts(dds)))
res_7vs1_med <- results(dds,c("Stages","S7","S1"), filter = rowMedians(counts(dds)))
res_8vs1_med <- results(dds,c("Stages","S8","S1"), filter = rowMedians(counts(dds)))
res_4vs2_med <- results(dds,c("Stages","S4","S2"), filter = rowMedians(counts(dds)))
res_5vs2_med <- results(dds,c("Stages","S5","S2"), filter = rowMedians(counts(dds)))
res_6vs2_med <- results(dds,c("Stages","S6","S2"), filter = rowMedians(counts(dds)))
res_7vs2_med <- results(dds,c("Stages","S7","S2"), filter = rowMedians(counts(dds)))
res_8vs2_med <- results(dds,c("Stages","S8","S2"), filter = rowMedians(counts(dds)))
res_5vs3_med <- results(dds,c("Stages","S5","S3"), filter = rowMedians(counts(dds)))
res_6vs3_med <- results(dds,c("Stages","S6","S3"), filter = rowMedians(counts(dds)))
res_7vs3_med <- results(dds,c("Stages","S7","S3"), filter = rowMedians(counts(dds)))
res_8vs3_med <- results(dds,c("Stages","S8","S3"), filter = rowMedians(counts(dds)))
res_6vs4_med <- results(dds,c("Stages","S6","S4"), filter = rowMedians(counts(dds)))
res_7vs4_med <- results(dds,c("Stages","S7","S4"), filter = rowMedians(counts(dds)))
res_8vs4_med <- results(dds,c("Stages","S8","S4"), filter = rowMedians(counts(dds)))
res_7vs5_med <- results(dds,c("Stages","S7","S5"), filter = rowMedians(counts(dds)))
res_8vs5_med <- results(dds,c("Stages","S8","S5"), filter = rowMedians(counts(dds)))
res_8vs6_med <- results(dds,c("Stages","S8","S6"), filter = rowMedians(counts(dds)))

#>>>>>>> 



#'  # Number of DE genes in different stages of SE
#'  Extract names of all significantly DE genes in the experiment
#' ##all genes
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
  if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
  if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
  if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
} 

#all
res_list <- list(res_2vs1_med, res_3vs2_med, res_4vs3_med, res_5vs4_med, res_6vs5_med, res_7vs6_med, res_8vs7_med)
names(res_list) <- c("res_2vs1_med", "res_3vs2_med", "res_4vs3_med", "res_5vs4_med", "res_6vs5_med", "res_7vs6_med", "res_8vs7_med")
res_sig_list <- lapply(res_list, sigDeg)

# export R object with all significantly DEGs
#save(res_sig_list, file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DEGs_DESeqResults_median_p001_lfc05.Rda")

all <- lapply(res_sig_list, function(x){
  rownames(x)
})

#up
up <- lapply(res_sig_list, function(x){
  ab <- sigDeg(x, genes="up")
  rownames(ab)
})

#down
down <- lapply(res_sig_list, function(x){
  ab <- sigDeg(x, genes="down")
  rownames(ab)
})

nr_DEgenes <- cbind(elementNROWS(up), elementNROWS(down))
colnames(nr_DEgenes) <- c("up-regulated", "down-regulated")
rownames(nr_DEgenes) <- c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8")

barplot(t(nr_DEgenes),
        beside = TRUE,
        legend.text = TRUE,
        xlab = "Stage",
        ylab = "Number of DE genes",
        ylim = c(0,19000),
        col = pal12[c(2,3)],
        args.legend = list(bty = "n", x = "top")
)

barplot2(nr_DEgenes, 
         beside = TRUE, 
         legend.text = TRUE,
         xlab = "Stage",
         ylab = "Number of DE genes",
         ylim = c(0,20000))

#' with bigger font size
barplot(t(nr_DEgenes), 
        beside = TRUE, 
        col = pal12[c(3,2)], 
        ylim = c(0, 12000),
        cex.axis = 1.4, 
        cex.names = 1.4)
mtext(side=1, line=3, "Stages", cex=1.6)
mtext(side=2, line=2.5, "Number of DE genes", cex=1.6)
mtext(side=3, line=2, "Number of up- and down-regulated genes", font=2, cex=1.8)

legend("top", bty = "n",
       fill = pal12[c(3,2)],
       legend=c("up-regulated", "down-regulated"), cex = 1.4)

#'  # List of DE genes
# create a folder first
dir.create("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list")
# filter genes by padj <= 0.01, abs(log2FoldChange) >= 0.5, padj is not NA, if they are too many, write only the first 5000.
rownames2vs1_med <- rownames(res_2vs1_med[res_2vs1_med$padj <= 0.01 & abs(res_2vs1_med$log2FoldChange) >= 0.5 & !is.na(res_2vs1_med$padj),])
list_DE2vs1_med <- sub("\\.1$","",rownames2vs1_med)
write(list_DE2vs1_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE2vs1_med")

#res_3vs2o_med <- res_3vs2_med[order(res_3vs2_med$padj),][1:5000,]
rownames3vs2_med <-rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & abs(res_3vs2_med$log2FoldChange) >= 0.5 & !is.na(res_3vs2_med$padj),])
list_DE3vs2_med <- sub("\\.1$","",rownames3vs2_med)
write(list_DE3vs2_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE3vs2_med")

rownames4vs3_med <- rownames(res_4vs3_med[res_4vs3_med$padj <= 0.01 & abs(res_4vs3_med$log2FoldChange) >= 0.5 & !is.na(res_4vs3_med$padj),])
list_DE4vs3_med <- sub("\\.1$","",rownames4vs3_med)
write(list_DE4vs3_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE4vs3_med")

rownames5vs4_med <- rownames(res_5vs4_med[res_5vs4_med$padj <= 0.01 & abs(res_5vs4_med$log2FoldChange) >= 0.5 & !is.na(res_5vs4_med$padj),])
list_DE5vs4_med <- sub("\\.1$","",rownames5vs4_med)
write(list_DE5vs4_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE5vs4_med")

#res_6vs5o <- res_6vs5[order(res_6vs5$padj),][1:5000,]
rownames6vs5_med <- rownames(res_6vs5_med[res_6vs5_med$padj <= 0.01 & abs(res_6vs5_med$log2FoldChange) >= 0.5 & !is.na(res_6vs5_med$padj),])
list_DE6vs5_med <- sub("\\.1$","",rownames6vs5_med)
write(list_DE6vs5_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE6vs5_med")

#res_7vs6o <- res_7vs6[order(res_7vs6$padj),][1:5000,]
rownames7vs6_med <- rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & abs(res_7vs6_med$log2FoldChange) >= 0.5 & !is.na(res_7vs6_med$padj),])
list_DE7vs6_med <- sub("\\.1$","",rownames7vs6_med)
write(list_DE7vs6_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE7vs6_med")

rownames8vs7_med <- rownames(res_8vs7_med[res_8vs7_med$padj <= 0.01 & abs(res_8vs7_med$log2FoldChange) >= 0.5 & !is.na(res_8vs7_med$padj),])
list_DE8vs7_med <- sub("\\.1$","",rownames8vs7_med)
write(list_DE8vs7_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE8vs7_med")

#listDE_both_up 
rownames3vs2_up_med <- rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & res_3vs2_med$log2FoldChange >= 0.5 & !is.na(res_3vs2_med$padj),])
list_DE3vs2_up_med <- sub("\\.1$","",rownames3vs2_up_med)
write(list_DE3vs2_up_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE3vs2_up_med")

rownames7vs6_up_med <- rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & res_7vs6_med$log2FoldChange >= 0.5 & !is.na(res_7vs6_med$padj),])
list_DE7vs6_up_med <- sub("\\.1$","",rownames7vs6_up_med)
write(list_DE7vs6_up_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE7vs6_up_med")

listDE_both_up <- intersect(list_DE3vs2_up_med, list_DE7vs6_up_med)
write(listDE_both_up,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/listDE_both_up")

#listDE_both_down
rownames3vs2_down_med <- rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & res_3vs2_med$log2FoldChange <= -0.5 & !is.na(res_3vs2_med$padj),])
list_DE3vs2_down_med <- sub("\\.1$","",rownames3vs2_down_med)
write(list_DE3vs2_down_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE3vs2_down_med")

rownames7vs6_down_med <- rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & res_7vs6_med$log2FoldChange <= -0.5 & !is.na(res_7vs6_med$padj),])
list_DE7vs6_down_med <- sub("\\.1$","",rownames7vs6_down_med)
write(list_DE7vs6_down_med,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE7vs6_down_med")

listDE_both_down <- intersect(list_DE3vs2_down_med, list_DE7vs6_down_med)
write(listDE_both_down,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/listDE_both_down")
#3up and 7down
listDE3up_7down <- intersect(list_DE3vs2_up_med, list_DE7vs6_down_med)
write(listDE3up_7down,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/listDE3up_7down")
#3down 7up
listDE3down_7up <- intersect(list_DE3vs2_down_med, list_DE7vs6_up_med)
write(listDE3down_7up,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/listDE3down_7up")

#' #Gopher enrichment

list_DE2vs1_med.enr <- gopher(genes=list_DE2vs1_med, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
list_DE3vs2_med.enr <- gopher(genes=list_DE3vs2_med, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
list_DE4vs3_med.enr <- gopher(genes=list_DE4vs3_med, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
list_DE5vs4_med.enr <- gopher(genes=list_DE5vs4_med, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
list_DE6vs5_med.enr <- gopher(genes=list_DE6vs5_med, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
list_DE7vs6_med.enr <- gopher(genes=list_DE7vs6_med, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
list_DE8vs7_med.enr <- gopher(genes=list_DE8vs7_med, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
listDE_both_up.enr <- gopher(genes=listDE_both_up, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
listDE_both_down.enr <- gopher(genes=listDE_both_down, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
listDE3up_7down.enr <- gopher(genes=listDE3up_7down, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
listDE3down_7up.enr <- gopher(genes=listDE3down_7up, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")

treemap(listDE3down_7up.enr$go,
        index = c('namespace', 'name'), vSize = ('nt'),
        type = "categorical", vColor = 'name', title="",
        inflate.labels = FALSE, lowerbound.cex.labels = 0,
        bg.labels = "#CCCCCC00", 
        position.legend = "none"
)

listDE3down_7up.enr$mapman$name <- gsub("\\."," ", listDE3down_7up.enr$mapman$name)

treemap(listDE3down_7up.enr$mapman,
        index =  'name', vSize = ('nt'),
        type = "categorical", vColor = 'name', title="",
        inflate.labels = TRUE, lowerbound.cex.labels = 0,
        bg.labels = "#CCCCCC00", 
        position.legend = "none"
)


########
enr_list<- list(list_DE2vs1_med.enr, list_DE3vs2_med.enr, list_DE4vs3_med.enr, list_DE5vs4_med.enr, list_DE6vs5_med.enr, list_DE7vs6_med.enr, list_DE8vs7_med.enr)
GOlist <- lapply(enr_list, function(x){
  x$go$name
})
MapManlist <- lapply(enr_list, function(x){
  x$mapman$name
})

lapply(GOlist, function(x){
  x %in% c("chromatin organization", "chromatin") 
   })

lapply(MapManlist, function(x){
  x %in% c("chromatin organization", "chromatin") 
})

GO_chromatin <- grepl("chromatin", GOlist)
GO_chromatin_g <- grep("chromatin", list_DE3vs2_med.enr$go$name)
GO_methylation <- grepl("methylation",GOlist)
GO_methylation_g <- grep("methylation",list_DE7vs6_med.enr$go$name)
GO_methylation2 <- grepl("methyl*",GOlist)
GO_methylation2_g <- grep("methyl*",list_DE4vs3_med.enr$go$name)
                                    #list_DE6vs5_med.enr$go$name)
                                    #list_DE7vs6_med.enr$go$name)
chromatin_mapman <- grepl("chromatin", MapManlist)
chromatin_mapman_g <- grep("chromatin",list_DE4vs3_med.enr$mapman$name)
                                       #list_DE5vs4_med.enr$mapman$name)
                                       #list_DE8vs7_med.enr$mapman$name)
methylation_mapman <- grepl("methylation",MapManlist)
methylation_mapman_g <- grep("methylation", list_DE2vs1_med.enr$mapman$name)
                                           #list_DE5vs4_med.enr$mapman$name)
                                           #list_DE8vs7_med.enr$mapman$name)

request_chromatin <- httr::POST(paste0("https://microasp.upsc.se",":","5432","/","pabies","/","term-to-gene"),
                                            body = list(
                                               target=c(
                                                 list(
                                                   name = "go",
                                                   terms= c("GO:0006325", "GO:0000785")
                                                 )
                                               )
                                             ),
                                             encode = "json")
parsed_chromatin <- jsonlite::fromJSON(httr::content(request_chromatin, as = "text",
                                           encoding = "UTF-8"))

request_methylation <- httr::POST(paste0("https://microasp.upsc.se",":","5432","/","pabies","/","term-to-gene"),
                                body = list(
                                  target=c(
                                    list(
                                      name = "go",
                                      terms= c("GO:0008170","GO:0042054","GO:0032259")
                                    )
                                  )
                                ),
                                encode = "json")

parsed_methylation <- jsonlite::fromJSON(httr::content(request_methylation, as = "text",
                                                     encoding = "UTF-8"))


#' #VennDiagram

Venn_2vs1_med_3vs2_med <- plot.new()
grid.draw(venn.diagram(list(
  stage2vs1_med=rownames(res_2vs1_med[res_2vs1_med$padj <= 0.01 & abs(res_2vs1_med$log2FoldChange) >= 0.5 & !is.na(res_2vs1_med$padj),]),
  stage3vs2_med=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & abs(res_3vs2_med$log2FoldChange) >= 0.5 & !is.na(res_3vs2_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_3vs2_med_4vs3_med <- plot.new()
grid.draw(venn.diagram(list(
  stage3vs2_med=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & abs(res_3vs2_med$log2FoldChange) >= 0.5 & !is.na(res_3vs2_med$padj),]),
  stage4vs3_med=rownames(res_4vs3_med[res_4vs3_med$padj <= 0.01 & abs(res_4vs3_med$log2FoldChange) >= 0.5 & !is.na(res_4vs3_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_4vs3_med_5vs4_med <- plot.new()
grid.draw(venn.diagram(list(
  stage4vs3_med=rownames(res_4vs3_med[res_3vs2_med$padj <= 0.01 & abs(res_4vs3_med$log2FoldChange) >= 0.5 & !is.na(res_4vs3_med$padj),]),
  stage5vs4_med=rownames(res_5vs4_med[res_5vs4_med$padj <= 0.01 & abs(res_5vs4_med$log2FoldChange) >= 0.5 & !is.na(res_5vs4_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_5vs4_med_6vs5_med <- plot.new()
grid.draw(venn.diagram(list(
  stage5vs4_med=rownames(res_5vs4_med[res_5vs4_med$padj <= 0.01 & abs(res_5vs4_med$log2FoldChange) >= 0.5 & !is.na(res_5vs4_med$padj),]),
  stage6vs5_med=rownames(res_6vs5_med[res_6vs5_med$padj <= 0.01 & abs(res_6vs5_med$log2FoldChange) >= 0.5 & !is.na(res_6vs5_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_6vs5_med_7vs6_med <- plot.new()
grid.draw(venn.diagram(list(
  stage6vs5_med=rownames(res_6vs5_med[res_6vs5_med$padj <= 0.01 & abs(res_6vs5_med$log2FoldChange) >= 0.5 & !is.na(res_6vs5_med$padj),]),
  stage7vs6_med=rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & abs(res_7vs6_med$log2FoldChange) >= 0.5 & !is.na(res_7vs6_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_7vs6_med_8vs7_med <- plot.new()
grid.draw(venn.diagram(list(
  stage7vs6_med=rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & abs(res_7vs6_med$log2FoldChange) >= 0.5 & !is.na(res_7vs6_med$padj),]),
  stage8vs7_med=rownames(res_8vs7_med[res_8vs7_med$padj <= 0.01 & abs(res_8vs7_med$log2FoldChange) >= 0.5 & !is.na(res_8vs7_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_3vs2_med_7vs6_med <- plot.new()
grid.draw(venn.diagram(list(
  stage3vs2_med=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & abs(res_3vs2_med$log2FoldChange) >= 0.5 & !is.na(res_3vs2_med$padj),]),
  stage7vs6_med=rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & abs(res_7vs6_med$log2FoldChange) >= 0.5 & !is.na(res_7vs6_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_3vs2_med_7vs6_med_upregulated <- plot.new()
grid.draw(venn.diagram(list(
  stage3vs2_med=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & res_3vs2_med$log2FoldChange >= 0.5 & !is.na(res_3vs2_med$padj),]),
  stage7vs6_med=rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & res_7vs6_med$log2FoldChange >= 0.5 & !is.na(res_7vs6_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_3vs2_med_7vs6_med_downregulated <- plot.new()
grid.draw(venn.diagram(list(
  stage3vs2_med=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & res_3vs2_med$log2FoldChange <= -0.5 & !is.na(res_3vs2_med$padj),]),
  stage7vs6_med=rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & res_7vs6_med$log2FoldChange <= -0.5 & !is.na(res_7vs6_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_3vs2_7vs6_u3d7_med <- plot.new
grid.draw(venn.diagram(list(
  stage3vs2_med=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & res_3vs2_med$log2FoldChange >= 0.5 & !is.na(res_3vs2_med$padj),]),
  stage7vs6_med=rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & res_7vs6_med$log2FoldChange <= -0.5 & !is.na(res_7vs6_med$padj),])),
  filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

Venn_3vs2_7vs6_d3u7_med <- plot.new()
grid.draw(venn.diagram(list(
  stage3vs2_med=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & res_3vs2_med$log2FoldChange <= -0.5 & !is.na(res_3vs2_med$padj),]),
  stage7vs6_med=rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & res_7vs6_med$log2FoldChange >= 0.5 & !is.na(res_7vs6_med$padj),])),
  filename= NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))


par(mar=c(0.1,0.1,0.1,0.1))
Venn_3vs2_med_7vs6_med_all4 <- plot.new()
grid.draw(venn.diagram(list(
  stage3vs2u=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & res_3vs2_med$log2FoldChange >= 0.5 & !is.na(res_3vs2_med$padj),]),
  stage3vs2d=rownames(res_3vs2_med[res_3vs2_med$padj <= 0.01 & res_3vs2_med$log2FoldChange <= -0.5 & !is.na(res_3vs2_med$padj),]),
  stage7vs6u=rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & res_7vs6_med$log2FoldChange >= 0.5 & !is.na(res_7vs6_med$padj),]),
  stage7vs6d=rownames(res_7vs6_med[res_7vs6_med$padj <= 0.01 & res_7vs6_med$log2FoldChange <= -0.5 & !is.na(res_7vs6_med$padj),])),
  filename= NULL, fill=pal[1:4])) 
#cex=1.5,height=1000, width=1500, resolution=500, imagetype="tiff", units="px"))

#'  # Plot the Median vs Average
DESeq2::plotMA(res_2vs1_med,0.01)
DESeq2::plotMA(res_3vs2_med,0.01)
DESeq2::plotMA(res_4vs3_med,0.01)
DESeq2::plotMA(res_5vs4_med,0.01)
DESeq2::plotMA(res_6vs5_med,0.01)
DESeq2::plotMA(res_7vs6_med,0.01)
DESeq2::plotMA(res_8vs7_med,0.01)

#' # The volcano plot shows the same results as the MA plot; a
# large number of genes show significant fold-changes
alpha=0.01
lfc=0.5
volcanoPlot(res_2vs1_med,alpha=alpha,lfc = lfc)
volcanoPlot(res_3vs2_med,alpha=alpha,lfc = lfc)
volcanoPlot(res_4vs3_med,alpha=alpha,lfc = lfc)
volcanoPlot(res_5vs4_med,alpha=alpha,lfc = lfc)
volcanoPlot(res_6vs5_med,alpha=alpha,lfc = lfc)
volcanoPlot(res_7vs6_med,alpha=alpha,lfc = lfc)
volcanoPlot(res_8vs7_med,alpha=alpha,lfc = lfc)

#' #### Soft clustering (Mfuzz)

#' #Create the eSet

vst.g <- read.table("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/vala/data/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection.tsv",
                    row.names=1)

alpha=0.01
lfc=0.5

#for all results


sel <- unique(unlist(all))
eset <- ExpressionSet(as.matrix(vst.g[rownames(vst.g) %in% sel,]))

#' Standardise
eset.s <- standardise(eset)

#' Average the replicates (mean)
samples <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3.csv")
eset.m <- ExpressionSet(
  sapply(split.data.frame(t(exprs(eset.s)),
                          samples$Stages),
         colMeans))
#' Estimate the fuzzification
m <- mestimate(eset.m)

#' Find the clusters (8 is based on the previous heatmap)
cl <- mfuzz(eset.m,centers=24,m=m)

#' There are a number of clusters that behave similarly
par(mar=c(0.1,0.1,0.1,0.1))
mfuzz.plot(eset.m,
           cl=cl,
           mfrow=c(6,4),
           time.labels = colnames(eset.m),
           new.window=FALSE)
par(mar = mar)


############
#' # Differential expression of all stages
 
#' Define pairs of stages to compare
names_stages <- unique(levels(dds$Stages))
# reverse the order of the stages, so you will always compare a later stage to an earlier stage in the experiment,
# when comparing expression of the DEgenes
comparisons <- combn(rev(names_stages), 2)

#' Extract DEgenes
res_med <- apply(comparisons, 2, function(f){
  r <- results(dds, contrast = c("Stages", f), filter = rowMedians(counts(dds)))
})

#' Filter them by log2fc and padj  
# The function to select the significantly DEGs (FDR < 0.01, |log2FC| => 0.5)
sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
  if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
  if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
  if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
} 

names(res_med) <- apply(comparisons, 2, paste, collapse="~")
res_sig_med <- lapply(res_med, sigDeg)

#' Construct a matrix with number of up- and down-regulated genes (comparisons of all the stages)  
# construct the matrix
de_mat <- matrix(NA, 
                 nrow = length(names_stages), 
                 ncol = length(names_stages))
# name the rows and columns
dimnames(de_mat) <- list(names_stages, names_stages)
# calculate number of up- and down-regulated genes
de_mat[t(comparisons)] <- unlist(lapply(res_sig_med, function(f){ - sum(f$log2FoldChange < 0, na.rm = TRUE) }))
de_mat[t(comparisons)[, c(2, 1)]] <- unlist(lapply(res_sig_med, function(f){ sum(f$log2FoldChange > 0, na.rm = TRUE) }))

# print the heatmap with number of DEGs between all the stages
pdf("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/nrDEgenes_AllvsAllStages.pdf", width = 12, height = 9)
pheatmap(de_mat, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 18)
dev.off()



#' ```{r empty,eval=FALSE,echo=FALSE}
#' ```
#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
