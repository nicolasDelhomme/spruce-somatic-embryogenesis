#' ---
#' title: "Spruce-somatic-embryogenesis Differential Expression"
#' author: "Nicolas Delhomme and Camilla Canovi"
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


#' Helper files
suppressMessages(source("~/Git/UPSCb/src/R/densityPlot.R"))
suppressMessages(source("~/Git/UPSCb/src/R/plotMA.R"))
suppressMessages(source("~/Git/UPSCb/src/R/volcanoPlot.R"))
suppressMessages(source("~/Git/UPSCb/src/R/plot.multidensity.R"))

#' Load DESeq object
load("DESeqDataSet-3.rda")

#' Annotation
# annot <- read.delim("/mnt/picea/storage/reference/Picea-abies/v1.1/annotation/Pabies-Pfam-description-annotation.tsv")
# colnames(annot) <- c("GeneID","Confidence","PfamDomains")

#' Setup graphics
pal=brewer.pal(8,"Dark2")
pal12 <- brewer.pal(12,"Paired")
hpal <- colorRampPalette(c("blue","white","red"))(100)
mar <- par("mar")

#' ```{r drop limma,echo=FALSE}
#' detach("package:limma")
#' ```
#' 
#' # Process

#' # Differential Expression
  dds <- DESeq(dds.counts)
  
  # Dispersion Estimation
  #
  # The dispersion estimation is adequate
  plotDispEsts(dds)
  
#'  ## Obtain the results 
   resultsNames(dds)
  
   res_2vs1 <- results(dds,c("Stages","2","1"))
   res_3vs2 <- results(dds,c("Stages","3","2"))
   res_4vs3 <- results(dds,c("Stages","4","3"))
   res_5vs4 <- results(dds,c("Stages","5","4"))
   res_6vs5 <- results(dds,c("Stages","6","5"))
   res_7vs6 <- results(dds,c("Stages","7","6"))
   res_8vs7 <- results(dds,c("Stages","8","7"))
  
#' ## Filter results  
#' by log2fc and FDR  

#'  # List of DE genes
# create a folder first
   dir.create("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list")
# filter genes by padj <= 0.01, abs(log2FoldChange) >= 0.5, padj is not NA, if they are too many, write only the first 5000.
  rownames2vs1 <- rownames(res_2vs1[res_2vs1$padj <= 0.01 & abs(res_2vs1$log2FoldChange) >= 0.5 & !is.na(res_2vs1$padj),])
  list_DE2vs1 <- sub("\\.1$",",",rownames2vs1)
  write(list_DE2vs1,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE2vs1")

  res_3vs2o <- res_3vs2[order(res_3vs2$padj),][1:5000,]
  rownames3vs2 <-rownames(res_3vs2o[res_3vs2o$padj <= 0.01 & abs(res_3vs2o$log2FoldChange) >= 0.5 & !is.na(res_3vs2o$padj),])
  list_DE3vs2 <- sub("\\.1$",",",rownames3vs2)
  write(list_DE3vs2,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE3vs2")

  rownames4vs3 <- rownames(res_4vs3[res_4vs3$padj <= 0.01 & abs(res_4vs3$log2FoldChange) >= 0.5 & !is.na(res_4vs3$padj),])
  list_DE4vs3 <- sub("\\.1$",",",rownames4vs3)
  write(list_DE4vs3,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE4vs3")

  rownames5vs4 <- rownames(res_5vs4[res_5vs4$padj <= 0.01 & abs(res_5vs4$log2FoldChange) >= 0.5 & !is.na(res_5vs4$padj),])
  list_DE5vs4 <- sub("\\.1$",",",rownames5vs4)
  write(list_DE5vs4,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE5vs4")

  res_6vs5o <- res_6vs5[order(res_6vs5$padj),][1:5000,]
  rownames6vs5 <- rownames(res_6vs5o[res_6vs5o$padj <= 0.01 & abs(res_6vs5o$log2FoldChange) >= 0.5 & !is.na(res_6vs5o$padj),])
  list_DE6vs5 <- sub("\\.1$",",",rownames6vs5)
  write(list_DE6vs5,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE6vs5")

  res_7vs6o <- res_7vs6[order(res_7vs6$padj),][1:5000,]
  rownames7vs6 <- rownames(res_7vs6o[res_7vs6o$padj <= 0.01 & abs(res_7vs6o$log2FoldChange) >= 0.5 & !is.na(res_7vs6o$padj),])
  list_DE7vs6 <- sub("\\.1$",",",rownames7vs6)
  write(list_DE7vs6,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE7vs6")

  rownames8vs7 <- rownames(res_8vs7[res_8vs7$padj <= 0.01 & abs(res_8vs7$log2FoldChange) >= 0.5 & !is.na(res_8vs7$padj),])
  list_DE8vs7 <- sub("\\.1$",",",rownames8vs7)
  write(list_DE8vs7,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE8vs7")


  #'  # Plot the Median vs Average
  DESeq2::plotMA(res_2vs1,0.01)
  DESeq2::plotMA(res_3vs2,0.01)
  DESeq2::plotMA(res_4vs3,0.01)
  DESeq2::plotMA(res_5vs4,0.01)
  DESeq2::plotMA(res_6vs5,0.01)
  DESeq2::plotMA(res_7vs6,0.01)
  DESeq2::plotMA(res_8vs7,0.01)

  #'  # Number of DE genes in different stages of SE
  #'  Extract names of all significantly DE genes in the experiment
  #' ##all genes
  alpha=0.01
  lfc=0.5

  allres <- lapply(2:8,function(i){
    res <- results(dds,c("Stages",as.character(i),as.character(i-1)))
    return(rownames(res[res$padj <= alpha & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) ,]))
  })

  #' ##upregulated genes
   allres_upregulated <- lapply(2:8,function(i){
    res <- results(dds,c("Stages",as.character(i),as.character(i-1)))
    return(rownames(res[res$padj <= alpha & res$log2FoldChange >= lfc & ! is.na(res$padj) ,]))
    })

   #' ##downregulated genes
  allres_downregulated <- lapply(2:8,function(i){
    res <- results(dds,c("Stages",as.character(i),as.character(i-1)))
    return(rownames(res[res$padj <= alpha & res$log2FoldChange <= -lfc & ! is.na(res$padj) ,]))
  })

  nr_DEgenes <- cbind(elementNROWS(allres_downregulated), elementNROWS(allres_upregulated))
  colnames(nr_DEgenes) <- c("down-regulated", "up-regulated")
  rownames(nr_DEgenes) <- c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8")

  barplot(t(nr_DEgenes),
          beside = TRUE,
          legend.text = TRUE,
          xlab = "Stage",
          ylab = "Number of DE genes",
          ylim = c(0,12000),
          col = pal12[c(2,3)],
          args.legend = list(x = "topleft")
          )

  barplot2(nr_DEgenes, beside = TRUE, legend.text = TRUE,xlab = "Stage", ylab = "Number of DE genes", ylim = c(0,20000))

  #' Prepare lists of differentially expressed genes with annotation

  #' Create function:
  #' Select significant genes (FDR < 0.01, |log2FC| => 0.5)
  sigDeg <- function(res, p = 0.01, log2fc = 0.5, genes = "all") {
    if(genes == "all") return(res[res$padj<p & !is.na(res$padj) & abs(res$log2FoldChange) >= log2fc,])
    if(genes == "up") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange >= log2fc,])
    if(genes == "down") return(res[res$padj<p & !is.na(res$padj) & res$log2FoldChange <= -log2fc,])
  }

  #' Load gene annotation
  gene_table <- read.delim("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/Pabies_gene_table_conf-PFAM-GO-bestAra.tsv",
                                stringsAsFactors = FALSE)

  annotated_DEGs <- lapply(list(res_2vs1, res_3vs2, res_4vs3, res_5vs4, res_6vs5, res_7vs6, res_8vs7), function(x) {
    sigx <- sigDeg(x, genes = "all")
    sigx <- sigx[order(sigx$padj), ]
    rownames(sigx) <- gsub(".1$", "", rownames(sigx))
    cbind(sigx[ , c("baseMean", "log2FoldChange", "padj")], gene_table[match(rownames(sigx), gene_table$Gene), -1])
  })

  names(annotated_DEGs) <- c("DEGs_2vs1", "DEGs_3vs2", "DEGs_4vs3", "DEGs_5vs4", "DEGs_6vs5", "DEGs_7vs6", "DEGs_8vs7")

  dev.null <- lapply(names(annotated_DEGs), function(x) {
    write.table(annotated_DEGs[[x]],
                file = file.path("doc/DE_list", paste0(x, "_one-percent-FDR-05-log2fc-cutoff_significant-genes-annotation.tsv")),
                col.names = NA,
                sep = "\t")
  })



  #'# TEs
  #' Load DESeq object
  load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/kallisto_LTR-TEs/DESeqDataSet-3_TEs_R_blindf.rda")

  # Differential Expression
  dds_TEs <- DESeq(dds.counts_TEs)

  # Dispersion Estimation
  #
  # The dispersion estimation is adequate
  plotDispEsts(dds_TEs)

  #'  ## Obtain the results
  resultsNames(dds_TEs)

  res_2vs1_TEs <- results(dds_TEs,c("Stages","2","1"))
  res_3vs2_TEs <- results(dds_TEs,c("Stages","3","2"))
  res_4vs3_TEs <- results(dds_TEs,c("Stages","4","3"))
  res_5vs4_TEs <- results(dds_TEs,c("Stages","5","4"))
  res_6vs5_TEs <- results(dds_TEs,c("Stages","6","5"))
  res_7vs6_TEs <- results(dds_TEs,c("Stages","7","6"))
  res_8vs7_TEs <- results(dds_TEs,c("Stages","8","7"))

  allres_TEs <- lapply(2:8,function(i){
    res_TEs <- results(dds_TEs,c("Stages",as.character(i),as.character(i-1)))
    return(rownames(res_TEs[res_TEs$padj <= alpha & abs(res_TEs$log2FoldChange) >= lfc & ! is.na(res_TEs$padj) ,]))
  })

  #' ##upregulated genes
  allres_upregulated_TEs <- lapply(2:8,function(i){
    res_TEs <- results(dds_TEs,c("Stages",as.character(i),as.character(i-1)))
    return(rownames(res_TEs[res_TEs$padj <= alpha & res_TEs$log2FoldChange >= lfc & ! is.na(res_TEs$padj) ,]))
  })

  #' ##downregulated genes
  allres_downregulated_TEs <- lapply(2:8,function(i){
    res_TEs <- results(dds_TEs,c("Stages",as.character(i),as.character(i-1)))
    return(rownames(res_TEs[res_TEs$padj <= alpha & res_TEs$log2FoldChange <= -lfc & ! is.na(res_TEs$padj) ,]))
  })

  nr_DEgenes <- cbind(elementNROWS(allres_downregulated), elementNROWS(allres_upregulated), elementNROWS(allres_downregulated_TEs), elementNROWS(allres_upregulated_TEs))
  colnames(nr_DEgenes) <- c("downregulated", "upregulated", "downregulated_TEs", "upregulated_TEs")
  rownames(nr_DEgenes) <- c("1-2", "2-3", "3-4", "4-5", "5-6", "6-7", "7-8")

  barplot2(t(nr_DEgenes), beside = TRUE, legend.text = TRUE, xlab = "Stage", ylab = "DE genes", log = "y")
  barplot2(nr_DEgenes, beside = TRUE, legend.text = TRUE,xlab = "Stage", ylab = "DE genes", log = "y")



  #' #VennDiagram
  library(VennDiagram)

  Venn_2vs1_3vs2 <- plot.new()
  grid.draw(venn.diagram(list(
    stage2vs1=rownames(res_2vs1[res_2vs1$padj <= 0.01 & abs(res_2vs1$log2FoldChange) >= 0.5 & !is.na(res_2vs1$padj),]),
    stage3vs2=rownames(res_3vs2[res_3vs2$padj <= 0.01 & abs(res_3vs2$log2FoldChange) >= 0.5 & !is.na(res_3vs2$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

  Venn_3vs2_4vs3 <- plot.new()
  grid.draw(venn.diagram(list(
    stage3vs2=rownames(res_3vs2[res_3vs2$padj <= 0.01 & abs(res_3vs2$log2FoldChange) >= 0.5 & !is.na(res_3vs2$padj),]),
    stage4vs3=rownames(res_4vs3[res_4vs3$padj <= 0.01 & abs(res_4vs3$log2FoldChange) >= 0.5 & !is.na(res_4vs3$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

  Venn_4vs3_5vs4 <- plot.new()
  grid.draw(venn.diagram(list(
    stage4vs3=rownames(res_4vs3[res_3vs2$padj <= 0.01 & abs(res_4vs3$log2FoldChange) >= 0.5 & !is.na(res_4vs3$padj),]),
    stage5vs4=rownames(res_5vs4[res_5vs4$padj <= 0.01 & abs(res_5vs4$log2FoldChange) >= 0.5 & !is.na(res_5vs4$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

  Venn_5vs4_6vs5 <- plot.new()
  grid.draw(venn.diagram(list(
    stage5vs4=rownames(res_5vs4[res_5vs4$padj <= 0.01 & abs(res_5vs4$log2FoldChange) >= 0.5 & !is.na(res_5vs4$padj),]),
    stage6vs5=rownames(res_6vs5[res_6vs5$padj <= 0.01 & abs(res_6vs5$log2FoldChange) >= 0.5 & !is.na(res_6vs5$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

  Venn_6vs5_7vs6 <- plot.new()
  grid.draw(venn.diagram(list(
    stage6vs5=rownames(res_6vs5[res_6vs5$padj <= 0.01 & abs(res_6vs5$log2FoldChange) >= 0.5 & !is.na(res_6vs5$padj),]),
    stage7vs6=rownames(res_7vs6[res_7vs6$padj <= 0.01 & abs(res_7vs6$log2FoldChange) >= 0.5 & !is.na(res_7vs6$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

  Venn_7vs6_8vs7 <- plot.new()
  grid.draw(venn.diagram(list(
    stage7vs6=rownames(res_7vs6[res_7vs6$padj <= 0.01 & abs(res_7vs6$log2FoldChange) >= 0.5 & !is.na(res_7vs6$padj),]),
    stage8vs7=rownames(res_8vs7[res_8vs7$padj <= 0.01 & abs(res_8vs7$log2FoldChange) >= 0.5 & !is.na(res_8vs7$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

  Venn_3vs2_7vs6 <- plot.new()
  grid.draw(venn.diagram(list(
    stage3vs2=rownames(res_3vs2[res_3vs2$padj <= 0.01 & abs(res_3vs2$log2FoldChange) >= 0.5 & !is.na(res_3vs2$padj),]),
    stage7vs6=rownames(res_7vs6[res_7vs6$padj <= 0.01 & abs(res_7vs6$log2FoldChange) >= 0.5 & !is.na(res_7vs6$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

  Venn_3vs2_7vs6_upregulated <- plot.new()
  grid.draw(venn.diagram(list(
    stage3vs2=rownames(res_3vs2[res_3vs2$padj <= 0.01 & res_3vs2$log2FoldChange >= 0.5 & !is.na(res_3vs2$padj),]),
    stage7vs6=rownames(res_7vs6[res_7vs6$padj <= 0.01 & res_7vs6$log2FoldChange >= 0.5 & !is.na(res_7vs6$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

   Venn_3vs2_7vs6_downregulated <- plot.new()
  grid.draw(venn.diagram(list(
    stage3vs2=rownames(res_3vs2[res_3vs2$padj <= 0.01 & res_3vs2$log2FoldChange <= -0.5 & !is.na(res_3vs2$padj),]),
    stage7vs6=rownames(res_7vs6[res_7vs6$padj <= 0.01 & res_7vs6$log2FoldChange <= -0.5 & !is.na(res_7vs6$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

  Venn_3vs2_7vs6_u3d7 <- plot.new
  grid.draw(venn.diagram(list(
    stage3vs2=rownames(res_3vs2[res_3vs2$padj <= 0.01 & res_3vs2$log2FoldChange >= 0.5 & !is.na(res_3vs2$padj),]),
    stage7vs6=rownames(res_7vs6[res_7vs6$padj <= 0.01 & res_7vs6$log2FoldChange <= -0.5 & !is.na(res_7vs6$padj),])),
    filename=NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))

  Venn_3vs2_7vs6_d3u7 <- plot.new()
  grid.draw(venn.diagram(list(
    stage3vs2=rownames(res_3vs2[res_3vs2$padj <= 0.01 & res_3vs2$log2FoldChange <= -0.5 & !is.na(res_3vs2$padj),]),
    stage7vs6=rownames(res_7vs6[res_7vs6$padj <= 0.01 & res_7vs6$log2FoldChange >= 0.5 & !is.na(res_7vs6$padj),])),
    filename= NULL, fill=pal[1:2])) #col = "red", label.col = "green", cat.col = "blue"))


  #rownames3vs2_7vs6_u3d7 <- rownames(res_3vs2[res_3vs2$padj <= 0.01 & res_3vs2$log2FoldChange >= 0.5 & !is.na(res_3vs2$padj)],
  #                          rownames(res_7vs6[res_7vs6$padj <= 0.01 & res_7vs6$log2FoldChange <= -0.5 & !is.na(res_7vs6$padj)]))

  #list_DE3vs2_7vs6_u3d7 <- sub("\\.1$",",",rownames3vs2_7vs6_u3d7)
  #write(list_DE2vs1,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DE_list/list_DE2vs1")


  # Plot the log10 odds (i.e. -log10 FDR) vs. log2 fold change
  #
  #' # The volcano plot shows the same results as the MA plot; a
  # large number of genes show significant fold-changes
  alpha=0.01
  lfc=0.5
  volcanoPlot(res_2vs1,alpha=alpha,lfc = lfc)
  volcanoPlot(res_3vs2,alpha=alpha,lfc = lfc)
  volcanoPlot(res_4vs3,alpha=alpha,lfc = lfc)
  volcanoPlot(res_5vs4,alpha=alpha,lfc = lfc)
  volcanoPlot(res_6vs5,alpha=alpha,lfc = lfc)
  volcanoPlot(res_7vs6,alpha=alpha,lfc = lfc)
  volcanoPlot(res_8vs7,alpha=alpha,lfc = lfc)

  #' #### Soft clustering (Mfuzz)
  #' !! Each time you run it, the results will be a bit different. Use stored data
  #' from the first experiment to do any additional analysis.
  suppressPackageStartupMessages(library(Mfuzz))
  suppressPackageStartupMessages(library(org.At.tair.db))

  #' #Create the eSet
  # Comment (Katja): not found
  # vst.g <- read.table("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection.tsv",
  #                    row.names=1)
  # the file with the same name exists here:
  vst.g <- read.table("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/vala/data/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection.tsv",
                    row.names=1)

  alpha=0.01
  lfc=0.5

  #for all results
  allres <- lapply(2:8,function(i){
    res <- results(dds,c("Stages",as.character(i),as.character(i-1)))
    return(rownames(res[res$padj <= alpha & abs(res$log2FoldChange) >= lfc & ! is.na(res$padj) ,]))
  })

  sel <- unique(unlist(allres))
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

  cluster11 <- names(cl$cluster[cl$cluster == 11])
  cluster12 <- names(cl$cluster[cl$cluster == 12])
  cluster18 <- names(cl$cluster[cl$cluster == 18])

  list_cluster11 <- sub("\\.1$",",",cluster11)
  write(list_cluster11,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster11")

  list_cluster12 <- sub("\\.1$",",",cluster12)
  write(list_cluster12,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster12")

  list_cluster18 <- sub("\\.1$",",",cluster18)
  write(list_cluster18,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster18")

  cluster1 <- names(cl$cluster[cl$cluster == 1])
  cluster2 <- names(cl$cluster[cl$cluster == 2])
  cluster3 <- names(cl$cluster[cl$cluster == 3])
  cluster4 <- names(cl$cluster[cl$cluster == 4])
  cluster5 <- names(cl$cluster[cl$cluster == 5])
  cluster6 <- names(cl$cluster[cl$cluster == 6])
  cluster7 <- names(cl$cluster[cl$cluster == 7])
  cluster8 <- names(cl$cluster[cl$cluster == 8])
  cluster9 <- names(cl$cluster[cl$cluster == 9])
  cluster10 <- names(cl$cluster[cl$cluster == 10])
  cluster13 <- names(cl$cluster[cl$cluster == 13])
  cluster14 <- names(cl$cluster[cl$cluster == 14])
  cluster15 <- names(cl$cluster[cl$cluster == 15])
  cluster16 <- names(cl$cluster[cl$cluster == 16])
  cluster17 <- names(cl$cluster[cl$cluster == 17])
  cluster19 <- names(cl$cluster[cl$cluster == 19])
  cluster20 <- names(cl$cluster[cl$cluster == 20])
  cluster21 <- names(cl$cluster[cl$cluster == 21])
  cluster22 <- names(cl$cluster[cl$cluster == 22])
  cluster23 <- names(cl$cluster[cl$cluster == 23])
  cluster24 <- names(cl$cluster[cl$cluster == 24])


  list_cluster1 <- sub("\\.1$",",",cluster1)
  write(list_cluster1,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster1")

  list_cluster2 <- sub("\\.1$",",",cluster2)
  write(list_cluster2,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster2")

  list_cluster3 <- sub("\\.1$",",",cluster3)
  write(list_cluster3,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster3")

  list_cluster4 <- sub("\\.1$",",",cluster4)
  write(list_cluster4,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster4")

  list_cluster5 <- sub("\\.1$",",",cluster5)
  write(list_cluster5,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster5")

  list_cluster6 <- sub("\\.1$",",",cluster6)
  write(list_cluster6,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster6")

  list_cluster7 <- sub("\\.1$",",",cluster7)
  write(list_cluster7,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster7")

  list_cluster8 <- sub("\\.1$",",",cluster8)
  write(list_cluster8,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster8")

  list_cluster9 <- sub("\\.1$",",",cluster9)
  write(list_cluster9,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster9")

  list_cluster10 <- sub("\\.1$",",",cluster10)
  write(list_cluster10,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster10")

  list_cluster13 <- sub("\\.1$",",",cluster13)
  write(list_cluster13,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster13")

  list_cluster14 <- sub("\\.1$",",",cluster14)
  write(list_cluster14,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster14")

  list_cluster15 <- sub("\\.1$",",",cluster15)
  write(list_cluster15,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster15")

  list_cluster16 <- sub("\\.1$",",",cluster16)
  write(list_cluster16,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster16")

  list_cluster17 <- sub("\\.1$",",",cluster17)
  write(list_cluster17,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster17")

  list_cluster19 <- sub("\\.1$",",",cluster19)
  write(list_cluster19,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster19")

  list_cluster20 <- sub("\\.1$",",",cluster20)
  write(list_cluster20,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster20")

  list_cluster21 <- sub("\\.1$",",",cluster21)
  write(list_cluster21,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster21")

  list_cluster22 <- sub("\\.1$",",",cluster22)
  write(list_cluster22,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster22")

  list_cluster23 <- sub("\\.1$",",",cluster23)
  write(list_cluster23,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster23")

  list_cluster24 <- sub("\\.1$",",",cluster24)
  write(list_cluster24,file = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster24")

  #' Size of the clusters
  clusters <- dir(path = "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/", full.names = TRUE, pattern = "list_cluster[0-9]+")
  cluster_list <- lapply(clusters, FUN = read.csv)
  cluster_size <- lapply(cluster_list, FUN = nrow)
  names(cluster_size) <- str_extract(clusters, "[0-9]+")
  cluster_size <- unlist(cluster_size)

  barplot(cluster_size[order(cluster_size)],
          col = "skyblue3",
          ylim = c(0, round(max(cluster_size), -2)+100),
          #ylim = c(0, 1400),
          ylab = "number of genes in the cluster",
          xlab = "clusters")

  abline(h = c(800, 1000, 1200), col = "grey")
  barplot(cluster_size[order(cluster_size)],
          col = "skyblue3",
          ylim = c(0, round(max(cluster_size), -2)+100),
          #ylim = c(0, 1400),
          ylab = "number of genes in the cluster",
          xlab = "clusters",
          add = TRUE)
