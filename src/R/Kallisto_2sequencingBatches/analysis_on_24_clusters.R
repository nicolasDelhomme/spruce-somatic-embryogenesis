
#' Set the working dir
setwd("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")
#' ```


library(VennDiagram)
suppressPackageStartupMessages(library(RColorBrewer))
pal=brewer.pal(8,"Dark2")

#' K: do they come from? ConGenIE list?
spruce.specific <- scan("analysis/spruce-specific/putative-specific-gene-IDs.txt",what="character")
#all_expressed_genes <- scan("analysis/vala/data/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection.tsv",what="character", sep = "\t")

all_expressed_genes_2 <- read.table("analysis/vala/data/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection.tsv", 
                                    sep = "\t", header = TRUE, row.names = 1)
all_expressed_genes_2_MA <- sub("\\.1$","", rownames(all_expressed_genes_2))

#all_expressed_genes_3 <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3.csv",
                              #      sep = "\t", header = TRUE, row.names = 1)
#rownames(all_expressed_genes_3)
spruce_specific_expressed <- plot.new()
grid.draw(venn.diagram(list(
  area1 =spruce.specific,
  area2= all_expressed_genes_2_MA),
  filename = NULL, fill=pal[1:2]))



#' K: When were these files created? Why is a number of genes different then in 
#' "~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/list_cluster/list_cluster1",
#' when they were created?
cluster1.genes <- scan("analysis/spruce-specific/list_cluster1",what="character")
cluster1.genes <- sub(",","", cluster1.genes)
cluster2.genes <- scan("analysis/spruce-specific/list_cluster2",what="character")
cluster2.genes <- sub(",","", cluster2.genes)
cluster3.genes <- scan("analysis/spruce-specific/list_cluster3",what="character")
cluster3.genes <- sub(",","", cluster3.genes)
cluster4.genes <- scan("analysis/spruce-specific/list_cluster4",what="character")
cluster4.genes <- sub(",","", cluster4.genes)
cluster5.genes <- scan("analysis/spruce-specific/list_cluster5",what="character")
cluster5.genes <- sub(",","", cluster5.genes)
cluster6.genes <- scan("analysis/spruce-specific/list_cluster6",what="character")
cluster6.genes <- sub(",","", cluster6.genes)
cluster7.genes <- scan("analysis/spruce-specific/list_cluster7",what="character")
cluster7.genes <- sub(",","", cluster7.genes)
cluster8.genes <- scan("analysis/spruce-specific/list_cluster8",what="character")
cluster8.genes <- sub(",","", cluster8.genes)
cluster9.genes <- scan("analysis/spruce-specific/list_cluster9",what="character")
cluster9.genes <- sub(",","", cluster9.genes)
cluster10.genes <- scan("analysis/spruce-specific/list_cluster10",what="character")
cluster10.genes <- sub(",","", cluster10.genes)
cluster11.genes <- scan("analysis/spruce-specific/list_cluster11",what="character")
cluster11.genes <- sub(",","", cluster11.genes)
cluster12.genes <- scan("analysis/spruce-specific/list_cluster12",what="character")
cluster12.genes <- sub(",","", cluster12.genes)
cluster13.genes <- scan("analysis/spruce-specific/list_cluster13",what="character")
cluster13.genes <- sub(",","", cluster13.genes)
cluster14.genes <- scan("analysis/spruce-specific/list_cluster14",what="character")
cluster14.genes <- sub(",","", cluster14.genes)
cluster15.genes <- scan("analysis/spruce-specific/list_cluster15",what="character")
cluster15.genes <- sub(",","", cluster15.genes)
cluster16.genes <- scan("analysis/spruce-specific/list_cluster16",what="character")
cluster16.genes <- sub(",","", cluster16.genes)
cluster17.genes <- scan("analysis/spruce-specific/list_cluster17",what="character")
cluster17.genes <- sub(",","", cluster17.genes)
cluster18.genes <- scan("analysis/spruce-specific/list_cluster18",what="character")
cluster18.genes <- sub(",","", cluster18.genes)
cluster19.genes <- scan("analysis/spruce-specific/list_cluster19",what="character")
cluster19.genes <- sub(",","", cluster19.genes)
cluster20.genes <- scan("analysis/spruce-specific/list_cluster20",what="character")
cluster20.genes <- sub(",","", cluster20.genes)
cluster21.genes <- scan("analysis/spruce-specific/list_cluster21",what="character")
cluster21.genes <- sub(",","", cluster21.genes)
cluster22.genes <- scan("analysis/spruce-specific/list_cluster22",what="character")
cluster22.genes <- sub(",","", cluster22.genes)
cluster23.genes <- scan("analysis/spruce-specific/list_cluster23",what="character")
cluster23.genes <- sub(",","", cluster23.genes)
cluster24.genes <- scan("analysis/spruce-specific/list_cluster24",what="character")
cluster24.genes <- sub(",","", cluster24.genes)

          
cluster.genes_list <- list(cluster1.genes, cluster2.genes, cluster3.genes, cluster4.genes, cluster5.genes, cluster6.genes, cluster7.genes, cluster8.genes,
                           cluster9.genes, cluster10.genes, cluster11.genes, cluster12.genes, cluster13.genes, cluster14.genes, cluster15.genes, cluster16.genes,
                           cluster17.genes, cluster18.genes, cluster19.genes, cluster20.genes, cluster21.genes, cluster22.genes, cluster23.genes, cluster24.genes)
names(cluster.genes_list) <- c("cluster1.genes", "cluster2.genes", "cluster3.genes", "cluster4.genes", "cluster5.genes", "cluster6.genes", "cluster7.genes", "cluster8.genes",
                               "cluster9.genes", "cluster10.genes", "cluster11.genes", "cluster12.genes", "cluster13.genes", "cluster14.genes", "cluster15.genes", "cluster16.genes",
                               "cluster17.genes", "cluster18.genes", "cluster19.genes", "cluster20.genes", "cluster21.genes", "cluster22.genes", "cluster23.genes", "cluster24.genes")

cluster.sprucespecificgenes_list <- lapply(cluster.genes_list, function(x) {
  x[x %in% spruce.specific]
})
  
#blabla <- "MA_3445307g0010" %in% unlist(cluster.sprucespecificgenes_list)
#lapply(cluster.sprucespecificgenes_list, function(x) {
#  "MA_344980" %in% x
#})
  
 # df$spruce.specific <- cluster1.genes[cluster1.genes %in% spruce.specific]

cluster1_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster1.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster2_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster2.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster3_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster3.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster4_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster4.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster5_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster5.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster6_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster6.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster7_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster7.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster8_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster8.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))
 
cluster9_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster9.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster10_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster10.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster11_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster11.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster12_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster12.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster13_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster13.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster14_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster14.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster15_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster15.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster16_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster16.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster17_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster17.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster18_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster18.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster19_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster19.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster20_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster20.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster21_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster21.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster22_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster22.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster23_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster23.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))

cluster24_spruce.specific <- plot.new()
grid.draw(venn.diagram(list(
  area1 =cluster24.genes,
  area2 =spruce.specific), 
  filename = NULL, fill=pal[1:2]))


#' #fisher test

fisher.test(matrix(c(90,586,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(97,801,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(101,554,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value
fisher.test(matrix(c(126,665,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value
fisher.test(matrix(c(175,1086,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(136,679,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value
fisher.test(matrix(c(179,878,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(164,681,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value
fisher.test(matrix(c(118,684,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(110,651,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(82,719,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(90,793,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value
fisher.test(matrix(c(116,743,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value
fisher.test(matrix(c(112,569,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(107,665,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(94,777,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value


fisher.test(matrix(c(148,790,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(105,855,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(134,730,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(86,749,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(101,645,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(203,743,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(147,661,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value

fisher.test(matrix(c(163,1062,
                     3476,24528),
                   ncol=2,byrow=TRUE))$p.value
clusters <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24)

percentages <- c(13,11,15,16,14,17,17,19,15,14,10,10,14,16,14,11,16,
                 11,16,10,14,21,18,13)

perc_of_clusters <- data.frame(clusters,percentages)

barplot(perc_of_clusters[order(perc_of_clusters$percentages), "percentages"], 
        names.arg = perc_of_clusters[order(perc_of_clusters$percentages), "clusters"], 
        xlab = "clusters", 
        ylab = "% of spruce-specific genes", 
        ylim = c(0,25), 
        col = "skyblue3")
abline(h = c(15,20), col = "grey")
barplot(perc_of_clusters[order(perc_of_clusters$percentages), "percentages"], 
        names.arg = perc_of_clusters[order(perc_of_clusters$percentages), "clusters"], 
        xlab = "clusters", 
        ylab = "% of spruce-specific genes", 
        ylim = c(0,25), 
        col = "skyblue3",
        ad = TRUE)
