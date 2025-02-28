#' ---
#' title: "Search of DE genes in InfoMap clusters (spruce SE)"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
#' # Setup  
#' Load libraries  
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(plyr))

# edit !!! suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/gopher.R"))

#' Import data
clusters <- read_tsv("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/backbone/infomapCluster.tsv")

#' # Analysis  
#' ## Gene co-expression clusters
#' 
#' How many clusters are there altogether and how many genes they include?
table(clusters$data.level1)

#' How many clusters and which ones contain more than 1 or at least 100/200 genes? ######cumulative number of clusters?????????
#' How many genes are included in these clusters alltogether?
sum(table(clusters$data.level1)>1)
clusters_1more <- table(clusters$data.level1)[table(clusters$data.level1)>1]
sum(clusters_1more)

sum(table(clusters$data.level1)>9)
clusters_10 <- table(clusters$data.level1)[table(clusters$data.level1)>9]
sum(clusters_10)

sum(table(clusters$data.level1)>99)
clusters_100 <- table(clusters$data.level1)[table(clusters$data.level1)>99]
sum(clusters_100)

sum(table(clusters$data.level1)>199)
clusters_200 <- table(clusters$data.level1)[table(clusters$data.level1)>199]
sum(clusters_200)

#' ## DEGs in the network  
#' 
#' In which clusters can we find differentially expressed genes?  
#' 
#' Import DEGs
load("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/DEGs_DESeqResults_median_p001_lfc05.Rda")

#' Prepare DE data for the network  
#' 
#' extract criteria needed to see in the network (lfc & padj) and return data frame
list_df <- lapply(res_sig_list, function(x){
    data.frame(x[ , c("log2FoldChange", "padj")])
})

df_DEGs <- list_df %>% 
    map(rownames_to_column, "gene") %>%
    purrr::reduce(full_join, by = "gene")

colnames(df_DEGs) <- c("gene",
                       "lfc_2v1", "padj_2v1", 
                       "lfc_3v2", "padj_3v2", 
                       "lfc_4v3", "padj_4v3",
                       "lfc_5v4", "padj_5v4", 
                       "lfc_6v5", "padj_6v5", 
                       "lfc_7v6", "padj_7v6", 
                       "lfc_8v7", "padj_8v7")

df_DEGs$gene <- sub(".1$", "", df_DEGs$gene)

#write_tsv(df_DEGs,"/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/DEGs4network.tsv", na = "")

#' How many genes do two datasets (InfoMap clusters, DEGs) have in common?
length(intersect(df_DEGs$gene, clusters$data.gene))

#' Let's focus in the further analyses all the clusters that have at least 100 genes
clusters_X <- clusters_100

#' In which InfoMap cluster are DEGs found?  
DEGs_clusters <- clusters[clusters$data.gene %in% df_DEGs$gene, c("data.gene", "data.level1")]
colnames(DEGs_clusters) <- sub("data.", "", colnames(DEGs_clusters))

table(DEGs_clusters$level1)

# in which clusters with at least X genes
table(DEGs_clusters$level1)[names(clusters_X)]

length(table(DEGs_clusters$level1)[names(clusters_X)]) == length(clusters_X)
# DEGs are found in all the clusters, that contain at least X genes  

#' How many DEGs would we expect in an InfoMap cluster by chance?   
# Number of all DEGs in the clusters with at least X genes
DEGs_in_clustersX <- intersect(df_DEGs$gene, clusters$data.gene[clusters$data.level1 %in% names(clusters_X)])
length(DEGs_in_clustersX)
# There is altogether this many genes in clusters with at least X genes
length(clusters$data.gene[clusters$data.level1 %in% names(clusters_X)])
# On average DEGs would represent this many % of all the genes in the cluster
average_DEGs_in_clusterX <- length(DEGs_in_clustersX)/length(clusters$data.gene[clusters$data.level1 %in% names(clusters_X)])*100
# Should the expectation be calculated on all the clusters, not only clusters with at least X genes?  

#' What percentage of clusters do DEGs represent?
percent_DEGs_in_clusters <- table(DEGs_clusters$level1)[1:70]/table(clusters$data.level1)[1:70]*100
boxplot(as.matrix(table(DEGs_clusters$level1)[1:70]/table(clusters$data.level1)[1:70]*100))  

# in clusters with at least X genes
sort(percent_DEGs_in_clusters[names(clusters_X)], decreasing = TRUE)
boxplot(as.matrix(percent_DEGs_in_clusters[names(clusters_X)]), 
        main = "Abundance of DEGs in the clusters",
        ylab = "% of DEGs in the cluster")

barplot(sort(percent_DEGs_in_clusters[names(clusters_X)], decreasing = TRUE),
        main = "Abundance of DEGs in the clusters",
        ylab = "% of DEGs in the cluster",
        xlab = "clusters",
        ylim = c(0, 100))
# add a line representing the average % of DEGs in clusters with at least X genes
abline(h = average_DEGs_in_clusterX,
       col = "red")

# How many clusters exceed the average content of DEGs?
table(percent_DEGs_in_clusters[names(clusters_X)]>average_DEGs_in_clusterX)

#' Are DEGs (DE between two stages) contained in one or more InfoMap clusters?
# merge info
df_DEGs_1level <- left_join(df_DEGs, DEGs_clusters, "gene")

# number of clusters in which DEGs from one comparison appear
lapply(colnames(df_DEGs_1level)[grep("lfc_", colnames(df_DEGs_1level))], function(x){
    length(table(df_DEGs_1level$level1[!is.na(df_DEGs_1level[ , x])]))
})

# for clusters with at least X genes
df_DEGs_1level_X <- df_DEGs_1level[df_DEGs_1level$level1 %in% names(clusters_X), ]

lapply(colnames(df_DEGs_1level_X)[grep("lfc_", colnames(df_DEGs_1level_X))], function(x){
    length(table(df_DEGs_1level_X$level1[!is.na(df_DEGs_1level_X[ , x])]))
})

# What percent of each cluster, where they are found, do DEGs from certain stage comparion represent?
# Are there any clusters with higher abundance of DEGs from certain stage comparison?
percent_DEGs_perStage_in_clusters <- lapply(colnames(df_DEGs_1level_X)[grep("lfc_", colnames(df_DEGs_1level_X))], function(x){
    nr <- table(df_DEGs_1level_X$level1[!is.na(df_DEGs_1level_X[ , x])])
    nr/clusters_X[names(nr)]*100
})
names(percent_DEGs_perStage_in_clusters) <- sub("lfc_", "", colnames(df_DEGs_1level_X)[grep("lfc_", colnames(df_DEGs_1level_X))])

# ??????how to merge it into dataframe, object names representing e.g. rows, names of clusters representing columns
prepare_matrix <- lapply(percent_DEGs_perStage_in_clusters, function(x){
    t(as.matrix(x))
})

percent_DEGs_perStage_in_clusters_df <- do.call(rbind.fill.matrix, prepare_matrix)
rownames(percent_DEGs_perStage_in_clusters_df) <- names(prepare_matrix)
percent_DEGs_perStage_in_clusters_df <- percent_DEGs_perStage_in_clusters_df[ , order(as.numeric(colnames(percent_DEGs_perStage_in_clusters_df)))]

barplot2(percent_DEGs_perStage_in_clusters_df)

# ! Not possible: It does not add up to 100 %, because DEGs between different stages can be the same  
#change NAs for 0 first, so all the clusters get plotted
percent_DEGs_perStage_in_clusters_df[is.na(percent_DEGs_perStage_in_clusters_df)] <- 0

apply(percent_DEGs_perStage_in_clusters_df, 1, function(x){
    barplot2(x, 
             col = "cyan3", 
             xlab = "clusters", 
             ylab = "% of DEGs", 
             plot.grid = TRUE) 
             # how to add colnames for the title?
    })

# it would be also interesting to see plots by cluster and decide if we can find cluters representing DEGs in certain stages  
lapply(seq_along(colnames(percent_DEGs_perStage_in_clusters_df)), function(x){
  barplot2(percent_DEGs_perStage_in_clusters_df[ , x], 
           col = "cyan3", 
           xlab = "stages of DE", 
           ylab = "% of DEGs", 
           ylim = c(0, 100),
           plot.grid = TRUE,
           main = paste("cluster", colnames(percent_DEGs_perStage_in_clusters_df)[x], sep = " "))
})

# write in pdf
pdf(file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/graphs_DEGsByStageInClusters.pdf", width=12, height=18)
par(mfrow=c(9,6), mar=c(2,2,2,2))
lapply(seq_along(colnames(percent_DEGs_perStage_in_clusters_df)), function(x){
  barplot2(percent_DEGs_perStage_in_clusters_df[ , x], 
           col = "cyan3", 
           xlab = "stages of DE", 
           ylab = "% of DEGs", 
           ylim = c(0, 100),
           plot.grid = TRUE,
           main = paste("cluster", colnames(percent_DEGs_perStage_in_clusters_df)[x], sep = " "))
})
dev.off()

# maybe a heatplot would be better  

#' ## miRNAs and their targets  
#' ### target prediction with default expectation value (>= 5)  
#' 
#' In which clusters are predicted targets of DE miRNAs?  
#' 
# import data
miRNA_details <- read.csv2("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_miRNA_results_all.csv", 
                          row.names = 1)
load("~/Git/UPSCb/projects/spruce-srna/doc/DEmiRNAs_timeline.rda")
DE_miRNA_dds <- res_sig_timeline
rm(res_sig_timeline)

# extract names of DE miRNAs
DE_miRNAs_names <- unique(unlist(lapply(DE_miRNA_dds, rownames)))

# extract targets of DE miRNAs
DE_miRNA_targets <- miRNA_details[DE_miRNAs_names, "targets"]
DE_miRNA_targets <- lapply(DE_miRNA_targets, function(x){
    y <- unlist(strsplit(as.character(x), split = ", "))
    sub(".1$", "", y)
})
names(DE_miRNA_targets) <- DE_miRNAs_names

# are predicted targets of each DE miRNA included in InfoMap clusters?
targets_clusters <- lapply(DE_miRNA_targets, function(x){
    x %in% clusters$data.gene
})

nr_targets_in_clusters <- sapply(targets_clusters, sum)
percent_targets_in_clusters <- nr_targets_in_clusters/lengths(targets_clusters)*100

barplot(percent_targets_in_clusters, 
        col = c("aquamarine2"))

# which predicted targets can be found in the clusters?
targets_clusters_genes <- lapply(seq_along(DE_miRNA_targets), function(x){
    DE_miRNA_targets[[x]][targets_clusters[[x]]]
})
names(targets_clusters_genes) <- names(DE_miRNA_targets)

# miRNA_20838-5p has more than one LTR-TE among predicted targets, present in the clusters. Quick Google search 
# shows that it was already shown to regulate LTR-TE. Check more on that later.  

#' How many miRNA targets are found in the clusters and are DE?  
targets_clusters_DE <- lapply(DE_miRNA_targets, function(x){
    x %in% DEGs_clusters$gene
})

nr_targets_in_clusters_DE <- sapply(targets_clusters_DE, sum)
nr_targets_in_clusters_notDE <- nr_targets_in_clusters-nr_targets_in_clusters_DE

nr_targets_final <- data.frame("DE_targets_in_clusters" = nr_targets_in_clusters_DE,
                               "notDE_targets_in_clusters" = as.integer(nr_targets_in_clusters_notDE),
                               "targets_not_in_clusters" = as.integer(lengths(targets_clusters)-nr_targets_in_clusters))

percent_targets_final <- nr_targets_final/rowSums(nr_targets_final)*100

barplot(t(percent_targets_final[,1:2 ]),
        col = c("cadetblue4", "aquamarine2"), 
        main = "miRNA targets in InfoMap clusters",
        ylim = c(0, 100),
        xlab = "miRNAs", 
        ylab = "% of predicted miRNA targets")
legend("top", legend = rownames(t(percent_targets_final[1:2])), fill = c("cadetblue4", "aquamarine2"), cex = 0.9)
boxplot(percent_targets_final)

#' Which targets are found in the clusters and are DE?
targets_clusters_DE_genes <- lapply(seq_along(DE_miRNA_targets), function(x){
    DE_miRNA_targets[[x]][targets_clusters_DE[[x]]]
})
names(targets_clusters_DE_genes) <- names(DE_miRNA_targets)

# number of unique DE targets in the InfoMap clusters
length(unique(unlist(targets_clusters_DE_genes)))

#' Check later: Do targets that are DE and included in the clusters have higher score from psRNATarget?  
#' 
#' In which clusters are DE targets of DE miRNA found?
clusters_targeted_by_DE_miRNAs <- lapply(targets_clusters_DE_genes, function(x){
    clusters[clusters$data.gene %in% x, c("data.gene", "data.level1")]
})

#' How many clusters do DE targets of DE miRNAs belong to?
nr_clusters_targeted_by_each_DE_miRNA <- sapply(clusters_targeted_by_DE_miRNAs, function(x){
    length(table(x[ , "data.level1"]))
})

barplot(table(nr_clusters_targeted_by_each_DE_miRNA),
     main = "DE miRNAs targeting DE genes in different number of gene clusters", 
     ylab = "number of miRNAs", 
     xlab = "number of clusters targeted by miRNA",
     col = "cadetblue3")

#' Group miRNAs by gene clusters they are targeting  
# get all the clusters containing DE genes, targeted by DE miRNAs
targeted_clusters_unique <- as.character(unique(unlist(lapply(clusters_targeted_by_DE_miRNAs, function(x){
    (x[, "data.level1"])
}))))

#' DE miRNAs target DE genes found in `r length(targeted_clusters_unique)` clusters.  
#'   
# prepare data frame
DEmiRNA_DEtarget_cluster_df <- ldply(clusters_targeted_by_DE_miRNAs, data.frame)
colnames(DEmiRNA_DEtarget_cluster_df) <- c("miRNA", "target", "cluster")

DEmiRNA_byCluster <- lapply(unique(DEmiRNA_DEtarget_cluster_df$cluster), function(x){
  unique(DEmiRNA_DEtarget_cluster_df$miRNA[DEmiRNA_DEtarget_cluster_df$cluster == x])
})
names(DEmiRNA_byCluster) <- unique(DEmiRNA_DEtarget_cluster_df$cluster)

lengths(DEmiRNA_byCluster)

barplot(sort(lengths(DEmiRNA_byCluster), decreasing = TRUE),
        ylim = c(0, 100),
        main = "Number of DE miRNAs targeting genes from each cluster",
        xlab = "Cluster",
        ylab = "Number of DE miRNAs")

#' Use only clusters with min 100 genes
DEmiRNA_byClusterX <- DEmiRNA_byCluster[names(clusters_X)]

barplot(sort(lengths(DEmiRNA_byClusterX), decreasing = TRUE),
        ylim = c(0, 100),
        main = "DE miRNAs targeting genes from each cluster",
        sub = "only clusters with at least 100 genes",
        xlab = "Cluster",
        ylab = "Number of DE miRNAs",
        col = "gold")

#' How many miRNA families target genes in each cluster?
DEmiRNA_family_byClusterX <- lapply(DEmiRNA_byClusterX, function(x){
  best_byClusterX <- miRNA_details$best.miRBase[rownames(miRNA_details) %in% x]
  unique(unlist(str_extract_all(best_byClusterX, "MIR[0-9]+")))
})

# NA means there was a novel miRNA, which does not have best hit in miRBase

lengths(DEmiRNA_family_byClusterX)

barplot(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE),
        ylim = c(0, 50),
        main = "Families of DE miRNAs targeting genes from each cluster",
        sub = "only clusters with at least 100 genes",
        xlab = "Cluster",
        ylab = "Number of miRNA families",
        col = "orange")

#' Plot size of these clusters and amount of DE targets in them (in the same order)
barplot(clusters_X[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        xlab = "Cluster",
        ylab = "Number of genes in the cluster",
        main = "Size of the clusters",
        col = "seagreen")

#' How many DE targets is in these clusters?
DEtargets_in_clustersX <- table(DEmiRNA_DEtarget_cluster_df$cluster[DEmiRNA_DEtarget_cluster_df$cluster %in% names(clusters_X)])

barplot(DEtargets_in_clustersX[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        ylim= c(0, 200),
        xlab = "Cluster",
        ylab = "DE targets",
        main = "Number of DE targets in each cluster",
        col = "lightblue")

# Put on the same plot - stacked bars
all(names(clusters_X) == names(DEtargets_in_clustersX))
DEtargets_and_clustersX_df <- data.frame(DEtargets_in_clustersX, clusters_X)
names(DEtargets_and_clustersX_df) <- c("cluster", "DEtargets", "clusterX", "allGenes")
all(DEtargets_and_clustersX_df$cluster == DEtargets_and_clustersX_df$clusterX)
rownames(DEtargets_and_clustersX_df) <- DEtargets_and_clustersX_df$cluster
DEtargets_and_clustersX_df <- DEtargets_and_clustersX_df[ , c(2,4)]
DEtargets_and_clustersX_df$otherGenes <- DEtargets_and_clustersX_df$allGenes - DEtargets_and_clustersX_df$DEtargets

barplot(t(DEtargets_and_clustersX_df[ , c("DEtargets", "otherGenes")]),
        col = c("lightblue", "seagreen"), 
        main = "Number of DE miRNA targets in InfoMap clusters",
        ylim = c(0, 1500),
        xlab = "Clusters", 
        ylab = "Number of genes")
legend("top", legend = rownames(t(DEtargets_and_clustersX_df[ , c("DEtargets", "otherGenes")])), fill = c("lightblue", "seagreen"), cex = 0.9)

#' Quick look at "number of DE miRNAs/miRNA families per gene/DEG"
barplot(sort(lengths(DEmiRNA_byClusterX), decreasing = TRUE)/clusters_X[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        main = "DE miRNA/gene",
        xlab = "Cluster")
barplot(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE)/clusters_X[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        main = "DE miRNA family/gene",
        xlab = "Cluster")

barplot(sort(lengths(DEmiRNA_byClusterX), decreasing = TRUE)/DEtargets_in_clustersX[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        main = "DE miRNA/DE target",
        xlab = "Cluster")
barplot(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE)/DEtargets_in_clustersX[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        main = "DE miRNA faily/DE target",
        xlab = "Cluster")

#' ### Target prediction with stricter settings (max E <= 4)  
#' As there are a lot of predicted miRNA targets, repeat the analysis with reduced set of targets  
#' 
# copy the code above and overwrite object in the beginning  

#' In which clusters are predicted targets of DE miRNAs?  
#' 
# import data
miRNA_details <- read.csv2("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_miRNA_results_all_maxE_4.csv", 
                           row.names = 1)
load("~/Git/UPSCb/projects/spruce-srna/doc/DEmiRNAs_timeline.rda")
DE_miRNA_dds <- res_sig_timeline
rm(res_sig_timeline)

# extract names of DE miRNAs
DE_miRNAs_names <- unique(unlist(lapply(DE_miRNA_dds, rownames)))

# extract targets of DE miRNAs
DE_miRNA_targets <- miRNA_details[DE_miRNAs_names, "targets"]

DE_miRNA_targets <- lapply(DE_miRNA_targets, function(x){
  y <- unlist(strsplit(as.character(x), split = ", "))
  sub(".1$", "", y)
})
names(DE_miRNA_targets) <- DE_miRNAs_names

# are there any DE miRNAs without predicted targets?
any(isEmpty(DE_miRNA_targets))
DE_miRNA_targets[isEmpty(DE_miRNA_targets)]
# if there are any, exclude them from further analysis about targets
DE_miRNA_targets <- DE_miRNA_targets[!isEmpty(DE_miRNA_targets)]

# are predicted targets of each DE miRNA included in InfoMap clusters?
targets_clusters <- lapply(DE_miRNA_targets, function(x){
  x %in% clusters$data.gene
})

nr_targets_in_clusters <- sapply(targets_clusters, sum)
percent_targets_in_clusters <- nr_targets_in_clusters/lengths(targets_clusters)*100

barplot(percent_targets_in_clusters, col = "aquamarine2")

# which predicted targets can be found in the clusters?
targets_clusters_genes <- lapply(seq_along(DE_miRNA_targets), function(x){
  DE_miRNA_targets[[x]][targets_clusters[[x]]]
})
names(targets_clusters_genes) <- names(DE_miRNA_targets)

# miRNA_20838-5p has more than one LTR-TE among predicted targets, present in the clusters. Quick Google search 
# shows that it was already shown to regulate LTR-TE. Check more on that later.

#' How many miRNA targets are found in the clusters and are DE?  
targets_clusters_DE <- lapply(DE_miRNA_targets, function(x){
  x %in% DEGs_clusters$gene
})

# !!! e.g. miRNA_41778-30 has 2 targets that are not DE; as there is only one value after 
# using table function (FALSE 2 and no TRUE value), cbind reuses value of 2?? 
# in nr_targets_in_clusters it says FALSE 2, TRUE 2. Also other values/miRNAs? are affected.

nr_targets_in_clusters_DE <- sapply(targets_clusters_DE, sum)
nr_targets_in_clusters_notDE <- nr_targets_in_clusters-nr_targets_in_clusters_DE

nr_targets_final <- data.frame("DE_targets_in_clusters" = nr_targets_in_clusters_DE,
                               "notDE_targets_in_clusters" = as.integer(nr_targets_in_clusters_notDE),
                               "targets_not_in_clusters" = as.integer(lengths(targets_clusters)-nr_targets_in_clusters))

percent_targets_final <- nr_targets_final/rowSums(nr_targets_final)*100

barplot(t(percent_targets_final[,1:2 ]),
        col = c("cadetblue4", "aquamarine2"), 
        main = "miRNA targets in InfoMap clusters",
        ylim = c(0, 100),
        xlab = "miRNAs", 
        ylab = "% of predicted miRNA targets")
legend("top", legend = rownames(t(percent_targets_final[1:2])), fill = c("cadetblue4", "aquamarine2"), cex = 0.9)
boxplot(percent_targets_final)

#' Which targets are found in the clusters and are DE?
targets_clusters_DE_genes <- lapply(seq_along(DE_miRNA_targets), function(x){
  DE_miRNA_targets[[x]][targets_clusters_DE[[x]]]
})
names(targets_clusters_DE_genes) <- names(DE_miRNA_targets)

# number of unique DE targets in the InfoMap clustes
length(unique(unlist(targets_clusters_DE_genes)))

#' Check later: Do targets that are DE and included in the clusters have higher score from psRNATarget?  
#' 
#' In which clusters are DE targets of DE miRNA found?
clusters_targeted_by_DE_miRNAs <- lapply(targets_clusters_DE_genes, function(x){
  clusters[clusters$data.gene %in% x, c("data.gene", "data.level1")]
})

#' How many clusters do DE targets of DE miRNAs belong to?
nr_clusters_targeted_by_each_DE_miRNA <- sapply(clusters_targeted_by_DE_miRNAs, function(x){
  length(table(x[ , "data.level1"]))
})

barplot(table(nr_clusters_targeted_by_each_DE_miRNA),
        ylim = c(0,25),
        main = "DE miRNAs targeting DE genes in different number of gene clusters", 
        ylab = "number of miRNAs", 
        xlab = "number of clusters targeted by miRNA",
        col = "cadetblue3")

#' Group miRNAs by gene clusters they are targeting  
# get all the clusters containing DE genes, targeted by DE miRNAs
targeted_clusters_unique <- as.character(unique(unlist(lapply(clusters_targeted_by_DE_miRNAs, function(x){
  (x[, "data.level1"])
}))))

#' DE miRNAs target DE genes found in `r length(targeted_clusters_unique)` clusters.  
#'   
# prepare data frame
DEmiRNA_DEtarget_cluster_df <- ldply(clusters_targeted_by_DE_miRNAs, data.frame)
colnames(DEmiRNA_DEtarget_cluster_df) <- c("miRNA", "target", "cluster")

DEmiRNA_byCluster <- lapply(unique(DEmiRNA_DEtarget_cluster_df$cluster), function(x){
  unique(DEmiRNA_DEtarget_cluster_df$miRNA[DEmiRNA_DEtarget_cluster_df$cluster == x])
})
names(DEmiRNA_byCluster) <- unique(DEmiRNA_DEtarget_cluster_df$cluster)

lengths(DEmiRNA_byCluster)

barplot(sort(lengths(DEmiRNA_byCluster), decreasing = TRUE),
        ylim = c(0, 40),
        main = "Number of DE miRNAs targeting genes from each cluster",
        xlab = "Cluster",
        ylab = "Number of DE miRNAs")

#' Use only clusters with min 100 genes  
#' 
#' Are DE targets found in all of the clusters with at least 100 genes?
table(names(clusters_X) %in% targeted_clusters_unique)
names(clusters_X)[(names(clusters_X) %in% targeted_clusters_unique) == FALSE]
# They can be found in all but one cluster (Cluster 40) with at least 100 genes.

DEmiRNA_byClusterX <- DEmiRNA_byCluster[names(clusters_X)]

barplot(sort(lengths(DEmiRNA_byClusterX), decreasing = TRUE),
        ylim = c(0, 40),
        main = "DE miRNAs targeting genes from each cluster",
        sub = "only clusters with at least 100 genes",
        xlab = "Cluster",
        ylab = "Number of DE miRNAs",
        col = "gold")

#' How many miRNA families target genes in each cluster?
DEmiRNA_family_byClusterX <- lapply(DEmiRNA_byClusterX, function(x){
  best_byClusterX <- miRNA_details$best.miRBase[rownames(miRNA_details) %in% x]
  unique(unlist(str_extract_all(best_byClusterX, "MIR[0-9]+")))
})

# NA means there was a novel miRNA, which does not have best hit in miRBase

lengths(DEmiRNA_family_byClusterX)

barplot(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE),
        ylim = c(0, 25),
        main = "Families of DE miRNAs targeting genes from each cluster",
        sub = "only clusters with at least 100 genes",
        xlab = "Cluster",
        ylab = "Number of miRNA families",
        col = "orange")

#' Plot size of these clusters and amount of DE targets in them (in the same order)
barplot(clusters_X[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        xlab = "Cluster",
        ylab = "Number of genes in the cluster",
        main = "Size of the clusters",
        col = "seagreen")

#' How many DE targets is in these clusters?

DEmiRNA_DEtarget_clusterX_df <- DEmiRNA_DEtarget_cluster_df[DEmiRNA_DEtarget_cluster_df$cluster %in% names(clusters_X), ]
DEtargets_in_clustersX <- table(DEmiRNA_DEtarget_clusterX_df$cluster)

barplot(DEtargets_in_clustersX[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        ylim= c(0, 50),
        xlab = "Cluster",
        ylab = "DE targets",
        main = "Number of DE targets in each cluster",
        col = "lightblue")

# # Put on the same plot - stacked bars
# all(names(clusters_X) == names(DEtargets_in_clustersX)) # cluster 40 is missing in the second object
# DEtargets_and_clustersX_df <- data.frame(DEtargets_in_clustersX, clusters_X)
# names(DEtargets_and_clustersX_df) <- c("cluster", "DEtargets", "clusterX", "allGenes")
# all(DEtargets_and_clustersX_df$cluster == DEtargets_and_clustersX_df$clusterX)
# rownames(DEtargets_and_clustersX_df) <- DEtargets_and_clustersX_df$cluster
# DEtargets_and_clustersX_df <- DEtargets_and_clustersX_df[ , c(2,4)]
# DEtargets_and_clustersX_df$otherGenes <- DEtargets_and_clustersX_df$allGenes - DEtargets_and_clustersX_df$DEtargets

barplot(t(DEtargets_and_clustersX_df[ , c("DEtargets", "otherGenes")]),
        col = c("lightblue", "seagreen"), 
        main = "Number of DE miRNA targets in InfoMap clusters",
        ylim = c(0, 1500),
        xlab = "Clusters", 
        ylab = "Number of genes")
legend("top", legend = rownames(t(DEtargets_and_clustersX_df[ , c("DEtargets", "otherGenes")])), fill = c("lightblue", "seagreen"), cex = 0.9)

#' Quick look at "number of DE miRNAs/miRNA families per gene/DEG"
barplot(sort(lengths(DEmiRNA_byClusterX), decreasing = TRUE)/clusters_X[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        main = "DE miRNA/gene",
        xlab = "Cluster")
barplot(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE)/clusters_X[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        main = "DE miRNA family/gene",
        xlab = "Cluster")

barplot(sort(lengths(DEmiRNA_byClusterX), decreasing = TRUE)/DEtargets_in_clustersX[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        main = "DE miRNA/DE target",
        xlab = "Cluster")
barplot(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE)/DEtargets_in_clustersX[names(sort(lengths(DEmiRNA_family_byClusterX), decreasing = TRUE))],
        main = "DE miRNA family/DE target",
        xlab = "Cluster")

#' ### miRNAs targeting max 5 clusters  
#' Extract those miRNAs
DE_miRNA_max5clusters <- nr_clusters_targeted_by_each_DE_miRNA[nr_clusters_targeted_by_each_DE_miRNA <= 5 & nr_clusters_targeted_by_each_DE_miRNA != 0]

#' Number of DE miRNAs targeting genes in max 5 clusters
length(DE_miRNA_max5clusters)

barplot(table(DE_miRNA_max5clusters),
        ylim = c(0,25),
        xlab = "number of clusters with targets of miRNA",
        ylab = "number of miRNAs")

#' How many known miRNA families do they belong to?
DE_miRNA_max5clusters_fam <- unique(unlist(str_extract_all(miRNA_details[names(DE_miRNA_max5clusters), "best.miRBase"], "MIR[0-9]+")))
sum(!is.na(DE_miRNA_max5clusters_fam))

#' Are there any novel miRNAs, without best hit in miRBase?
any(is.na(DE_miRNA_max5clusters_fam))
DE_miRNA_max5clusters[is.na(DE_miRNA_max5clusters_fam)]

#' How many different targets do they have and how many clusters targets belong to?  
# nr of all the targets
length(DEmiRNA_DEtarget_cluster_df[DEmiRNA_DEtarget_cluster_df$miRNA %in% names(DE_miRNA_max5clusters), "target"])
# nr of unique targets
length(unique(DEmiRNA_DEtarget_cluster_df[DEmiRNA_DEtarget_cluster_df$miRNA %in% names(DE_miRNA_max5clusters), "target"]))
# nr of unique clusters
length(unique(DEmiRNA_DEtarget_cluster_df[DEmiRNA_DEtarget_cluster_df$miRNA %in% names(DE_miRNA_max5clusters), "cluster"]))

# the same in clusters with at least 100 genes  

# nr of all the targets
length(DEmiRNA_DEtarget_clusterX_df[DEmiRNA_DEtarget_clusterX_df$miRNA %in% names(DE_miRNA_max5clusters), "target"])
# nr of unique targets
length(unique(DEmiRNA_DEtarget_clusterX_df[DEmiRNA_DEtarget_clusterX_df$miRNA %in% names(DE_miRNA_max5clusters), "target"]))
# nr of unique clusters
length(unique(DEmiRNA_DEtarget_clusterX_df[DEmiRNA_DEtarget_clusterX_df$miRNA %in% names(DE_miRNA_max5clusters), "cluster"]))

#' Group by miRNA family
DEmiRNA_DEtarget_clusterX_df$miRNA_fam <- str_extract_all(miRNA_details[DEmiRNA_DEtarget_clusterX_df$miRNA, "best.miRBase"], "MIR[0-9]+")
DEmiRNA_DEtarget_clusterX_df$miRNA_fam <- lapply(DEmiRNA_DEtarget_clusterX_df$miRNA_fam, unique)

#!! There is miRNA, which has best hit to two different miRNA families! 
any(lapply(DEmiRNA_DEtarget_clusterX_df$miRNA_fam, function(x){length(x) > 1}))
DEmiRNA_DEtarget_clusterX_df[which(unlist(lapply(DEmiRNA_DEtarget_clusterX_df$miRNA_fam, function(x){length(x) > 1}))), ]
# Search by sequence on miRBase site shows that MIR156 is a better hit, so I will use this one and discard hit to MIR529.
DEmiRNA_DEtarget_clusterX_df$miRNA_fam[which(unlist(lapply(DEmiRNA_DEtarget_clusterX_df$miRNA_fam, function(x){length(x) > 1})))] <- "MIR156"

DEmiRNA_DEtarget_clusterX_df$miRNA_fam <- unlist(DEmiRNA_DEtarget_clusterX_df$miRNA_fam)

DEmiRNA_DEtarget_clusterX_df_max5cl <- DEmiRNA_DEtarget_clusterX_df[DEmiRNA_DEtarget_clusterX_df$miRNA %in% names(DE_miRNA_max5clusters), ]

# nr of unique targets/fam
tapply(DEmiRNA_DEtarget_clusterX_df_max5cl$target, DEmiRNA_DEtarget_clusterX_df_max5cl$miRNA_fam, function(x){length(unique(x))})
# nr of clusters/fam
tapply(DEmiRNA_DEtarget_clusterX_df_max5cl$cluster, DEmiRNA_DEtarget_clusterX_df_max5cl$miRNA_fam, function(x){length(unique(x))})

# check miRNAs without hit to known families  
DEmiRNA_DEtarget_clusterX_df_max5cl_novel <- DEmiRNA_DEtarget_clusterX_df_max5cl[is.na(DEmiRNA_DEtarget_clusterX_df_max5cl$miRNA_fam),]

# how many are they?
length(unique(DEmiRNA_DEtarget_clusterX_df_max5cl_novel$miRNA))
# Their RNA sequences are distinct
miRNA_details[DEmiRNA_DEtarget_clusterX_df_max5cl_novel$miRNA,]

# nr of all the targets
length(DEmiRNA_DEtarget_clusterX_df_max5cl_novel[DEmiRNA_DEtarget_clusterX_df_max5cl_novel$miRNA %in% names(DE_miRNA_max5clusters), "target"])
# nr of unique targets
length(unique(DEmiRNA_DEtarget_clusterX_df_max5cl_novel[DEmiRNA_DEtarget_clusterX_df_max5cl_novel$miRNA %in% names(DE_miRNA_max5clusters), "target"]))
# nr of unique clusters
length(unique(DEmiRNA_DEtarget_clusterX_df_max5cl_novel[DEmiRNA_DEtarget_clusterX_df_max5cl_novel$miRNA %in% names(DE_miRNA_max5clusters), "cluster"]))

#' When are expressed DE miRNA and their DE targets?  
#' Since sampling points are weeks apart, we can only check miRNAs and targets, which are DE at the same stage 
# which miRNAs are expressed at which stage?
DEmiRNAs_stage_max5cl <- sapply(DE_miRNA_dds, function(x){
  names(DE_miRNA_max5clusters) %in% rownames(x)
  })

rownames(DEmiRNAs_stage_max5cl) <- names(DE_miRNA_max5clusters)

DEmiRNAs_stage_max5cl_names <- apply(DEmiRNAs_stage_max5cl, 2, function(x){
  names(which(x))
})

#' Majority of miRNAs are DE in multiple stages! Percents of all miRNAs targeting max 5 clusters DE at different stages:
colSums(DEmiRNAs_stage_max5cl)/nrow(DEmiRNAs_stage_max5cl)*100
#' Suggestion: take stage with the biggest change in expression (lfc) or take higher treshold for filtering?  
#' 
#' additional filtering on dds
# lfc 2? padj value?
#' 
# extract all of their DE unique targets (which are also in the clusters >= 100)  
DEmiRNAs_stage_max5cl_targets <- lapply(DEmiRNAs_stage_max5cl_names, function(x){
  unique(DEmiRNA_DEtarget_clusterX_df_max5cl[DEmiRNA_DEtarget_clusterX_df_max5cl$miRNA %in% x, "target"])
})

# in which stages are these targets expressed?
DEmiRNAs_stage_max5cl_targets_stage <- lapply(DEmiRNAs_stage_max5cl_targets, function(x){
  uniq_targets_stages <- df_DEGs_1level_X[df_DEGs_1level_X$gene %in% x, grepl("lfc", colnames(df_DEGs_1level_X))]
  apply(uniq_targets_stages, 2, function(x){
    nr_targets <- sum(!is.na(x)) # how many targets is DE in each stage comparison
    median_lfc <- median(abs(x), na.rm = TRUE) # what is median of lfc (absolute value) in each stage comparison
    ret.df <- t(data.frame(nr_targets, median_lfc))
    rownames(ret.df) <- c("nr_targets", "median_abs_lfc")
    return(ret.df)
  })
})

# stacked barplot?
# lapply(names(DEmiRNAs_stage_max5cl_targets_stage), function(x){
#   barplot2(DEmiRNAs_stage_max5cl_targets_stage[x][1,])
# })

#' ## miRNA-target pairs with an opposite expression pattern  
#' 
# find pairs of miRNA-target DE at the same stage with the opposite lfc; for all 
# targets (max E <=4), not just those included in the network (filter later)
#
# prepare df for miRNAs (similar as for DEGs for the network) 
# extract lfc & padj and return a data frame
list_df_miRNA <- lapply(DE_miRNA_dds, function(x){
  data.frame(x[ , c("log2FoldChange", "padj")])
})

df_miRNAs <- list_df_miRNA %>% 
  map(rownames_to_column, "miRNA") %>%
  purrr::reduce(full_join, by = "miRNA")

colnames(list_df) <- c("miRNA",
                       "lfc_2v1", "padj_2v1", 
                       "lfc_3v2", "padj_3v2", 
                       "lfc_4v3", "padj_4v3",
                       "lfc_5v4", "padj_5v4", 
                       "lfc_6v5", "padj_6v5", 
                       "lfc_7v6", "padj_7v6", 
                       "lfc_8v7", "padj_8v7")

stage_names <- c("2v1", "3v2", "4v3", "5v4", "6v5", "7v6", "8v7")

DEmiRNA_DEtarget_sameStage <- lapply(stage_names, function(x){
  # prepare
  collfc <- paste0("lfc_", x)
  colpadj <- paste0("padj_", x)
  stage_miRNAs <- df_miRNAs$miRNA[!is.na(df_miRNAs[, collfc])]
  stage_DEGs <- df_DEGs$gene[!is.na(df_DEGs[, collfc])]
  # find DE targets of each DE miRNA in the same stage
  stage_comb <- lapply(DE_miRNA_targets[stage_miRNAs], function(y){
    y[y %in% stage_DEGs]
  })
  # transform list to data frame
  stage_comb_df <- ldply(stage_comb, data.frame)
  colnames(stage_comb_df) <- c("DE_miRNA", "DE_target")
  # add info from DE analysis: lfc, padj
  stage_comb_anot_df <- data.frame(stage_comb_df,
                                   lfc_miRNA = df_miRNAs[match(stage_comb_df$DE_miRNA, df_miRNAs$miRNA), collfc],
                                   padj_miRNA = df_miRNAs[match(stage_comb_df$DE_miRNA, df_miRNAs$miRNA), colpadj],
                                   lfc_target = df_DEGs[match(stage_comb_df$DE_target, df_DEGs$gene), collfc],
                                   padj_target = df_DEGs[match(stage_comb_df$DE_target, df_DEGs$gene), colpadj])
  # add note about DE pattern of miRNA and its target
  stage_comb_anot_df$DE_pattern <- apply(stage_comb_anot_df, 1, function(z){
    if (as.numeric(z["lfc_miRNA"]) > 0 & as.numeric(z["lfc_target"]) > 0){
      paste("up")
    } else if (as.numeric(z["lfc_miRNA"]) < 0 & as.numeric(z["lfc_target"]) < 0){
      paste("down")
    } else {
      paste("opposite")
    }
    })
  return(stage_comb_anot_df)
})
names(DEmiRNA_DEtarget_sameStage) <- stage_names

DEmiRNA_DEtarget_sameStage_opposite <- lapply(DEmiRNA_DEtarget_sameStage, function(x){
  x[x["DE_pattern"] == "opposite", ]
})

#' Some statistics about miRNA-target pairs with the opposite DE pattern  
#' 
#' How many miRNA-target pairs have an opposite pattern of DE in each stage?
stats_opp_pairs <- sapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
  nr_pairs <- nrow(x)
  nr_unique_pairs <- nrow(unique(x[c("DE_miRNA", "DE_target")]))
  nr_miRNAs <- nrow(unique(x["DE_miRNA"]))
  nr_targets <- nrow(unique(x["DE_target"]))
  df <- rbind(nr_pairs, nr_unique_pairs, nr_miRNAs, nr_targets)
  # why the code below does not work?
  # rownames(df) <- c("nr_pairs", "nr_unique_pairs", "nr_miRNAs", "nr_targets")
  return(df)
  })
rownames(stats_opp_pairs) <- c("nr_pairs", "nr_unique_pairs", "nr_miRNAs", "nr_targets")

#' Calculate the same stats for the whole experiment  
# How many miRNAs are represented in these pairs in the whole experiment?
nr_miRNAs_exp <- length(unique(unlist(lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
  x["DE_miRNA"]
}))))

# How many DEGs are represented in these pairs in the whole experiment?
nr_targets_exp <- length(unique(unlist(lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
  x["DE_target"]
}))))

# How many miRNA-target pairs have the opposite expression pattern in the whole experiment?
nr_pairs_exp <- nrow(do.call(rbind,
                                    lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
                                      x[c("DE_miRNA", "DE_target")]
                                      })
                                    ))

nr_pairs_unique_exp <- nrow(unique(do.call(rbind,
                                    lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
                                      x[c("DE_miRNA", "DE_target")]
                                      })
                                    )))


# Add to df
stats_opp_pairs <- cbind(stats_opp_pairs, all = c(nr_pairs_exp, nr_unique_pairs_exp, nr_miRNAs_exp, nr_targets_exp))
stats_opp_pairs

# How to check for all the targets with the opposite DE pattern of one miRNA? E.g.:
lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
  x[x["DE_miRNA"] == "miRNA_10861-5p", ]
})

# How to check for all the miRNAs with the opposite DE pattern to one gene? E.g.:
lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
  x[x["DE_target"] == "MA_181770g0010", ]
})

#' Extract names of the targets in each stage  
DEtargets_opp <- lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
  unique(x["DE_target"])
  })
DEtargets_opp_all <- unique(unlist(DEtargets_opp))

# write them in the file to be able to explore their function
lapply(names(DEtargets_opp), function(x){
  write.table(DEtargets_opp[x],
             file = paste0("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/miRNAtargets_maxE4_oppDE_", x, ".csv"), 
             row.names = FALSE, 
             col.names = FALSE,
             quote = FALSE)
  })

write.table(DEtargets_opp_all,
            file = "/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/miRNAtargets_maxE4_oppDE_all.csv", 
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)
            
#' In how many InfoMap clusters can these targets be found?
DEtargets_opp_all_clusters <- clusters[clusters$data.gene %in% DEtargets_opp_all, c("data.gene", "data.level1")]
colnames(DEtargets_opp_all_clusters) <- sub("data.", "", colnames(DEtargets_opp_all_clusters))
table(DEtargets_opp_all_clusters$level1)

# in which clusters with at least X genes
table(DEtargets_opp_all_clusters$level1)[names(clusters_X)]

#' There are 1-23 targets in almost all of the clusters with >= 100 genes.

#' Functional enrichment ?? here or congenie? example below:
# enrichment_DEtargets_opp <- lapply(DEtargets_opp, function(x){
#   enr <- gopher(x, url="pabies", task = list('go', 'pfam', 'kegg', 'mapman'), alpha = 0.05)
# })
# names(enrichment_DEtargets_opp) <- names(DEtargets_opp)
# 
# temp$go <- temp$go[temp$go$namespace != 'CC',]
# temp$go$newp <- abs(log10(temp$go$padj))
# clusters$`Cluster 34`
# treemap(temp$go,
#         index = c('namespace', 'name'), vSize ="newp",
#         type = "categorical", vColor = 'name', title="",
#         inflate.labels = FALSE, lowerbound.cex.labels = 0,
#         bg.labels = "#CCCCCC00",
#         position.legend = "none"

#' ## LTR-TEs  
#' 
#' In which clusters can be LTR-TEs found?
clusters[grepl("MA_", clusters$data.gene) == FALSE, ]

# There is only 20 of them and they belong to 5 different clusters. 12 LTR-TEs belong to the cluster 1.  
  
#' NOTE: for the SE paper exclude TEs from the pool of targets!?

#' Session info
sessionInfo()


