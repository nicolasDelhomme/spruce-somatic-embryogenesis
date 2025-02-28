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

# edit !!! suppressPackageStartupMessages(source("~/Git/UPSCb/src/R/gopher.R"))

#' Import data
clusters <- read_tsv("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/genes_miRNAs/infomapClusters.tsv")

#' # Analysis  
#' ## Gene co-expression clusters
#' 
#' How many clusters are there altogether and how many genes they include?  
#' Number of clusters with at least 3 genes:
length(table(clusters$cluster))

table(clusters$cluster)

#' How many clusters and which ones contain more than 1 or at least 100/200 genes? ######cumulative number of clusters?????????
#' How many genes are included in these clusters alltogether?
sum(table(clusters$cluster)>9)
clusters_10 <- table(clusters$cluster)[table(clusters$cluster)>9]
sum(clusters_10)

sum(table(clusters$cluster)>99)
clusters_100 <- table(clusters$cluster)[table(clusters$cluster)>99]
sum(clusters_100)

sum(table(clusters$cluster)>199)
clusters_200 <- table(clusters$cluster)[table(clusters$cluster)>199]
sum(clusters_200)

#' ## DEGs in the network  
#' 
#' In which clusters can we find differentially expressed genes?  
#' 
#' Import DEGs
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genes_padj001_lfc05.rda")

#' Prepare DE data for the network  
#' 
#' extract criteria needed to see in the network (lfc & padj) and return data frame
list_df <- lapply(res_sig_genes, function(x){
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

#' How many genes do two datasets (InfoMap clusters, DEGs) have in common?
length(intersect(df_DEGs$gene, clusters$gene))

#' How many % of all DE genes can be found in the clusters?
length(intersect(df_DEGs$gene, clusters$gene))/length(df_DEGs$gene)*100

#' Let's focus in the further analyses all the clusters that have at least 100 genes
clusters_X <- clusters_10

#' In which InfoMap cluster are DEGs found?  
DEGs_clusters <- clusters[clusters$gene %in% df_DEGs$gene, c("gene", "cluster")]

table(DEGs_clusters$cluster)

# in which clusters with at least X genes
table(DEGs_clusters$cluster)[names(clusters_X)]

length(table(DEGs_clusters$cluster)[names(clusters_X)]) == length(clusters_X)
# DEGs are found in all the clusters, that contain at least X genes  

#' How many DEGs would we expect in an InfoMap cluster by chance?   
# Number of all DEGs in the clusters with at least X genes
DEGs_in_clustersX <- intersect(df_DEGs$gene, clusters$gene[clusters$cluster %in% names(clusters_X)])
length(DEGs_in_clustersX)
# There is altogether this many genes in clusters with at least X genes
length(clusters$gene[clusters$cluster %in% names(clusters_X)])
# On average DEGs would represent this many % of all the genes in the cluster
average_DEGs_in_clusterX <- length(DEGs_in_clustersX)/length(clusters$gene[clusters$cluster %in% names(clusters_X)])*100
# Should the expectation be calculated on all the clusters, not only clusters with at least X genes?  

#' What percentage of clusters do DEGs represent?
percent_DEGs_in_clusters <- table(DEGs_clusters$cluster)[1:70]/table(clusters$cluster)[1:70]*100
boxplot(as.matrix(table(DEGs_clusters$cluster)[1:70]/table(clusters$cluster)[1:70]*100))  

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
    length(table(df_DEGs_1level$cluster[!is.na(df_DEGs_1level[ , x])]))
})

# for clusters with at least X genes
df_DEGs_1level_X <- df_DEGs_1level[df_DEGs_1level$cluster %in% names(clusters_X), ]

lapply(colnames(df_DEGs_1level_X)[grep("lfc_", colnames(df_DEGs_1level_X))], function(x){
    length(table(df_DEGs_1level_X$cluster[!is.na(df_DEGs_1level_X[ , x])]))
})

# What percent of each cluster, where they are found, do DEGs from certain stage comparion represent?
# Are there any clusters with higher abundance of DEGs from certain stage comparison?
percent_DEGs_perStage_in_clusters <- lapply(colnames(df_DEGs_1level_X)[grep("lfc_", colnames(df_DEGs_1level_X))], function(x){
    nr <- table(df_DEGs_1level_X$cluster[!is.na(df_DEGs_1level_X[ , x])])
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
           main = colnames(percent_DEGs_perStage_in_clusters_df)[x])
})

# write in pdf
pdf(file="~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/graphs_DEGsByStageInClusters_only2ndBatch.pdf", width=12, height=18)
par(mfrow=c(9,6), mar=c(2,2,2,2))
lapply(seq_along(colnames(percent_DEGs_perStage_in_clusters_df)), function(x){
  barplot2(percent_DEGs_perStage_in_clusters_df[ , x], 
           col = "cyan3", 
           xlab = "stages of DE", 
           ylab = "% of DEGs", 
           ylim = c(0, 100),
           plot.grid = TRUE,
           main = colnames(percent_DEGs_perStage_in_clusters_df)[x])
})
dev.off()

# maybe a heatplot would be better  

#' ## miRNAs and their targets  
#' ### target prediction with expectation value <= 4  
#' 
#' In which clusters are predicted targets of DE miRNAs?  
#' 
# import data
miRNA_details <- read.csv2("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_miRNA_results_all_maxE_4.csv", 
                          row.names = 1)
load("~/Git/UPSCb/projects/spruce-srna/doc/DEmiRNAs_timeline_2ndBatch.rda")
DE_miRNA_dds <- res_sig_timeline_W
rm(res_sig_timeline_W)

# extract names of DE miRNAs
DE_miRNAs_names <- unique(unlist(lapply(DE_miRNA_dds, rownames)))

# extract targets of DE miRNAs
DE_miRNA_targets <- miRNA_details[DE_miRNAs_names, "targets"]
DE_miRNA_targets <- lapply(DE_miRNA_targets, function(x){
    y <- unlist(strsplit(as.character(x), split = ", "))
    sub(".1$", "", y)
})
names(DE_miRNA_targets) <- DE_miRNAs_names

#' How many DE miRNAs were found in the clusters?
sum(DE_miRNAs_names %in% clusters$gene)

#' In which clusters and how many miRNAs?  
#' all predicted miRNAs:
clusters_miRNAs <- clusters[grepl("miRNA", clusters$gene), ]
table(clusters_miRNAs$cluster)

#' DE miRNAs:
clusters_DEmiRNAs <- clusters[clusters$gene %in% DE_miRNAs_names, ]
table(clusters_DEmiRNAs$cluster)

#' Which miRNAs are in the same cluster?  
# get best miRBase hit
clusters_DEmiRNAs$best.miRBase <- miRNA_details[clusters_DEmiRNAs$gene, "best.miRBase"]
# extract miRNA family
clusters_DEmiRNAs$MIRfam <- lapply(clusters_DEmiRNAs$best.miRBase, function(x){
  unique(unlist(str_extract_all(x, "[0-9]+")))
})
# check if there is always only one family per mature miRNA sequence
any(lengths(clusters_DEmiRNAs$MIRfam) != 1)

# How many and which families of DE miRNAs are in the same cluster?
miRfam_in_cluster <- sapply(unique(clusters_DEmiRNAs$cluster), function(cl){
  fam <- unlist(clusters_DEmiRNAs$MIRfam[clusters_DEmiRNAs$cluster == cl])
  fam[is.na(fam)] <- "novel"
  table(fam)
})

#' df with info about miRNA in df
clmiR <- unique(clusters_DEmiRNAs$cluster)[order(unique(clusters_DEmiRNAs$cluster))]
cl_with_miRNA_df <- data.frame(nr_miRNAs = table(clusters_DEmiRNAs$cluster)[clmiR],
                               nr_miRfam = lengths(miRfam_in_cluster)[clmiR])
                               
miRfam_in_cluster[clmiR]


#' Targets of DE miRNAs  
# are predicted targets of each DE miRNA included in InfoMap clusters?
targets_clusters <- lapply(DE_miRNA_targets, function(x){
    x %in% clusters$gene
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

#' How many DE MIRNAs do have DE targets in the clusters?
DEmiRNAs_wDEtargetsInClusters <- unlist(lapply(targets_clusters_DE, function(x){
  any(x)
}))
table(DEmiRNAs_wDEtargetsInClusters)

#' How many of those miRNAs are also contained in the clusters?
length(intersect(DE_miRNAs_names[DE_miRNAs_names %in% clusters$gene], names(which(DEmiRNAs_wDEtargetsInClusters))))

#' Check later: Do targets that are DE and included in the clusters have higher score from psRNATarget?  
#' 
#' In which clusters are DE targets of DE miRNA found?
clusters_targeted_by_DE_miRNAs <- lapply(targets_clusters_DE_genes, function(x){
    clusters[clusters$gene %in% x, c("gene", "cluster")]
})

#' How many clusters do DE targets of DE miRNAs belong to?
nr_clusters_targeted_by_each_DE_miRNA <- sapply(clusters_targeted_by_DE_miRNAs, function(x){
    length(table(x[ , "cluster"]))
})

barplot(table(nr_clusters_targeted_by_each_DE_miRNA),
     main = "DE miRNAs targeting DE genes in different number of gene clusters", 
     ylab = "number of miRNAs", 
     xlab = "number of clusters targeted by miRNA",
     ylim = c(0, 50),
     col = "cadetblue3")

#' Group miRNAs by gene clusters they are targeting  
# get all the clusters containing DE genes, targeted by DE miRNAs
targeted_clusters_unique <- as.character(unique(unlist(lapply(clusters_targeted_by_DE_miRNAs, function(x){
    (x[, "cluster"])
}))))

#' DE miRNAs target DE genes found in `r length(targeted_clusters_unique)` clusters.  
#'   
# prepare data frame
DEmiRNA_DEtarget_cluster_df <- ldply(clusters_targeted_by_DE_miRNAs, data.frame)
colnames(DEmiRNA_DEtarget_cluster_df) <- c("miRNA", "target", "gene_cluster")

DEmiRNA_DEtarget_cluster_df$miRNA_cluster <- lapply(DEmiRNA_DEtarget_cluster_df$miRNA, function(x){
  paste(unlist(clusters_DEmiRNAs[clusters_DEmiRNAs$gene == x, "cluster"]))
})
DEmiRNA_DEtarget_cluster_df$miRNA_cluster[which(unlist(lapply(DEmiRNA_DEtarget_cluster_df$miRNA_cluster, is_empty)))] <- NA

DEmiRNA_DEtarget_cluster_df$MIRfam <- lapply(DEmiRNA_DEtarget_cluster_df$miRNA, function(x){
  unique(unlist(str_extract_all(miRNA_details[x, "best.miRBase"], "MIR[0-9]+")))
})


#!! There is miRNA, which has best hit to two different miRNA families! 
any(lapply(DEmiRNA_DEtarget_cluster_df$MIRfam, function(x){length(x) > 1}))
DEmiRNA_DEtarget_cluster_df[which(unlist(lapply(DEmiRNA_DEtarget_cluster_df$MIRfam, function(x){length(x) > 1}))), ]
# Search by sequence on miRBase site shows that MIR156 is a better hit, so I will use this one and discard hit to MIR529.
DEmiRNA_DEtarget_cluster_df$MIRfam[which(unlist(lapply(DEmiRNA_DEtarget_cluster_df$MIRfam, function(x){length(x) > 1})))] <- "MIR156"



DEmiRNA_byCluster <- lapply(unique(DEmiRNA_DEtarget_cluster_df$gene_cluster), function(x){
  unique(DEmiRNA_DEtarget_cluster_df$miRNA[DEmiRNA_DEtarget_cluster_df$gene_cluster == x])
})
names(DEmiRNA_byCluster) <- unique(DEmiRNA_DEtarget_cluster_df$gene_cluster)

lengths(DEmiRNA_byCluster)

barplot(sort(lengths(DEmiRNA_byCluster), decreasing = TRUE),
        ylim = c(0, 50),
        main = "Number of DE miRNAs targeting genes from each cluster",
        xlab = "Cluster",
        ylab = "Number of DE miRNAs")

barplot(table(lengths(DEmiRNA_byCluster)))

# What is the size of these clusters?
table(clusters$cluster)[names(DEmiRNA_byCluster)]

#' How many of the DE miRNA and their DE targets with the opposite expression pattern can be found in the network clusters?  

#' get all (unique) pairs with opposite DE pattern in any of the stages
pairs_opp_sameStage <- read.table("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/pairs_maxE4_oppDE_all.csv")

#' A) Criteria: target needs to be in any of the clusters (with min 3 nodes)  
# how to compare rows of unequal dfs?
pairs_clusters <- paste0(DEmiRNA_DEtarget_cluster_df$miRNA, "_", DEmiRNA_DEtarget_cluster_df$target)
pairs_opp <- paste0(pairs_opp_sameStage$V1, "_", pairs_opp_sameStage$V2)

sum(pairs_opp %in% pairs_clusters)

#' B) Criteria: target and miRNA need to be in the clusters (with min 3 nodes)  
# how to compare rows of unequal dfs?
both_in_the_clusters <- DEmiRNA_DEtarget_cluster_df[DEmiRNA_DEtarget_cluster_df$miRNA_cluster != "NA" & DEmiRNA_DEtarget_cluster_df$gene_cluster != "NA", ]
pairs_both_clusters <- paste0(both_in_the_clusters$miRNA, "_", both_in_the_clusters$target)

sum(pairs_opp %in% pairs_both_clusters)

#'
#' Write out interesting info:  
# get rid of the lists from the df
DEmiRNA_DEtarget_cluster_df$miRNA_cluster <- vapply(DEmiRNA_DEtarget_cluster_df$miRNA_cluster, paste, collapse = ", ", character(1L))
DEmiRNA_DEtarget_cluster_df$MIRfam <- vapply(DEmiRNA_DEtarget_cluster_df$MIRfam, paste, collapse = ", ", character(1L))

write.csv(DEmiRNA_DEtarget_cluster_df, 
          file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/network/DEmiRNA_DEtarget_cluster_MIRfam.csv", 
          quote = FALSE, row.names = FALSE)

#' 
#' Session info
sessionInfo()


