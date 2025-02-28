#' Aim: Cluster genes from the spruce SE with InfoMap  
#' 
#' 
#' Load libraries
library(RLinuxModules)
library(tidyverse)
library(treemap)

module("load bioinfo-tools seidr-devel")
module("load bioinfo-tools InfoMap")
source("~/Git/UPSCb/src/R/gopher.R")
<<<<<<< HEAD


#' Set working directory
setwd("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")
=======
setwd("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/backbone/")
>>>>>>> 80b17dcdbe2a3bac56c16d34fdb42ffed7ba4b76


#' Import data
workFolder <- paste0("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/backbone")
aggregatedFile <- paste0(workFolder, "/", "backbone-2-percent.sf")
indexFile <- paste0(workFolder, "/", "outIndex.txt")
treeFile <- paste0(workFolder, "/", "outIndex.tree")

#' Analysis
markovTime <- 0.5
system(paste0("Infomap ", indexFile," -z --markov-time ", markovTime," ", workFolder))

#system(paste0("/mnt/picea/home/bastian/Git/seidr-devel/build/seidr resolve -s ", aggregatedFile, " ", treeFile, " > resolve2.txt") )
#infomapRes <- system(paste0("/mnt/picea/home/bastian/Git/seidr-devel/build/seidr resolve -s ", aggregatedFile, " ", treeFile), intern=TRUE)
#infomapTable <-as_tibble(data.frame(do.call(rbind, strsplit(infomapRes, "\t"))))

data <- read.table("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/backbone/resolve2.txt",
                   sep="\t", header=F)
colnames(data) <- c("path","L1","L2","L3","L4","L5", "number","index","gene")

data$level1 <- data$L1
data$level2 <- paste0(data$L1,":",data$L2)
data$level3 <- paste0(data$L1,":",data$L2,":",data$L3)

# Check if 50-60% genes belong in the 20 biggest clusters

freq1 <- as.data.frame(table(data$level1))
freq1$Var1 <- names(table(data$level1) )
counts <-freq1[order(freq1$Freq, decreasing = T),]
sum(counts[1:20,2]) #how many genes in the first 20 clusters

# Check the size of the first 50 cluster
counts[1:50,]
minClusterSize = 200

# How many clusters do we keep with a minimal size of 200?
numberOfClusters <- length(counts[counts$Freq>minClusterSize, 2])

clusters <- lapply(counts$Var1[1:numberOfClusters], function(x){
    as.character(data$gene[data$level1 == x])
})

# We name the clusters with 'Cluster prefix'
names(clusters) <- paste0(rep("Cluster",length(clusters)), counts$Var1[1:length(clusters)])

# Check the densest cluster in the network: our guys
our_guys.enrichment <- gopher(clusters$Cluster4, alpha = 0.05, task=list("go", "mapman", "kegg", "pfam"), url="pabies", endpoint = "enrichment")

our_guys.enrichment$go <- our_guys.enrichment$go[our_guys.enrichment$go$namespace != 'CC',]
our_guys.enrichment$go$newp <- abs(log10(our_guys.enrichment$go$padj))
treemap(our_guys.enrichment$go,
        index = c('namespace', 'name'), vSize ="newp",
        type = "categorical", vColor = 'name', title="",
        inflate.labels = FALSE, lowerbound.cex.labels = 0,
        bg.labels = "#CCCCCC00",
        position.legend = "none"
)
our_guys.enrichment$mapman$name <- gsub("\\."," ", our_guys.enrichment$mapman$name)
our_guys.enrichment$mapman$newp <- abs(log10(our_guys.enrichment$mapman$padj))
treemap(our_guys.enrichment$mapman,
        index =  'name', vSize ="newp",
        type = "categorical", vColor = 'name', title="",
        inflate.labels = FALSE, lowerbound.cex.labels = 0,
        bg.labels = "#CCCCCC00",
        position.legend = "none"
)
ourguys <- clusters$Cluster4

library(ggplot2)
library(DESeq2)
data <- read.table("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/vala/data/library-size-corrected_blindFALSE-variance-stabilised_gene-expression_data-3_selection.tsv",
                   row.names=1, header = TRUE)
meta <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/technical_samples-3_TEs_R.csv", row.names = 1, header = TRUE)
data <- t(data) 
rownames(data) <- meta$Stages
data <- data[,colMads(data) > 0]
plotEigen <- function(data, genes, inverse = FALSE){
    require(ggplot2)
    
    d <- data[, genes]
    
    expr <- NA
    stage <- NA
    
    if(length(genes) == 1 ){
        expr <- d
        stage <- as.integer(substr(names(d),2,2))
        
    } else if (length(genes > 1)) {
        pca <- prcomp(scale(d))
        pc1 <- pca$x[, 1]
        
        
        if(sum(sign(cor(pca$x[,1,drop = FALSE], d))) < 0) {
            pc1 <- pc1 * -1
        }
        if(inverse)
            pc1 <- pc1 * -1
        
        stage <- as.integer(substr(rownames(d),2,2))
        expr <- pc1
    } else {
        stop("this function needs at least one gene")
    }
    
    
    
    ggplot( data.frame(x = stage, y = expr),
            aes(x = x, y = y)) +
        stat_summary(fun.data = mean_se, geom = "ribbon", fill = "grey", alpha = 0.75) +
        stat_summary(fun.data = mean_se, geom = "line",  aes(col= "red"), lwd = 2) +
        xlab("Stage") +
        ylab("Expression") + 
        ggtitle(paste("Eigengene plot.",  "Genes:", length(colnames(d)))) +
        theme(legend.position = "none")
}

#plotEigen(data, colnames(data)[1:20])
plotEigen(data, paste)



# We enrich each cluster and store its results,
# It will take some time to execute
clusters.enriched <- lapply(clusters, function(x){
    enr <- gopher(x, url="pabies", task = list('go', 'pfam', 'kegg', 'mapman'), alpha = 0.05)
})
names(clusters.enriched) <- names(clusters)

temp$go <- temp$go[temp$go$namespace != 'CC',]
temp$go$newp <- abs(log10(temp$go$padj))
clusters$`Cluster 34`
treemap(temp$go,
        index = c('namespace', 'name'), vSize ="newp",
        type = "categorical", vColor = 'name', title="",
        inflate.labels = FALSE, lowerbound.cex.labels = 0,
        bg.labels = "#CCCCCC00",
        position.legend = "none"
)

#level2 
freq2 <- as.data.frame(table(data$level2))
freq2$Var1 <- names(table(data$level2) )
counts2 <-freq2[order(freq2$Freq, decreasing = T),]
sum(counts2[1:20,2]) #how many genes in the first 20 clusters



#with markov-time = 0.5
#we get 59 % of the genes in the top 20 clusters



# data export for cytoscape
exportData <- data.frame(data$gene, data$level1, data$level2, data$level3)
write.table(exportData, file="infomapCluster.tsv", sep='\t', row.names = F, col.names=TRUE, quote=FALSE)

