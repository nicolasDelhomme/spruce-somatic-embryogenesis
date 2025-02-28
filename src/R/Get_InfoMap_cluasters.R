#' Aim: Get InfoMap clusters  

library(data.table)

infomapTable <- read.delim("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/genes_miRNAs/reducedIndexResolve.txt", header=FALSE, stringsAsFactors=FALSE)

# Table was subsetted from a bigger network, including lncRNAs. It still contains names that are not relevant for this project,
# therefore remove all the rows that contain "TRINITY"
infomapTable <- infomapTable[!grepl("TRINITY", infomapTable$V9), ]
# There are 29048 genes and miRNAs left in the co-expression network

# calculate the number of levels, total columns -3 info columns -1 path column
numberLevels <- length(colnames(infomapTable)) - 4
# create names for the levels
levels <- sapply(1:numberLevels, function(x){
    paste0("Level", x)
})
# add the names to the data frame
colnames(infomapTable) <- c('Path', levels, 'Flow', 'Index', 'Gene')

# create a new results data frame
res <- data.frame(Gene = infomapTable[,length(colnames(infomapTable))],
                  Path= infomapTable[,1],
                  stringsAsFactors = F)

# Data always have level 1 and none of the data is NA
res$Level1 <- infomapTable[,2]

# loop though the rest of the levels and attach the name from the previous ones
for (level in 2:numberLevels) {
    currentLevel <- paste0("Level",level)
    prevLevel <- paste0("Level",(level-1))
    # join names
    res[[currentLevel]] <- paste0(res[[prevLevel]], ":", infomapTable[[currentLevel]])
    # if there is an NA inside the current name, that gene doesn't belong to a cluster in that level, it turns into NA
    res[[currentLevel]] <- ifelse(res[[currentLevel]] %like% "NA", NA, res[[currentLevel]])  
}
df <- res
level='Level1'
freq1 <- as.data.frame(table(df[level]) )
freq1$Var1 <- names(table(df[level]) )
counts <-freq1[order(freq1$Freq, decreasing = T),]

#how many genes in the first 20 clusters
clusters <- length(unique(freq1$Var1))
print(paste(level,"Clusters:", clusters))
genes20 <- 0
if (clusters < 20) {
    genes20 <- sum(counts[1:clusters,2])
} else {
    genes20 <- sum(counts[1:20,2])  
}
genesTotal <- length(rownames(df))

print(paste("Genes in the top 20 clusters", genes20))
print(paste("Genes in the network",genesTotal))
print(round(genes20/genesTotal, digits = 4)*100)

getCountsPerCluster <- function(df, level='Level1'){
    freq1 <- as.data.frame(table(df[level]) )
    freq1$Var1 <- names(table(df[level]) )
    counts <- freq1[order(freq1$Freq, decreasing = T),]
    return(counts)
}


counts <- getCountsPerCluster(df, level=level)
hist(counts$Freq[counts$Freq > 9])

# number of clusters with more than NrGenes:
NrGenes <- c(1, 9, 99 ,199)
sapply(NrGenes, function(x){
    nrClusters <- sum(counts$Freq > x)
    sumOfGenes <- sum(counts$Freq[counts$Freq > x])
    return <- c(nrClusters, sumOfGenes)
})

numberOfGenesInCluster = 3
numberOfClusters=sum(counts$Freq >= numberOfGenesInCluster)
clusters <- lapply(counts$Var1[1:numberOfClusters], function(x){
    as.character(df$Gene[df[level] == x])
})
names(clusters) <- paste0(rep("Cluster",length(clusters)), counts$Var1[1:length(clusters)])
file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr/genes_miRNAs/infomapClusters.tsv"
fileConn<-file(file,"w")
writeLines(paste0("gene","\t","cluster"), fileConn)
clusterNames <- names(clusters)
for (i in 1:length(clusterNames)) {
    thisCluster <- clusterNames[i]
    print (thisCluster)
    for (gene in clusters[i]) {
        writeLines(paste0(gene,"\t",thisCluster), fileConn)
    }
}
close(fileConn)
