#' Read in file with information about predicted miRNA targets
miRNA_details <- read.csv2("/mnt/picea/projects/functional-genomics-course/spring2019/group08/data/Pabies_SE_miRNA_results_all.csv",
                           row.names = 1)


#' Import candidate genes
TFs <- read.table("~/AFG2019/DE_TF_PABIES_Hoja1.tsv", header = TRUE, sep = "\t")

#' Check if candidate genes are among predicted targets of miRNAs
miRNA_details_TFs <- lapply(TFs$GENE, function(x){
    miRNA_details[grep(pattern = x, miRNA_details$targets), ]
})

names(miRNA_details_TFs) <- TFs$GENE

#' How many miRNA are potential regulators of candidate genes?
sapply(miRNA_details_TFs, nrow)


#' Plot expression profile of miRNA
library(ggplot2)
data <- read.csv("/mnt/picea/projects/functional-genomics-course/spring2019/group08/data/Pabies_SE_miRNA_filtered_exp1rep3tp1_noOutlier_zinbNormalised_counts.csv", row.names = 1, header = TRUE)
meta <- read.csv("/mnt/picea/projects/functional-genomics-course/spring2019/group08/data/Pabies_SE_miRNA_sampleInfo_noOutlier.csv", row.names = 1, header = TRUE)
data <- t(data) 
rownames(data) <- meta$stage
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
        ggtitle(paste("Eigengene plot.",  "miRNA:", genes)) +
        theme(legend.position = "none")
}


plotEigen(data, "miRNA_20838-5p")
plotEigen()

