#' ---
#' title: "Somatic Embryogenesis Network data preparation"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Set the working dir
setwd("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")
#' ```

suppressPackageStartupMessages(library(stringr))

#' Helper
source("~/Git/UPSCb/src/R/featureSelection.R")

#' # Data  
#' ## Genes
dat <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/library-size-corrected_blindFALSE-rlog_gene-expression_data-3_genes_and_TEs.csv",
                row.names=1)

colnames(dat) <- sub(".*_P","P",colnames(dat))

samples <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/samples-3.csv",as.is=TRUE)

samples$ScilifeID <- sub(".*_P","P",sub("_L00[0-9]$","",samples$ScilifeID))

samples <- samples[match(colnames(dat),samples$ScilifeID),]

#' ### Filter
sels <- lapply(1:18,function(i){
  featureSelect(counts=as.matrix(dat),conditions=as.factor(samples$Stages),exp=i,nrep=3)
})

plot(sapply(sels,sum),main="Number of selected genes by cutoff values",type="b",ylab="Number of genes",xlab="cutoff")
abline(v=8,lty=2)

dat_sel <- dat[sels[[8]],]

#' ## miRNAs  
miRNA <- read.csv("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_miRNA_filtered_exp1rep3tp1_noOutlier_zinbNormalised_counts.csv",
         row.names = 1)

sample_info_miRNA <- read.csv("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_miRNA_sampleInfo_noOutlier.csv",
                              row.names = 1)

#' ### Filter  
#' Raw reads of miRNAs were filtered before the normalisation; expression value had to be 
#' at least 1 in at least 3 replicates of at least one time point of the experiment.  
#' 

#' ## Join gene and miRNA expression tables  
#' 
#' Check how many replicates there are in each dataset (miRNAs, genes). 
#' Discard possible replicates in miRNA dataset, which are not present in the gene dataset.
#' Add possible replicates that are not present in miRNA data, but they are in gene dataset, and fill them with zeros?  
#'  
# Add missing information about submitted ID to miRNA data
submittedID <- samples$SubmittedID
# submitted ID for sRNA data has S instead of T in the name, but for now keep the names the same
#submittedID <- sub("T", "S", submittedID)

# Remove name of the replicate that is not present in miRNA dataset
submittedID <- submittedID[submittedID != "#PabK14-03RO_TR2"]

# Order in the same way and concatenate
sample_info_miRNA <- sample_info_miRNA[order(rownames(sample_info_miRNA), decreasing = FALSE), ]
submittedID[1:7] <- submittedID[1:7][order(submittedID[1:7], decreasing = FALSE)]
sample_info_miRNA$SubmittedID <- submittedID

# There is one more sample in the gene dataset
length(sample_info_miRNA$SubmittedID) == length(samples$SubmittedID)

table(samples$Stages)
table(sample_info_miRNA$stage)

#' Add sample to the miRNA dataset and fill it with zeros. We do this because we are adding 
#' miRNA information to already existing network, otherwise it would be better to remove the sample from
#' the gene dataset.
miRNA$P7854_120 <- 0

sample_info_miRNA <- rbind(sample_info_miRNA, c("S7", "B2", "R7-4", "roots", "#PabK14-03RO_TR2"))
levels(sample_info_miRNA$replicate) <- c(levels(sample_info_miRNA$replicate), "R7-4")
rownames(sample_info_miRNA)[31] <- "P7854_120"
sample_info_miRNA["P7854_120", "replicate"] <- "R7-4"

# Order sample info again by stages
sample_info_miRNA <- sample_info_miRNA[order(sample_info_miRNA$stage), ]

#' Join datasets  
# Order miRNA dataset in the same way as gene dataset and bind them
colnames(dat_sel) <- samples[match(colnames(dat_sel), samples$ScilifeID), "SubmittedID"]
colnames(miRNA) <- sample_info_miRNA[colnames(miRNA), "SubmittedID"]

length(colnames(dat_sel)) == length(colnames(miRNA))
colnames(dat_sel) == colnames(miRNA)
colnames(dat_sel) %in% colnames(miRNA)

miRNA <- miRNA[, colnames(dat_sel)]
colnames(dat_sel) == colnames(miRNA)

# Bind dataframes
df_merged <- rbind(dat_sel, miRNA)

#' # Export
dir.create("seidr_geneTEmiRNA",showWarnings=FALSE)

#' * gene by column, without names matrix
write.table(t(df_merged),
            file="seidr_geneTEmiRNA/headless.tsv",
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' * gene names, one row
write.table(t(sub("\\.1$","",rownames(df_merged))),
            file="seidr_geneTEmiRNA/genes.tsv",
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' * only miRNA names, without genes and TEs
write.table(rownames(df_merged)[grepl("miRNA", rownames(df_merged))],
            file="seidr_geneTEmiRNA/miRNAs.tsv",
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)


#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
