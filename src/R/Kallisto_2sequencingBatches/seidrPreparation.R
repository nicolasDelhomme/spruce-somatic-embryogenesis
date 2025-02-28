#' ---
#' title: "Somatic Embryogenesis Network data preparation"
#' author: "Iryna Shutava and Nicolas Delhomme"
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

#' Helper
source("~/Git/UPSCb/src/R/featureSelection.R")

#' # Data
dat <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/library-size-corrected_blindFALSE-rlog_gene-expression_data-3_genes_and_TEs.csv",
                row.names=1)

colnames(dat) <- sub(".*_P","P",colnames(dat))

samples <- read.csv("~/Git/UPSCb/projects/spruce-somatic-embryogenesis/doc/samples-3.csv",as.is=TRUE)

samples$ScilifeID <- sub(".*_P","P",sub("_L00[0-9]$","",samples$ScilifeID))

samples <- samples[match(colnames(dat),samples$ScilifeID),]

#' # Filter
sels <- lapply(1:18,function(i){
  featureSelect(counts=as.matrix(dat),conditions=as.factor(samples$Stages),exp=i,nrep=3)
})

plot(sapply(sels,sum),main="Number of selected genes by cutoff values",type="b",ylab="Number of genes",xlab="cutoff")
abline(v=8,lty=2)

#' # Export
dir.create("seidr",showWarnings=FALSE)

#' * gene by column, without names matrix
write.table(t(dat[sels[[8]],]),
            file="seidr/headless.tsv",
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' * gene names, one row
write.table(t(sub("\\.1$","",rownames(dat)[sels[[8]]])),
            file="seidr/genes.tsv",
            col.names=FALSE,
            row.names=FALSE,
            sep="\t",quote=FALSE)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
