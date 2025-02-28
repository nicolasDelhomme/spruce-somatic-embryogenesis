#' ---
#' title: "Genes of interest - DE, miRNA targets"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
#' # Aim  
#' Check list of genes if they are DE in any stage of SE and if they are predicted targets of (DE) miRNAs.  
#' 
#' # Setup  

#' ## Load libraries
packageStartupMessage(library("readxl"))
packageStartupMessage(library("stringr"))

#' ## Import data  
#' Genes of interest  
genesOfInterest <- read_xlsx("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/MA_list_known_SE_genes_transfer.xlsx")
# change names of the columns - substitute space for underscore
colnames(genesOfInterest) <- sub(" ", "_", colnames(genesOfInterest))

#' DE genes  
DEGs_byStage <- read.delim("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DEGs_padj_lfc_byStage.tsv")

#' DE miRNAs  
load("~/Git/UPSCb/projects/spruce-srna/doc/DEmiRNAs_timeline_2ndBatch.rda")
DE_miRNA_dds <- res_sig_timeline_W
rm(res_sig_timeline_W)

miRNA_details <- read.csv2("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_miRNA_results_all_maxE_4.csv", 
                           row.names = 1)
miRNA_details_maxE5 <- read.csv2("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_miRNA_results_all.csv", 
                           row.names = 1)
miRNAs_byStage <- read.delim("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_DEmiRNA_lfc_padj_byStage.tsv")


#' DE miRNAs & DE targets in the same stage
load("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/DEmiRNA_DEtarget_sameStage.rda")  


#' # Analysis  
#' 
#' ## DEGs  
#' 
#' Which genes from the list are DE expressed and in which stages?  
DEGsOfInterest <- DEGs_byStage[DEGs_byStage$gene %in% genesOfInterest$MA_ID, ]
nrow(DEGsOfInterest)
nrow(genesOfInterest)

#' 50 out of 68 genes of interest are differentially expressed in at least one of the stages of SE  

# write them out
write.table(DEGsOfInterest, 
            "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/GenesOfInterest_DEG.tsv", 
            quote = FALSE, 
            sep = "\t")


#' ## miRNA targets  
#' 
#' We are interested in DE genes that are predicted targets of DE miRNAs. 
#' Possibly they have an opposite expression in the same stage of the experiment. 
#' Genes that would satisfy above criteria would make stronger candidates for being true miRNA targets.  
#' 
#' ### max E < 4  
#' 
#' Which DE genes are predicted targets of DE miRNAs?  

# extract names of DE miRNAs
DE_miRNAs_names <- unique(unlist(lapply(DE_miRNA_dds, rownames)))

# extract targets of DE miRNAs
DE_miRNA_targets <- miRNA_details[DE_miRNAs_names, "targets"]
# change the form of the object
DE_miRNA_targets <- lapply(DE_miRNA_targets, function(x){
    y <- unlist(strsplit(as.character(x), split = ", "))
    # modify gene names
    sub(".1$", "", y)
})
names(DE_miRNA_targets) <- DE_miRNAs_names

# how many DE genes of interest is among predicted targets of DE mRNAs?
DEtargets_DEmiRNAs <- sapply(DEGsOfInterest$gene, function(x){
    lapply(DE_miRNA_targets, function(y){
        x %in% y
    })
    })

apply(DEtargets_DEmiRNAs, 2, any)

#' No DE gene of interest is predicted to be a target of DE miRNA.  


# positive control to test if the code is working properly
lapply(DE_miRNA_targets, function(y){
    "MA_10079279g0010" %in% y
})

test_positive <- sapply(c(DEGsOfInterest$gene, "MA_10079279g0010"), function(x){
    lapply(DE_miRNA_targets, function(y){
        x %in% y
    })
})

any(test_positive)
apply(test_positive, 2, any)

#' The code indeed reports a positive result, if there is one.  

#' Is any gene of interest, DE or not, predicted target of identified miRNA, whether it is DE or not?  

all_miRNA_targets <- lapply(miRNA_details$targets, function(x){
    y <- unlist(strsplit(as.character(x), split = ", "))
    # modify gene names
    sub(".1$", "", y)
})
names(all_miRNA_targets) <- rownames(miRNA_details)

test_all <- sapply(genesOfInterest$MA_ID, function(x) {
    lapply(all_miRNA_targets, function(y){
        x %in% y
    })
})

any(test_all)
targets_miRNAs_true_genes <- which(apply(test_all, 2, any))
# get miRNA-target pairs
targets_miRNAs_true_comb <- lapply(names(targets_miRNAs_true_genes), function(x){
    which(unlist(test_all[ , x]))
})
names(targets_miRNAs_true_comb) <- names(targets_miRNAs_true_genes)
# print them
targets_miRNAs_true_comb

#' There are four miRNA-target pairs that include genes of interest, but gene and/or miRNA are not DE  
# check if genes are differentially expressed  
any(names(targets_miRNAs_true_comb) %in% DEGsOfInterest)
#' None of the genes is differentially expressed.  

# check if miRNAs are diferentially expressed
str_extract(names(unlist(targets_miRNAs_true_comb)), "miRNA_.+") %in% DE_miRNAs_names

#' Only one miRNA is differentially expressed: miRNA_10861-5p.  
#' 

#' ### max E < 5  
#' 
#' Are genes of interest among predicted targets, if criteria for target prediction is more relaxed (max E < 5, default setting in psRNATarget)

# extract targets of DE miRNAs
DE_miRNA_targets_maxE5 <- miRNA_details_maxE5[DE_miRNAs_names, "targets"]

# change the form of the object
DE_miRNA_targets_maxE5 <- lapply(DE_miRNA_targets_maxE5, function(x){
    y <- unlist(strsplit(as.character(x), split = ", "))
    # modify gene names
    sub(".1$", "", y)
})
names(DE_miRNA_targets_maxE5) <- DE_miRNAs_names

# how many DE genes of interest is among predicted targets of DE mRNAs?
DEtargets_DEmiRNAs_maxE5 <- sapply(DEGsOfInterest$gene, function(x){
    lapply(DE_miRNA_targets_maxE5, function(y){
        x %in% y
    })
})

# are there any?
any(DEtargets_DEmiRNAs_maxE5)
# get names of the genes
DEtargets_DEmiRNAs_maxE5_true_genes <- which(apply(DEtargets_DEmiRNAs_maxE5, 2, any))

# get miRNA-target pairs
DEtargets_DEmiRNAs_maxE5_true_comb <- lapply(names(DEtargets_DEmiRNAs_maxE5_true_genes), function(x){
    which(unlist(DEtargets_DEmiRNAs_maxE5[ , x]))
})
names(DEtargets_DEmiRNAs_maxE5_true_comb) <- names(DEtargets_DEmiRNAs_maxE5_true_genes)
# print them
DEtargets_DEmiRNAs_maxE5_true_comb

#' When criteria for predicting miRNA targets is more relaxed, 3 DE genes of interest are predicted targets of DE miRNAs.  

#' Are these conserved or novel miRNAs?  
miRNA_details_maxE5[c("miRNA_31495-3p", "miRNA_16605-3p", "miRNA_1029-3p"), ]

#' One of them is novel miRNA, not found in miRBase, the other two are most similar to MIR169 and MIR159. 
#' It is often the case that novel/non-conserved miRNAs have targets that do not follow strict pairing rules, 
#' like targets of known, more conserved miRNAs.  

#' More information about these miRNA-target pairs:  
DEGsOfInterest[DEGsOfInterest$gene %in% names(DEtargets_DEmiRNAs_maxE5_true_comb), ]
names_DEmiRNAs_maxE5 <- unique(unlist(str_extract_all(DEtargets_DEmiRNAs_maxE5_true_comb, "miRNA_[0-9]+-[35]p")))
miRNAs_byStage[miRNAs_byStage$miRNA %in% names_DEmiRNAs_maxE5, ]

# expression pattern
load("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_miRNA_filtered_exp1rep2tp1_only2ndBatch_zinbNormalised.rda", verbose = TRUE)
dds_miRNA <- DESeqDataSet(miRNA_filt_zinb, ~stage)

load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/salmon_dds_genes_TEs/dds_genes.rda")

lapply(names_DEmiRNAs_maxE5, function(x){
    plotCounts(dds_miRNA, gene = x, intgroup = "stage")
})

lapply(names(DEtargets_DEmiRNAs_maxE5_true_comb), function(x){
    plotCounts(dds_genes, gene = paste0(x, ".1"), intgroup = "Stages")
})

#' Is any gene of interest, DE or not, predicted target of identified miRNA, whether it is DE or not?  

all_miRNA_targets_maxE5 <- lapply(miRNA_details_maxE5$targets, function(x){
    y <- unlist(strsplit(as.character(x), split = ", "))
    # modify gene names
    sub(".1$", "", y)
})
names(all_miRNA_targets_maxE5) <- rownames(miRNA_details_maxE5)

test_all_maxE5 <- sapply(genesOfInterest$MA_ID, function(x) {
    lapply(all_miRNA_targets_maxE5, function(y){
        x %in% y
    })
})

any(test_all_maxE5)
targets_miRNAs_true_genes_maxE5 <- which(apply(test_all_maxE5, 2, any))
# get miRNA-target pairs
targets_miRNAs_true_comb_maxE5 <- lapply(names(targets_miRNAs_true_genes_maxE5), function(x){
    which(unlist(test_all_maxE5[ , x]))
})
names(targets_miRNAs_true_comb_maxE5) <- names(targets_miRNAs_true_genes_maxE5)
# print them
targets_miRNAs_true_comb_maxE5
