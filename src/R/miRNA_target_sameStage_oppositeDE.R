#' ---
#' title: "miRNA-target pairs DE in the same stage (spruce SE)"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
#' # Aim  
#' To find pairs of miRNAs and their targets that are differentially expressed at the same stage of experiment.
#' In addition explore pairs that have the opposite lfc in the same stage of the experiment. 
#' Use predicted targets when E ≤ 4. 
#
#' # Setup  
#' Load libraries
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))

#' Get helper files
suppressMessages(source("~/Git/UPSCb/UPSCb-common/src/R/gopher.R"))

#' Load data  
# miRNAs
load("~/Git/UPSCb/projects/spruce-srna/doc/DEmiRNAs_timeline_2ndBatch.rda")
DE_miRNA_dds <- res_sig_timeline_W
rm(res_sig_timeline_W)

miRNA_details <- read.csv2("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/Pabies_SE_miRNA_results_all_maxE_4.csv", 
                           row.names = 1)

# DEGs
load("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DE_genes_padj001_lfc05.rda")



#' # Analysis  
#' ## Transform data  
#' 
# DEGs
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

write_tsv(df_DEGs, 
          path = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/DEGs_padj_lfc_byStage.tsv")

# miRNAs
list_df_miRNA <- lapply(DE_miRNA_dds, function(x){
    data.frame(x[ , c("log2FoldChange", "padj")])
})

df_miRNAs <- list_df_miRNA %>% 
    map(rownames_to_column, "miRNA") %>%
    purrr::reduce(full_join, by = "miRNA")

colnames(df_miRNAs) <- c("miRNA",
                       "lfc_2v1", "padj_2v1", 
                       "lfc_3v2", "padj_3v2", 
                       "lfc_4v3", "padj_4v3",
                       "lfc_5v4", "padj_5v4", 
                       "lfc_6v5", "padj_6v5", 
                       "lfc_7v6", "padj_7v6", 
                       "lfc_8v7", "padj_8v7")

stage_names <- c("2v1", "3v2", "4v3", "5v4", "6v5", "7v6", "8v7")

#' ## Extract targets of DE miRNAs

# extract names of DE miRNAs
DE_miRNAs_names <- unique(unlist(lapply(DE_miRNA_dds, rownames)))

# extract targets of DE miRNAs
DE_miRNA_targets <- miRNA_details[DE_miRNAs_names, "targets"]

# exclude TEs, as they are not relevant for this project and modify gene names
DE_miRNA_targets <- lapply(DE_miRNA_targets, function(x){
    y <- unlist(strsplit(as.character(x), split = ", "))
    # keep only genes, not TEs
    y <- y[grepl("^MA_", y)]
    # modify gene names
    sub(".1$", "", y)
})
names(DE_miRNA_targets) <- DE_miRNAs_names

#' ## miRNAs and genes DE in the same stage

DEmiRNA_DEtarget_sameStage <- lapply(stage_names, function(x){
    # prepare
    collfc <- paste0("lfc_", x)
    colpadj <- paste0("padj_", x)
    stage_miRNAs <- df_miRNAs$miRNA[!is.na(df_miRNAs[, collfc])]
    stage_DEGs <- df_DEGs$gene[!is.na(df_DEGs[, collfc])]
    # find DE targets of each DE miRNA in the same stage
    if(!isEmpty(stage_miRNAs)){
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
        }
})

names(DEmiRNA_DEtarget_sameStage) <- stage_names

DEmiRNA_DEtarget_sameStage_opposite <- lapply(DEmiRNA_DEtarget_sameStage, function(x){
    x[x["DE_pattern"] == "opposite", ]
})

# save object
save(DEmiRNA_DEtarget_sameStage, file = "/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/DEmiRNA_DEtarget_sameStage.rda")

# all DE pairs from the same stage
lapply(names(DEmiRNA_DEtarget_sameStage), function(x){
    write_tsv(data.frame(DEmiRNA_DEtarget_sameStage[[x]]), 
              path = paste0("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/DEmiRNA_DEtarget_sameStage_padj_lfc_", x, ".tsv"), 
              col_names = TRUE)
})

# only DE pairs with the opposite expression pattern
lapply(names(DEmiRNA_DEtarget_sameStage_opposite), function(x){
    write_tsv(data.frame(DEmiRNA_DEtarget_sameStage_opposite[[x]]), 
              path = paste0("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/DEmiRNA_DEtarget_sameStage_padj_lfc_", x, ".tsv"), 
              col_names = TRUE)
})
#' ## miRNA-target pairs with the opposite pattern  
#' 
#' Some statistics about miRNA-target pairs with the opposite DE pattern  
#' 
#' How many miRNA-target pairs have an opposite pattern of DE in each stage?
stats_opp_pairs <- sapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
    nr_pairs <- nrow(x)
    nr_miRNAs <- nrow(unique(x["DE_miRNA"]))
    nr_targets <- nrow(unique(x["DE_target"]))
    df <- rbind(nr_pairs, nr_miRNAs, nr_targets)
    # why the code below does not work?
    # rownames(df) <- c("nr_pairs", "nr_miRNAs", "nr_targets")
    if(!is.null(df)){
        return(df)
    }else{
        df0 <- rbind("nr_pairs" = 0, "nr_miRNAs" = 0, "nr_targets" = 0)
        return(df0)
    }
})
rownames(stats_opp_pairs) <- c("nr_pairs", "nr_miRNAs", "nr_targets")

#' Calculate the same stats for the whole experiment  
# How many miRNAs are represented in these pairs in the whole experiment?
nr_miRNAs_exp <- length(unique(unlist(lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
    x["DE_miRNA"]
}))))

# How many DEGs are represented in these pairs in the whole experiment?
nr_targets_exp <- length(unique(unlist(lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
    x["DE_target"]
}))))

# How many miRNA-target pairs have the opposite expression pattern in the whole experiment? (count unique, non redundant)
nr_pairs_exp <- nrow(unique(do.call(rbind,
                                           lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
                                               x[c("DE_miRNA", "DE_target")]
                                           })
)))

# Add to df
stats_opp_pairs <- cbind(stats_opp_pairs, all = c(nr_pairs_exp, nr_miRNAs_exp, nr_targets_exp))
stats_opp_pairs

# How to check for all the targets with the opposite DE pattern of a miRNA? E.g.:
lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
    x[x["DE_miRNA"] == "miRNA_10861-5p", ]
})

# How to check for all the targets with the opposite DE pattern of a miRNA based on best hit in miRBase:  
# create function
getOppDETarg <- function(mirna){
miRNAID <- rownames(miRNA_details[grep(mirna, miRNA_details$best.miRBase), ])
miRNAID_DE <- miRNAID[miRNAID %in% DE_miRNAs_names]
miRNAID_opptarg <- lapply(miRNAID_DE, function(mID){
    opptarg <- lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
            x[x["DE_miRNA"] == mID, ]
    })
})
names(miRNAID_opptarg) <- miRNAID_DE
return(miRNAID_opptarg)
}

getOppDETarg("MIR390")

# How to check for all the miRNAs with the opposite DE pattern to a gene? E.g.:
lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
    x[x["DE_target"] == "MA_181770g0010", ]
})

#' Extract names of the targets in each stage  
DEtargets_opp <- lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
    unique(x["DE_target"])
})
DEtargets_opp_all <- unique(unlist(DEtargets_opp))

# write them in the file
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

#' Extract names of the miRNAs in each stage  
DEmiRNAs_opp <- lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
    unique(x["DE_miRNA"])
})
DEmiRNAs_opp_all <- unique(unlist(DEmiRNAs_opp))

# write them in the file
lapply(names(DEmiRNAs_opp), function(x){
    write.table(DEmiRNAs_opp[x],
                file = paste0("/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/miRNAs_maxE4_oppDE_", x, ".csv"), 
                row.names = FALSE, 
                col.names = FALSE,
                quote = FALSE)
})

write.table(DEmiRNAs_opp_all,
            file = "/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/miRNAs_maxE4_oppDE_all.csv", 
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)

#' Extract all the pairs with the opposite expression pattern regardless of the stage
DEpairs_opp_all <- unique(do.call(rbind,
                         lapply(DEmiRNA_DEtarget_sameStage_opposite, function(x){
                             x[c("DE_miRNA", "DE_target")]
                         })
))

write.table(DEpairs_opp_all,
            file = "/mnt/picea/projects/spruce/nstreet/sRNA/SE/ShortStack_genome/target_prediction/targets_with_opposite_sameStageDE_maxE4/pairs_maxE4_oppDE_all.csv", 
            row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)


#' Functional enrichment/function of targets with the opposite expression pattern from miRNAs  
#' There should be enough genes to perform functional enrichment analysis, therefore check only stages, where number of genes ≥ 30.
#' For stages with less than 30 genes, only check function of all the genes
enr_opp <- lapply(DEtargets_opp, function(x){
        if(!is.null(x)){
            genes <- x$DE_target
            if(length(genes >= 30)){
                gopher(genes=genes, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
 #           }else{
 #               gopher(genes=genes, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "gene-to-term")
            }
        }
    })

enr_opp_all <- gopher(genes=DEtargets_opp_all, 
                      alpha = 0.05, 
                      task=list("go", "mapman"), 
                      url="pabies", 
                      endpoint = "enrichment")

save(enr_opp, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/DE_miRNAtargets_sameStage_oppDE_stages_gopher_enrichment.rda")
save(enr_opp_all, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/DE_miRNAtargets_sameStage_oppDE_all_gopher_enrichment.rda")

#' plot treemaps  
#' Below is code for plotting treemaps from Alonso:
clusterTreemapColors <- rep(c("#9E1F63","#662D91","#2B3990","#1B75BC","#27AAE1",
                              "#2AB592","#035E31","#009444","#65BC46","#A5CE42",
                              "#F9ED32","#FBB040","#F15A29","#EF4136","#BE1E2D"),15)
clusterTreemapText <- rep(c("white","white","white","white","white",
                            "white","white","white","white","white",
                            "black","black","white","white","white"),15)

plotEnrichedTreemap <- function(x, enrichment = c('go','mapman', 'kegg', 'pfam', 'ko_pathway', 'ko', 'kog','cog'), 
                                namespace = c('none', 'BP', 'MF', 'CC'), 
                                title = "", 
                                de = c("none", "up", "down"), 
                                clusterColor = "#9E1F63", clusterText='black',
                                nameCol="name",
                                namespaceCol="namespace",
                                sizeCol="padj",
                                colorCol="padj", 
                                convertSize=TRUE,legend=TRUE) {
    
    require(treemap)
    require(RColorBrewer)
    
    enrichment <- match.arg(enrichment)
    namespace <- match.arg(namespace)
    de <- match.arg(de)
    
    enrData <- if (is.data.frame(x)){
        x 
    } else if (is.list(x)) {
        x[[enrichment]]
    } else {
        stop("x object must be either a list or a data.frame")
    }
    
    # calculate the size based on padj
    
    enrData$size <- if(convertSize) {
        abs(log10(enrData[[colorCol]]))
    } else {
        enrData[[colorCol]]
    }
    
    #default treemap
    index = nameCol
    fontcolor.labels=clusterText
    fontface.labels=c(2)
    fontsize.labels=c(12)
    inflate.labels = TRUE
    align.labels=list(c("center", "center"))
    position.legend <- if (!legend) {"none"} else {"bottom"}
    #vColor = 'name'
    border.col='black'
    border.lwds=c(4,2)
    palette <- colorRampPalette(c("white", clusterColor))
    palette <- palette(10)
    vColor = c("size")
    if (sizeCol=="padj") {
        vSize = vColor
    } else {
        vSize = sizeCol
    }
    type = "value"
    title.legend="abs(log10(pAdj))"
    bg.labels= 0
    
    if(enrichment =='go' ){
        if(namespace=='none') {
            index = c('namespace', 'name')
            palette = "Set1"
            
            fontcolor.labels=c("#FF000000", "black")
            fontface.labels=c(1, 2)
            fontsize.labels=c(1, 20)
            inflate.labels=TRUE
            align.labels=list(c("left", "top"), c("center", "center") )
            type = "index"
            border.col=c('black', 'white')
            title.legend="GO namespace"
        } else {
            enrData <- enrData[enrData[[namespaceCol]]==namespace,]
        }
    }
    
    # mapman name fix for better visualization
    if(enrichment =='mapman') {
        enrData[[nameCol]] <- gsub("\\."," ",enrData[[nameCol]])
    } 
    
    # Paint it properly if it is up or down regulated
    if(de !='none') {
        if(de=='up') {
            palette = "OrRd"
        }
        if(de=='down') {
            palette = "GnBu"
        }
    } 
    
    # generate treemap
    treemap(enrData, 
            index = index,
            vSize = vSize, 
            palette = palette,
            type = type, 
            vColor = vColor,
            title=title, 
            fontcolor.labels=fontcolor.labels, 
            fontface.labels=fontface.labels,     
            #fontsize.labels=fontsize.labels,
            bg.labels=bg.labels, 
            inflate.labels = inflate.labels ,
            lowerbound.cex.labels = 0, 
            position.legend = position.legend,
            border.col=border.col,         
            border.lwds=border.lwds,
            align.labels=align.labels,
            title.legend=title.legend,
            overlap.labels = 1
    )
}

#' Use the function above to plot treemaps for miRNA targets with the opposite expression as miRNA in the same stage  
# save them in the directory below
enr_path = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/"

#' per stage  
#' 
# gene ontology, all namespace categories
lapply(names(enr_opp), function(x){
    if(!is.null(enr_opp[[x]][["go"]])) {
        png(paste0(enr_path, "treemap_go_", x, "_oppmiRNA.png"), units="px", width=1920, height=1080, res=300)
        plotEnrichedTreemap(x = enr_opp[[x]], 
                            enrichment = "go",
                            de = "none", 
                            title = x)
        dev.off()
    }
})

# mapman
lapply(names(enr_opp), function(x){
    if(!is.null(enr_opp[[x]][["mapman"]])) {
        png(paste0(enr_path, "treemap_mapman_", x, "_oppmiRNA.png"), units="px", width=1920, height=1080, res=300)
        plotEnrichedTreemap(x = enr_opp[[x]], 
                            enrichment = "mapman",
                            de = "none", 
                            title = x)
        dev.off()
    }
})

#' all targets with the opposite expression pattern in the experiment  
#' 
# gene ontology, all namespace categories
png(paste0(enr_path, "treemap_go_all_oppmiRNA.png"), units="px", width=1920, height=1080, res=300)
plotEnrichedTreemap(x = enr_opp_all,
                    enrichment = "go",
                    de = "none", 
                    title = "all targets with opposite DE to miRNAs")
dev.off()

# mapman
png(paste0(enr_path, "treemap_mapman_all_oppmiRNA.png"), units="px", width=1920, height=1080, res=300)
plotEnrichedTreemap(x = enr_opp_all,
                    enrichment = "mapman",
                    de = "none", 
                    title = "all targets with opposite DE to miRNAs")
dev.off()

#' Export enrichment results  
#' 
lapply(names(enr_opp), function(x){
if(!is.null(enr_opp[[x]][["go"]])) {
    write.table(enr_opp[[x]][["go"]], 
                file = paste0(enr_path, "enr_details_go_", x, "_oppmiRNA.tsv"), 
                sep = "\t",
                quote = FALSE, 
                row.names = FALSE)
    }
})

lapply(names(enr_opp), function(x){
    if(!is.null(enr_opp[[x]][["mapman"]])) {
        write.table(enr_opp[[x]][["mapman"]], 
                    file = paste0(enr_path, "enr_details_mapman_", x, "_oppmiRNA.tsv"), 
                    sep = "\t",
                    quote = FALSE,
                    row.names = FALSE)
    }
})

write.table(enr_opp_all["go"],
            file = paste0(enr_path, "enr_details_go_all_oppmiRNA.tsv"), 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(enr_opp_all["mapman"],
            file = paste0(enr_path, "enr_details_mapman_all_oppmiRNA.tsv"), 
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

#' CHECK In how many InfoMap clusters can these targets be found?
#' There are XX targets in XX of the clusters with >= 100 genes.


#' Session info
sessionInfo()
