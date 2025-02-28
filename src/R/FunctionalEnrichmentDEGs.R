#' ---
#' title: "Functional enrichment analysis of DEGs in spruce SE"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---  

#' # Setup  
#' 
#' Load data: DEGs between consecutive stages of somatic embryogenesis
load(file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/DE_genes_padj001_lfc05.rda")

#' Get helper files
suppressMessages(source("~/Git/UPSCb/UPSCb-common/src/R/gopher.R"))

#' # Get functional enrichment  
#'  
#' Run only the first time:  
#' all DEGs  
# enr <- lapply(res_sig_genes, function(x){
#     genes <- sub("\\.[0-9]", "", rownames(x))
#     gopher(genes=genes, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
# })
# 
# save(enr, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/DE_genes_gopher_enrichment.rda")

#' up-regulated genes  

# enr_up <- lapply(res_sig_genes, function(x){
#     x_up <- x[x$log2FoldChange > 0, ]
#     genes <- sub("\\.[0-9]", "", rownames(x_up))
#     gopher(genes=genes, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
# })

# save(enr_up, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/up_genes_gopher_enrichment.rda")

#' down-regulated genes  
# 
# enr_down <- lapply(res_sig_genes, function(x){
#     x_down <- x[x$log2FoldChange < 0, ]
#     genes <- sub("\\.[0-9]", "", rownames(x_down))
#     gopher(genes=genes, alpha = 0.05, task=list("go", "mapman"), url="pabies", endpoint = "enrichment")
# })
# 
# save(enr_down, file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/down_genes_gopher_enrichment.rda")

#' Load the results instead of running gopher again
load(file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/DE_genes_gopher_enrichment.rda")
load(file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/up_genes_gopher_enrichment.rda")
load(file = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/down_genes_gopher_enrichment.rda")


#' # Plot treemaps  
#' 
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

#' Use the function above to plot treemaps for DEGs  
# save them in the directory below
enr_path = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/enrichment/"

# gene ontology, all namespace categories
lapply(names(enr_down), function(x){
    if(!is.null(enr_down[[x]][["go"]])) {
        png(paste0(enr_path, "treemap_go_", x, "_downreg.png"), units="px", width=1920, height=1080, res=300)
        plotEnrichedTreemap(x = enr_down[[x]], 
                            enrichment = "go",
                            de = "down", 
                            title = x)
        dev.off()
    }
})

lapply(names(enr_up), function(x){
    if(!is.null(enr_up[[x]][["go"]])) {
        png(paste0(enr_path, "treemap_go_", x, "_upreg.png"), units="px", width=1920, height=1080, res=300)
        plotEnrichedTreemap(x = enr_up[[x]], 
                            enrichment = "go",
                            de = "up", 
                            title = x)
        dev.off()
    }
})


# # gene ontology, biological process
# lapply(names(enr_down), function(x){
#     if(any(enr_down[[x]][["go"]][["namespace"]]== "BP")) {
#         png(paste0(enr_path, "treemap_goBP_", x, "_downreg.png"), units="px", width=1920, height=1080, res=300)
#         plotEnrichedTreemap(x = enr_down[[x]], 
#                             enrichment = "go", 
#                             namespace = "BP",
#                             de = "down", 
#                             title = x)
#         dev.off()
#     }
# })
# 
# lapply(names(enr_up), function(x){
#     if(any(enr_up[[x]][["go"]][["namespace"]]== "BP")) {
#         png(paste0(enr_path, "treemap_goBP_", x, "_upreg.png"), units="px", width=1920, height=1080, res=300)
#         plotEnrichedTreemap(x = enr_up[[x]], 
#                             enrichment = "go", 
#                             namespace = "BP",
#                             de = "up", 
#                             title = x)
#         dev.off()
#     }
# })
  

# mapman

lapply(names(enr_down), function(x){
    if(!is.null(enr_down[[x]]$mapman)) {
        png(paste0(enr_path, "treemap_mapman_", x, "_downreg.png"), units="px", width=1920, height=1080, res=300)
        plotEnrichedTreemap(x = enr_down[[x]], 
                            enrichment = "mapman",
                            de = "down", 
                            title = x)
        dev.off()
    }
})

lapply(names(enr_up), function(x){
    if(!is.null(enr_up[[x]]$mapman)) {
        png(paste0(enr_path, "treemap_mapman_", x, "_upreg.png"), units="px", width=1920, height=1080, res=300)
        plotEnrichedTreemap(x = enr_up[[x]], 
                            enrichment = "mapman",
                            de = "up", 
                            title = x)
        dev.off()
    }
})

#' Export details of enrichment results 
lapply(names(enr_down), function(x){
     if(!is.null(enr_down[[x]][["mapman"]])) {
       write.table(enr_down[[x]][["mapman"]], 
       file = paste0(enr_path, "enr_details_mapman_", x, "_downreg.tsv"), 
       sep = "\t",
       quote = FALSE,
       row.names = FALSE)
     }
})

lapply(names(enr_up), function(x){
  if(!is.null(enr_up[[x]][["mapman"]])) {
    write.table(enr_up[[x]][["mapman"]], 
                file = paste0(enr_path, "enr_details_mapman_", x, "_upreg.tsv"), 
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
  }
})

lapply(names(enr_down), function(x){
  if(!is.null(enr_down[[x]][["go"]])) {
    write.table(enr_down[[x]][["go"]], 
                file = paste0(enr_path, "enr_details_go_", x, "_downreg.tsv"), 
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
  }
})

lapply(names(enr_up), function(x){
  if(!is.null(enr_up[[x]][["go"]])) {
    write.table(enr_up[[x]][["go"]], 
                file = paste0(enr_path, "enr_details_go_", x, "_upreg.tsv"), 
                sep = "\t",
                quote = FALSE,
                row.names = FALSE)
  }
})
