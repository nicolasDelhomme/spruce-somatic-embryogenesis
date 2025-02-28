#' # Aim: add description of DE genes from TAIR database  

#' # Load data  
#' Annotation

araport11_functional_descriptions <- read.delim("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/Araport11_functional_descriptions_20190930.txt", 
                                                stringsAsFactors=FALSE, quote ="" )
araport11_functional_descriptions[grepl("\\.1",araport11_functional_descriptions$name),]
gene_aliases <- read.delim("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2/gene_aliases_20190930.txt", 
                           quote="", stringsAsFactors=FALSE)

#' List of DEGs
files <- dir(path = "/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis_salmon/DESeq2", 
             pattern = "*_padj001_lfc05_annotated.tsv", 
             full.names = TRUE)
          
#' # Add TAIR info  
#' Write in the files with new name (*_TAIR.tsv)
sapply(files, function(file){
    print(file)
    theFile <- read.delim(file, row.names=1, stringsAsFactors=FALSE)
    tempAT <- theFile$Best.BLAST.Arabidopsis
    tempAT <- gsub("\\..*","",tempAT)
    theFile$otherNames <- sapply(tempAT,  function(gene){
        #print(gene)
        x <- gene_aliases[gene_aliases$name == gene,]$full_name
        y <- gene_aliases[gene_aliases$name == gene,]$symbol
        z <- union(x, y)
        z <- z[z != ""]
        z <- paste0(z, collapse="|")
        #print(x)
    })
    daRows <- rownames(theFile) 
    #genedescription, araport11 has unique names
    theFile <- left_join(theFile, araport11_functional_descriptions, by = c("Best.BLAST.Arabidopsis" = "name"))
    rownames(theFile) <-  daRows
    filename <- gsub("\\.tsv", "_TAIR\\.tsv", file)
    write.table(theFile, file=filename, sep='\t', col.names = NA)
})
