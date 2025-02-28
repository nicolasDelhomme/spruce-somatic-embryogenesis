#' ---
#' title: "Sruce somatic embryogenesis"
#' author: "Nicolas Delhomme & Camilla Canovi"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' 
#' # Setup
#' ## Environment
#' Set the working dir
setwd("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")
#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")
#' ```

#' Libraries
suppressPackageStartupMessages(library(igraph))

#' # Process
#' Read the data
mat <- read.delim("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/vala/aggregate/aggregate_threshold-0.99998.txt",header=FALSE)

#' Create the undirected graph
sel <- mat[,4]=="Undirected"
u.graf <- graph.edgelist(as.matrix(mat[sel,1:2]),directed = FALSE)
u.graf <- set_edge_attr(u.graf,name="rank",value=mat[sel,3])
u.graf <- set_edge_attr(u.graf,name="scores",value=mat[sel,5])

#' Create the directed graph
d.graf <- graph.edgelist(as.matrix(mat[!sel,1:2]),directed = TRUE)
d.graf <- set_edge_attr(d.graf,name="rank",value=mat[!sel,3])
d.graf <- set_edge_attr(d.graf,name="scores",value=mat[!sel,5])

#' Combine the graphs
graf <- d.graf+as.directed(u.graf,mode="mutual")

#' Just display the graph info
clusters(graf)

#' # Export
write_graph(graf,file="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/analysis/vala/aggregate/aggregate_threshold-0.99998.graphml",format="graphml")

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
