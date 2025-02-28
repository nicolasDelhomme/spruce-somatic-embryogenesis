#' ---
#' title: "seidr gold standard AUC"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---
#' # Setup
#' Environment
#' Set the working dir
setwd("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr")

#' ```{r set up, echo=FALSE}
#' knitr::opts_knit$set(root.dir="/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project/seidr")
#' ```

#' Libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(LSD))
suppressPackageStartupMessages(library(pander))
suppressPackageStartupMessages(library(pracma))
suppressPackageStartupMessages(library(RColorBrewer))

#' Palette
hpal <- colorRampPalette(c("blue","white","red"))(100)
cols <- brewer.pal(12,"Paired")
mar <- par("mar")

#' # Data
AUCs <- sapply(dir("roc",
                   pattern=".*_roc.tsv$",full.names=TRUE),
               function(f){
                 
                 algo <- sub("[^_]+_","",sub("_roc\\.tsv","",f))
                 
                 dat <- read.delim(f,header=FALSE,skip=1,col.names = c("TP","FP","PR"))
                 
                 head <- scan(f,nmax = 3,what = "character",sep="\t")
                 
                 message(sprintf("Found %s GS edges out of %s edges",head[2],sum(as.integer(head[2:3]))))
                 
                 auc <- round(trapz(dat[,2],dat[,1]),digits=3)
                 
                 plot(dat[,2],dat[,1],type="l",main=sprintf("%s (AUC = %s)",algo,auc),
                      xlab="False Positive Rate",ylab="True Positive Rate",
                      sub=sprintf("%s Gold Standard edges out of %s edges\n%s",head[2],sum(as.integer(head[2:3])),
                                  basename(f)))
                 
                 abline(0,1,lty=2)
                 
                 message(sprintf("The AUC is %s",auc))
                 
                 return(c(sub("_.*","",basename(f)),algo,auc))
               })

tab <- split.data.frame(t(AUCs[2:3,]),AUCs[1,])

stopifnot(nrow(unique(t(sapply(tab,function(ta){ta[,1]})))) == 1)

res <- sapply(tab,function(ta){as.numeric(ta[,2])})
rownames(res) <- tab[[1]][,1]
res <- res[,c(2:ncol(res),1)]

heatmap.2(res,
          trace="none",
          margins = c(10.1,7.1),
          col = hpal)

heatmap.2(res,Colv = FALSE,
          dendrogram = "row",
          trace="none",
          margins = c(10.1,7.1),
          col = hpal)

par(mar=c(10.1,3.1,3.1,0.1))
linesplot(res,cols = colorRampPalette(c("blue","red"))(20),addboxes = TRUE,las=2,main="AUC distribution")
par(mar=mar)

par(mar=c(10.1,4.1,3.1,2.1))
plot(0,0,ylim=c(0,1),
     xlim=c(1,ncol(res)),type="n",
     xaxt="n",xlab="",ylab="AUC",main="AUC evolution")
axis(1,1:ncol(res),las=2,labels=colnames(res))
sapply(1:nrow(res),function(i){lines(res[i,],col=cols[i],lwd=2,lty=2)})
legend("bottomright",legend = rownames(res),col=cols[1:nrow(res)],lty=1,lwd=2)
par(mar=mar)

pander(res)
write.csv(res,file="roc/AUCs.csv",quote=FALSE)

#' # Session Info
#' ```{r session info, echo=FALSE}
#' sessionInfo()
#' ```
