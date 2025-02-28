cl <- cluster.genes_list[[1]]

setwd("/mnt/picea/projects/spruce/uegertsdotter/22_Somatic_Embryogenesis_Project")

imap <- read.delim("analysis/vala/aggregate/Subset/aggregate_threshold-0.99998.txt.final.h.txt",as.is = TRUE)

str(imap)
imap$ID <- sub("\\.1$","", imap$ID)
table1 <- table(imap$C2[imap$ID %in% cl])

round(table1/sum(table1)*100,digits=2)

tab <- table(imap$C2)

round(tab/sum(tab)*100,digits=2)[1]



table2 <- table(paste(imap$C2,imap$C3,sep=".")[imap$ID %in% cl])

tab2 <- round(table2/sum(table2)*100,digits=2)

sort(tab2,decreasing = TRUE)[1:6]

tab <- table(paste(imap$C2,imap$C3,sep="."))

round(tab/sum(tab)*100,digits=2)["3.1"]

round(tab/sum(tab)*100,digits=2)["1.9"]

fisher.test(matrix(ncol=2,c(5.96,0.2,94.04,99.8)))


pop <- table(paste(imap$C2,imap$C3,sep="."))

checkCluster <- function(cl,pop,fdr=.01){
  set <- table(paste(imap$C2,imap$C3,sep=".")[imap$ID %in% cl])
  
  local.pop <- pop[names(set)]
  
  df <- as.data.frame(t(sapply(as.character(names(local.pop)),function(n){
    message(sprintf("processing %s",n))
    ft <- fisher.test(matrix(ncol=2,
                             c(set[n],
                               pop[n],
                               sum(set)-set[n],
                               sum(pop)-pop[n]
                               )))
    data.frame(clID=n,
               cluster.set=set[n],
               cluster.pop=pop[n],
               set=sum(set),
               pop=sum(pop),
               p.val=ft$p.value,
               low.conf.int=ft$conf.int[1],
               up.conf.int=ft$conf.int[2],
               estimate=ft$estimate,
               stringsAsFactors = FALSE
               )
  })))

  df$adj <- p.adjust(df$p.val,method="BH",nrow(df))
    
  df[df$adj<= fdr,]
  
}

all_clusters <- lapply(cluster.genes_list,checkCluster,pop)
