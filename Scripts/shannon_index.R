#!/usr/bin/env Rscript


library(vegan)
#library(phyloseq)

data <- read.table("../Data/pathways_species2.tsv" , header = TRUE , sep = "\t")

rutas <- unique(data$Pathway)


dist <- list()


for (i in 1:length(rutas)){
  ruta_i <- which(data$Pathway == rutas[i])
  data_i <- data[ruta_i , ]
  shannon_i <- c()
  for (j in 3:dim(data)[2]){
    s_j <- diversity(data_i[,j])
    shannon_i <- c(shannon_i , s_j)
  }
  dist[[i]] <- shannon_i
  #print(length(shannon_i))
}

distribuciones <- matrix(dist[[1]] , ncol = length(dist[[1]]))
for (i in 2:length(rutas)){
  distribuciones <- rbind(distribuciones , dist[[i]])
}

distribuciones <- as.data.frame(distribuciones)
colnames(distribuciones) <- colnames(data)[3:dim(data)[2]]
row.names(distribuciones) <- rutas

write.table( distribuciones , file = "pathway_shannon_index.tsv" , sep = "\t" )

