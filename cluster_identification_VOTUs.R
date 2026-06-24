#Import libraries 
library(readr)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(stringi)
library(dplyr)
library(tidyr)
library(stringr)
library(igraph)
library(Biostrings)
library(cluster)
library(tidyverse)


#VOTU cluster_identification 

setwd("C:/Users/Lawal/Desktop/Data/predicted_virus_genomes/vOTU")

cluster_table=read_tsv("clusterRes_viromes_cluster.tsv", col_names = c("x1", "x2"))
graph1=graph_from_data_frame(cluster_table, directed = F)
cluster1=igraph::components(graph1)
check1=as.data.frame(cluster1$membership)%>%rownames_to_column(var="Contig")
clusters_final=check1
clusters_final$VOTUs=paste0("vOTU", clusters_final$`cluster1$membership`)

cluster_ids=gsub("^.*(GC[AF]_.*)$", "\\1", clusters_final$Contig)
Genome_ids=gsub("_genomic_virus.*$", "", viral_genome$Genome)
short_genome_ids=gsub("^(GC[AF]_[0-9]+\\.[0-9]+)_.*$", "\\1", Genome_ids)
vOTU_Genome=match(cluster_ids, short_genome_ids)

clusters_final$Genome2=Genome_ids[vOTU_Genome]

vOTU_table=full_join(Quality_table,clusters_final, by = "Genome2")

vOTU_tablee=filter(vOTU_table)%>%dcast(VOTUs+Genome2~., value.var = "Genome2", length)%>%
  select(Viral_Count=".", Genome_Count=".", VOTUs, Genome2)
vOTU_tablee$VOTUs[is.na(vOTU_tablee$VOTUs)]=0