#PLASMID-ANTIPHAGE SYSTEMS
setwd("C:/Users/Lawal/Desktop/Data/plasmid-antiphage-systems")

list13=list.files(path = ".", pattern = "genomic_plasmid_proteins_defense_finder_systems.tsv", recursive = T)
plasmid_defense_finder_tables=list()
for (i in 1:length(list13)){
  file1= read_table(paste0(list13[i]))
  file1$Genome=paste0(list13[i])
  plasmid_defense_finder_tables=rbind(plasmid_defense_finder_tables, file1)  
  print(paste0(i))
}
df=as.data.frame(plasmid_defense_finder_tables)
plasmid_defense_finder_tables$Genome= gsub("_genomic.*$", "",plasmid_defense_finder_tables$Genome)
plasmid_defense_finder_tables$Genome2= gsub("_genomic.*$", "",plasmid_defense_finder_tables$Genome)

plasmid_defense_systems_prime=dcast(plasmid_defense_finder_tables, Genome~type, value.var = "type", length)%>%gather(-Genome, value = "Count", key="Type")

systems=c("Cas","Dnd", "Mok_Hok_Sok", "Prometheus", "Retron", "RM")

freq_dataIII=plasmid_defense_systems_prime%>%
  filter(Type %in% systems) %>%
  group_by(Type) %>%
  summarise(Total_Count=sum(Count > 0))

ggplot(freq_dataIII, aes(x = Type, y = Total_Count, fill = Type)) + geom_bar(stat = "identity", color = "black", width = 0.7) + geom_text(aes(label = Total_Count), vjust = -0.5, fontface = "bold") + theme_classic() +
  scale_fill_brewer(palette = "Set3") + labs(x = "Defense System Type", y = "Plasmid Proteins") +theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),legend.position = "none")