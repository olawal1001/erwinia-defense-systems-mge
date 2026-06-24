#PHAGE-ANTIPHAGE SYSTEM 
setwd("C:/Users/Lawal/Desktop/Data/checkv_defensefinder")
phage_defense_finder_table=list()
listC= list.files(path = ".", pattern = "viruses_defense_finder_systems.tsv", recursive = T)
phage_defense_finder_tables=list()
for (i in 1:length(listC)){
  file1= read_table(paste0(listC[i]))
  file1$Genome=paste0(listC[i])
  phage_defense_finder_table=rbind(phage_defense_finder_table, file1)  
  print(paste0(i))
}

df=as.data.frame(phage_defense_finder_table)
phage_defense_finder_table$Genome2= gsub("*viruses_defense_finder_systems.*","" , gsub("_genomic_virus/.*$", "", gsub("DefenseFinder.", "",phage_defense_finder_table$Genome)))

phage_defense_system_prime=dcast(phage_defense_finder_table, Genome2~type, value.var = "type", length)%>%gather(-Genome2, value = "Count", key="Type")
systems=c("Cas","Dnd", "Mok_Hok_Sok", "Prometheus", "Retron", "RM")

freq_dataII=phage_defense_system_prime%>%
  filter(Type %in% systems) %>%
  group_by(Type) %>%
  summarise(Total_Count=sum(Count > 0))
ggplot(freq_dataII, aes(x = Type, y = Total_Count, fill = Type)) +geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = Total_Count), vjust = -0.5, fontface = "bold") + theme_classic() + scale_fill_brewer(palette = "Set2") + labs(
    x = "Defense System Type",y = "Phages present in genomes") +theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
                                                                      legend.position = "none")

