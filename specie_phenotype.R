#PATHOGENICITY PHENOTYPES
setwd("C:/Users/Lawal/Desktop/Data/mge-plasmidII")
list11=list.files(path = ".", pattern = "plasmid_proteins.faa.tsv", recursive = F)
col_names=c("qseqid","sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
plasmid_tables=list()
for (i in 1:length(list11)){
  file1= read_tsv(list11[i],  col_names = col_names)
  file1$Genome=paste0(list11[i])
  plasmid_tables=rbind(plasmid_tables, file1)  
  print(paste0(i))
}

df=as.data.frame(plasmid_tables)
plasmid_tables$Genome = gsub("plasmid_proteins.faa.tsv","",plasmid_tables$Genome)
plasmid_tables$Genome = gsub("_genomic_","",plasmid_tables$Genome)

plasmid_tables=plasmid_tables%>%
  separate(`sseqid`,
           into = c("mobileOG ID","Gene","Accession Number","Category","Initiation","Copy"),
           sep = "\\|", fill = "right", extra = "merge")

setwd("C:/Users/Lawal/Desktop/Data/phi-plasmidII")
list12=list.files(path = ".", pattern = "plasmid_proteins.faa.tsv", recursive = F)
col_names=c("qseqid","sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
phi_plasmid_tables=list()
for (i in 1:length(list12)){
  file1= read_tsv(list12[i],  col_names = col_names)
  file1$Genome=paste0(list12[i])
  phi_plasmid_tables=rbind(phi_plasmid_tables, file1)  
  print(paste0(i))
}

df=as.data.frame(phi_plasmid_tables)
phi_plasmid_tables$Genome = gsub("_genomic_plasmid_proteins.faa.tsv","",phi_plasmid_tables$Genome)
phi_plasmid_tables=phi_plasmid_tables%>%
  separate(`sseqid`,
           into = c("UniprotID","Accession","Gene","Hostcode","Specie","Phenotype"),
           sep = "#", fill = "right", extra = "merge")

phi_plasmid_tables$contig=gsub(".1_.*",".1",
                               phi_plasmid_tables$qseqid)


plasmid_tables$contig=gsub(".1_.*",".1",
                            plasmid_tables$qseqid)
pathogen_table= phi_plasmid_tables%>%
  filter (contig%in%plasmid_tables$contig&Genome%in%plasmid_tables$Genome)
  #select(c(Genome, Specie, Phenotype))



#To filter the phenotype for easier classification 
filtered_plasmids=pathogen_table%>%
  filter(
    str_detect(Phenotype, "reduced_virulence") | 
      str_detect(Phenotype, "loss_of_pathogenicity") | 
      str_detect(Phenotype, "plant_avirulence_determinant") 
  )

#Classification of the phenotypes 
filtered_plasmids=filtered_plasmids%>%
  mutate(pathogenic_traits = case_when(
    str_detect(Phenotype, "loss_of_pathogenicity") ~ "Essential Pathogenicity",
    str_detect(Phenotype, "reduced_virulence") ~ "Virulence Enhancer",
    str_detect(Phenotype, "plant_avirulence_determinant") ~ "Immune System Effector",
    TRUE ~ "Other"
  ))

trait_counts=filtered_plasmids %>%
  group_by(pathogenic_traits)%>%
  summarise(Count = n())


v_plot=ggplot(trait_counts, aes(x =pathogenic_traits, y = Count,  fill = pathogenic_traits)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") + 
  labs(
    y = "Number of Associated Genes", 
    x = "Functional Category"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


bubble_summary=filtered_plasmids %>%
  group_by(Specie, pathogenic_traits) %>%
  summarise(Gene_Count = n(), .groups = 'drop')

pathogenic_genes=ggplot(bubble_summary, aes(x = Specie, y = pathogenic_traits, fill = log10(Gene_Count))) +
  geom_point(size=10, shape = 21, color = "black", alpha = 0.8) +
  scale_fill_viridis_c(option = "plasma") +
  theme_bw() +
  labs(
    x = "Species where originally the gene was found in PHI-database",
    y = "",
    size = "Number of genes"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, size= 15, hjust = 1, face = "italic"),
    axis.text.y = element_text(size =15, hjust =1, face =NULL), 
    #panel.grid.major = element_line(linetype = "dashed", color = "black"),
    legend.position = "right"
  )
ggsave("pathogenic_genes.png", plot =  pathogenic_genes, width = 7, height = 5, dpi = 600)
dev.off()

#DEFENSE VS PATHOGENIC TRAITS

phenotypes =phi_tables%>%
  mutate(pathogenic_traits = case_when(
    str_detect(Phenotype, "loss_of_pathogenicity") ~ "Essential Pathogenicity",
    str_detect(Phenotype, "reduced_virulence") ~ "Virulence Enhancer",
    str_detect(Phenotype, "plant_avirulence_determinant") ~ "Immune System Effector",
    TRUE ~ "Other"))

phenotypes1=phenotypes%>% 
  filter(pathogenic_traits != "Other")
phenotypes2=dcast(phenotypes1, Genome+pathogenic_traits~., value.var = "Genome", length)

phenotypes3=full_join(select(phenotypes2, Genome2=Genome, everything()), defense_system_prime)
phenotypes4=filter(phenotypes3,DF_system_pr!="NA" )

#Correlation

path_correlation = ggplot(phenotypes4, aes(y = `.`, x =DF_system_Count, colour = pathogenic_traits)) + 
  geom_point() + geom_smooth(method = "lm")+
  theme_bw() + 
  labs(
    x = "traits", 
    y = "Defense System number per Genome"
  ) + 
  facet_wrap(~paste(type,"vs",pathogenic_traits), scales="free")+
  theme(
    strip.background = element_rect(fill = "white", colour="white"), 
    strip.text.y = element_text(angle = 0, face = "bold", size = 10), 
    strip.placement = "outside", 
    legend.position = "right", 
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10))+geom_smooth(method = "lm")+ stat_cor(method="spearman")
ggsave(file = "patho_correlation.png",plot = path_correlation, width =  7, height = 5, dpi = 600)
dev.off()

#Boxplot
path_boxplot= ggplot(phenotypes4, aes(y = `.`, x = pathogenic_traits, fill=DF_system_pr)) + 
  geom_boxplot()+
  theme_bw() + 
  labs(
    x = "traits", 
    y = "Defense System number per Genome"
  ) + 
  facet_wrap(~paste(type), scales="free")+
  theme(
    strip.background = element_rect(fill = "white", colour="white"), 
    strip.text.y = element_text(angle = 0, face = "bold", size = 10), 
    strip.placement = "outside", 
    legend.position = "right", 
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10))+geom_smooth(method = "lm")+ stat_compare_means(method = "wilcox.test")

ggsave(file = "path_boxplot.png",plot = path_boxplot, width =  7, height = 5, dpi = 600 )
dev.off()


