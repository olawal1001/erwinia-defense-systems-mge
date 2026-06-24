#PATHOGENICITY PHENOTYPES 

library(treemapify)

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

combined_plasmid_table=plasmid_tables%>% 
  full_join(phi_plasmid_tables, by = c("Genome"))


filtered_plasmids=combined_plasmid_table%>%
  filter(
    str_detect(Phenotype, "reduced_virulence") | 
      str_detect(Phenotype, "loss_of_pathogenicity") | 
      str_detect(Phenotype, "plant_avirulence_determinant") 
  )

filtered_plasmids=filtered_plasmids%>%
  mutate(Summary = case_when(
    str_detect(Phenotype, "loss_of_pathogenicity") ~ "Essential Pathogenicity",
    str_detect(Phenotype, "reduced_virulence") ~ "Virulence Enhancer",
    str_detect(Phenotype, "plant_avirulence_determinant") ~ "Immune System Effector",
    TRUE ~ "Other"
  ))


summary_counts=filtered_plasmids %>%
  group_by(Summary)%>%
  summarise(Count = n())

#chi_test=chisq.test(summary_counts$Count)
#p_val_text=paste("p-value:", format.pval(chi_test$p.value))


v_plot <- ggplot(summary_counts, aes(x =Summary, y = Count,  fill = Summary)) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2") + 
  labs(
    #title = "Plasmid Pathogenicity Phenotypes",
    #subtitle = p_val_text,
    y = "Number of Genes", 
    x = "Functional Category"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Grouping and counting
bubble_summary <- filtered_plasmids %>%
  group_by(Specie, Summary) %>%
  summarise(Gene_Count = n(), .groups = 'drop')

# The Plot
bubble_plot <- ggplot(bubble_summary, aes(x = Specie, y = Summary)) +
  # Use size for the count and color for the category
  geom_point(aes(size = Gene_Count, fill = Summary), shape = 21, color = "black", alpha = 0.8) +
  # Customizing the bubble sizes
  scale_size_continuous(range = c(2, 12), breaks = c(1, 5, 10, 20, 50)) +
  scale_fill_viridis_d(option = "plasma") + # Professional color palette
  theme_light() +
  labs(
    #title = "Plasmid Virulence Reservoir",
    #subtitle = "Size indicates number of genes; Colors represent functional categories",
    x = "Species (PHI-base prediction)",
    y = "Predicted Phenotype",
    size = "Gene Count"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    panel.grid.major = element_line(linetype = "dashed", color = "black"),
    legend.position = "right"
  )

#DEFENSE VS PATHOGENIC TRAITS
df_patho = full_join(filtered_plasmids, defense_system_prime, by =c("Genome" ="Genome2"))%>%
  select(Summary, type, DF_system_Count, Viral_Count)

correlation=full_join(filtered_plasmids, defense_system_prime, by = c("Genome" = "Genome2"))%>%
  select(Genome, Summary, type) %>%
  filter(!is.na(Summary), !is.na(type)) %>%
  group_by(Genome, Summary, type) %>%
  summarize(Count = n(), .groups = "drop") %>% 
  group_by(Summary, type) %>%
  summarize(
    rho = cor(Count, as.numeric(factor(Genome)), method = "spearman", use = "complete.obs"),
    p_val = cor.test(Count, as.numeric(factor(Genome)), method = "spearman", exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    stars = case_when(
      p_val < 0.001 ~ "***", 
      p_val < 0.01  ~ "**", 
      p_val < 0.05  ~ "*", 
      TRUE          ~ "ns" 
    ),
    label_text = paste0("R = ", round(rho, 2), " (", stars, ")")
  )

print(correlation)
correlation$Summary <- factor(correlation$Summary, 
                              levels = c("Immune System Effector", "Virulence Enhancer", "Essential Pathogenicity"))

# 2. Plotting Script with Asterisks positioned above points
ggplot(correlation, aes(x = Summary, y = rho)) +
  # Draw the individual calculated correlation points
  geom_point(aes(fill = Summary), size = 4.5, shape = 21, color = "black", show.legend = FALSE) + 
  
  # Connect the points with a trend trajectory line within each panel
  geom_smooth(method = "lm", se = FALSE, aes(group = 1), #color = "Trend Line"), 
              linewidth = 0.8, linetype = "solid", alpha = 0.7, show.legend = TRUE) + 
  
  # 2. Plots your points cleanly without breaking the grid grouping
  geom_point(aes(fill = Summary), size = 4.5, shape = 21, color = "black", show.legend = TRUE) +
  # FIXED: Places the asterisks (stars) directly above each individual point
  geom_text(aes(label = stars), vjust = -0.6, fontface = "plain", size = 5.5, color = "black") +
  
  
  # Separate by Defense System type
  facet_wrap(~type, ncol = 3) + 
  
  # Colors and styling
  scale_fill_manual(values = c("Immune System Effector" = "maroon", 
                               "Virulence Enhancer" = "#99D8C8", 
                               "Essential Pathogenicity" = "#FFFDC9")) +
  # Expand Y limits slightly so top asterisks don't hit the box roof
  scale_y_continuous(limits = c(-0.42, -0.15)) +
  theme_bw(base_size = 12) +
  labs(
    x = "Pathogenicity Traits",
    y = "Spearman Correlation Coefficient",
    #title = "Immune Strain Trends Across Virulence Classes"
  ) +
  theme(
    plot.title = element_text(face = "plain", hjust = 0.5, size = 13, margin = margin(b=12)),
       axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, face = "plain", color = "black", margin = margin(t = 5)),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "#E0ECF4", color = "black"),
    strip.text = element_text(face = "plain", size = 11),
    panel.grid.minor = element_blank(),
    plot.margin = margin(t = 10, r = 10, b = 25, l = 10, unit = "pt")
  )