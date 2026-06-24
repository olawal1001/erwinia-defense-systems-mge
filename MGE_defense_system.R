#MOBILE-GENETIC-ELEMENTS
setwd("C:/Users/Lawal/Desktop/Data/mobile_genetic_elements")
listP=list.files(path = ".", pattern = "genomic_prots.faa.tsv", recursive = F)
col_names=c("qseqid","sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
mge_tables=list()
for (i in 1:length(listP)){
  file1= read_tsv(listP[i],  col_names = col_names)
  file1$Genome=paste0(listP[i])
  mge_tables=rbind(mge_tables, file1)  
  print(paste0(i))
}
df=as.data.frame(mge_tables)
mge_tables=mge_tables%>%
  separate(`sseqid`,
           into = c("mobileOG ID","Gene","Accession Number","Category","Initiation","Copy"),
           sep = "\\|", fill = "right", extra = "merge")
df_mge=mge_tables%>%
  mutate(Class=case_when(
    grepl("tnp|transpos|nt|insertion sequence|IS", Gene,ignore.case = TRUE)~ "Transposons",
    grepl("Plasmid", Copy, ignore.case = TRUE )~ "Plasmids",
    grepl("integrase|intI|recombinase|xer|excision", Category, ignore.case = TRUE )~ "Integrases", 
    TRUE~"Other Group"
  ))
systems=c("Cas","Dnd", "Mok_Hok_Sok", "Prometheus", "Retron", "RM")
classes=c("Plasmids", "Transposons", "Integrases")
df_mge=mge_tables%>%
  mutate(Class=case_when(
    grepl("tnp|transpos|nt|insertion sequence|IS", Gene,ignore.case = TRUE)~ "Transposons",
    grepl("Plasmid", Copy, ignore.case = TRUE )~ "Plasmids",
    grepl("integrase|intI|recombinase|xer|excision", Category, ignore.case = TRUE )~ "Integrases", 
    TRUE~"Other Group"
  ))
systems=c("Cas","Dnd", "Mok_Hok_Sok", "Prometheus", "Retron", "RM")
classes=c("Plasmids", "Transposons", "Integrases")
df_counts= defense_finder_tables %>%
  group_by(Genome2) %>%
  summarise(
    Cas = sum(type == "Cas", na.rm = TRUE),
    RM  = sum(type == "RM", na.rm = TRUE),
    Retron = sum(type == "Retron", na.rm = TRUE),
    Dnd = sum(type == "Dnd", na.rm = TRUE),
    Mok_Hok_Sok = sum(type == "Mok_Hok_Sok", na.rm = TRUE),
    Prometheus = sum(type == "Prometheus", na.rm = TRUE)
  )

mge_counts = df_mge %>%
  group_by(Genome2) %>%
  summarise(
    Plasmids = sum(Class == "Plasmids", na.rm = TRUE),
    Transposons = sum(Class == "Transposons", na.rm = TRUE),
    Integrases = sum(Class == "Integrases", na.rm = TRUE)
  )
mge_counts=inner_join(df_counts, mge_counts, by = "Genome2")

stats_list = list()

for (mge in classes) {
  for (sys in systems) {
    
    x = mge_counts[[mge]]
    y = mge_counts[[sys]]
    
    group = ifelse(x > 0, "Present", "Absent")
    
    if (length(unique(group)) < 2) {
      message(paste("Skipping", mge, "-", sys, ": All genomes are either all Present or all Absent"))
      next
      }
    group_f = as.factor(group)
    testI = shapiro.test(resid(lm(y ~ group_f)))
    if (testI$p.value < 0.05) {
      # Non-normal: Use Rank-based T-test
      testII = t.test(base::rank(y) ~ group_f)
    } else {
      # Data is normal: Check log transform (using y+1 to avoid log(0))
      test_log = shapiro.test(resid(lm(log10(y + 1) ~ group_f)))
      
      if (test_log$p.value > 0.05) {
        testII = t.test(log10(y + 1) ~ group_f)
      } else {
        testII = t.test(y ~ group_f)
      }
    }
    stats_list[[paste(mge, sys, sep="_")]] <- data.frame(
      MGE_Class = mge,
      Defense_System = sys,
      P_Value = testII$p.value
    )
  }
}


stats_results <- do.call(rbind, stats_list)
stats_results$p_adj = p.adjust(stats_results$P_Value, method = "BH")
stats_results$Significance = ifelse(stats_results$p_adj < 0.05, "Significant", "NS")

print(stats_results)

mge_correlation = function(df, x, y) { 
  res <- expand.grid(MGE = x, Defense_Systems = y, stringsAsFactors = FALSE)
  stats <- mapply(function(m, d) {
    test <- cor.test(df[[m]], df[[d]], method = "spearman", exact = FALSE)
    return(c(R = test$estimate, p = test$p.value))
  }, res$MGE, res$Defense_Systems)
  res$R <- stats["R.rho", ]
  res$p_val <- stats["p", ]
  res %>%
    mutate(
      p_adj = p.adjust(p_val, method = "BH"),
      stars = cut(p_adj, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), 
                  labels = c("***", "**", "*", ""))
    )
}
data_mge = mge_correlation(mge_counts, classes, systems)

ggplot(data_mge, aes(x = MGE, y = Defense_Systems, fill = R)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = paste0(round(R, 2), "\n", stars)), 
            color = "black", size = 3.5, fontface = "bold") +
  scale_fill_gradient2(low = "#99D8C9", 
                       mid = "#FFFFCC", 
                       high = "orange", 
                       midpoint = 0, 
                       limit = c(-1, 1),
                       name = "Spearman") +
  labs(
       x = "Mobile Genetic Elements",
       y = "Defense Systems") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

data_mge%>% 
  select(MGE, Defense_Systems, R, p_val, p_adj, stars) %>%
  arrange(p_val)

mge_boxplot=mge_counts%>% 
  pivot_longer(cols = c(Cas, RM, Retron, Dnd, Mok_Hok_Sok, Prometheus), 
          names_to = "Defense_System",
          values_to = "Defense_Count")%>% 
  mutate(Status = ifelse(Defense_Count > 0, "Present", "Absent")) %>%
  mutate(Status = factor(Status, levels = c("Absent", "Present")))
#1. Integrases 
ggplot(mge_boxplot, aes(x = Status, y = log10(Integrases+1), fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1, color = "grey30") +

  facet_wrap(~Defense_System, scales = "free_y") +
  scale_fill_manual(values = c("Absent" = "#FFFFCC", "Present" = "#99D8C9")) +

  stat_compare_means(label = "p.signif", method = "wilcox.test", label.x = 1.5) +
  labs(
       x = "Defense System Status",
       y = "log10(Integrases)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

#2. Plasmids
ggplot(mge_boxplot, aes(x = Status, y = log10(Plasmids+1), fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1, color = "grey30") +
  
  facet_wrap(~Defense_System, scales = "free_y") +
  scale_fill_manual(values = c("Absent" = "#FFFFCC", "Present" = "#99D8C9")) +
  
  stat_compare_means(label = "p.signif", method = "wilcox.test", label.x = 1.5) +
  labs(
    x = "Defense System Status",
    y = "log10(Plasmids)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

#3. Transposons
ggplot(mge_boxplot, aes(x = Status, y = log10(Transposons+1), fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1, color = "grey30") +
  
  facet_wrap(~Defense_System, scales = "free_y") +
  scale_fill_manual(values = c("Absent" = "#FFFFCC", "Present" = "#99D8C9")) +
  
  stat_compare_means(label = "p.signif", method = "wilcox.test", label.x = 1.5) +
  labs(
    x = "Defense System Status",
    y = "log10(Transposons)") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )


#Multi-faceted Plot

#This option was not selected 
mge_boxplot %>%
  pivot_longer(
    cols = c(Integrases, Plasmids, Transposons),
    names_to = "MGE_Type",
    values_to = "Count"
  ) %>%
  ggplot(aes(x = Status, y = log10(Count + 1), fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1, color = "grey30") +
  facet_grid(MGE_Type ~ Defense_System, scales = "free_y") +
  scale_fill_manual(values = c("Absent" = "#FFFFCC", "Present" = "#99D8C9")) +
  stat_compare_means(label = "p.signif", method = "wilcox.test", label.x = 1.5) +
  
  labs(
    x = "Defense System Status",
    y = "log10(MGE Abundance + 1)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", color = "grey80"),
    panel.grid.minor = element_blank()
  )
  

#This option was selected-  pre-review
mge_boxplot %>%
  pivot_longer(
    cols = c(Integrases, Plasmids, Transposons),
    names_to = "MGE_Type",
    values_to = "Count"
  ) %>%
  ggplot(aes(x = MGE_Type, y = log10(Count + 1), fill = Status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6, position = position_dodge(0.75)) +
  geom_jitter(aes(group = Status), alpha = 0.3, size = 1, 
              color = "grey30", position = position_dodge(0.75)) +
  facet_wrap(~Defense_System, scales = "free_y") +
  
  scale_fill_manual(values = c("Absent" = "#FFFFCC", "Present" = "#99D8C9")) +
  stat_compare_means(aes(group = Status), label = "p.signif", method = "wilcox.test") +
  
  labs(
    x = "Mobile Genetic Element (MGE) Type",
    y = "log10(MGE_count + 1)",
    fill = "Defense Systems"
  ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

