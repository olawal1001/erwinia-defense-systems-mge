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

df_counts=dcast(defense_system_prime, Genome2~type, value.var = "DF_system_Count" )
  
  
mge_counts = df_mge %>%
  group_by(Genome) %>%
  summarise(
    Plasmids = sum(Class == "Plasmids", na.rm = TRUE),
    Transposons = sum(Class == "Transposons", na.rm = TRUE),
    Integrases = sum(Class == "Integrases", na.rm = TRUE)
  )
mge_counts=full_join(df_counts,select( mge_counts, Genome2 ="Genome", everything()))%>%
  filter(Genome2%in%c(df_counts$Genome2))

stats_list = list()

for (mge in classes) {
  for (sys in systems) {
    
    x = mge_counts[[mge]]
    y = mge_counts[[sys]]
    
    group = ifelse(y > 0, "Present", "Absent")
    
    if (length(unique(group)) < 2) {
      message(paste("Skipping", mge, "-", sys, ": All genomes are either all Present or all Absent"))
      next
      }
    group_f = as.factor(group)
    testI = shapiro.test(resid(lm(x ~ group_f)))
    if (testI$p.value < 0.05) {
      # Non-normal: Use Rank-based T-test
      testII = t.test(base::rank(x) ~ group_f)
    } else {
      # Data is normal: Check log transform (using y+1 to avoid log(0))
      test_log = shapiro.test(resid(lm(log10(x + 1) ~ group_f)))
      
      if (test_log$p.value > 0.05) {
        testII = t.test(log10(x + 1) ~ group_f)
      } else {
        testII = t.test(x ~ group_f)
      }
    }
    stats_list[[paste(mge, sys, sep="_")]] <- data.frame(
      MGE_Class = mge,
      Defense_System = sys,
      P_Value = testII$p.value
    )
  }
}

stats_results$p_adj = p.adjust(stats_results$P_Value, method = "BH")
stats_results$Significance = ifelse(stats_results$p_adj < 0.05, "Significant", "NS")

print(stats_results)

stats_results$s = as.numeric(stats_results$p_adj)

stats_results$s[stats_results$p_adj>0.05]=""
stats_results$s[stats_results$p_adj<0.05]="*"
stats_results$s[stats_results$p_adj<0.01]="**"
stats_results$s[stats_results$p_adj<0.001]="***"
stats_results$s[stats_results$p_adj<0.0001]="****"


label_data=stats_results %>%
  filter(Defense_System %in% c("Cas", "Dnd", "Mok_Hok_Sok", "Prometheus", "Retron", "RM"))

mge_counts2=mge_counts%>%
  gather(key="Defense_System", value = "Defense_System_counts",
         -c(paste0(classes), "Genome2"))%>%
  gather(key="MGE_Class", value = "HGT_Marker_counts",
         -c("Defense_System", 
                   "Defense_System_counts", "Genome2"))

#The boxplot
boxplot=full_join(mge_counts2, label_data)
boxplot$Defense_System_pr=boxplot$Defense_System_counts
boxplot$Defense_System_pr[ boxplot$Defense_System_counts>0]="Presence "
boxplot$Defense_System_pr[ boxplot$Defense_System_counts==0]="Absence "

boxplots1=ggplot(boxplot, aes(y = Defense_System, x =HGT_Marker_counts , fill = Defense_System_pr)) + 
  geom_boxplot( width = 0.6, alpha = 0.6) + 
  geom_jitter(position = position_jitterdodge(jitter.width =0.2,dodge.width = 0.6), alpha = 0.75, size = 2, shape=21) +
  scale_fill_manual(values=c( "#FFFFC8", "#99D8C9")) +
  geom_text( aes(label = s, x =60), size=7)+
  theme_bw() + facet_wrap(~paste0( MGE_Class,"-associated ORFs"), scales="free")+
  labs(
    x = expression((Predicted ~ HGT-Marker~ ORFs ~ per ~ Genome)), 
    y = "Defense System Type", 
    fill = "Status"
  ) + 
  theme(
    strip.background = element_rect(fill = "white", colour="white"), 
    strip.text.y = element_text(angle = 0, face = "bold", size = 0), 
    strip.placement = "outside", 
    legend.position = "right", 
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10))
  
ggsave(file = "mge_boxplot.png",plot = boxplots1, width =  10, height = 5, dpi = 600 )
dev.off()
#correlation  
mge_correlation=ggplot(boxplot, aes(y = HGT_Marker_counts, x =Defense_System_counts)) + 
  geom_point() + geom_smooth(method = "lm")+
  theme_bw() + 
  labs(
    x = "HGT-markers", 
    y = "Defense System number per Genome"
  ) + 
  facet_wrap(~paste(MGE_Class,"vs",Defense_System), scales="free")+
  theme(
    strip.background = element_rect(fill = "white", colour="white"), 
    strip.text.y = element_text(angle = 0, face = "bold", size = 10), 
    strip.placement = "outside", 
    legend.position = "right", 
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10))+geom_smooth(method = "lm")+ stat_cor(method="spearman")

ggsave(file = "mge_correlation.png",plot = mge_correlation, width =  7, height = 5, dpi = 600 )
dev.off()






