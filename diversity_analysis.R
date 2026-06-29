library(vegan)
library(ggplot2)
library(corrplot)
#Richness
matrix=vOTU_tablee%>%dcast(VOTUs~Genome2, value.var = "Genome2", length)%>%
  filter(VOTUs!="NA")%>%
  column_to_rownames(var="VOTUs")


richness=data.frame(Richness=apply(matrix, 2, function(x) sum(x[x>0])))%>%
  rownames_to_column(var ="Genome2")
richness2=full_join(defense_system_prime, richness, by = "Genome2")%>%
  filter(!is.na(Viral_Count))





#Stats 
all_stats_results = data.frame()
for (sys in systems){ 
  test_data1=richness2%>%
    filter(type==sys)%>%
    filter(!is.na(Richness))
  test_data1$Richness=as.numeric(test_data1$Richness)
  
  test1=shapiro.test(resid(lm((test_data1$Richness+1)~test_data1$DF_system_pr)))
  if (test1$p.value<0.05){
    test2=t.test(test_data1$Richness~test_data1$DF_system_pr)
  } else {
    test1=shapiro.test(resid(lm(log10(test_data1$Richness+1)~test_data1$DF_system_pr)))
    if (test1$p.value>0.05){
      test2=t.test(log10(test_data1$Richness+1)~test_data1$DF_system_pr)
    }  else{
      test2=t.test(base::rank(test_data1$Richness)~test_data1$DF_system_pr)
    }
  }
  
  m1=data.frame(
    Group = paste0(sys),      
    t=paste0(test1$p.value),
    s=paste0(test2$p.value)
  )
  print(m1)
  all_stats_results=rbind(all_stats_results, m1)
  
}

all_stats_results$p_adj = p.adjust(all_stats_results$s, method = "BH")
all_stats_results$Significance = ifelse(all_stats_results$p_adj < 0.05, "Significant", "NS")
print(all_stats_results)

all_stats_results$s = as.numeric(all_stats_results$p_adj)

all_stats_results$s[all_stats_results$p_adj>0.05]=""
all_stats_results$s[all_stats_results$p_adj<0.05]="*"
all_stats_results$s[all_stats_results$p_adj<0.01]="**"
all_stats_results$s[all_stats_results$p_adj<0.001]="***"
all_stats_results$s[all_stats_results$p_adj<0.0001]="****"


label_data=all_stats_results %>%
  select(type=Group, everything()) %>%
  filter(type %in% c("Cas", "Dnd", "Mok_Hok_Sok", "Prometheus", "Retron", "RM"))


#The Boxplot
richness3=full_join(richness2, label_data)
richness_plot=ggplot(richness3, aes(y = type, x =Richness , fill = DF_system_pr)) + 
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.6) + 
  geom_jitter(position = position_jitterdodge(jitter.width =0.2,dodge.width = 0.6), alpha = 0.75, size = 2, shape=21) +
  scale_fill_manual(values = c("Presence" = "#FFFFC8", "Absence" = "#99D8C9")) + stat_compare_means(method = "wilcox.test") +
  geom_text( aes(label = s, x = 7), size=7)+
  theme_bw() + 
  labs(
    x = "vOTU richness per Genome", 
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
    axis.text = element_text(size=10)
    
  )
ggsave(file = "vOTu_richness.png",plot = richness_plot, width =  10, height = 7, dpi = 600 )
dev.off()
#Correlation

vOTU_correlation=ggplot(richness2, aes(y = DF_system_Count, x =Richness)) + 
#  geom_point( size = 2, shape=21) +
  theme_bw() + 
  labs(
    x = "VoTUs per Genome", 
    y = "Defense System number per Genome"
  ) + 
  facet_wrap(~type, scales="free")+
  theme(
    strip.background = element_rect(fill = "white", colour="white"), 
    strip.text.y = element_text(angle = 0, face = "bold", size = 10), 
    strip.placement = "outside", 
    legend.position = "right", 
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size=10))+geom_smooth(method = "lm")+ 
  stat_cor(method="spearman")+ylim(-0.25, 6)+
  geom_jitter(position = position_jitterdodge(jitter.height = 0.1, jitter.width = 0.25 ), size = 2, shape=21) 
ggsave(file = "vOTu_richness_correlation.png",plot = vOTU_correlation, width =  10, height = 7, dpi = 600 )
dev.off()
