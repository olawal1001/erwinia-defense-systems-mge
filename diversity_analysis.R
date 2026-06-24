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

richness_plot=ggplot(richness2, aes(x=type, y=Richness, fill=DF_system_pr))+ geom_boxplot() + geom_jitter(position=position_jitterdodge(dodge.width = 1))+stat_compare_means(data = NULL)
richness_plot2=ggplot(richness2, aes(x=type, y=Richness, fill=DF_system_pr))+geom_point(stat="summary") + geom_smooth(method="lm") + facet_wrap(~type, scales="free") + stat_cor(method ="spearman")
richness_plot3=ggplot(richness2, aes(x=type, y=Richness, fill=DF_system_pr))+ geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(position=position_jitter(height = 0.1), alpha = 0.3, size = 0.8)+stat_compare_means(data = NULL)+ 
  scale_fill_manual(values = c("Presence" = "#FFFFCC", "Absence" = "#99D8C9"))
  
print(richness_plot2)

