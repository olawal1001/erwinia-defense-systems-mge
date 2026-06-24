library(readr)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(stringi)
library(dplyr)
library(tidyr)
library(stringr)
library(gridExtra)
library(corrplot)

setwd("C:/Users/Lawal/Desktop/Data/CheckV")

list1= list.dirs(recursive = F)
list1=list.files(path = ".", pattern = "quality_summary.tsv", recursive = T)

checkv_tables=list()
for (i in 1:length(list1)){
  file1= read_table(paste0(list1[i]))
  file1$Genome=paste0(list1[i])
  checkv_tables=rbind(checkv_tables, file1)  
  print(paste0(i))
}
checkv_tables$Genome2 = gsub("_genomic_virus*.*", "",
                             gsub("CheckV.","",checkv_tables$Genome))

checkv_tables2=filter(checkv_tables, completeness>=95)%>%dcast(Genome2~., value.var = "Genome2", length)%>%
  select(Viral_Count=".", Genome2)

setwd("C:/Users/Lawal/Desktop/Data/geNomad")

#GENOMAD 
list4= list.files(path = ".", pattern = "_genomic_virus_summary.tsv", recursive = T)

Genomad_tables=list()
for (i in 1:length(list4)){
  file1= read_table(paste0(list4[i]))
  file1$Genome=paste0(list4[i])
  Genomad_tables=rbind(Genomad_tables, file1)  
  print(paste0(i))
}

df= as.data.frame(Genomad_tables)
Genomad_tables$Genome2 = gsub("_genomic_summary*.*","" , gsub("geNomad.", "",Genomad_tables$Genome))
Genomad_tables=Genomad_tables%>%
  separate(Kingdom,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";", fill = "right", extra = "merge")


prediction_table = full_join(x = Genomad_tables,
                             y = checkv_tables,
                             by=c("Genome2")) 

setwd("C:/Users/Lawal/Desktop/Data/CheckM2")

#CHECKM2
list2=list.files(path = ".", pattern = "quality_report.tsv", recursive = T)
checkm2_tables=list()
for (i in 1:length(list2)){
  file1=read_table(paste0(list2[i]))
  file1$Genome=paste0(list2[i])
  checkm2_tables=rbind(checkm2_tables, file1)
  print(paste0(i))
}
df=as.data.frame(checkm2_tables)


df$Genome= gsub("1_.*.*", "1", gsub("checkm.", "",df$Genome))
checkm2_tables$Genome2=  gsub("*_genomic/quality_report.*","" , gsub("CheckM2.", "",checkm2_tables$Genome))


write.csv(df, file = 'checkmresults.csv')

setwd("C:/Users/Lawal/Desktop/Data/DefenseFinder")

#DEFENSE_FINDER 
list3= list.files(path = ".", pattern = "genomic_defense_finder_systems.tsv", recursive = T)
defense_finder_tables=list()
for (i in 1:length(list3)){
  file1= read_table(paste0(list3[i]))
  file1$Genome=paste0(list3[i])
  defense_finder_tables=rbind(defense_finder_tables, file1)  
  print(paste0(i))
}

df=as.data.frame(defense_finder_tables)
df$Genome= gsub("1_.*.*","1" , gsub("defense_finder.", "",df$Genome))

defense_finder_tables$Genome2= gsub("*_genomic_defense_finder_systems.*","" , gsub("_genomic/.*$", "", gsub("DefenseFinder.", "",defense_finder_tables$Genome)))

write.csv(df, file = 'defense_finder.csv')


defense_finder_table2=dcast(defense_finder_tables, type+Genome2~., value.var = "Genome2", length)%>%
  select(DF_system_Count=".", type, Genome2)

prediction_table = full_join(x = Genomad_tables,
                             y = checkv_tables,
                             by=c("Genome2")) 

prediction_table$virus_presence=prediction_table$completeness
for (i in 1:nrow(prediction_table)){ 
  if (is.na(prediction_table[i, "completeness"])){
    prediction_table[i, "virus_presence"]=0} 
  else if 
  (prediction_table[i, "completeness"]<100){
    prediction_table[i, "virus_presence"]=0
  }
  else{
    prediction_table[i, "virus_presence"]=1
  } 
}

prediction_table$virus_presence=prediction_table$completeness
prediction_table$virus_presence[is.na(prediction_table$completeness)]=0
prediction_table$virus_presence[(prediction_table$completeness)<100]=0
prediction_table$virus_presence[(prediction_table$completeness)!=0]=1

virus_presence_absence=dcast(prediction_table,
                             Genome2~., value.var ="virus_presence", sum)%>%
  select(phage_pr=".", Genome2)
Quality_table=filter(checkm2_tables, Completeness>=95&Contamination<=5) 


plot_data = full_join(virus_presence_absence, Quality_table , by = "Genome2")%>%
  filter( Completeness>=90&Contamination<=5) %>%filter(!str_detect(Genome2,"MAG"))

checkv_tables2=filter(checkv_tables, completeness>=95)%>%dcast(Genome~., value.var = "Genome", length)%>%
  select(Viral_Count=".", Genome2=Genome)

defense_finder_table2=dcast(defense_finder_tables, type+Genome2~., value.var = "Genome2", length)%>%
  select(DF_system_Count=".", type, Genome2)



merged_table= full_join(checkv_tables2, defense_finder_table2)%>%
  full_join(select( checkm2_tables, Genome2, everything()))%>%
  filter(Completeness>=95&Contamination<=5)



merged_table$Viral_Count[is.na(merged_table$Viral_Count)]=0
merged_table=as.data.frame(merged_table)
merged_table_2=dcast(data=merged_table, Genome2+Viral_Count~type, value.var = "DF_system_Count", sum, fill = 0)%>%
  gather(-c("Genome2", "Viral_Count"), key="type", value="DF_system_Count")
merged_table_2$DF_system_pr=merged_table_2$DF_system_Count
merged_table_2$DF_system_pr[merged_table_2$DF_system_Count>0]="Presence"
merged_table_2$DF_system_pr[merged_table_2$DF_system_Count==0]="Absence"

defense_system_prime=merged_table_2%>%
  filter(type %in% c("Cas","Dnd", "Mok_Hok_Sok", "Prometheus", "Retron","RM"))


all_stats_results = data.frame()
for (sys in systems){ 
  test_data=defense_system_prime%>%
    filter(type==sys)%>%
    filter(!is.na(Viral_Count))
  test1=shapiro.test(resid(lm((test_data$Viral_Count+1)~test_data$DF_system_pr)))
  if (test1$p.value<0.05){
    test2=t.test((test_data$Viral_Count+1)~test_data$DF_system_pr)
  } else {
    test1=shapiro.test(resid(lm(log10(test_data$Viral_Count+1)~test_data$DF_system_pr)))
    if (test1$p.value>0.05){
      test2=t.test(log10(test_data$Viral_Count+1)~test_data$DF_system_pr)
    }  else{
      test2=t.test(base::rank(test_data$Viral_Count)~test_data$DF_system_pr)
    }
  }
  m1=data.frame(
    Group = paste0(sys),
    t=paste0(test1$p.value),
    s=paste0(test2$p.value)
  )
  print(m1)
  all_stats_results=rbind(all_stats_results, m1)
  
  plot=ggplot(defense_system_prime, aes(x=DF_system_pr, y=Viral_Count + 1, fill=DF_system_pr))+geom_boxplot(outlier.shape = NA, width = 0.6) + geom_jitter(position = position_jitter(height = 0.1), alpha = 0.3, size = 0.8) + labs(title = paste("Defense system:", sys),subtitle = paste("p-value:", round(test2$p.value, 6)))
}
print(plot)
all_stats_results$s = as.numeric(all_stats_results$s)
all_stats_results$p_adj = p.adjust(all_stats_results$s, method = "BH")
all_stats_results$Significance = ifelse(all_stats_results$p_adj < 0.05, "Significant", "NS")
print(all_stats_results)

test_data2 = defense_system_prime %>%
  filter(!is.na(Viral_Count)) %>%
  filter(type %in% systems)


stats_for_plot=all_stats_results %>%
  mutate(s = as.numeric(s)) %>%
  mutate(p_stars = case_when(
    s = 0.0001 ~ "****",
    s = 0.001 ~ "***",
    s = 0.01  ~ "**",
    s = 0.05  ~ "*",
    TRUE      ~ "ns"
  ))

label_data=stats_for_plot %>%
  rename(Type=Group) %>%
  filter(Type %in% c("Cas", "Dnd", "Mok_Hok_Sok", "Prometheus", "Retron", "RM"))
ggplot(test_data2, aes(x = Viral_Count, y =DF_system_pr, fill = DF_system_pr)) + 
  geom_boxplot(outlier.shape = NA, width = 0.6) + 
  geom_jitter(position = position_jitter(height = 0.1), alpha = 0.3, size = 0.8) +
  facet_wrap(~Type, ncol = 1, strip.position = "left") + 
  scale_fill_manual(values = c("Presence" = "maroon", "Absence" = "#99D8C9")) +
  geom_text(data = label_data, aes(label = p_stars), x = Inf, y = Inf,  hjust = 1.5, vjust = 1.5, size = 4, fontface = "bold", inherit.aes = FALSE) +
  theme_bw() + theme(strip.background = element_rect(fill = "white"), strip.text.y.right  = element_text(angle = 0, face = "bold", size = 10), strip.placement = "outside", legend.position = "right", axis.title = element_text(face = "bold")
  ) + labs( x = expression((Predicted~ Phage ~ Number ~ per ~ Genome)), y = "Defense System Type", fill = "Status")


freq_data=defense_system_prime %>%
  filter(type %in% systems) %>%
  group_by(type) %>%
  summarise(Total_Count = sum(DF_system_Count > 0))
ggplot(freq_data, aes(x = type, y = Total_Count, fill = type)) + geom_bar(stat = "identity", color = "black", width = 0.7) + geom_text(aes(label = Total_Count), vjust = -0.5, fontface = "bold") +
  theme_classic() + scale_fill_brewer(palette = "Set3") + labs(x = "Defense System Type", y = "Number of Genomes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), legend.position = "none")