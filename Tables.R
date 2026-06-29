#Import libraries 
library(readr)
library(ggpubr)
library(ggplot2)
library(reshape2)
library(stringi)
library(dplyr)
library(tidyr)
library(stringr)

#set working directory to folder that contained the data 
setwd("C:/Users/Lawal/Desktop/Data/CheckV")

#make a list from the file in the  working  directory
#The recursive function, if TRUE list all directories including all subdirectories; if false, returns only the main directory in the specified working directory/path. 
#list1= list.dirs(recursive = F)
list1=list.files(path = ".", pattern = "quality_summary.tsv", recursive = T)

#CHECKV

#creating a table from list
#paste0() is a function used to concatenate/combine strings without any separators (" ", ;) between them 

checkv_tables=list()
for (i in 1:length(list1)){
  file1= read_table(paste0(list1[i]))
  file1$Genome=paste0(list1[i])
  checkv_tables=rbind(checkv_tables, file1)  
  print(paste0(i))
}
#Generating a data frame from a table 

df=as.data.frame(checkv_tables)
#data cleaning, removing separators 
df$Genome= gsub("_genomic_virus*.*","", df$Genome)

checkv_tables$Genome = gsub("_genomic_virus*.*", "",
                            gsub("CheckV.","",checkv_tables$Genome))

df$contig_fasta=gsub("\\|.*","", df$contig_id )
df$Genomefasta=gsub("\\.1.*",".1", gsub( "\\.2.*",".2",df$Genome))

#save csv file locally, writing data to local computer 
write.csv(df, file = "CheckVresultsII.csv")


count=sum(checkv_tables$completeness ==100, na.rm =TRUE)

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


write.csv(df, file = 'checkm2.csv')

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

Genomad_tables$Kingdom = Genomad_tables$taxonomy

names(Genomad_tables)

Genomad_tables <- Genomad_tables%>%
  separate(Kingdom,
           into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";", fill = "right", extra = "merge")
names(Genomad_tables)

write.csv(df, file = 'Genomadresults.csv')

checkv_tables= select(checkv_tables, completeness, seq_name=contig_id, Genome2)

prediction_table = full_join(x = Genomad_tables,
                             y = checkv_tables,
                             by=c("Genome2", "seq_name")) 

checkv_tables2=filter(checkv_tables, completeness>=95)%>%dcast(Genome~., value.var = "Genome", length)%>%
  select(Viral_Count=".", Genome2=Genome)

defense_finder_table2=dcast(defense_finder_tables, type+Genome2~., value.var = "Genome2", length)%>%
  select(DF_system_Count=".", type, Genome2)
defense_finder_table2$Genome =defense_finder_table2$Genome2
defense_finder_table2= full_join(checkv_tables2, defense_finder_table2)%>%
full_join(select( checkm2_tables, Genome2, everything()))%>%filter(
  Completeness>=95&Contamination<=5)



merged_table$Viral_Count[is.na(merged_table$Viral_Count)]=0
merged_table=as.data.frame(merged_table)
merged_table_2=dcast(data=merged_table, Genome2+Viral_Count~type, value.var = "DF_system_Count", sum, fill = 0)%>%
  gather(-c("Genome2", "Viral_Count"), key="type", value="DF_system_Count" )
merged_table_3=full_join(merged_table_2, mge_tables, by=Genome)


merged_table_2$DF_system_pr=merged_table_2$DF_system_Count
merged_table_2$DF_system_pr[merged_table_2$DF_system_Count>0]="Presence"
merged_table_2$DF_system_pr[merged_table_2$DF_system_Count==0]="Absence"

defense_system_prime=merged_table_2%>%
  filter(type %in% c("Cas","Dnd", "Mok_Hok_Sok", "Prometheus", "Retron","RM"))

#filtering the quality table by completeness >=95% and contamination <=5%
Quality_table=filter(checkm2_tables, Completeness>=95&Contamination<=5) 

#making sure of no NA values
Quality_table$Frequency[is.na(Quality_table$Frequency)]=0
View(Quality_table)


#DIAMOND

#creating a table from list
#paste0() is a function used to concatenate/combine strings without any separators (" ", ;) between them 
setwd("C:/Users/Lawal/Desktop/Diamond/Diamond")

listP=list.files(path = ".", pattern = "genomic_prots.tsv", recursive = F)

col_names=c("qseqid","sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

diamond_tables=list()
for (i in 1:length(listP)){
  file1= read_tsv(listP[i], col_names = col_names)
  file1$Genome=paste0(listP[i])
  diamond_tables=rbind(diamond_tables, file1)  
  print(paste0(i))
}
#Generating a data frame from a table 

df=as.data.frame(diamond_tables)
#data cleaning, removing separators 
diamond_tables$Genome = gsub("_genomic_prots.tsv","",diamond_tables$Genome)
write.csv(df, file = 'diamond.csv')
#MOBILE GENETIC ELEMENTS (MGE)
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
#data cleaning, removing separators 
mge_tables$Genome = gsub("_genomic_prots.faa.tsv","",mge_tables$Genome)
View(mge_tables)

mge_tables <- mge_tables%>%
  separate(`sseqid`,
           into = c("mobileOG ID","Gene","Accession Number","Category","Initiation","Copy"),
           sep = "\\|", fill = "right", extra = "merge")
write.csv(df, file = 'mge_tables.csv')
#PATHOGEN-HOST INTERACTIONS
setwd("C:/Users/Lawal/Desktop/Data/PHI")
listP=list.files(path = ".", pattern = "genomic_prots.faa.tsv", recursive = F)
col_names=c("qseqid","sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")

phi_tables=list()
for (i in 1:length(listP)){
  file1= read_tsv(listP[i],  col_names = col_names)
  file1$Genome=paste0(listP[i])
  phi_tables=rbind(phi_tables, file1)  
  print(paste0(i))
}

df=as.data.frame(phi_tables)
#data cleaning, removing separators 
phi_tables$Genome = gsub("_genomic_prots.faa.tsv","",phi_tables$Genome)
View(phi_tables)

phi_tables <- phi_tables%>%
  separate(`sseqid`,
           into = c("UniprotID","Accession","Gene","Hostcode","Specie","Phenotype"),
           sep = "#", fill = "right", extra = "merge")
write.csv(df, file = 'phi_tables.csv')

#PHAGE-DEFENSE-SYSTEM TABLE 
setwd("C:/Users/Lawal/Desktop/Data/checkv_defensefinder")
phage_defense_finder_table=list()
listC= list.files(path = ".", pattern = "viruses_defense_finder_systems.tsv", recursive = T)
defense_finder_tables=list()
for (i in 1:length(listC)){
  file1= read_table(paste0(listC[i]))
  file1$Genome=paste0(listC[i])
  phage_defense_finder_table=rbind(phage_defense_finder_table, file1)  
  print(paste0(i))
}

df=as.data.frame(phage_defense_finder_table)
#df$Genome= gsub("_genomic_viruses/*.*" , gsub("defense_finder.", "",df$Genome))

phage_defense_finder_table$Genome2= gsub("*viruses_defense_finder_systems.*","" , gsub("_genomic_virus/.*$", "", gsub("DefenseFinder.", "",phage_defense_finder_table$Genome)))
#GCA_001078075.1_ASM107807v1_genomic_virus/viruses_defense_finder_systems.tsv


write.csv(df, file = 'phage_defense_finder.csv')


#ANTIBIOTIC RESISTANCE GENES
setwd("C:/Users/Lawal/Desktop/FINAL_CONCATENATED_RESULTS/ARGs")

list10=list.files(path = ".", pattern = "genomic_prots.faa.tsv", recursive = F)
col_names=c("qseqid","sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
arg_tables=list()
for (i in 1:length(list10)){
  file1= read_tsv(list10[i],  col_names = col_names)
  file1$Genome=paste0(list10[i])
  arg_tables=rbind(arg_tables, file1)  
  print(paste0(i))
}

df=as.data.frame(arg_tables)
#data cleaning, removing separators 
arg_tables$Genome = gsub("_genomic_prots.faa.tsv","",arg_tables$Genome)
arg_tables=arg_tables %>%
  mutate(sseqid = str_trim(sseqid)) %>% 
  separate_wider_delim(
    cols = sseqid,
    delim = "|",
    names = c("db_prefix", "accession", "aro_id", "gene_symbol"),
    too_many = "merge",
    too_few = "align_start"
  )
filepath= "C:/Users/Lawal/Downloads/card-ontology/aro.tsv"
aro_metadata=read_tsv(filepath)


#Predicted plasmid sequences from geNomad Validation 
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

write.csv(plasmid_tables, file = 'plasmid_tables.csv')

#plasmid anti-phage systems 
setwd("C:/Users/Lawal/Desktop/Data/plasmid-antiphage-systems")

list13=list.files(path = ".", pattern = "genomic_plasmid_proteins_defense_finder_systems.tsv", recursive = T)
#col_names=c("qseqid","sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
plasmid_defense_finder_tables=list()
for (i in 1:length(list13)){
  file1= read_table(paste0(list13[i]))
  file1$Genome=paste0(list13[i])
  plasmid_defense_finder_tables=rbind(plasmid_defense_finder_tables, file1)  
  print(paste0(i))
}

df=as.data.frame(plasmid_defense_finder_tables)
plasmid_defense_finder_tables$Genome = gsub("plasmid_proteins_defense_finder_systems.tsv","",plasmid_defense_finder_tables$Genome)
plasmid_defense_finder_tables$Genome = gsub("_genomic_","",plasmid_defense_finder_tables$Genome)
plasmid_defense_finder_tables$Genome= gsub("1_.*.*","1" , gsub("defense_finder.", "",df$Genome))

plasmid_defense_finder_tables$Genome= gsub("*_genomic_plasmid_proteins_defense_finder_systems.*","",plasmid_defense_finder_tables$Genome)
plasmid_defense_finder_tables$Genome2= gsub("*gemomic_plasmid_proteins_defense_finder_systems.*", "", gsub("DefenseFinder.", "",plasmid_defense_finder_tables$Genome))


write.csv(plasmid_defense_finder_tables, file = 'plasmid_defense_finder_tables.csv')

#Predicted plasmid sequences pathogenicity inference 
setwd("C:/Users/Lawal/Desktop/phi-plasmidII")
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
write.csv(phi_plasmid_tables, file = 'phi-plasmid_tables.csv')




