#viral_genomes for VOTU prediction 

setwd("C:/Users/Lawal/Desktop/Data/predicted_virus_genomes/vOTU")
viral_genome=checkv_tables%>% 
  filter(completeness >=95.0)
viral_genome$Genome2=viral_genome$Genome_name
fasta=readDNAStringSet("viromes.combined.fa")

names1=gsub("\\.1 .*", ".1", names(fasta))
names1=gsub("\\.2 .*", ".2", names1)
names1=gsub(".GC.*", "", names1)


names2=gsub("*.*.GC", "GC", names(fasta))

names2=gsub("\\.1.*",".1", gsub( "\\.2.*",".2",names2))
names(fasta)=paste0( names1,names2)

checkv_ids=paste0(viral_genome$contig_id, viral_genome$Genomefasta)

quality_fasta = fasta[names(fasta) %in% checkv_ids]
writeXStringSet(quality_fasta, file ="quality_viromes.fa", format ="fasta")
