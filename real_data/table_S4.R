setwd("~/mbDriver/real_data/UC")
library(pheatmap)
library(magrittr)
library(dplyr)
library(tidyr)
library(purrr)
library(readr)

abundance_example <- read.delim("KEGG.PathwayL2.raw.txt", row.names = 1, sep = '\t', check.names = FALSE)
abundance_example <- t(abundance_example)
map <- read.delim("Group.txt", sep = '\t', check.names = FALSE)
G_abundance <- as.data.frame(abundance_example)
G_abundance$ID <- rownames(abundance_example)
G_metadata <- map
merged_df <- merge(G_metadata, G_abundance, by = "ID")
merged_df<- merged_df[,-1]

KO_pathway <- merged_df %>%
  group_by(Group) %>%
  summarize(across(everything(), mean))
KO_pathway <- KO_pathway[,-1]
rownames(KO_pathway) <- c("H1","H2","UC1","UC2")

KO_pathway_write <- t(KO_pathway)
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}
KO_pathway_write <- adjustdata(KO_pathway_write)
#write.table(KO_pathway_write,file ="Supplementary_table3.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###KW test
dt_ <- read_delim("KEGG.PathwayL2.raw.txt",delim="\t")
mdata <- read_delim("Group.txt", delim="\t")

dt_long <- dt_ %>% 
  pivot_longer(-dplyr::contains("PathwayL2"), names_to = "Sample", values_to = "abundance") %>%
  inner_join(mdata, by=c("Sample"="ID")) %>% mutate(Group=factor(.$Group))

kw_model <- dt_long %>%
  group_by(PathwayL2) %>%
  nest() %>%
  mutate(mod=map(data, function(dt_) kruskal.test(abundance ~ Group, dt_))) %>%
  mutate(pvalue = map_dbl(mod, function(md) md$p.value))

sig_pathway <- kw_model %>% filter(pvalue <=0.05)
nosig_pathway <- kw_model %>% filter(pvalue >=0.05)
sig_pathway_p <- sig_pathway[,c(1,4)]

#tableS3
tableS3 <- sig_pathway_p[sig_pathway_p$PathwayL2 %in% c("Lipid metabolism", "Amino acid metabolism"), ]
write.table(tableS3,file ="Supplementary Table 3.txt",row.names = F,col.names = T, sep = "\t",quote = F)

