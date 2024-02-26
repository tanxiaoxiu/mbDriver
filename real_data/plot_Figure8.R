library(tidyverse)
library(readxl)
library(ggpubr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(tidyr)
library(readr)
library(dplyr)
library(ggrepel)
library(forcats)
library(ggpubr)
library(pheatmap)

setwd("~/mbDriver/real_data/UC")

###Figure8B
MAX_N_OTU =10
l_tb <- read.delim("UC_driver_summary.txt",  sep = '\t', check.names = FALSE)
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "Score")

w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "Score",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")

p1 <- pheatmap(w_tb, cluster_cols = F,angle_col=0,fontsize=14)
ggsave(filename="Figure8B.png",plot=p1,device="png",dpi=600,units="in",width=6,height=7)

###Figure8C
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

pathway <-  t(read.delim("KEGG.PathwayL2.raw.txt",  sep = '\t', header = FALSE, row.names = 1, check.names = FALSE))
colnames(pathway)[1] <- "SampleID"
data_UC <- read.delim("data_UC.txt",  sep = '\t', check.names = FALSE)
Group <- data_UC[,1:4]

colnames(Group) <- c("SampleID","subject","time","group")

pathway_group <- merge(Group, pathway, by = "SampleID", all = T)

species_top10 <- data_UC[,c(1,5:14)]
colnames(species_top10)[1] <- "SampleID"
pathway_group <- na.omit(pathway_group)
pathway_species <- merge(pathway_group, species_top10, by = "SampleID", all = F)

#Lipid metabolism /Energy metabolism /Amino acid metabolism
#########################################################################H1
pathway_species_H1 <- pathway_species[pathway_species$group == "H1", ]
Driver_H1 <- pathway_species_H1[,c("time","Lipid metabolism","Amino acid metabolism",
                                   "Porphyromonadaceae","Akkermansiaceae","Desulfovibrionaceae",
                                   "Bacteroidaceae","Prevotellaceae",
                                   "group","subject")]

colnames(Driver_H1) <- c("time","Lipid_metabolism","Amino_acid_metabolism",
                         "Porphyromonadaceae","Akkermansiaceae","Desulfovibrionaceae",
                         "Bacteroidaceae","Prevotellaceae",
                         "group","subject")

Driver_H1_scale <- Driver_H1
Driver_H1_scale[,2] <- as.numeric(Driver_H1_scale[,2])
Driver_H1_scale[,3] <- as.numeric(Driver_H1_scale[,3])
Driver_H1_scale[,2:8] <- log10(Driver_H1_scale[,2:8]+1)
S <- Driver_H1_scale$subject

###Lipid_metabolism
model.f1 <- lmer(Lipid_metabolism ~ Porphyromonadaceae + Akkermansiaceae + Bacteroidaceae +
                   Prevotellaceae + Desulfovibrionaceae  + (1|S), data = Driver_H1_scale, REML=FALSE)
summary <- summary(model.f1)
Fixed_effects_H1_Lipid_metabolism <- adjustdata(summary$coefficients)
Fixed_effects_H1_Lipid_metabolism <- cbind(class = "H1:Lipid metabolism",Fixed_effects_H1_Lipid_metabolism)
colnames(Fixed_effects_H1_Lipid_metabolism)[2] <- "var"

###Amino_acid_metabolism
model.f3 <- lmer(Amino_acid_metabolism ~ Porphyromonadaceae + Akkermansiaceae + Bacteroidaceae +
                   Prevotellaceae + Desulfovibrionaceae + (1|S), data = Driver_H1_scale, REML=FALSE)
summary <- summary(model.f3)
Fixed_effects_H1_Amino_acid <- adjustdata(summary$coefficients)
Fixed_effects_H1_Amino_acid <- cbind(class = "H1:Amino acid metabolism",Fixed_effects_H1_Amino_acid)
colnames(Fixed_effects_H1_Amino_acid)[2] <- "var"

#########################################################################H2
pathway_species_H2 <- pathway_species[pathway_species$group == "H2", ]
Driver_H2 <- pathway_species_H2[,c("time","Lipid metabolism","Amino acid metabolism",
                                   "Ruminococcaceae", "Akkermansiaceae","Acidaminococcaceae",
                                   "Porphyromonadaceae","Lachnospiraceae",
                                   "group","subject")]

colnames(Driver_H2) <- c("time","Lipid_metabolism","Amino_acid_metabolism",
                         "Ruminococcaceae","Akkermansiaceae","Acidaminococcaceae",
                         "Porphyromonadaceae","Lachnospiraceae",
                         "group","subject")

Driver_H2_scale <- Driver_H2
Driver_H2_scale[,2] <- as.numeric(Driver_H2_scale[,2])
Driver_H2_scale[,3] <- as.numeric(Driver_H2_scale[,3])
Driver_H2_scale[,2:8] <- log10(Driver_H2_scale[,2:8]+1)
S <- Driver_H2_scale$subject

###Lipid_metabolism
model.f1 <- lmer(Lipid_metabolism ~ Ruminococcaceae + Akkermansiaceae + Acidaminococcaceae +
                   Porphyromonadaceae + Lachnospiraceae  + (1|S), data = Driver_H2_scale, REML=FALSE)
summary <- summary(model.f1)
Fixed_effects_H2_Lipid_metabolism <- adjustdata(summary$coefficients)
Fixed_effects_H2_Lipid_metabolism <- cbind(class = "H2:Lipid metabolism",Fixed_effects_H2_Lipid_metabolism)
colnames(Fixed_effects_H2_Lipid_metabolism)[2] <- "var"

###Amino_acid_metabolism
model.f3 <- lmer(Amino_acid_metabolism ~ Ruminococcaceae + Akkermansiaceae + Acidaminococcaceae +
                   Porphyromonadaceae + Lachnospiraceae  + (1|S), data = Driver_H2_scale, REML=FALSE)
summary <- summary(model.f3)
Fixed_effects_H2_Amino_acid <- adjustdata(summary$coefficients)
Fixed_effects_H2_Amino_acid <- cbind(class = "H2:Amino acid metabolism",Fixed_effects_H2_Amino_acid)
colnames(Fixed_effects_H2_Amino_acid)[2] <- "var"


##############################################################################UC1
pathway_species_UC1 <- pathway_species[pathway_species$group == "UC1", ]

Driver_UC1 <- pathway_species_UC1[,c("time","Lipid metabolism","Amino acid metabolism",
                                     "Lachnospiraceae","Porphyromonadaceae",
                                     "Desulfovibrionaceae","Ruminococcaceae","Bacteroidaceae",
                                     "group","subject")]


colnames(Driver_UC1) <- c("time","Lipid_metabolism","Amino_acid_metabolism",
                          "Lachnospiraceae","Porphyromonadaceae",
                          "Desulfovibrionaceae","Ruminococcaceae","Bacteroidaceae",
                          "group","subject")

Driver_UC1_scale <- Driver_UC1
Driver_UC1_scale[,2] <- as.numeric(Driver_UC1_scale[,2])
Driver_UC1_scale[,3] <- as.numeric(Driver_UC1_scale[,3])
Driver_UC1_scale[,2:8] <- log10(Driver_UC1_scale[,2:8]+1)
S <- Driver_UC1_scale$subject

###Lipid_metabolism
model.f1 <- lmer(Lipid_metabolism ~ Lachnospiraceae + Porphyromonadaceae + Desulfovibrionaceae +
                   Ruminococcaceae + Bacteroidaceae  + (1|S), data = Driver_UC1_scale, REML=FALSE)
summary <- summary(model.f1)
Fixed_effects_UC1_Lipid_metabolism <- adjustdata(summary$coefficients)
Fixed_effects_UC1_Lipid_metabolism <- cbind(class = "UC1:Lipid metabolism",Fixed_effects_UC1_Lipid_metabolism)
colnames(Fixed_effects_UC1_Lipid_metabolism)[2] <- "var"

###Amino_acid_metabolism
model.f3 <- lmer(Amino_acid_metabolism ~ Lachnospiraceae + Porphyromonadaceae + Desulfovibrionaceae +
                   Ruminococcaceae + Bacteroidaceae + (1|S), data = Driver_UC1_scale, REML=FALSE)
summary <- summary(model.f3)
Fixed_effects_UC1_Amino_acid <- adjustdata(summary$coefficients)
Fixed_effects_UC1_Amino_acid <- cbind(class = "UC1:Amino acid metabolism",Fixed_effects_UC1_Amino_acid)
colnames(Fixed_effects_UC1_Amino_acid)[2] <- "var"

##############################################################################UC2
pathway_species_UC2 <- pathway_species[pathway_species$group == "UC2", ]

Driver_UC2 <- pathway_species_UC2[,c("time","Lipid metabolism","Amino acid metabolism",
                                     "Enterobacteriaceae","Acidaminococcaceae",
                                     "Ruminococcaceae","Prevotellaceae","Desulfovibrionaceae",
                                     "group","subject")]

colnames(Driver_UC2) <- c("time","Lipid_metabolism","Amino_acid_metabolism",
                          "Enterobacteriaceae","Acidaminococcaceae",
                          "Ruminococcaceae","Prevotellaceae","Desulfovibrionaceae",
                          "group","subject")

Driver_UC2_scale <- Driver_UC2
Driver_UC2_scale[,2] <- as.numeric(Driver_UC2_scale[,2])
Driver_UC2_scale[,3] <- as.numeric(Driver_UC2_scale[,3])
Driver_UC2_scale[,2:8] <- log10(Driver_UC2_scale[,2:8]+1)
S <- Driver_UC2_scale$subject

###Lipid_metabolism
model.f1 <- lmer(Lipid_metabolism ~ Enterobacteriaceae + Acidaminococcaceae + Ruminococcaceae +
                   Prevotellaceae + Desulfovibrionaceae  + (1|S), data = Driver_UC2_scale, REML=FALSE)
summary <- summary(model.f1)
Fixed_effects_UC2_Lipid_metabolism <- adjustdata(summary$coefficients)
Fixed_effects_UC2_Lipid_metabolism <- cbind(class = "UC2:Lipid metabolism",Fixed_effects_UC2_Lipid_metabolism)
colnames(Fixed_effects_UC2_Lipid_metabolism)[2] <- "var"

###Amino_acid_metabolism
model.f3 <- lmer(Amino_acid_metabolism ~ Enterobacteriaceae + Acidaminococcaceae + Ruminococcaceae +
                   Prevotellaceae + Desulfovibrionaceae + (1|S), data = Driver_UC2_scale, REML=FALSE)
summary <- summary(model.f3)
Fixed_effects_UC2_Amino_acid <- adjustdata(summary$coefficients)
Fixed_effects_UC2_Amino_acid <- cbind(class = "UC2:Amino acid metabolism",Fixed_effects_UC2_Amino_acid)
colnames(Fixed_effects_UC2_Amino_acid)[2] <- "var"

UC_scatter_data <- rbind(Fixed_effects_H1_Lipid_metabolism,Fixed_effects_H1_Amino_acid,
                         Fixed_effects_H2_Lipid_metabolism,Fixed_effects_H2_Amino_acid,
                         Fixed_effects_UC1_Lipid_metabolism,Fixed_effects_UC1_Amino_acid,
                         Fixed_effects_UC2_Lipid_metabolism,Fixed_effects_UC2_Amino_acid)

colnames(UC_scatter_data) <- c("class","var","Estimate","Std.Error","df","t value","P_value")
write.table(UC_scatter_data,file ="UC_scatter_data.txt",row.names = F,col.names = T, sep = "\t",quote = F)

##plot
##Figure8C
infile <- "UC_scatter_data.txt"
scatter_data <- read_delim(infile)
dt4plot <- scatter_data %>% 
  filter(!var == "(Intercept)") %>%
  mutate(sig=ifelse(.$P_value <= 0.05, "yes", "no"))
dt4plot <- dt4plot[, c("class", "var", "Estimate", "P_value", "sig")]
dt4label <- dt4plot %>% filter(sig == "yes")

p2 <- ggplot(dt4plot, aes(x=Estimate, y=-log10(P_value), color=sig, group=class)) +
  geom_point() +
  geom_text_repel(
    aes(x=Estimate, y=-log10(P_value), label=var), data=dt4label, inherit.aes = FALSE) +
  facet_wrap(. ~ class, ncol=2) +
  theme_bw() +
  theme(legend.position="none", panel.grid = element_blank())
p2 <- p2 + labs(x = "Coefficient")
p2 <- p2 + labs(y = "-log10(P-value)")
p2 <- p2 + scale_color_manual(
  values = c("yes" = "#fc6c64", "no" = "#08bdc5"), 
  breaks = c("yes", "no"),                   
)

p2 <- p2 + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")
p2 <- p2 + geom_vline(xintercept = 0, linetype = "dashed", color = "grey")
ggsave(filename = "Figure9C.png", plot=p2, width=8, height = 8, dpi=600)







