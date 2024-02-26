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

adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}


setwd("~/mbDriver/real_data/UC")
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
##############################################################################UC1
pathway_species_UC1 <- pathway_species[pathway_species$group == "UC1", ]

Driver_UC1 <- pathway_species_UC1[,c("time","Lipid metabolism","Amino acid metabolism",
                                   "Lachnospiraceae","Porphyromonadaceae",
                                   "Desulfovibrionaceae","Prevotellaceae","Ruminococcaceae",
                                   "group","subject")]


colnames(Driver_UC1) <- c("time","Lipid_metabolism","Amino_acid_metabolism",
                          "Lachnospiraceae","Porphyromonadaceae",
                          "Desulfovibrionaceae","Prevotellaceae","Ruminococcaceae",
                          "group","subject")

Driver_UC1_scale <- Driver_UC1
Driver_UC1_scale[,2] <- as.numeric(Driver_UC1_scale[,2])
Driver_UC1_scale[,3] <- as.numeric(Driver_UC1_scale[,3])
Driver_UC1_scale[,2:8] <- log10(Driver_UC1_scale[,2:8]+1)
S <- Driver_UC1_scale$subject

###Lipid_metabolism
model.f1 <- lmer(Lipid_metabolism ~ Lachnospiraceae + Porphyromonadaceae + Desulfovibrionaceae +
                   Prevotellaceae + Ruminococcaceae  + (1|S), data = Driver_UC1_scale, REML=FALSE)
summary <- summary(model.f1)
Fixed_effects_UC1_Lipid_metabolism <- adjustdata(summary$coefficients)
Fixed_effects_UC1_Lipid_metabolism <- cbind(class = "UC1:Lipid metabolism",Fixed_effects_UC1_Lipid_metabolism)
colnames(Fixed_effects_UC1_Lipid_metabolism)[2] <- "var"

###Amino_acid_metabolism
model.f3 <- lmer(Amino_acid_metabolism ~ Lachnospiraceae + Porphyromonadaceae + Desulfovibrionaceae +
                   Prevotellaceae + Ruminococcaceae + (1|S), data = Driver_UC1_scale, REML=FALSE)
summary <- summary(model.f3)
Fixed_effects_UC1_Amino_acid <- adjustdata(summary$coefficients)
Fixed_effects_UC1_Amino_acid <- cbind(class = "UC1:Amino acid metabolism",Fixed_effects_UC1_Amino_acid)
colnames(Fixed_effects_UC1_Amino_acid)[2] <- "var"


UC_scatter_data <- rbind(Fixed_effects_UC1_Lipid_metabolism,Fixed_effects_UC1_Amino_acid)
                 
colnames(UC_scatter_data) <- c("class","var","Estimate","Std.Error","df","t value","P_value")
write.table(UC_scatter_data,file ="UC_scatter_data_UC1.txt",row.names = F,col.names = T, sep = "\t",quote = F)

##plot
##Figure9
infile <- "UC_scatter_data_UC1.txt"
scatter_data <- read_delim(infile)
dt4plot <- scatter_data %>% 
  filter(!var == "(Intercept)") %>%
  mutate(sig=ifelse(.$P_value <= 0.05, "yes", "no"))
dt4plot <- dt4plot[, c("class", "var", "Estimate", "P_value", "sig")]
dt4label <- dt4plot %>% filter(sig == "yes")

p1 <- ggplot(dt4plot, aes(x=Estimate, y=-log10(P_value), color=sig, group=class)) +
  geom_point() +
  geom_text_repel(
    aes(x=Estimate, y=-log10(P_value), label=var), data=dt4label, inherit.aes = FALSE) +
  facet_wrap(. ~ class, ncol=2) +
  theme_bw() +
  theme(legend.position="none", panel.grid = element_blank())
p1 <- p1 + labs(x = "Coefficient")
p1 <- p1 + labs(y = "-log10(P-value)")
p1 <- p1 + scale_color_manual(
  values = c("yes" = "#fc6c64", "no" = "#08bdc5"), 
  breaks = c("yes", "no"),                   
)

p1 <- p1 + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")
p1 <- p1 + geom_vline(xintercept = 0, linetype = "dashed", color = "grey")
ggsave(filename = "FigureS11.png", plot=p1, width=8, height = 2, dpi=600)






