library(tidyverse)
library(readxl)
library(ggpubr)
library(pheatmap)
library(ggplot2)
library(ggprism)
setwd("~/mbDriver/real_data/Fiber_diet")

###Figure6B
MAX_N_OTU =10
l_tb <- read.delim("Fiber_driver_summary.txt",  sep = '\t', check.names = FALSE)
l_tb$Driver <- gsub("Parabacteroides-goldsteinii", "P. goldsteinii", l_tb$Driver)
l_tb$Driver <- gsub("Lachnospiraceae-NK4A136-group", "L. NK4A136-group", l_tb$Driver)
l_tb$Driver <- gsub("Bacteroides-acidifaciens", "B. acidifaciens", l_tb$Driver)
l_tb$Driver <- gsub("Akkermansia-muciniphila", "A. muciniphila", l_tb$Driver)

l_tb %>% pivot_wider(names_from = "Group" ,values_from = "Score")
w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "Score",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")
p1 <- pheatmap(w_tb, cluster_cols = F,angle_col=0,fontsize=14)
ggsave(filename="Figure6B.png",plot=p1,device="png",dpi=600,units="in",width=5,height=7)

###Figure6C
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

metabolite <-  read.delim("SCFA.csv",  sep = ',',  check.names = FALSE)
data_fiber_diet <- read.delim("data_fiber.txt",  sep = '\t', check.names = FALSE)
Group <- data_fiber_diet[,c(1:5)]
colnames(Group)[1] <- "SampleID"
Group <- Group[-nrow(Group), ]

metabolite_group <- merge(metabolite, Group, by = "SampleID", all = T)
metabolite_group <- metabolite_group[order(metabolite_group$Time), ]
metabolite_group$Time <- factor(metabolite_group$Time, level=unique(metabolite_group$Time))

species_top10 <- data_fiber_diet[,c(1,6:15)]
colnames(species_top10)[1] <- "SampleID"
metabolite_species <- merge(metabolite_group, species_top10, by = "SampleID", all = F)
selected_Times <- c(0, 1, 3, 5, 8, 13, 19, 25, 31)
metabolite_species <- metabolite_species[metabolite_species$Time %in% selected_Times, ]
#metabolite_group <- na.omit(metabolite_group)

metabolite_species[,13:22] <- log(metabolite_species[,13:22] +1)

metabolite_species$Group <- gsub("Control", "Con", metabolite_species$Group)
metabolite_species$Group <- gsub("Resistant starch", "Rs", metabolite_species$Group)
metabolite_species$Group <- gsub("Inulin", "In", metabolite_species$Group)
metabolite_species$Group <- factor(metabolite_species$Group, levels = c("Con", "Rs", "In"))
#Faecalibaculum
color=c("#1597A5","#FEB3AE","#FFC24B")
p2 <- ggplot(metabolite_species,aes(x=Group,y=Faecalibaculum))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("log(Abandance)")+
  ggtitle('Faecalibaculum') +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","Rs"),
                                 c("Con","In")),
              map_signif_level = T, 
              step_increase = 0.08,
              test = t.test, 
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),
              size=0.7,color="black")+
  theme_prism(
    base_fontface = "plain", 
    #base_family = "serif", 
    base_family = "sans",
    base_size = 20,  
    base_line_size = 0.8, 
    axis_text_angle = 0)+ 
  theme(legend.position = 'none')
ggsave(filename="Figure6C_1.png",plot=p2,device="png",dpi=600,units="in",width=6,height=5)

#Butyrate
p3 <- ggplot(metabolite_species,aes(x=Group,y=Butyrate))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("Abandance")+
  ggtitle('Butyrate') +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","Rs"),
                                 c("Con","In")),
              map_signif_level = T, 
              step_increase = 0.08,
              test = t.test, 
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),
              size=0.7,color="black")+
  theme_prism(
    base_fontface = "plain", 
    #base_family = "serif", 
    base_family = "sans",
    base_size = 20,  
    base_line_size = 0.8, 
    axis_text_angle = 0)+ 
  theme(legend.position = 'none')
ggsave(filename="Figure6C_2.png",plot=p3,device="png",dpi=600,units="in",width=6,height=5)
