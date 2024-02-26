library(tidyverse)
library(readxl)
library(ggpubr)
library(ggplot2)
setwd("~/mbDriver/real_data/Fiber_diet")

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
metabolite_species[,13:22] <- log(metabolite_species[,13:22] +1)
metabolite_species$Group <- gsub("Control", "Con", metabolite_species$Group)
metabolite_species$Group <- gsub("Resistant starch", "Rs", metabolite_species$Group)
metabolite_species$Group <- gsub("Inulin", "In", metabolite_species$Group)


##FigS6A
#Rs
subset_df <- metabolite_species[metabolite_species$Group %in% c("Con", "Rs"), ]

#Faecalibaculum
color=c("#1597A5","#FEB3AE")
p1 <- ggplot(subset_df,aes(x=Group,y=Faecalibaculum))+
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
  geom_signif(comparisons = list(c("Con","Rs")),
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
ggsave(filename="Rs_Faecalibaculum.png",plot=p1,device="png",dpi=600,units="in",width=6,height=5)


#Bacteroides acidifaciens
color=c("#1597A5","#FEB3AE")
p2 <- ggplot(subset_df,aes(x = Group,y = subset_df$`Bacteroides-acidifaciens`))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("log(Abandance)")+
  ggtitle('Bacteroides acidifaciens') +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","Rs")),
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
ggsave(filename="Rs_Bacteroides acidifaciens.png",plot=p2,device="png",dpi=600,units="in",width=6,height=5)

##Parabacteroides goldsteinii
color=c("#1597A5","#FEB3AE")
p3 <- ggplot(subset_df,aes(x = Group,y = subset_df$`Parabacteroides-goldsteinii`))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("log(Abandance)")+
  ggtitle('Parabacteroides goldsteinii') +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","Rs")),
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
ggsave(filename="Rs_Parabacteroides goldsteinii.png",plot=p3,device="png",dpi=600,units="in",width=6,height=5)

# Parasutterella
color=c("#1597A5","#FEB3AE")
p4 <- ggplot(subset_df,aes(x = Group,y = subset_df$Parasutterella))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("log(Abandance)")+
  ggtitle("Parasutterella") +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","Rs")),
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
ggsave(filename="Rs_Parasutterella.png",plot=p4,device="png",dpi=600,units="in",width=6,height=5)

#Lachnospiraceae
color=c("#1597A5","#FEB3AE")
p5 <- ggplot(subset_df,aes(x = Group,y = subset_df$Lachnospiraceae))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("log(Abandance)")+
  ggtitle("Lachnospiraceae") +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","Rs")),
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
ggsave(filename="Rs_Lachnospiraceae.png",plot=p5,device="png",dpi=600,units="in",width=6,height=5)


##FigS6B
###In
#Faecalibaculum
subset_df <- metabolite_species[metabolite_species$Group %in% c("Con", "In"), ]

#Faecalibaculum
color=c("#1597A5","#FFC24B")
p1 <- ggplot(subset_df,aes(x=Group,y=Faecalibaculum))+
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
  geom_signif(comparisons = list(c("Con","In")),
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
ggsave(filename="In_Faecalibaculum.png",plot=p1,device="png",dpi=600,units="in",width=6,height=5)

#Bacteroides acidifaciens
p2 <- ggplot(subset_df,aes(x = Group,y = subset_df$`Bacteroides-acidifaciens`))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("log(Abandance)")+
  ggtitle('Bacteroides acidifaciens') +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","In")),
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
ggsave(filename="In_Bacteroides acidifaciens.png",plot=p2,device="png",dpi=600,units="in",width=6,height=5)

#Akkermansia-muciniphila
p3 <- ggplot(subset_df,aes(x = Group,y = subset_df$`Akkermansia-muciniphila`))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("log(Abandance)")+
  ggtitle('Akkermansia muciniphila') +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","In")),
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
ggsave(filename="In_Akkermansia-muciniphila.png",plot=p3,device="png",dpi=600,units="in",width=6,height=5)

#Muribaculaceae
p4 <- ggplot(subset_df,aes(x = Group,y = subset_df$Muribaculaceae))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("log(Abandance)")+
  ggtitle('Muribaculaceae') +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","In")),
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
ggsave(filename="In_Muribaculaceae.png",plot=p4,device="png",dpi=600,units="in",width=6,height=5)

#Alloprevotella
p5 <- ggplot(subset_df,aes(x = Group,y = subset_df$Alloprevotella))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  ylab("log(Abandance)")+
  ggtitle('Alloprevotella') +
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("Con","In")),
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
ggsave(filename="In_Alloprevotella.png",plot=p5,device="png",dpi=600,units="in",width=6,height=5)

