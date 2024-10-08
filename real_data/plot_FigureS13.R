library(tidyverse)
library(readxl)
library(ggpubr)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(ggprism)
library(vegan)
library(picante)
library(dplyr)
library(RColorBrewer)

setwd("~/mbDriver/real_data/UC")

data_UC <- read.delim("data_UC.txt",  sep = '\t', check.names = FALSE)

Group <- data_UC[,c(1,4)]
colnames(Group) <- c("ID","Group")

otu <- data_UC[,-c(2:4)]
row.names(otu) <- otu[,1]
OTU <- otu[,-1]

#Alpha-diversity
df <- t(OTU)
Shannon <- diversity(OTU, index = "shannon", MARGIN = 1, base = exp(1))
Simpson <- diversity(OTU, index = "simpson", MARGIN = 1, base =  exp(1))
Richness <- specnumber(df, MARGIN = 2)
index <- as.data.frame(cbind(Shannon, Simpson, Richness))
tdf <- t(df)
tdf<-ceiling(as.data.frame(t(df)))
obs_chao_ace <- t(estimateR(tdf))
obs_chao_ace <- obs_chao_ace[rownames(index),]
index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]
index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(df ==1) / colSums(df)
#write.table(cbind(sample=c(rownames(index)),index),'diversity.index.txt', row.names = F, sep = '\t', quote = F)
index$samples <- rownames(index)
groups <- Group
colnames(groups)[1:2] <- c('samples','Group')
df2 <- merge(index,groups,by = 'samples')

###Figure_S13A
#Shannon
color=c("#1597A5","#6f6da1","#FFC24B","#FEB3AE")
p1 <- ggplot(df2,aes(x=Group,y=Shannon))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("H1","H2"),
                                 c("H1","UC1"),
                                 c("UC1","UC2"),
                                 c("H2","UC2")),
              map_signif_level = T, 
              step_increase = 0.08,
              test = wilcox.test, 
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),
              size=0.7,color="black")+
  theme_prism(
    base_fontface = "plain", 
    base_family = "serif", 
    base_size = 20,  
    base_line_size = 0.8, 
    axis_text_angle = 0)+ 
  
  theme(legend.position = 'none')
p1
ggsave(filename="FigureS13A.png",plot=p1,device="png",dpi=600,units="in",width=6,height=5)

###Figure_S13B
#Simpson
color=c("#1597A5","#6f6da1","#FFC24B","#FEB3AE")
p2 <- ggplot(df2,aes(x=Group,y=Simpson))+
  stat_boxplot(geom = "errorbar", width=0.1,size=0.8)+
  geom_boxplot(aes(fill=Group), 
               outlier.colour="white",size=0.8)+
  theme(panel.background =element_blank(), 
        axis.line=element_line(),
        plot.title = element_text(size=14))+
  scale_fill_manual(values=color)+
  geom_jitter(width = 0.2)+
  geom_signif(comparisons = list(c("H1","H2"),
                                 c("H1","UC1"),
                                 c("UC1","UC2"),
                                 c("H2","UC2")),
              map_signif_level = T, 
              step_increase = 0.08,
              test = wilcox.test,
              tip_length = c(c(0,0),
                             c(0,0),
                             c(0,0)),
              size=0.7,color="black")+
  theme_prism(
    base_fontface = "plain", 
    base_family = "serif", 
    base_size = 20, 
    base_line_size = 0.8, 
    axis_text_angle = 0)+ 
  theme(legend.position = 'none')
p2
ggsave(filename="FigureS13B.png",plot=p2,device="png",dpi=600,units="in",width=6,height=5)

###Figure_S13C
#Beta-diversity
#PCoA/bray-curtis
OTU_BC <- OTU
log_OTU <- log(OTU_BC + 1)
dist <- vegdist(log_OTU, method="bray", binary=F)
pcoa <- cmdscale(dist, k=(nrow(log_OTU) - 1), eig=T)
pcoa_points <- as.data.frame(pcoa$points)
sum_eig <- sum(pcoa$eig)
eig_percent <- round(pcoa$eig/sum_eig*100,1)
colnames(pcoa_points) <- paste0("PCoA", 1:3)
pcoa_points$ID <- rownames(pcoa_points)
pcoa_result <- merge(pcoa_points, Group, by = 'ID', all.x = TRUE)

dune.div <- adonis2(log_OTU ~ Group, data = Group, permutations = 999, method="bray")
dune_adonis2_2 <- paste0("R²=",round(dune.div$R2,2), ", P-value=", dune.div$`Pr(>F)`)

##PCoA
color=c("#1597A5","#6f6da1","#FFC24B","#FEB3AE")
p3 <- ggplot(pcoa_result,aes(x=PCoA1,y=PCoA2,
                             color= Group,shape=Group))+
  theme_bw()+
  geom_point(size=1.8)+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste0("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title= dune_adonis2_2
  )+
  stat_ellipse(data=pcoa_result,
               geom = "polygon",level=0.95,
               linetype = 2,size=0.5,
               aes(fill=Group),
               alpha=0.2,
               show.legend = T)+
  scale_color_manual(values = color) +
  scale_fill_manual(values = c("#1597A5","#6f6da1","#FFC24B","#FEB3AE"))+
  theme(axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14,angle=90),
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=12),
        panel.grid=element_blank())
p3
ggsave(filename="FigureS13C.png",plot=p3,device="png",dpi=600,units="in",width=6,height=5)





