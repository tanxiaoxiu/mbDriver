library(tidyverse)
library(readxl)
library(ggpubr)
library(ggplot2)

###plot FigureS10A
setwd("~/mbDriver/real_data/Fiber_diet/p10")
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

###FigureS10A
#Faecalibaculum

species_top10 <- data_fiber_diet[,c(1,6:15)]
colnames(species_top10)[1] <- "SampleID"
metabolite_species <- merge(metabolite_group, species_top10, by = "SampleID", all = F)
selected_days <- c(0, 1, 3, 5, 8, 13, 19, 25, 31)
metabolite_species <- metabolite_species[metabolite_species$Time %in% selected_days, ]
#metabolite_group <- na.omit(metabolite_group)
metabolite_species_Con <- metabolite_species[metabolite_species$Group == "Control", ]

Faecalibaculum_Con <- metabolite_species_Con[,c(14,11,12)]
Faecalibaculum_Con_mean <- Faecalibaculum_Con %>%
  group_by(Time) %>%
  summarize(Faecalibaculum = mean(Faecalibaculum),
  )
Faecalibaculum_Con_mean$group <- "Con"

metabolite_species_Rs <- metabolite_species[metabolite_species$Group == "Resistant starch", ]
Faecalibaculum_Rs <- metabolite_species_Rs[,c(14,11,12)]
Faecalibaculum_Rs_mean <- Faecalibaculum_Rs %>%
  group_by(Time) %>%
  summarize(Faecalibaculum = mean(Faecalibaculum),
  )
Faecalibaculum_Rs_mean$group <- "Rs"

metabolite_species_In <- metabolite_species[metabolite_species$Group == "Inulin", ]
Faecalibaculum_In <- metabolite_species_In[,c(14,11,12)]
Faecalibaculum_In_mean <- Faecalibaculum_In %>%
  group_by(Time) %>%
  summarize(Faecalibaculum = mean(Faecalibaculum),
  )
Faecalibaculum_In_mean$group <- "In"


Faecalibaculum_mean_group <- rbind(Faecalibaculum_Con_mean,Faecalibaculum_Rs_mean,Faecalibaculum_In_mean) 
Faecalibaculum_mean_group$Faecalibaculum <- log(Faecalibaculum_mean_group$Faecalibaculum)
selected_days <- c(0, 1, 3, 5, 8, 13, 19, 25, 31)
Faecalibaculum_mean_group <- Faecalibaculum_mean_group[Faecalibaculum_mean_group$Time %in% selected_days, ]

p1 = ggplot(data = Faecalibaculum_mean_group,aes(x=Time,y=Faecalibaculum,group = group,color = group))+
  geom_point(size=1.2)+ 
  geom_line(size=1.2)+
  xlab("Time")+
  ylab("log(Abandance)")+
  theme_bw()+ 
  ggtitle('Faecalibaculum') +
  scale_color_manual(values = c("#1597A5", "#FFC24B", "#FEB3AE")) +
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=16))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(color="black"))

p1 <- p1 + guides(color = FALSE)
p1 <- p1 + annotate("text", x = 1.7, y = 17.4, label = "P-value=6.24e-30")
p1 <- p1 + annotate("text", x = 1.7, y = 17.1, label = "P-value=9.67e-26")
ggsave(filename="FigureS10A.png",plot=p1,device="png",dpi=600,width=6,height=3)


###FigureS10B
#Butyrate
selected_days <- c(0, 1, 3, 5, 8, 13, 19, 25, 31)
metabolite_group <- metabolite_group[metabolite_group$Time %in% selected_days, ]
metabolite_group <- na.omit(metabolite_group)

metabolite_Con <- metabolite_group[metabolite_group$Group == "Control", ]
metabolite_Rs <- metabolite_group[metabolite_group$Group == "Resistant starch", ]
metabolite_In <- metabolite_group[metabolite_group$Group == "Inulin", ]

##
Butyrate_Con_Rs <- t.test(metabolite_Con$Butyrate, metabolite_Rs$Butyrate)
Butyrate_Con_In <- t.test(metabolite_Con$Butyrate, metabolite_In$Butyrate)

dune_adonis2_2 <- "P-value<0.05"

#In
metabolite_In <- na.omit(metabolite_In)
Butyrate_In <- metabolite_In[,c(4,8,11,12)]
Butyrate_In_mean <- Butyrate_In %>%
  group_by(Time) %>%
  summarize(Butyrate = mean(Butyrate),
            Total = mean(Total),
  )
Butyrate_In_mean$group <- "In"

#Con
metabolite_Con <- na.omit(metabolite_Con)
Butyrate_Con <- metabolite_Con[,c(4,8,11,12)]
Butyrate_Con_mean <- Butyrate_Con %>%
  group_by(Time) %>%
  summarize(Butyrate = mean(Butyrate),
            Total = mean(Total),
  )
Butyrate_Con_mean$group <- "Con"

#Rs
metabolite_Rs <- na.omit(metabolite_Rs)
Butyrate_Rs <- metabolite_Rs[,c(4,8,11,12)]
Butyrate_Rs_mean <- Butyrate_Rs %>%
  group_by(Time) %>%
  summarize(Butyrate = mean(Butyrate),
            Total = mean(Total),
  )
Butyrate_Rs_mean$group <- "Rs"

Butyrate_mean_group <- rbind(Butyrate_Con_mean,Butyrate_Rs_mean,Butyrate_In_mean) 
selected_days <- c(0, 1, 3, 5, 8, 13, 19, 25, 31)
Butyrate_mean_group <- Butyrate_mean_group[Butyrate_mean_group$Time %in% selected_days, ]

Butyrate_mean_group$group <- as.factor(Butyrate_mean_group$group)

p2 = ggplot(data = Butyrate_mean_group,aes(x=Time,y= Butyrate,
                                           group = group,color = group))+
  geom_point(size=1.2)+ 
  geom_line(size=1.2)+
  xlab("Time")+
  ylab("Abandance")+
  theme_bw()+ 
  ggtitle('Butyrate') +
  scale_color_manual(values = c("#1597A5", "#FFC24B", "#FEB3AE")) +
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=16))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        legend.box.background = element_rect(color="black")
  )

p2 <- p2 + guides(color = FALSE)
p2 <- p2 + annotate("text", x = 1.7, y = 8, label = "P-value=2.15e-29")
p2 <- p2 + annotate("text", x = 1.7, y = 7.3, label = "P-value=2.35e-30")
ggsave(filename="FigureS10B.png",plot=p2,device="png",dpi=600,width=6,height=3)

