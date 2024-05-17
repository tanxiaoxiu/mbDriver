setwd("~/mbDriver/real_data/Fiber_diet")

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

metabolite <-  read.delim("SCFA.csv",  sep = ',',  check.names = FALSE)
data_fiber_diet <- read.delim("data_fiber.txt",  sep = '\t', check.names = FALSE)
Group <- data_fiber_diet[,c(1:5)]
colnames(Group)[1] <- "SampleID"
Group <- Group[-nrow(Group), ]
species_top10 <- data_fiber_diet[,c(1,6:15)]
colnames(species_top10)[1] <- "SampleID"

metabolite_group <- merge(metabolite, Group, by = "SampleID", all = T)
metabolite_group <- metabolite_group[order(metabolite_group$Time), ]
metabolite_group$Time <- factor(metabolite_group$Time, level=unique(metabolite_group$Time))

metabolite_group <- na.omit(metabolite_group)
metabolite_species <- merge(metabolite_group, species_top10, by = "SampleID", all = F)

#################################Rs group
metabolite_species_Rs <- metabolite_species[metabolite_species$Group == "Resistant starch", ]

Driver_Rs <- metabolite_species_Rs[,c("Time","Total","Acetate","Propionate","Butyrate","Iso-butyrate","Iso-valerate","Valerate",
                                      "Parabacteroides-goldsteinii","Bacteroides-acidifaciens","Parasutterella",
                                      "Faecalibaculum","Lachnospiraceae",
                                      "Subject")]

colnames(Driver_Rs) <- c("Time","Total","Acetate","Propionate","Butyrate","Iso_butyrate","Iso_valerate","Valerate",
                         "P.goldsteinii","B.acidifaciens","Parasutterella",
                         "Faecalibaculum","Lachnospiraceae",
                         "Subject")


#Total
Driver_Rs_scale <- Driver_Rs
Driver_Rs_scale[,9:13] <- log(Driver_Rs_scale[,9:13]+1)

S <- Driver_Rs_scale$Subject
model.f <- lmer(Total ~ P.goldsteinii + Parasutterella + Faecalibaculum +
                  B.acidifaciens + Lachnospiraceae + (1|S), data = Driver_Rs_scale, REML=FALSE)
summary <- summary(model.f)
Fixed_effects_Rs_Total <- adjustdata(summary$coefficients)
Fixed_effects_Rs_Total <- cbind(class = "Total",Fixed_effects_Rs_Total)
colnames(Fixed_effects_Rs_Total)[2] <- "var"
#write.table(Fixed_effects_Rs_Total,file ="Fixed_effects_Rs_Total.txt",row.names = F,col.names = T, sep = "\t",quote = F)


###Butyrate
model.f2 <- lmer(Butyrate ~ P.goldsteinii + Parasutterella + Faecalibaculum +
                   B.acidifaciens + Lachnospiraceae + (1|S), data = Driver_Rs_scale, REML=FALSE)
summary <- summary(model.f2)
Fixed_effects_Rs_Butyrate <- adjustdata(summary$coefficients)
Fixed_effects_Rs_Butyrate <- cbind(class = "Butyrate",Fixed_effects_Rs_Butyrate)
colnames(Fixed_effects_Rs_Butyrate)[2] <- "var"
#write.table(Fixed_effects_Rs_Butyrate,file ="Fixed_effects_Rs_Butyrate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Iso_Butyrate
model.f3 <- lmer(Iso_butyrate ~ P.goldsteinii + Parasutterella + Faecalibaculum +
                   B.acidifaciens + Lachnospiraceae + (1|S), data = Driver_Rs_scale, REML=FALSE)

summary <- summary(model.f3)
Fixed_effects_Rs_Iso_butyrate <- adjustdata(summary$coefficients)
Fixed_effects_Rs_Iso_butyrate <- cbind(class = "Iso-butyrate",Fixed_effects_Rs_Iso_butyrate)
colnames(Fixed_effects_Rs_Iso_butyrate)[2] <- "var"
#write.table(Fixed_effects_Rs_Iso_butyrate,file ="Fixed_effects_Rs_Iso_butyrate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Acetate
model.f4 <- lmer(Acetate ~ P.goldsteinii + Parasutterella + Faecalibaculum +
                   B.acidifaciens + Lachnospiraceae + (1|S), data = Driver_Rs_scale, REML=FALSE)

summary <- summary(model.f4)
Fixed_effects_Rs_Acetate <- adjustdata(summary$coefficients)
Fixed_effects_Rs_Acetate <- cbind(class = "Acetate",Fixed_effects_Rs_Acetate)
colnames(Fixed_effects_Rs_Acetate)[2] <- "var"
#write.table(Fixed_effects_Rs_Acetate,file ="Fixed_effects_Rs_Acetate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Propionate
model.f5 <- lmer(Propionate ~ P.goldsteinii + Parasutterella + Faecalibaculum +
                   B.acidifaciens + Lachnospiraceae + (1|S), data = Driver_Rs_scale, REML=FALSE)

summary <- summary(model.f5)
Fixed_effects_Rs_Propionate <- adjustdata(summary$coefficients)
Fixed_effects_Rs_Propionate <- cbind(class = "Propionate",Fixed_effects_Rs_Propionate)
colnames(Fixed_effects_Rs_Propionate)[2] <- "var"
#write.table(Fixed_effects_Rs_Propionate,file ="Fixed_effects_Rs_Propionate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Valerate
model.f6 <- lmer(Valerate ~ P.goldsteinii + Parasutterella + Faecalibaculum +
                   B.acidifaciens + Lachnospiraceae + (1|S), data = Driver_Rs_scale, REML=FALSE)

summary <- summary(model.f6)
Fixed_effects_Rs_Valerate <- adjustdata(summary$coefficients)
Fixed_effects_Rs_Valerate <- cbind(class = "Valerate",Fixed_effects_Rs_Valerate)
colnames(Fixed_effects_Rs_Valerate)[2] <- "var"
#write.table(Fixed_effects_Rs_Valerate,file ="Fixed_effects_Rs_Valerate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Iso-valerate
model.f7 <- lmer(Iso_valerate ~ P.goldsteinii + Parasutterella + Faecalibaculum +
                   B.acidifaciens + Lachnospiraceae + (1|S), data = Driver_Rs_scale, REML=FALSE)

summary <- summary(model.f7)
Fixed_effects_Rs_Iso_valerate <- adjustdata(summary$coefficients)
Fixed_effects_Rs_Iso_valerate <- cbind(class = "Iso-valerate",Fixed_effects_Rs_Iso_valerate)
colnames(Fixed_effects_Rs_Iso_valerate)[2] <- "var"
#write.table(Fixed_effects_Rs_Iso_valerate,file ="Fixed_effects_Rs_Iso_valerate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

RS_scatter_data <- rbind(Fixed_effects_Rs_Total,Fixed_effects_Rs_Butyrate,Fixed_effects_Rs_Iso_butyrate,
                         Fixed_effects_Rs_Acetate,Fixed_effects_Rs_Propionate,
                         Fixed_effects_Rs_Valerate,Fixed_effects_Rs_Iso_valerate)
colnames(RS_scatter_data) <- c("class","var","Estimate","Std.Error","df","t value","P_value")
write.table(RS_scatter_data,file ="RS_scatter_data.txt",row.names = F,col.names = T, sep = "\t",quote = F)



#################################In group
metabolite_species_In <- metabolite_species[metabolite_species$Group == "Inulin", ]

Driver_In <- metabolite_species_In[,c("Time","Total","Acetate","Propionate","Butyrate","Iso-butyrate","Iso-valerate","Valerate",
                                      "Akkermansia-muciniphila","Muribaculaceae",
                                      "Faecalibaculum","Bacteroides-acidifaciens","Alloprevotella",
                                      "Subject")]

colnames(Driver_In) <- c("Time","Total","Acetate","Propionate","Butyrate","Iso_butyrate","Iso_valerate","Valerate",
                         "A.muciniphila","Muribaculaceae",
                         "Faecalibaculum","B.acidifaciens","Alloprevotella",
                         "Subject")


Driver_In_scale <- Driver_In
Driver_In_scale[,9:13] <- log(Driver_In_scale[,9:13]+1)

S <- Driver_In_scale$Subject
model.f <- lmer(Total ~ A.muciniphila + Muribaculaceae + Faecalibaculum +
                  B.acidifaciens + Alloprevotella  + (1|S), data = Driver_In_scale, REML=FALSE)
summary <- summary(model.f)
Fixed_effects_In_Total <- adjustdata(summary$coefficients)
Fixed_effects_In_Total <- cbind(class = "Total",Fixed_effects_In_Total)
colnames(Fixed_effects_In_Total)[2] <- "var"
#write.table(Fixed_effects_In_Total,file ="Fixed_effects_In_Total.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Butyrate
model.f2 <- lmer(Butyrate ~ A.muciniphila + Muribaculaceae + Faecalibaculum +
                   B.acidifaciens + Alloprevotella  + (1|S), data = Driver_In_scale, REML=FALSE)
summary <- summary(model.f2)
Fixed_effects_In_Butyrate <- adjustdata(summary$coefficients)
Fixed_effects_In_Butyrate <- cbind(class = "Butyrate",Fixed_effects_In_Butyrate)
colnames(Fixed_effects_In_Butyrate)[2] <- "var"
#write.table(Fixed_effects_In_Butyrate,file ="Fixed_effects_In_Butyrate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Iso_Butyrate
model.f3 <- lmer(Iso_butyrate ~ A.muciniphila + Muribaculaceae + Faecalibaculum +
                   B.acidifaciens + Alloprevotella  + (1|S), data = Driver_In_scale, REML=FALSE)
summary <- summary(model.f3)
Fixed_effects_In_Iso_butyrate <- adjustdata(summary$coefficients)
Fixed_effects_In_Iso_butyrate <- cbind(class = "Iso-butyrate",Fixed_effects_In_Iso_butyrate)
colnames(Fixed_effects_In_Iso_butyrate)[2] <- "var"
#write.table(Fixed_effects_In_Iso_butyrate,file ="Fixed_effects_In_Iso_butyrate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Acetate
model.f4 <- lmer(Acetate ~ A.muciniphila + Muribaculaceae + Faecalibaculum +
                   B.acidifaciens + Alloprevotella  + (1|S), data = Driver_In_scale, REML=FALSE)
summary <- summary(model.f4)
Fixed_effects_In_Acetate <- adjustdata(summary$coefficients)
Fixed_effects_In_Acetate <- cbind(class = "Acetate",Fixed_effects_In_Acetate)
colnames(Fixed_effects_In_Acetate)[2] <- "var"
#write.table(Fixed_effects_In_Acetate,file ="Fixed_effects_In_Acetate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Propionate
model.f5 <- lmer(Propionate ~ A.muciniphila + Muribaculaceae + Faecalibaculum +
                   B.acidifaciens + Alloprevotella  + (1|S), data = Driver_In_scale, REML=FALSE)
summary <- summary(model.f5)
Fixed_effects_In_Propionate <- adjustdata(summary$coefficients)
Fixed_effects_In_Propionate <- cbind(class = "Propionate",Fixed_effects_In_Propionate)
colnames(Fixed_effects_In_Propionate)[2] <- "var"
#write.table(Fixed_effects_In_Propionate,file ="Fixed_effects_In_Propionate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Valerate
model.f6 <- lmer(Valerate ~ A.muciniphila + Muribaculaceae + Faecalibaculum +
                   B.acidifaciens + Alloprevotella  + (1|S), data = Driver_In_scale, REML=FALSE)
summary <- summary(model.f6)
Fixed_effects_In_Valerate <- adjustdata(summary$coefficients)
Fixed_effects_In_Valerate <- cbind(class = "Valerate",Fixed_effects_In_Valerate)
colnames(Fixed_effects_In_Valerate)[2] <- "var"
#write.table(Fixed_effects_In_Valerate,file ="Fixed_effects_In_Valerate.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Iso-valerate
model.f7 <- lmer(Iso_valerate ~ A.muciniphila + Muribaculaceae + Faecalibaculum +
                   B.acidifaciens + Alloprevotella  + (1|S), data = Driver_In_scale, REML=FALSE)
summary <- summary(model.f7)
Fixed_effects_In_Iso_valerate <- adjustdata(summary$coefficients)
Fixed_effects_In_Iso_valerate <- cbind(class = "Iso-valerate",Fixed_effects_In_Iso_valerate)
colnames(Fixed_effects_In_Valerate)[2] <- "var"
#write.table(Fixed_effects_In_Iso_valerate,file ="Fixed_effects_In_Iso_valerate.txt",row.names = F,col.names = T, sep = "\t",quote = F)
In_scatter_data <- rbind(Fixed_effects_In_Total,Fixed_effects_In_Butyrate,Fixed_effects_In_Iso_butyrate,
                         Fixed_effects_In_Acetate,Fixed_effects_In_Propionate,
                         Fixed_effects_In_Valerate,Fixed_effects_In_Iso_valerate)
colnames(In_scatter_data) <- c("class","var","Estimate","Std.Error","df","t value","P_value")
write.table(In_scatter_data,file ="In_scatter_data.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###Figure7A
#RS group
infile <- "RS_scatter_data.txt"
scatter_data <- read_delim(infile)
dt4plot <- scatter_data %>% 
  filter(!class == "Total", !var == "(Intercept)") %>%
  mutate(sig=ifelse(.$P_value <= 0.05, "yes", "no"))

dt4plot <- dt4plot[, c("class", "var", "Estimate", "P_value", "sig")]

dt4plot_factor <- dt4plot %>% 
  mutate(class=fct_relevel(class,c("Acetate","Butyrate","Propionate","Valerate","Iso-butyrate","Iso-valerate")))%>% 
  arrange(class)

dt4label <- dt4plot_factor %>% filter(sig == "yes")

p1 <- ggplot(dt4plot_factor, aes(x=Estimate, y=-log10(P_value), color=sig, group=class)) +
  #geom_jitter(width = 0, height = 0) + 
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
ggsave(filename = "Figure7A.png", plot=p1, width=8, height = 6, dpi=600)

###Figure7B
###In group
infile <- "In_scatter_data.txt"
scatter_data <- read_delim(infile)
dt4plot <- scatter_data %>% 
  filter(!class == "Total", !var == "(Intercept)") %>%
  mutate(sig=ifelse(.$P_value <= 0.05, "yes", "no")) 

dt4plot <- dt4plot[, c("class", "var", "Estimate", "P_value", "sig")]

dt4plot_factor <- dt4plot %>% 
  mutate(class=fct_relevel(class,c("Acetate","Butyrate","Propionate","Valerate","Iso-butyrate","Iso-valerate")))%>% 
  arrange(class)
dt4label <- dt4plot_factor %>% filter(sig == "yes")

p2 <- ggplot(dt4plot_factor, aes(x=Estimate, y=-log10(P_value), color=sig, group=class)) +
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
ggsave(filename = "Figure7B.png", plot=p2, width=8, height = 6, dpi=600)


