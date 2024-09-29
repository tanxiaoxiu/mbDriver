###FigureS4
library(purrr)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(gridExtra)
library(cowplot)

setwd("~/mbDriver/simulation/parameter_estimation/supplementary/other_methods")
main_result_baseDir = "~/mbDriver/simulation/parameter_estimation/supplementary/other_methods"
main_result_timer= c("p10","p15","n15")
regulation_ID = c("t13/theta3")
#helper
read_data <- function(tm, reg, file){
  dt <- read.delim(
    file.path(main_result_baseDir, tm, reg, file), sep="\t", check.names=FALSE)
  return(dt)
}

#helper
read_data <- function(tm, reg, file){
  dt <- read.delim(
    file.path(main_result_baseDir, tm, reg, file), sep="\t", check.names=FALSE)
  return(dt)
}

other_methods_dt_list = map(
  main_result_timer, function(x) {df=read_data(
    x, regulation_ID, "mdsine_mean_re_rmse_tables.txt"); 
  colnames(df)[1] <- "Method"; 
  return(df)})
names(other_methods_dt_list) = main_result_timer


##Ridge based spline
setwd("~/mbDriver/simulation/parameter_estimation/supplementary")
main_result_baseDir = "~/mbDriver/simulation/parameter_estimation/supplementary"
main_result_timer= c("p10","p15","n15")
regulation_ID = c("t13/theta3")

sp_dt_list = map(
  main_result_timer, function(x) {df=read_data(
    x, regulation_ID, "sp_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Method"; 
    df$Method <- gsub("lasso", "Lasso", df$Method);
    df$Method <- gsub("ridge", "Ridge", df$Method);
    df$Method <- gsub("elastic", "Elastic net", df$Method);
    df <- df[df$Method == "Ridge", ]
    return(df)})
names(sp_dt_list) = main_result_timer
dt4plt_list = map2(sp_dt_list, other_methods_dt_list, function(x, y) rbind(x, y))

#p10
#RMSE_Interaction
plot_p10 <- dt4plt_list$p10
plot_p10$Method <- factor(plot_p10$Method, level=unique(plot_p10$Method))
g <- ggplot(plot_p10, aes(Method))
plot1.1 <- g + geom_col(aes(fill = Method,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81"))+
  ylab("Relative RMSE")+
  theme_bw()+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
#ggsave(filename="RMSE_Interaction_p10.png",plot=plot1.1,device="png",dpi=600,units="in",width=8,height=8)

#RMSE_Growth
g <- ggplot(plot_p10, aes(Method))
plot1.2 <- g + geom_col(aes(fill = Method,y= RMSE_r), position = position_dodge())+
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81"))+
  ylab("Relative RMSE")+
  theme_bw()+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
#ggsave(filename="RMSE_Growth_p10.png",plot=plot1.2,device="png",dpi=600,units="in",width=8,height=8)

#p15
#RMSE_Interaction
plot_p15 <- dt4plt_list$p15
plot_p15$Method <- factor(plot_p15$Method, level=unique(plot_p15$Method))
g <- ggplot(plot_p15, aes(Method))
plot2.1 <- g + geom_col(aes(fill = Method,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81"))+
  ylab("Relative RMSE")+
  theme_bw()+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
#ggsave(filename="RMSE_Interaction_p15.png",plot=plot2.1,device="png",dpi=600,units="in",width=8,height=8)

#RMSE_Growth
g <- ggplot(plot_p15, aes(Method))
plot2.2 <- g + geom_col(aes(fill = Method,y= RMSE_r), position = position_dodge())+
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81"))+
  ylab("Relative RMSE")+
  theme_bw()+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
#ggsave(filename="RMSE_Growth_p15.png",plot=plot2.2,device="png",dpi=600,units="in",width=8,height=8)

#n15
#RMSE_Interaction
plot_n15 <- dt4plt_list$n15
plot_n15$Method <- factor(plot_n15$Method, level=unique(plot_n15$Method))
g <- ggplot(plot_n15, aes(Method))
plot3.1 <- g + geom_col(aes(fill = Method,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81"))+
  ylab("Relative RMSE")+
  theme_bw()+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
#ggsave(filename="RMSE_Interaction_n15.png",plot=plot3.1,device="png",dpi=600,units="in",width=8,height=8)

#RMSE_Growth
g <- ggplot(plot_n15, aes(Method))
plot3.2 <- g + geom_col(aes(fill = Method,y= RMSE_r), position = position_dodge())+
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81"))+
  ylab("Relative RMSE")+
  theme_bw()+ 
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

#ggsave(filename="RMSE_Growth_n15.png",plot=plot3.2,device="png",dpi=600,units="in",width=8,height=8)

plot_FigureS4 <- plot_grid(plot1.1,plot2.1,plot3.1,plot1.2,plot2.2,plot3.2,ncol=3, nrow=2,labels = c("A","B","C","D","E","F"))
ggsave(filename="plot_FigureS4.png",plot=plot_FigureS4,device="png",dpi=600,units="in",width=20,height=16)

