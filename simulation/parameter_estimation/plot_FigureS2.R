###FigureS2
library(purrr)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(gridExtra)
library(cowplot)


setwd("~/mbDriver/simulation/parameter_estimation/supplementary/p15")
main_result_baseDir = "~/mbDriver/simulation/parameter_estimation/supplementary/p15"

###Figure S2A
main_result_timer= c("t13")
regulation_ID = c("theta1","theta3","theta5")#c("theta1", "theta3", "theta5")

#helper
read_data <- function(tm, reg, file){
  dt <- read.delim(
    file.path(main_result_baseDir, tm, reg, file), sep="\t", check.names=FALSE)
  return(dt)
}

di_dt_list = map(
  regulation_ID, function(x) {df = read_data(
    main_result_timer, x, "di_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Methods"; 
    df$Methods <- gsub("lasso", "Lasso", df$Methods);
    df$Methods <- gsub("ridge", "Ridge", df$Methods);
    df$Methods <- gsub("elastic", "Elastic net", df$Methods);
    df$Group = "Difference"; return(df)})
names(di_dt_list) = regulation_ID

sp_dt_list = map(
  regulation_ID, function(x) {df=read_data(
    main_result_timer, x, "sp_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Methods"; 
    df$Methods <- gsub("lasso", "Lasso", df$Methods);
    df$Methods <- gsub("ridge", "Ridge", df$Methods);
    df$Methods <- gsub("elastic", "Elastic net", df$Methods);
    df$Group= "Spline";return(df)})
names(sp_dt_list) = regulation_ID

dt4plt_list = map2(di_dt_list, sp_dt_list, function(x, y) rbind(x, y))

##theta1
plot_theta1 <- dt4plt_list$theta1
plot_theta1 <- plot_theta1[plot_theta1$Methods != "LSE", ]
plot_theta1$Methods <- factor(plot_theta1$Methods, level=unique(plot_theta1$Methods))
g <- ggplot(plot_theta1, aes(Methods)) 
plot1 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c( "#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic("\U03C6") == 1)) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##theta3
plot_theta3 <- dt4plt_list$theta3
plot_theta3 <- plot_theta3[plot_theta3$Methods != "LSE", ]
plot_theta3$Methods <- factor(plot_theta3$Methods, level=unique(plot_theta3$Methods))

g <- ggplot(plot_theta3, aes(Methods)) 
plot2 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c( "#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic("\U03C6") == 3)) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##theta5
plot_theta5 <- dt4plt_list$theta5
plot_theta5 <- plot_theta5[plot_theta5$Methods != "LSE", ]
plot_theta5$Methods <- factor(plot_theta5$Methods, level=unique(plot_theta5$Methods))
g <- ggplot(plot_theta5, aes(Methods)) 
plot3 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c( "#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic("\U03C6") == 5)) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
plot_FigureS2A <- plot_grid(plot1,plot2,plot3, ncol=3, nrow=1)
ggsave(filename="plot_FigureS2A.png",plot=plot_FigureS2A,device="png",dpi=600,units="in",width=18,height=6)


###Figure S2B
main_result_timer= c("t18")
regulation_ID = c("theta1","theta3","theta5")

#helper
read_data <- function(tm, reg, file){
  dt <- read.delim(
    file.path(main_result_baseDir, tm, reg, file), sep="\t", check.names=FALSE)
  return(dt)
}

di_dt_list = map(
  regulation_ID, function(x) {df = read_data(
    main_result_timer, x, "di_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Methods"; 
    df$Methods <- gsub("lasso", "Lasso", df$Methods);
    df$Methods <- gsub("ridge", "Ridge", df$Methods);
    df$Methods <- gsub("elastic", "Elastic net", df$Methods);
    df$Group = "Difference"; return(df)})
names(di_dt_list) = regulation_ID

sp_dt_list = map(
  regulation_ID, function(x) {df=read_data(
    main_result_timer, x, "sp_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Methods"; 
    df$Methods <- gsub("lasso", "Lasso", df$Methods);
    df$Methods <- gsub("ridge", "Ridge", df$Methods);
    df$Methods <- gsub("elastic", "Elastic net", df$Methods);
    df$Group= "Spline";return(df)})
names(sp_dt_list) = regulation_ID

dt4plt_list = map2(di_dt_list, sp_dt_list, function(x, y) rbind(x, y))

##theta1
plot_theta1 <- dt4plt_list$theta1
plot_theta1 <- plot_theta1[plot_theta1$Methods != "LSE", ]
plot_theta1$Methods <- factor(plot_theta1$Methods, level=unique(plot_theta1$Methods))
g <- ggplot(plot_theta1, aes(Methods)) 
plot4 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c( "#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic("\U03C6") == 1)) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##theta3
plot_theta3 <- dt4plt_list$theta3
plot_theta3 <- plot_theta3[plot_theta3$Methods != "LSE", ]
plot_theta3$Methods <- factor(plot_theta3$Methods, level=unique(plot_theta3$Methods))

g <- ggplot(plot_theta3, aes(Methods)) 
plot5 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c( "#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic("\U03C6") == 3)) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##theta5
plot_theta5 <- dt4plt_list$theta5
plot_theta5 <- plot_theta5[plot_theta5$Methods != "LSE", ]
plot_theta5$Methods <- factor(plot_theta5$Methods, level=unique(plot_theta5$Methods))
g <- ggplot(plot_theta5, aes(Methods)) 
plot6 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c( "#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic("\U03C6") == 5)) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())
plot_FigureS2B <- plot_grid(plot4,plot5,plot6, ncol=3, nrow=1)
ggsave(filename="plot_FigureS2B.png",plot=plot_FigureS2B,device="png",dpi=600,units="in",width=18,height=6)
