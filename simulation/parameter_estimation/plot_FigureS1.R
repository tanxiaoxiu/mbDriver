library(purrr)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(gridExtra)
library(cowplot)

###Figure S1A
setwd("~/mbDriver/simulation/parameter_estimation/sparse")
main_result_baseDir = "~/mbDriver/simulation/parameter_estimation/sparse/"
main_result_timer= c("t8", "t13", "t18", "t25")
regulation_ID = "theta3"#c("theta1", "theta3", "theta5")

#helper
read_data <- function(tm, reg, file){
  dt <- read.delim(
    file.path(main_result_baseDir, tm, reg, file), sep="\t", check.names=FALSE)
  return(dt)
}

di_dt_list = map(
  main_result_timer, function(x) {df = read_data(
    x, regulation_ID, "di_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Methods"; 
    df$Methods <- gsub("lasso", "Lasso", df$Methods);
    df$Methods <- gsub("ridge", "Ridge", df$Methods);
    df$Methods <- gsub("elastic", "Elastic net", df$Methods);
    df$Group = "Difference"; return(df)})
names(di_dt_list) = main_result_timer

sp_dt_list = map(
    main_result_timer, function(x) {df=read_data(
      x, regulation_ID, "sp_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Methods"; 
      df$Methods <- gsub("lasso", "Lasso", df$Methods);
      df$Methods <- gsub("ridge", "Ridge", df$Methods);
      df$Methods <- gsub("elastic", "Elastic net", df$Methods);
      df$Group= "Spline";return(df)})
names(sp_dt_list) = main_result_timer
dt4plt_list = map2(di_dt_list, sp_dt_list, function(x, y) rbind(x, y))

plot_t8 <- dt4plt_list$t8
plot_t8 <- plot_t8[plot_t8$Methods != "LSE", ]
plot_t8$Methods <- factor(plot_t8$Methods, level=unique(plot_t8$Methods))

g <- ggplot(plot_t8, aes(Methods))
plot1 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic(m) * " = 8")) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##t13
plot_t13 <- dt4plt_list$t13
plot_t13 <- plot_t13[plot_t13$Methods != "LSE", ]
plot_t13$Methods <- factor(plot_t13$Methods, level=unique(plot_t13$Methods))

g <- ggplot(plot_t13, aes(Methods)) 
plot2 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c( "#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic(m) * " = 13")) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##t18
plot_t18 <- dt4plt_list$t18
plot_t18 <- plot_t18[plot_t18$Methods != "LSE", ]
plot_t18$Methods <- factor(plot_t18$Methods, level=unique(plot_t18$Methods))

g <- ggplot(plot_t18, aes(Methods))
plot3 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic(m) * " = 18")) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##t25
plot_t25 <- dt4plt_list$t25
plot_t25 <- plot_t25[plot_t25$Methods != "LSE", ]
plot_t25$Methods <- factor(plot_t25$Methods, level=unique(plot_t25$Methods))

g <- ggplot(plot_t25, aes(Methods))
plot4 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic(m) * " = 25")) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

plot_Figure_S1A <- plot_grid(plot1,plot2,plot3,plot4, ncol=2, nrow=2)
ggsave(filename="plot_Figure_S1A.png",plot=plot_Figure_S1A,device="png",dpi=600,units="in",width=12,height=12)

###Figure S1B
setwd("~/mbDriver/simulation/parameter_estimation/sparse")
main_result_baseDir = "~/mbDriver/simulation/parameter_estimation/sparse/"
main_result_timer= c("t8", "t13", "t18", "t25")
regulation_ID = "theta5"#c("theta1", "theta3", "theta5")

#helper
read_data <- function(tm, reg, file){
  dt <- read.delim(
    file.path(main_result_baseDir, tm, reg, file), sep="\t", check.names=FALSE)
  return(dt)
}

di_dt_list = map(
  main_result_timer, function(x) {df = read_data(
    x, regulation_ID, "di_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Methods"; 
    df$Methods <- gsub("lasso", "Lasso", df$Methods);
    df$Methods <- gsub("ridge", "Ridge", df$Methods);
    df$Methods <- gsub("elastic", "Elastic net", df$Methods);
    df$Group = "Difference"; return(df)})
names(di_dt_list) = main_result_timer

sp_dt_list = map(
  main_result_timer, function(x) {df=read_data(
    x, regulation_ID, "sp_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Methods"; 
    df$Methods <- gsub("lasso", "Lasso", df$Methods);
    df$Methods <- gsub("ridge", "Ridge", df$Methods);
    df$Methods <- gsub("elastic", "Elastic net", df$Methods);
    df$Group= "Spline";return(df)})
names(sp_dt_list) = main_result_timer
dt4plt_list = map2(di_dt_list, sp_dt_list, function(x, y) rbind(x, y))

plot_t8 <- dt4plt_list$t8
plot_t8 <- plot_t8[plot_t8$Methods != "LSE", ]
plot_t8$Methods <- factor(plot_t8$Methods, level=unique(plot_t8$Methods))

g <- ggplot(plot_t8, aes(Methods))
plot1 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic(m) * " = 8")) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##t13
plot_t13 <- dt4plt_list$t13
plot_t13 <- plot_t13[plot_t13$Methods != "LSE", ]
plot_t13$Methods <- factor(plot_t13$Methods, level=unique(plot_t13$Methods))

g <- ggplot(plot_t13, aes(Methods)) 
plot2 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c( "#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic(m) * " = 13")) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##t18
plot_t18 <- dt4plt_list$t18
plot_t18 <- plot_t18[plot_t18$Methods != "LSE", ]
plot_t18$Methods <- factor(plot_t18$Methods, level=unique(plot_t18$Methods))

g <- ggplot(plot_t18, aes(Methods))
plot3 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic(m) * " = 18")) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

##t25
plot_t25 <- dt4plt_list$t25
plot_t25 <- plot_t25[plot_t25$Methods != "LSE", ]
plot_t25$Methods <- factor(plot_t25$Methods, level=unique(plot_t25$Methods))

g <- ggplot(plot_t25, aes(Methods))
plot4 <- g + geom_col(aes(fill = Group,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#06a088","#df7676"))+
  xlab("")+
  ylab("Relative RMSE")+
  theme_bw()+ 
  ggtitle(expression(italic(m) * " = 25")) +
  theme(legend.position = "none")+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

plot_Figure_S1B <- plot_grid(plot1,plot2,plot3,plot4, ncol=2, nrow=2)
ggsave(filename="plot_Figure_S1B.png",plot=plot_Figure_S1B,device="png",dpi=600,units="in",width=12,height=12)
