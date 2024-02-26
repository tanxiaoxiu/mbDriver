#Figure5
library(purrr)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)

##plot Interaction,growth,time
###Figure5A
setwd("~/mbDriver/simulation/parameter_estimation/other_methods")
main_result_baseDir = "~/mbDriver/simulation/parameter_estimation/other_methods"
main_result_timer= c("t8","t13","t18","t25")
regulation_ID = "theta1"#c("theta1", "theta3", "theta5", "theta10")
#helper
read_data <- function(tm, reg, file){
  dt <- read.delim(
    file.path(main_result_baseDir, tm, reg, file), sep="\t", check.names=FALSE)
  return(dt)
}
mdsine_dt_list = map(
  main_result_timer, function(x) {df=read_data(
    x, regulation_ID, "mdsine_mean_re_rmse_tables.txt"); 
  colnames(df)[1] <- "Method"; 
  return(df)})
names(mdsine_dt_list) = main_result_timer

##Ridge based spline
setwd("~/mbDriver/simulation/parameter_estimation/sparse")
main_result_baseDir = "~/mbDriver/simulation/parameter_estimation/sparse"
main_result_timer= c("t8","t13","t18","t25")
regulation_ID = "theta1"#c("theta1", "theta3", "theta5", "theta10")
sp_dt_list = map(
  main_result_timer, function(x) {df=read_data(
    x, regulation_ID, "sp_mean_re_rmse_tables.txt"); colnames(df)[1] <- "Method"; 
    df$Method <- gsub("lasso", "Lasso", df$Method);
    df$Method <- gsub("ridge", "Ridge", df$Method);
    df$Method <- gsub("elastic", "Elastic net", df$Method);
    df <- df[df$Method == "Ridge", ]
    return(df)})
names(sp_dt_list) = main_result_timer
dt4plt_list = map2(sp_dt_list, mdsine_dt_list, function(x, y) rbind(x, y))
t8 <- dt4plt_list$t8
t8$time <- "8"
t13 <- dt4plt_list$t13
t13$time <- "13"
t18 <- dt4plt_list$t18
t18$time <- "18"
t25 <- dt4plt_list$t25
t25$time <- "25"
plot_summary <- rbind(t8,t13,t18,t25)
plot_summary$time <- factor(plot_summary$time, level=unique(plot_summary$time))
plot_summary$Method <- factor(plot_summary$Method, level=unique(plot_summary$Method))

setwd("~/mbDriver/simulation/parameter_estimation/")
g <- ggplot(plot_summary, aes(time))
plot1 <- g + geom_col(aes(fill = Method,y= RMSE_A), position = position_dodge())+
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81"))+
  xlab(expression(italic(m)))+
  ylab("Relative RMSE")+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

p1 <-  plot1 + theme(legend.spacing.y = unit(0.5, 'cm')) + guides(fill = guide_legend(byrow = TRUE))
ggsave(filename="Figure5A.png",plot=p1,device="png",dpi=600,units="in",width=12,height= 8)


###Figure5B
##plot growth
g <- ggplot(plot_summary, aes(time))
plot2 <- g + geom_col(aes(fill = Method,y= RMSE_r), position = position_dodge())+
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81"))+
  xlab(expression(italic(m)))+
  ylab("Relative RMSE")+
  theme_bw()+ 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

p2 <-  plot2 + theme(legend.spacing.y = unit(0.5, 'cm')) + guides(fill = guide_legend(byrow = TRUE))
ggsave(filename="Figure5B.png",plot=p2,device="png",dpi=600,units="in",width=12,height= 8)


###Figure5C
setwd("~/mbDriver/simulation/parameter_estimation/other_methods")
main_result_baseDir = "~/mbDriver/simulation/parameter_estimation/other_methods"
main_result_timer= c("t8","t13","t18","t25")
regulation_ID = "theta1"#c("theta1", "theta3", "theta5", "theta10")
#helper
read_data <- function(tm, reg, file){
  dt <- read.delim(
    file.path(main_result_baseDir, tm, reg, file), sep="\t", check.names=FALSE)
  return(dt)
}
mdsine_dt_list = map(
  main_result_timer, function(x) {df=read_data(
    x, regulation_ID, "mean_time_list.txt"); 
  colnames(df)[1] <- "Method"; 
  return(df)})
names(mdsine_dt_list) = main_result_timer

##Ridge based spline
setwd("~/mbDriver/simulation/parameter_estimation/sparse")
main_result_baseDir = "~/mbDriver/simulation/parameter_estimation/sparse"
main_result_timer= c("t8","t13","t18","t25")
regulation_ID = "theta1"#c("theta1", "theta3", "theta5", "theta10")
sp_dt_list = map(
  main_result_timer, function(x) {df=read_data(
    x, regulation_ID, "sp_mean_run_time_tables.txt"); colnames(df)[1] <- "Method"; 
    df$Method <- gsub("lasso", "Lasso", df$Method);
    df$Method <- gsub("ridge", "Ridge", df$Method);
    df$Method <- gsub("elastic", "Elastic net", df$Method);
    df <- df[df$Method == "Ridge", ]
    return(df)})
names(sp_dt_list) = main_result_timer

dt4plt_list = map2(sp_dt_list, mdsine_dt_list, function(x, y) rbind(x, y))

t8 <- dt4plt_list$t8
t8$group <- "8"
t13 <- dt4plt_list$t13
t13$group <- "13"
t18 <- dt4plt_list$t18
t18$group <- "18"
t25 <- dt4plt_list$t25
t25$group <- "25"

plot_summary <- rbind(t8,t13,t18,t25)
plot_summary$group <- factor(plot_summary$group, level=unique(plot_summary$group))
plot_summary$Method <- factor(plot_summary$Method, level=unique(plot_summary$Method))
plot_summary$time <- log(plot_summary$time + 1)

setwd("~/mbDriver/simulation/parameter_estimation/")
g <- ggplot(plot_summary, aes(group))
plot3 <- g + geom_col(aes(fill = Method,y= time), position = position_dodge())+
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81"))+
  xlab(expression(italic(m)))+
  ylab("log(Seconds + 1)")+
  theme_bw()+ 
  #ggtitle('Time') +
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))+
  theme(panel.grid.major=element_line(colour=NA),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank())

p3 <-  plot3 + theme(legend.spacing.y = unit(0.5, 'cm')) + guides(fill = guide_legend(byrow = TRUE))
ggsave(filename="Figure5C.png",plot=p3,device="png",dpi=600,units="in",width=12,height=8)


