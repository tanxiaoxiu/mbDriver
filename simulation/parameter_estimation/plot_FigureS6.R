library(tidyr)
library(ggplot2)
library(cowplot)
setwd("~/mbDriver/simulation/parameter_estimation/trajectory/theta3")

###p10
##Relative RMSE
load("p10/r_rmse_p10.Rdata")
data_df <- do.call(rbind, r_rmse_list)
data_long <- pivot_longer(data_df, cols = names(data_df), names_to = "Method", values_to = "Result")
data_long$Method <- factor(data_long$Method, levels = c("Ridge", "MLRR", "MLCRR", "BAL", "BVS"))

plot1.1 <- ggplot(data_long, aes(x = Method, y = Result, fill = Method)) +
  geom_boxplot(na.rm = TRUE) + 
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81")) + 
  labs(x = "Method", y = "Relative RMSE") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size=30),
        axis.title = element_text(size=30),
        axis.text = element_text(size=25),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color="black", size=0.5)) 
#ggsave(filename="r_rmse_p10.png",plot=plot1.1,device="png",dpi=600,units="in",width=8,height=8)

##Correlation
load("p10/corr_p10.Rdata")
data_df <- do.call(rbind, corr_list)
data_long <- pivot_longer(data_df, cols = names(data_df), names_to = "Method", values_to = "Result")
data_long$Method <- factor(data_long$Method, levels = c("Ridge", "MLRR", "MLCRR", "BAL", "BVS"))

plot1.2 <- ggplot(data_long, aes(x = Method, y = Result, fill = Method)) +
  geom_boxplot(na.rm = TRUE) + 
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81")) + 
  labs(x = "Method", y = "Correlation") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size=30),
        axis.title = element_text(size=30),
        axis.text = element_text(size=25),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color="black", size=0.5)) 
#ggsave(filename="corr_p10.png",plot=plot1.2,device="png",dpi=600,units="in",width=8,height=8)

###p15
##Relative RMSE
load("p15/r_rmse_p15.Rdata")
data_df <- do.call(rbind, r_rmse_list)
data_long <- pivot_longer(data_df, cols = names(data_df), names_to = "Method", values_to = "Result")
data_long$Method <- factor(data_long$Method, levels = c("Ridge", "MLRR", "MLCRR", "BAL", "BVS"))

plot2.1 <- ggplot(data_long, aes(x = Method, y = Result, fill = Method)) +
  geom_boxplot(na.rm = TRUE) + 
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81")) + 
  labs(x = "Method", y = "Relative RMSE") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size=30),
        axis.title = element_text(size=30),
        axis.text = element_text(size=25),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color="black", size=0.5)) 
#ggsave(filename="r_rmse_p15.png",plot=plot2.1,device="png",dpi=600,units="in",width=8,height=8)

##Correlation
load("p15/corr_p15.Rdata")
data_df <- do.call(rbind, corr_list)
data_long <- pivot_longer(data_df, cols = names(data_df), names_to = "Method", values_to = "Result")
data_long$Method <- factor(data_long$Method, levels = c("Ridge", "MLRR", "MLCRR", "BAL", "BVS"))

plot2.2 <- ggplot(data_long, aes(x = Method, y = Result, fill = Method)) +
  geom_boxplot(na.rm = TRUE) + 
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81")) + 
  labs(x = "Method", y = "Correlation") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size=30),
        axis.title = element_text(size=30),
        axis.text = element_text(size=25),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color="black", size=0.5)) 
#ggsave(filename="corr_p15.png",plot=plot2.2,device="png",dpi=600,units="in",width=8,height=8)


###n15
##Relative RMSE
load("n15/r_rmse_n15.Rdata")
data_df <- do.call(rbind, r_rmse_list)
data_long <- pivot_longer(data_df, cols = names(data_df), names_to = "Method", values_to = "Result")
data_long$Method <- factor(data_long$Method, levels = c("Ridge", "MLRR", "MLCRR", "BAL", "BVS"))

plot3.1 <- ggplot(data_long, aes(x = Method, y = Result, fill = Method)) +
  geom_boxplot(na.rm = TRUE) + 
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81")) + 
  labs(x = "Method", y = "Relative RMSE") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size=30),
        axis.title = element_text(size=30),
        axis.text = element_text(size=25),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color="black", size=0.5)) 
#ggsave(filename="r_rmse_n15.png",plot=plot3.1,device="png",dpi=600,units="in",width=8,height=8)

##Correlation
load("n15/corr_n15.Rdata")
data_df <- do.call(rbind, corr_list)
data_long <- pivot_longer(data_df, cols = names(data_df), names_to = "Method", values_to = "Result")
data_long$Method <- factor(data_long$Method, levels = c("Ridge", "MLRR", "MLCRR", "BAL", "BVS"))

plot3.2 <- ggplot(data_long, aes(x = Method, y = Result, fill = Method)) +
  geom_boxplot(na.rm = TRUE) + 
  scale_fill_manual(values = c("#df7676","#79add2", "#13679e", "#3c5587","#f09b81")) + 
  labs(x = "Method", y = "Correlation") +
  theme_minimal() +
  theme(legend.position = "none",
        text = element_text(size=30),
        axis.title = element_text(size=30),
        axis.text = element_text(size=25),
        panel.background = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(color="black", size=0.5)) 
#ggsave(filename="corr_n15.png",plot=plot3.2,device="png",dpi=600,units="in",width=8,height=8)

plot_FigureS6 <- plot_grid(plot1.1,plot2.1,plot3.1,plot1.2,plot2.2,plot3.2,ncol=3, nrow=2,labels = c("A","B","C","D","E","F"),label_size = 26)
ggsave(filename="plot_FigureS6.png",plot=plot_FigureS6,device="png",dpi=600,units="in",width=22,height=14)

