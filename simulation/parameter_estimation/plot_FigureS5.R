library(tidyr)
library(ggplot2)
setwd("~/mbDriver/simulation/trajectory/theta1")

###t8
##Relative RMSE
load("t8/r_rmse_t8.Rdata")
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
plot1.1
#ggsave(filename="r_rmse_t8.png",plot=plot1.1,device="png",dpi=600,units="in",width=8,height=8)

##Correlation
load("t8/corr_t8.Rdata")
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
#ggsave(filename="corr_t8.png",plot=plot1.2,device="png",dpi=600,units="in",width=8,height=8)

###t13
##Relative RMSE
load("t13/r_rmse_t13.Rdata")
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
#ggsave(filename="r_rmse_t13.png",plot=plot2.1,device="png",dpi=600,units="in",width=8,height=8)

##Correlation
load("t13/corr_t13.Rdata")
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
#ggsave(filename="corr_t13.png",plot=plot2.2,device="png",dpi=600,units="in",width=8,height=8)

###t18
##Relative RMSE
load("t18/r_rmse_t18.Rdata")
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
#ggsave(filename="r_rmse_t18.png",plot=plot3.1,device="png",dpi=600,units="in",width=8,height=8)

##Correlation
load("t18/corr_t18.Rdata")
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
#ggsave(filename="corr_t18.png",plot=plot3.2,device="png",dpi=600,units="in",width=8,height=8)

###t25
##Relative RMSE
load("t25/r_rmse_t25.Rdata")
data_df <- do.call(rbind, r_rmse_list)
data_long <- pivot_longer(data_df, cols = names(data_df), names_to = "Method", values_to = "Result")
data_long$Method <- factor(data_long$Method, levels = c("Ridge", "MLRR", "MLCRR", "BAL", "BVS"))
plot4.1 <- ggplot(data_long, aes(x = Method, y = Result, fill = Method)) +
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
#ggsave(filename="r_rmse_t25.png",plot=plot4.1,device="png",dpi=600,units="in",width=8,height=8)

##Correlation
load("t25/corr_t25.Rdata")
data_df <- do.call(rbind, corr_list)
data_long <- pivot_longer(data_df, cols = names(data_df), names_to = "Method", values_to = "Result")
data_long$Method <- factor(data_long$Method, levels = c("Ridge", "MLRR", "MLCRR", "BAL", "BVS"))
plot4.2 <- ggplot(data_long, aes(x = Method, y = Result, fill = Method)) +
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
#ggsave(filename="corr_t25.png",plot=plot4.2,device="png",dpi=600,units="in",width=8,height=8)

plot_FigureS5 <- plot_grid(plot1.1,plot2.1,plot3.1,plot4.1,plot1.2,plot2.2,plot3.2,plot4.2,ncol=4, nrow=2,labels = c("A","B","C","D","E","F","G","H"),label_size = 26)
ggsave(filename="plot_FigureS5.png",plot=plot_FigureS5,device="png",dpi=600,units="in",width=30,height=14)
