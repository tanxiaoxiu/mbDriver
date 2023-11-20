library(tidyverse)
library(readxl)
library(ggpubr)
library(pheatmap)

setwd("~/mbDriver/real_data_new/UC")

###Figure8B
MAX_N_OTU =10
l_tb <- read.delim("UC_driver_summary.txt",  sep = '\t', check.names = FALSE)
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "Score")

w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "Score",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")

p1 <- pheatmap(w_tb, cluster_cols = F,angle_col=0,fontsize=14)
ggsave(filename="Figure8B.png",plot=p1,device="png",dpi=600,units="in",width=6,height=7)







