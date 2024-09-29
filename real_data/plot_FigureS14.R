library(tidyverse)
library(readxl)
library(ggpubr)
library(pheatmap)
library(ggplot2)
library(ggprism)


###FigureS14A
setwd("~/mbDriver/real_data/UC/p15")
MAX_N_OTU =15
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
ggsave(filename="FigureS14A.png",plot=p1,device="png",dpi=600,units="in",width=7,height=6)

###FigureS14B
setwd("~/mbDriver/real_data")
UC_mbDriver <- read.delim("UC/p10/UC_driver_summary.txt",  sep = '\t', check.names = FALSE)
UC_MDSINE <- read.delim("other/p10/UC/UC_driver_MDSINE.txt",  sep = '\t', check.names = FALSE)

draw_venn <- function(venn_data, group_name, label) {
  venn.plot <- venn.diagram(
    x = venn_data,
    category.names = c("mbDriver", "MDSINE"), 
    fill = c("#F59797", "#5686AD"),  
    alpha = 0.5,  
    cex = 1.5,  
    cat.cex = 1.5, 
    cat.col = c("black", "black"),  
    cat.pos = c(-30, 30),  
    cat.dist = c(0.05, 0.05),  
    lwd = 0, 
    lty = "blank",  
    filename = NULL,  
    output = TRUE,
    print.mode = "percent",
    force.unique = TRUE,
    scale = FALSE  
  )
  venn_plot_labeled <- grobTree(venn.plot, textGrob(label, x = 0.5, y = 0.88, gp = gpar(fontsize = 20, fontface = "bold")))
  return(venn_plot_labeled)
}
H1_mbDriver <- UC_mbDriver %>% filter(Group == "H1") %>% pull(Driver) %>% na.omit()
H1_MDSINE <- UC_MDSINE %>% filter(Group == "H1") %>% pull(Driver) %>% na.omit()
H2_mbDriver <- UC_mbDriver %>% filter(Group == "H2") %>% pull(Driver) %>% na.omit()
H2_MDSINE <- UC_MDSINE %>% filter(Group == "H2") %>% pull(Driver) %>% na.omit()
UC1_mbDriver <- UC_mbDriver %>% filter(Group == "UC1") %>% pull(Driver) %>% na.omit()
UC1_MDSINE <- UC_MDSINE %>% filter(Group == "UC1") %>% pull(Driver) %>% na.omit()
UC2_mbDriver <- UC_mbDriver %>% filter(Group == "UC2") %>% pull(Driver) %>% na.omit()
UC2_MDSINE <- UC_MDSINE %>% filter(Group == "UC2") %>% pull(Driver) %>% na.omit()
venn_data_H1 <- list(UC_mbDriver = H1_mbDriver, UC_MDSINE = H1_MDSINE)
venn_data_H2 <- list(UC_mbDriver = H2_mbDriver, UC_MDSINE = H2_MDSINE)
venn_data_UC1 <- list(UC_mbDriver = UC1_mbDriver, UC_MDSINE = UC1_MDSINE)
venn_data_UC2 <- list(UC_mbDriver = UC2_mbDriver, UC_MDSINE = UC2_MDSINE)
venn_H1 <- draw_venn(venn_data_H1, "H1", "H1")
venn_H2 <- draw_venn(venn_data_H2, "H2", "H2")
venn_UC1 <- draw_venn(venn_data_UC1, "UC1", "UC1")
venn_UC2 <- draw_venn(venn_data_UC2, "UC2", "UC2")
png("FigureS14B.png", width = 16, height = 4, units = "in", res = 600)
grid.arrange(venn_H1, venn_H2, venn_UC1,venn_UC2,ncol = 4)
dev.off()

###FigureS14C
setwd("~/mbDriver/real_data")
UC_mbDriver <- read.delim("UC/p15/UC_driver_summary.txt",  sep = '\t', check.names = FALSE)
UC_MDSINE <- read.delim("other/p15/UC/UC_driver_MDSINE.txt",  sep = '\t', check.names = FALSE)

draw_venn <- function(venn_data, group_name, label) {
  venn.plot <- venn.diagram(
    x = venn_data,
    category.names = c("mbDriver", "MDSINE"),  
    fill = c("#F59797", "#5686AD"),  
    alpha = 0.5,  
    cex = 1.5,  
    cat.cex = 1.5,  
    cat.col = c("black", "black"),  
    cat.pos = c(-30, 30),  
    cat.dist = c(0.05, 0.05),  
    lwd = 0,  
    lty = "blank",  
    filename = NULL,  
    output = TRUE,
    print.mode = "percent",
    force.unique = TRUE,
    scale = FALSE  
  )
  venn_plot_labeled <- grobTree(venn.plot, textGrob(label, x = 0.5, y = 0.88, gp = gpar(fontsize = 20, fontface = "bold")))
  return(venn_plot_labeled)
}

H1_mbDriver <- UC_mbDriver %>% filter(Group == "H1") %>% pull(Driver) %>% na.omit()
H1_MDSINE <- UC_MDSINE %>% filter(Group == "H1") %>% pull(Driver) %>% na.omit()
H2_mbDriver <- UC_mbDriver %>% filter(Group == "H2") %>% pull(Driver) %>% na.omit()
H2_MDSINE <- UC_MDSINE %>% filter(Group == "H2") %>% pull(Driver) %>% na.omit()
UC1_mbDriver <- UC_mbDriver %>% filter(Group == "UC1") %>% pull(Driver) %>% na.omit()
UC1_MDSINE <- UC_MDSINE %>% filter(Group == "UC1") %>% pull(Driver) %>% na.omit()
UC2_mbDriver <- UC_mbDriver %>% filter(Group == "UC2") %>% pull(Driver) %>% na.omit()
UC2_MDSINE <- UC_MDSINE %>% filter(Group == "UC2") %>% pull(Driver) %>% na.omit()

venn_data_H1 <- list(UC_mbDriver = H1_mbDriver, UC_MDSINE = H1_MDSINE)
venn_data_H2 <- list(UC_mbDriver = H2_mbDriver, UC_MDSINE = H2_MDSINE)
venn_data_UC1 <- list(UC_mbDriver = UC1_mbDriver, UC_MDSINE = UC1_MDSINE)
venn_data_UC2 <- list(UC_mbDriver = UC2_mbDriver, UC_MDSINE = UC2_MDSINE)

venn_H1 <- draw_venn(venn_data_H1, "H1", "H1")
venn_H2 <- draw_venn(venn_data_H2, "H2", "H2")
venn_UC1 <- draw_venn(venn_data_UC1, "UC1", "UC1")
venn_UC2 <- draw_venn(venn_data_UC2, "UC2", "UC2")

png("FigureS14C.png", width = 16, height = 4, units = "in", res = 600)
grid.arrange(venn_H1, venn_H2, venn_UC1,venn_UC2,ncol = 4)
dev.off()
