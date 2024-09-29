library(tidyverse)
library(readxl)
library(ggpubr)
library(pheatmap)
library(ggplot2)
library(ggprism)

###FigureS11A
setwd("~/mbDriver/real_data/Fiber_diet/p15")
MAX_N_OTU =15
l_tb <- read.delim("Fiber_driver_summary.txt",  sep = '\t', check.names = FALSE)
l_tb$Driver <- gsub("Parabacteroides-goldsteinii", "P. goldsteinii", l_tb$Driver)
l_tb$Driver <- gsub("Lachnospiraceae-NK4A136-group", "L. NK4A136-group", l_tb$Driver)
l_tb$Driver <- gsub("Bacteroides-acidifaciens", "B. acidifaciens", l_tb$Driver)
l_tb$Driver <- gsub("Akkermansia-muciniphila", "A. muciniphila", l_tb$Driver)

l_tb$Driver <- gsub("Lachnospiraceae-bacterium-28-4", "L. bacterium-28-4", l_tb$Driver)
l_tb$Driver <- gsub("Rikenellaceae-RC9-gut-group", "Rikenellaceae", l_tb$Driver)

l_tb %>% pivot_wider(names_from = "Group" ,values_from = "Score")
w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "Score",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")
p1 <- pheatmap(w_tb, cluster_cols = F,angle_col=0,fontsize=14)
ggsave(filename="FigureS11A.png",plot=p1,device="png",dpi=600,units="in",width=6,height=6)


###FigureS11B
setwd("~/mbDriver/real_data/Fiber_diet/p20")
MAX_N_OTU =20
l_tb <- read.delim("Fiber_driver_summary.txt",  sep = '\t', check.names = FALSE)
l_tb$Driver <- gsub("Parabacteroides-goldsteinii", "P. goldsteinii", l_tb$Driver)
l_tb$Driver <- gsub("Lachnospiraceae-NK4A136-group", "L. NK4A136-group", l_tb$Driver)
l_tb$Driver <- gsub("Bacteroides-acidifaciens", "B. acidifaciens", l_tb$Driver)
l_tb$Driver <- gsub("Akkermansia-muciniphila", "A. muciniphila", l_tb$Driver)

l_tb$Driver <- gsub("Lachnospiraceae-bacterium-28-4", "L. bacterium-28-4", l_tb$Driver)
l_tb$Driver <- gsub("Rikenellaceae-RC9-gut-group", "Rikenellaceae", l_tb$Driver)

l_tb %>% pivot_wider(names_from = "Group" ,values_from = "Score")
w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "Score",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")
p1 <- pheatmap(w_tb, cluster_cols = F,angle_col=0,fontsize=14)
ggsave(filename="FigureS11B.png",plot=p1,device="png",dpi=600,units="in",width=6,height=6)


###FigureS11C
setwd("~/mbDriver/real_data")
###Fiber
Fiber_mbDriver <- read.delim("Fiber_diet/p10/Fiber_driver_summary.txt",  sep = '\t', check.names = FALSE)
Fiber_MDSINE <- read.delim("other/p10/Fiber/Fiber_driver_MDSINE.txt",  sep = '\t', check.names = FALSE)

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

con_mbDriver <- Fiber_mbDriver %>% filter(Group == "Con") %>% pull(Driver) %>% na.omit()
con_MDSINE <- Fiber_MDSINE %>% filter(Group == "Con") %>% pull(Driver) %>% na.omit()
rs_mbDriver <- Fiber_mbDriver %>% filter(Group == "Rs") %>% pull(Driver) %>% na.omit()
rs_MDSINE <- Fiber_MDSINE %>% filter(Group == "Rs") %>% pull(Driver) %>% na.omit()
in_mbDriver <- Fiber_mbDriver %>% filter(Group == "In") %>% pull(Driver) %>% na.omit()
in_MDSINE <- Fiber_MDSINE %>% filter(Group == "In") %>% pull(Driver) %>% na.omit()
venn_data_con <- list(Fiber_mbDriver = con_mbDriver, Fiber_MDSINE = con_MDSINE)
venn_data_rs <- list(Fiber_mbDriver = rs_mbDriver, Fiber_MDSINE = rs_MDSINE)
venn_data_in <- list(Fiber_mbDriver = in_mbDriver, Fiber_MDSINE = in_MDSINE)
venn_con <- draw_venn(venn_data_con, "Con", "Con")
venn_rs <- draw_venn(venn_data_rs, "Rs", "Rs")
venn_in <- draw_venn(venn_data_in, "In", "In")
png("FigureS11C.png", width = 15, height = 5, units = "in", res = 600)
grid.arrange(venn_con, venn_rs, venn_in, ncol = 3)
dev.off()


###FigureS11D
setwd("~/mbDriver/real_data")
###Fiber
Fiber_mbDriver <- read.delim("Fiber_diet/p15/Fiber_driver_summary.txt",  sep = '\t', check.names = FALSE)
Fiber_MDSINE <- read.delim("other/p15/Fiber/Fiber_driver_MDSINE.txt",  sep = '\t', check.names = FALSE)

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

con_mbDriver <- Fiber_mbDriver %>% filter(Group == "Con") %>% pull(Driver) %>% na.omit()
con_MDSINE <- Fiber_MDSINE %>% filter(Group == "Con") %>% pull(Driver) %>% na.omit()
rs_mbDriver <- Fiber_mbDriver %>% filter(Group == "Rs") %>% pull(Driver) %>% na.omit()
rs_MDSINE <- Fiber_MDSINE %>% filter(Group == "Rs") %>% pull(Driver) %>% na.omit()
in_mbDriver <- Fiber_mbDriver %>% filter(Group == "In") %>% pull(Driver) %>% na.omit()
in_MDSINE <- Fiber_MDSINE %>% filter(Group == "In") %>% pull(Driver) %>% na.omit()
venn_data_con <- list(Fiber_mbDriver = con_mbDriver, Fiber_MDSINE = con_MDSINE)
venn_data_rs <- list(Fiber_mbDriver = rs_mbDriver, Fiber_MDSINE = rs_MDSINE)
venn_data_in <- list(Fiber_mbDriver = in_mbDriver, Fiber_MDSINE = in_MDSINE)
venn_con <- draw_venn(venn_data_con, "Con", "Con")
venn_rs <- draw_venn(venn_data_rs, "Rs", "Rs")
venn_in <- draw_venn(venn_data_in, "In", "In")
png("FigureS11D.png", width = 15, height = 5, units = "in", res = 600)
grid.arrange(venn_con, venn_rs, venn_in, ncol = 3)
dev.off()




