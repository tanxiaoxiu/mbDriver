setwd("~/mbDriver/simulation/toy_example")
library(tidyverse)
library(readxl)
library(ggpubr)
library(pheatmap)
library(igraph)

####toy example
###D1: D1_Score
D1_Score <- function(A,r){
  x_star <- -solve(A) %*% r
  p <- length(r)
  D_square <- as.data.frame(matrix(nrow=p,ncol=1))
  for (i in 1:p){
    Ai <- A
    Ai[-i,i] <- 0
    Ai[i,-i] <- 0
    ri <- r
    z_star <- -solve(Ai) %*% ri
    di <- x_star - z_star
    D_square[i,] <- sum(as.numeric(di*di))
  }
  rownames(D_square) <- row.names(A)
  D_square <- cbind(rownames(D_square),D_square)
  return(D_square)
}


###D2: D2_Score
D2_Score <- function(A,r){
  x_star <- -solve(A) %*% r
  p <- length(r)
  D_square <- as.data.frame(matrix(nrow=p,ncol=1))
  for (i in 1:p){
    Ai <- A
    Ai[i,-i] <- 0
    ri <- r
    z_star <- -solve(Ai) %*% ri
    di <- x_star - z_star
    D_square[i,] <- sum(as.numeric(di*di))
  }
  rownames(D_square) <- row.names(A)
  D_square <- cbind(rownames(D_square),D_square)
  return(D_square)
}

###D3: D3_Score 
D3_Score <- function(A,r){
  x_star <- -solve(A) %*% r
  p <- length(r)
  D_square <- as.data.frame(matrix(nrow=p,ncol=1))
  for (i in 1:p){
    Ai <- A
    Ai[-i,i] <- 0
    ri <- r
    z_star <- -solve(Ai) %*% ri
    di <- x_star - z_star
    D_square[i,] <- sum(as.numeric(di*di))
  }
  rownames(D_square) <- row.names(A)
  D_square <- cbind(rownames(D_square),D_square)
  return(D_square)
}

### write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

#############Fig3A
#############Causal graph structure:fork
p = 3 #number of species
g <- graph(edges=c("A", "B", "A", "C"), directed=TRUE)
adj_matrix <- as_adjacency_matrix(g, sparse=FALSE)
A0 <- t(adj_matrix)
diag(A0) <- rep(-1,times=p)
### generate A,r,x0
r0 <- rep(1,times=p)
r0 <- as.matrix(r0)
rownames(r0) <- c("A","B","C")
colname2 <- "r"

#Index:D1,D2,D3
D1_fork <- D1_Score(A0,r0)
D2_fork <- D2_Score(A0,r0)
D3_fork <- D3_Score(A0,r0)
compare_D_fork <- rbind(D1_fork,D2_fork,D3_fork)
compare_D_fork$method <- "D"
compare_D_fork$Group <- c(rep("D1", p), rep("D2", p), rep("D3", p))
colnames(compare_D_fork) <- c("Driver","D","Method","Group")
write.table(compare_D_fork,file ="compare_D_fork.txt",row.names = F,col.names = T, sep = "\t",quote = F)
######Figure3A left
MAX_N_OTU =p
l_tb <-  compare_D_fork
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "D")
w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "D",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")

p1 <- pheatmap(w_tb, cluster_cols = F,cluster_rows = F,angle_col=0,fontsize=14)
ggsave(filename="p1.png",plot=p1,device="png",dpi=600,units="in",width=4,height=4)


#############Causal graph structure:collider
p = 3 #number of species
### generate A,r,x0
edges <- c("B", "A", "C", "A")
g <- graph(edges, directed = TRUE)
adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)
adj_matrix <- adj_matrix[c(2, 1, 3), c(2, 1, 3)]
A0 <- t(adj_matrix)
diag(A0) <- rep(-1,times=p)

r0 <- rep(1,times=p)
r0 <- as.matrix(r0)
rownames(r0) <- c("A","B","C")
colname2 <- "r"

#Index:D1,D2,D3
D1_collider <- D1_Score(A0,r0)
D2_collider <- D2_Score(A0,r0)
D3_collider <- D3_Score(A0,r0)
compare_D_collider <- rbind(D1_collider,D2_collider,D3_collider)
compare_D_collider$method <- "D"
compare_D_collider$Group <- c(rep("D1", p), rep("D2", p), rep("D3", p))
colnames(compare_D_collider) <- c("Driver","D","Method","Group")
write.table(compare_D_collider,file ="compare_D_collider.txt",row.names = F,col.names = T, sep = "\t",quote = F)
######Figure3A center 
MAX_N_OTU = p
l_tb <-  compare_D_collider
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "D")
w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "D",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")

p2 <- pheatmap(w_tb, cluster_cols = F,cluster_rows = F,angle_col=0,fontsize=14)
ggsave(filename="p2.png",plot=p2,device="png",dpi=600,units="in",width=4,height=4)


#############Causal graph structure:chain
p = 3 #number of species
### generate A,r,x0
edges <- c("A", "B", "B", "C")
g <- graph(edges, directed = TRUE)
adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)
A0 <- t(adj_matrix)
diag(A0) <- rep(-1,times=p)

r0 <- rep(1,times=p)
r0 <- as.matrix(r0)
rownames(r0) <- c("A","B","C")
colname2 <- "r"

#Index:D1,D2,D3
D1_chain <- D1_Score(A0,r0)
D2_chain <- D2_Score(A0,r0)
D3_chain <- D3_Score(A0,r0)
compare_D_chain <- rbind(D1_chain,D2_chain,D3_chain)
compare_D_chain$method <- "D"
compare_D_chain$Group <- c(rep("D1", p), rep("D2", p), rep("D3", p))
colnames(compare_D_chain) <- c("Driver","D","Method","Group")
write.table(compare_D_chain,file ="compare_D_chain.txt",row.names = F,col.names = T, sep = "\t",quote = F)
######Figure3A right
MAX_N_OTU =p
l_tb <-  compare_D_chain
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "D")
w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "D",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")

p3 <- pheatmap(w_tb, cluster_cols = F,cluster_rows = F,angle_col=0,fontsize=14)
ggsave(filename="p3.png",plot=p3,device="png",dpi=600,units="in",width=4,height=4)


#############Figure3B
#############Causal graph structure:collider
p = 3 #number of species
### generate A,r,x0
g <- graph(edges=c("A", "B", "A", "C"), directed=TRUE)
adj_matrix <- as_adjacency_matrix(g, sparse=FALSE)
adj_matrix[adj_matrix != 0] <- c(2,3)
A0 <- t(adj_matrix)
diag(A0) <- rep(-1,times=p)
r0 <- rep(1,times=p)
r0 <- as.matrix(r0)
rownames(r0) <- c("A","B","C")
colname2 <- "r"

#Index:D1,D2,D3
D1_fork <- D1_Score(A0,r0)
D2_fork <- D2_Score(A0,r0)
D3_fork <- D3_Score(A0,r0)
compare_D_fork <- rbind(D1_fork,D2_fork,D3_fork)
compare_D_fork$method <- "D"
compare_D_fork$Group <- c(rep("D1", p), rep("D2", p), rep("D3", p))
colnames(compare_D_fork) <- c("Driver","D","Method","Group")
write.table(compare_D_fork,file ="compare_D_fork_interaction.txt",row.names = F,col.names = T, sep = "\t",quote = F)
######Figure3B left
MAX_N_OTU =3
l_tb <-  compare_D_fork
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "D")
w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "D",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")

p4 <- pheatmap(w_tb, cluster_cols = F,cluster_rows = F,angle_col=0,fontsize=14)
ggsave(filename="p4.png",plot=p4,device="png",dpi=600,units="in",width=4,height=4)


#############Causal graph structure:fork
p = 3 #number of species
### generate A,r,x0
edges <- c("B", "A", "C", "A")
g <- graph(edges, directed = TRUE)
adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)
adj_matrix <- adj_matrix[c(2, 1, 3), c(2, 1, 3)]
adj_matrix[adj_matrix != 0] <- c(2,3)
A0 <- t(adj_matrix)
diag(A0) <- rep(-1,times=p)
r0 <- rep(1,times=p)
r0 <- as.matrix(r0)
rownames(r0) <- c("A","B","C")
colname2 <- "r"

#Index:D1,D2,D3
D1_collider <- D1_Score(A0,r0)
D2_collider <- D2_Score(A0,r0)
D3_collider <- D3_Score(A0,r0)
compare_D_collider <- rbind(D1_collider,D2_collider,D3_collider)
compare_D_collider$method <- "D"
compare_D_collider$Group <- c(rep("D1", p), rep("D2", p), rep("D3", p))
colnames(compare_D_collider) <- c("Driver","D","Method","Group")
write.table(compare_D_collider,file ="compare_D_collider_interaction.txt",row.names = F,col.names = T, sep = "\t",quote = F)
######Figure3B center 
MAX_N_OTU = p
l_tb <-  compare_D_collider
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "D")
w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "D",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")
p5 <- pheatmap(w_tb, cluster_cols = F,cluster_rows = F,angle_col=0,fontsize=14)
ggsave(filename="p5.png",plot=p5,device="png",dpi=600,units="in",width=4,height=4)


#############Causal graph structure:chain
p = 3 #number of species
edges <- c("A", "B", "B", "C")
g <- graph(edges, directed = TRUE)
adj_matrix <- as_adjacency_matrix(g, sparse = FALSE)
adj_matrix[adj_matrix != 0] <- c(2,3)
A0 <- t(adj_matrix)
diag(A0) <- rep(-1,times=p)

r0 <- rep(1,times=p)
r0 <- as.matrix(r0)
rownames(r0) <- c("A","B","C")
colname2 <- "r"

#Index:D1,D2,D3
D1_chain <- D1_Score(A0,r0)
D2_chain <- D2_Score(A0,r0)
D3_chain <- D3_Score(A0,r0)
compare_D_chain <- rbind(D1_chain,D2_chain,D3_chain)
compare_D_chain$method <- "D"
compare_D_chain$Group <- c(rep("D1", p), rep("D2", p), rep("D3", p))
colnames(compare_D_chain) <- c("Driver","D","Method","Group")
write.table(compare_D_chain,file ="compare_D_chain_interaction.txt",row.names = F,col.names = T, sep = "\t",quote = F)
######Figure3B right
MAX_N_OTU = p
l_tb <-  compare_D_chain
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "D")
w_tb <- l_tb %>% 
  pivot_wider(id_cols = "Driver", 
              names_from = "Group", 
              values_from = "D",
              values_fill=0) %>% 
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")
p6 <- pheatmap(w_tb, cluster_cols = F,cluster_rows = F,angle_col=0,fontsize=14)
ggsave(filename="p6.png",plot=p6,device="png",dpi=600,units="in",width=4,height=4)
