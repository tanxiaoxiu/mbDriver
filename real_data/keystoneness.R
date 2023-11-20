rm(list = ls(all = TRUE)) 
library(MASS)
library(glmnet)
library(Rfast)
library(mgcv)

###D1_Score:keystoneness
D1_Score <- function(A,r){
  x_star <- -solve(A) %*% r
  x_star[x_star < 0] <- 0
  p <- length(r)
  D_square <- as.data.frame(matrix(nrow=p,ncol=1))
  for (i in 1:p){
    Ai <- A[-i,-i]
    ri <- r[-i]
    z_star <- -solve(Ai) %*% ri
    z_star[z_star < 0] <- 0
    di <- x_star[-i] - z_star
    D_square[i,] <- sum(as.numeric(di*di))
  }
  rownames(D_square) <- row.names(A)
  D_square <- cbind(rownames(D_square),D_square)
  D_square_sort <-D_square[order(-D_square[,2]),]
  return(D_square_sort)
}

### write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

#######################Fiber_diet case
setwd("~/mbDriver/real_data/Fiber_diet")
####Con group
A_ridge <- read.delim("A_Con.txt",  
                      row.names = 1, sep = '\t', check.names = FALSE)
r_ridge <- read.delim("r_Con.txt",
                      row.names = 1, sep = '\t', check.names = FALSE)
A_ridge <- as.matrix(A_ridge)
r_ridge <- as.matrix(r_ridge)
##Driver
Con_D1_Score <- D1_Score(A_ridge,r_ridge)
Con_D1_driver <- Con_D1_Score[c(1:5),]
colnames(Con_D1_driver) <- c("Driver","Score")
Con_D1_driver$Group <- "Con"

#####Rs group
A_ridge <- read.delim("A_Rs.txt",  
                      row.names = 1, sep = '\t', check.names = FALSE)
r_ridge <- read.delim("r_Rs.txt",
                      row.names = 1, sep = '\t', check.names = FALSE)
A_ridge <- as.matrix(A_ridge)
r_ridge <- as.matrix(r_ridge)
##Driver
Rs_D1_Score <- D1_Score(A_ridge,r_ridge)
Rs_D1_driver <- Rs_D1_Score[c(1:5),]
colnames(Rs_D1_driver) <- c("Driver","Score")
Rs_D1_driver$Group <- "Rs"

#####In group
A_ridge <- read.delim("A_In.txt",  
                      row.names = 1, sep = '\t', check.names = FALSE)
r_ridge <- read.delim("r_In.txt",
                      row.names = 1, sep = '\t', check.names = FALSE)
A_ridge <- as.matrix(A_ridge)
r_ridge <- as.matrix(r_ridge)
##Driver
In_D1_Score <- D1_Score(A_ridge,r_ridge)
In_D1_driver <- In_D1_Score[c(1:5),]
colnames(In_D1_driver) <- c("Driver","Score")
In_D1_driver$Group <- "In"
summary <- rbind(Con_D1_driver,Rs_D1_driver,In_D1_driver)

#Supplementary Table 2
write.table(summary,file ="keystoneness_summary.txt",row.names = F,col.names = T, sep = "\t",quote = F)

#######################UC case
setwd("~/mbDriver/real_data/UC")
####H1 group
A_ridge <- read.delim("A_H1.txt",  
                      row.names = 1, sep = '\t', check.names = FALSE)
r_ridge <- read.delim("r_H1.txt",
                      row.names = 1, sep = '\t', check.names = FALSE)
A_ridge <- as.matrix(A_ridge)
r_ridge <- as.matrix(r_ridge)
##Driver
H1_D1_Score <- D1_Score(A_ridge,r_ridge)
H1_D1_driver <- H1_D1_Score[c(1:5),]
colnames(H1_D1_driver) <- c("Driver","Score")
H1_D1_driver$Group <- "H1"

####H2 group
A_ridge <- read.delim("A_H2.txt",  
                      row.names = 1, sep = '\t', check.names = FALSE)
r_ridge <- read.delim("r_H2.txt",
                      row.names = 1, sep = '\t', check.names = FALSE)
A_ridge <- as.matrix(A_ridge)
r_ridge <- as.matrix(r_ridge)
##Driver
H2_D1_Score <- D1_Score(A_ridge,r_ridge)
H2_D1_driver <- H2_D1_Score[c(1:5),]
colnames(H2_D1_driver) <- c("Driver","Score")
H2_D1_driver$Group <- "H2"

####UC1 group
A_ridge <- read.delim("A_UC2.txt",  
                      row.names = 1, sep = '\t', check.names = FALSE)
r_ridge <- read.delim("r_UC2.txt",
                      row.names = 1, sep = '\t', check.names = FALSE)
A_ridge <- as.matrix(A_ridge)
r_ridge <- as.matrix(r_ridge)
##Driver
UC1_D1_Score <- D1_Score(A_ridge,r_ridge)
UC1_D1_driver <- UC1_D1_Score[c(1:5),]
colnames(UC1_D1_driver) <- c("Driver","Score")
UC1_D1_driver$Group <- "H2"

#####UC2 group
A_ridge <- read.delim("A_UC2.txt",  
                      row.names = 1, sep = '\t', check.names = FALSE)
r_ridge <- read.delim("r_UC2.txt",
                      row.names = 1, sep = '\t', check.names = FALSE)
A_ridge <- as.matrix(A_ridge)
r_ridge <- as.matrix(r_ridge)
##Driver
UC2_D1_Score <- D1_Score(A_ridge,r_ridge)
UC2_D1_driver <- UC2_D1_Score[c(1:5),]
colnames(UC2_D1_driver) <- c("Driver","Score")
UC2_D1_driver$Group <- "H2"
summary <- rbind(H1_D1_driver,H2_D1_driver,UC1_D1_driver,UC2_D1_driver)
#Supplementary Table 4
write.table(summary,file ="keystoneness_summary.txt",row.names = F,col.names = T, sep = "\t",quote = F)

