rm(list = ls(all = TRUE)) 
library(MASS)
library(glmnet)
library(Rfast)
library(mgcv)

###calculate f
f <- function(X){ 
  Y<-c()
  X <- as.data.frame(X)
  fmla <- formula("X~s(t, m=2, k=5,bs = 'cr')")
  dat <- cbind(X,t)
  mod_nb <- gam(fmla, data = dat, family = nb(), method="REML")
  Y <- fitted.values(mod_nb)
  return(Y)
} 

###calculate f deriv
f_deriv <- function(X){
  Y<-c()
  X <- as.data.frame(X)
  fmla <- formula("X~s(t, m=2, k=5,bs = 'cr')")
  dat <- cbind(X,t)
  mod_nb <- gam(fmla, data = dat, family = nb(), method="REML")
  fitted.values <- fitted.values(mod_nb)
  newdat <- as.data.frame(t)
  X0 <- predict(mod_nb,newdata = newdat,type="lpmatrix")
  eps <- 1e-7 
  t <- t + eps 
  newdat <- as.data.frame(t)
  X1 <- predict(mod_nb,newdata = newdat,type="lpmatrix")
  Xp <- (X1-X0)/eps
  df <- Xp %*% coef(mod_nb)
  Y <- fitted.values * df
  return(Y)
}

###calculate Driver Score
Driver_Score <- function(A,r){
  x_star <- -solve(A) %*% r
  x_star[x_star < 0] <- 0
  p <- length(r)
  D_square <- as.data.frame(matrix(nrow=p,ncol=1))
  for (i in 1:p){
    Ai <- A
    Ai[i,-i] <- 0
    ri <- r
    z_star <- -solve(Ai) %*% ri
    z_star[z_star < 0] <- 0
    di <- x_star - z_star
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

setwd("~/mbDriver/real_data_new/UC")

#Input data
abundance <- read.delim("absolute_abundance.txt",row.names = 1, sep = '\t', check.names = FALSE)
metadata <- read.delim("metadata.txt", row.names = 1, sep = '\t', check.names = FALSE)

#preprocessing
data_abundance <- rbind(abundance,Total=colSums(abundance))
data_abundance_sort <- data_abundance[,order(-as.numeric(data_abundance["Total",]))]
data_UC <- merge(metadata,data_abundance_sort,by = 'row.names', all = F)
write.table(data_UC,file = "data_UC.txt",row.names = F,col.names = T, sep = "\t",quote = F)

p = 10 
n1 = 5
n2 = 10
group1 = "H1"
group2 = "H2"
group3 = "UC1"
group4 = "UC2"

#####################H1 group
###Data preprocessing
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
for (i in 2:n1){
  subject_i <- data_UC[data_UC$Subject==i,]
  subject_i <- subject_i[subject_i$Group == group1,]
  subject_i <- subject_i[,-c(1,2,4)]
  subject_i_top <- subject_i[,c(1:(p+1))]
  otu_table <- subject_i_top[order(subject_i_top[,"Time"]),] 
  t <- otu_table[,1]
  s_f <- apply(otu_table[,-1],2,f)
  s_f_deriv <- apply(otu_table[,-1],2,f_deriv)
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
}

F_d <- as.matrix(F_d)
F_dt <- as.matrix(F_dt)
F_F_dt <- F_dt/F_d

###Parameter estimation
X <- as.matrix(F_d)
Y <- as.matrix(F_F_dt)
pe_r <- as.data.frame(matrix(nrow=(p+1)))
for (i in 1:p){
  set.seed(12345)
  Yr <- as.matrix(Y[,i])
  ridge_fit <- glmnet(X, Yr,family = "gaussian",alpha=0)
  r_cv_fit <- cv.glmnet(X, Yr,family = "gaussian",type.measure="mse",alpha=0)
  ridge_best <- glmnet(X, Yr,family = "gaussian",lambda = r_cv_fit$lambda.min,alpha=0)
  coefficients <- coef(ridge_best,r_cv_fit$lambda.min)
  pe_r[,i] <- as.matrix(coefficients)
}
pe_r <- t(pe_r)
r_ridge <- as.matrix(pe_r[,1])
A_ridge <- pe_r[,2:(p+1)]
colnames(A_ridge) <- colnames(X)
row.names(A_ridge) <- colnames(A_ridge)
A_ridge_write <- adjustdata(A_ridge)
filename_A <- paste0("A_",group1,".txt")
write.table(A_ridge_write,file = filename_A,row.names = F,col.names = T, sep = "\t",quote = F)
filename_r <- paste0("r_",group1,".txt")
write.table(r_ridge,file = filename_r,row.names = T,col.names = T, sep = "\t",quote = F)

###Driver prediction
group1_Driver_Score <- Driver_Score(A_ridge,r_ridge)
filename_Score <- paste0(group1,"_Driver_Score.txt")
write.table(group1_Driver_Score,file = filename_Score,row.names = F,col.names = T, sep = "\t",quote = F)
group1_driver <- group1_Driver_Score[c(1:5),]
colnames(group1_driver) <- c("Driver","Score")
group1_driver$Group <- group1
filename_driver <- paste0(group1,"_driver.txt")
write.table(group1_driver,file = filename_driver,row.names = F,col.names = T, sep = "\t",quote = F)

#####################H2 group
###Data preprocessing
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
for (i in 2:n1){
  subject_i <- data_UC[data_UC$Subject==i,]
  subject_i <- subject_i[subject_i$Group == group2,]
  subject_i <- subject_i[,-c(1,2,4)]
  subject_i_top <- subject_i[,c(1:(p+1))]
  otu_table <- subject_i_top[order(subject_i_top[,"Time"]),] 
  t <- otu_table[,1]
  s_f <- apply(otu_table[,-1],2,f)
  s_f_deriv <- apply(otu_table[,-1],2,f_deriv)
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
}

F_d <- as.matrix(F_d)
F_dt <- as.matrix(F_dt)
F_F_dt <- F_dt/F_d


###Parameter estimation
X <- as.matrix(F_d)
Y <- as.matrix(F_F_dt)
pe_r <- as.data.frame(matrix(nrow=(p+1)))
for (i in 1:p){
  set.seed(12345)
  Yr <- as.matrix(Y[,i])
  ridge_fit <- glmnet(X, Yr,family = "gaussian",alpha=0)
  r_cv_fit <- cv.glmnet(X, Yr,family = "gaussian",type.measure="mse",alpha=0)
  ridge_best <- glmnet(X, Yr,family = "gaussian",lambda = r_cv_fit$lambda.min,alpha=0)
  coefficients <- coef(ridge_best,r_cv_fit$lambda.min)
  pe_r[,i] <- as.matrix(coefficients)
}
pe_r <- t(pe_r)
r_ridge <- as.matrix(pe_r[,1])
A_ridge <- pe_r[,2:(p+1)]
colnames(A_ridge) <- colnames(X)
row.names(A_ridge) <- colnames(A_ridge)
A_ridge_write <- adjustdata(A_ridge)
filename_A <- paste0("A_",group2,".txt")
write.table(A_ridge_write,file = filename_A,row.names = F,col.names = T, sep = "\t",quote = F)
filename_r <- paste0("r_",group2,".txt")
write.table(r_ridge,file = filename_r,row.names = T,col.names = T, sep = "\t",quote = F)

###Driver prediction
group2_Driver_Score <- Driver_Score(A_ridge,r_ridge)
filename_Score <- paste0(group2,"_Driver_Score.txt")
write.table(group2_Driver_Score,file = filename_Score,row.names = F,col.names = T, sep = "\t",quote = F)

group2_driver <- group2_Driver_Score[c(1:5),]
colnames(group2_driver) <- c("Driver","Score")
group2_driver$Group <- group2
filename_driver <- paste0(group2,"_driver.txt")
write.table(group2_driver,file = filename_driver,row.names = F,col.names = T, sep = "\t",quote = F)


#####################UC1 group
data_UC <- read.delim("data_UC.txt",  sep = '\t', check.names = FALSE)

###Data preprocessing
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
for (i in 6:n2){
  subject_i <- data_UC[data_UC$Subject==i,]
  subject_i <- subject_i[subject_i$Group == group3,]
  subject_i <- subject_i[,-c(1,2,4)]
  subject_i_top <- subject_i[,c(1:(p+1))]
  otu_table <- subject_i_top[order(subject_i_top[,"Time"]),] 
  t <- otu_table[,1]
  s_f <- apply(otu_table[,-1],2,f)
  s_f_deriv <- apply(otu_table[,-1],2,f_deriv)
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
}

F_d <- as.matrix(F_d)
F_dt <- as.matrix(F_dt)
F_F_dt <- F_dt/F_d

###Parameter estimation
X <- as.matrix(F_d)
Y <- as.matrix(F_F_dt)
pe_r <- as.data.frame(matrix(nrow=(p+1)))
for (i in 1:p){
  set.seed(12345)
  Yr <- as.matrix(Y[,i])
  ridge_fit <- glmnet(X, Yr,family = "gaussian",alpha=0)
  r_cv_fit <- cv.glmnet(X, Yr,family = "gaussian",type.measure="mse",alpha=0)
  ridge_best <- glmnet(X, Yr,family = "gaussian",lambda = r_cv_fit$lambda.min,alpha=0)
  coefficients <- coef(ridge_best,r_cv_fit$lambda.min)
  pe_r[,i] <- as.matrix(coefficients)
}
pe_r <- t(pe_r)
r_ridge <- as.matrix(pe_r[,1])
A_ridge <- pe_r[,2:(p+1)]
colnames(A_ridge) <- colnames(X)
row.names(A_ridge) <- colnames(A_ridge)
A_ridge_write <- adjustdata(A_ridge)
filename_A <- paste0("A_",group3,".txt")
write.table(A_ridge_write,file = filename_A,row.names = F,col.names = T, sep = "\t",quote = F)
filename_r <- paste0("r_",group3,".txt")
write.table(r_ridge,file = filename_r,row.names = T,col.names = T, sep = "\t",quote = F)

###Driver prediction
group3_Driver_Score <- Driver_Score(A_ridge,r_ridge)
filename_Score <- paste0(group3,"_Driver_Score.txt")
write.table(group3_Driver_Score,file = filename_Score,row.names = F,col.names = T, sep = "\t",quote = F)
group3_driver <- group3_Driver_Score[c(1:5),]
colnames(group3_driver) <- c("Driver","Score")
group3_driver$Group <- group3
filename_driver <- paste0(group3,"_driver.txt")
write.table(group3_driver,file = filename_driver,row.names = F,col.names = T, sep = "\t",quote = F)

#####################UC2 group
data_UC <- read.delim("data_UC.txt",  sep = '\t', check.names = FALSE)

###Data preprocessing
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
for (i in 6:n2){
  subject_i <- data_UC[data_UC$Subject==i,]
  subject_i <- subject_i[subject_i$Group == group4,]
  subject_i <- subject_i[,-c(1,2,4)]
  subject_i_top <- subject_i[,c(1:(p+1))]
  otu_table <- subject_i_top[order(subject_i_top[,"Time"]),]
  t <- otu_table[,1]
  s_f <- apply(otu_table[,-1],2,f)
  s_f_deriv <- apply(otu_table[,-1],2,f_deriv)
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
}
F_d <- as.matrix(F_d)
F_dt <- as.matrix(F_dt)
F_F_dt <- F_dt/F_d

###Parameter estimation
X <- as.matrix(F_d)
Y <- as.matrix(F_F_dt)
pe_r <- as.data.frame(matrix(nrow=(p+1)))
for (i in 1:p){
  set.seed(12345)
  Yr <- as.matrix(Y[,i])
  ridge_fit <- glmnet(X, Yr,family = "gaussian",alpha=0)
  r_cv_fit <- cv.glmnet(X, Yr,family = "gaussian",type.measure="mse",alpha=0)
  ridge_best <- glmnet(X, Yr,family = "gaussian",lambda = r_cv_fit$lambda.min,alpha=0)
  coefficients <- coef(ridge_best,r_cv_fit$lambda.min)
  pe_r[,i] <- as.matrix(coefficients)
}
pe_r <- t(pe_r)
r_ridge <- as.matrix(pe_r[,1])
A_ridge <- pe_r[,2:(p+1)]
colnames(A_ridge) <- colnames(X)
row.names(A_ridge) <- colnames(A_ridge)
A_ridge_write <- adjustdata(A_ridge)
filename_A <- paste0("A_",group4,".txt")
write.table(A_ridge_write,file = filename_A,row.names = F,col.names = T, sep = "\t",quote = F)
filename_r <- paste0("r_",group4,".txt")
write.table(r_ridge,file = filename_r,row.names = T,col.names = T, sep = "\t",quote = F)

###Driver prediction
group4_Driver_Score <- Driver_Score(A_ridge,r_ridge)
filename_Score <- paste0(group4,"_Driver_Score.txt")
write.table(group4_Driver_Score,file = filename_Score,row.names = F,col.names = T, sep = "\t",quote = F)
group4_driver <- group4_Driver_Score[c(1:5),]
colnames(group4_driver) <- c("Driver","Score")
group4_driver$Group <- group4
filename_driver <- paste0(group4,"_driver.txt")
write.table(group4_driver,file = filename_driver,row.names = F,col.names = T, sep = "\t",quote = F)

###Summary
summary <- rbind(group1_driver,group2_driver,group3_driver,group4_driver)
write.table(summary,file ="UC_driver_summary.txt",row.names = F,col.names = T, sep = "\t",quote = F)
