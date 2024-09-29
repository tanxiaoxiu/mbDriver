rm(list = ls(all = TRUE)) 
library(MASS)
library(glmnet)
library(Rfast)
library(mgcv)
library(tidyverse)
library(dplyr)
library(mbDriver)

setwd("~/mbDriver/real_data/Fiber_diet")

#Input data
abundance <- read.delim("absolute_abundance.txt", row.names = 1, sep = '\t', check.names = FALSE)
metadata <- read.delim("metadata.txt", row.names = 1, sep = '\t', check.names = FALSE)
 
#preprocessing
data_abundance <- rbind(abundance,Total=colSums(abundance))
data_abundance_sort <- data_abundance[,order(-as.numeric(data_abundance["Total",]))]
data_fiber_diet <- merge(metadata,data_abundance_sort,by = 'row.names', all = F)
#write.table(data_fiber_diet,file = "data_fiber.txt",row.names = F,col.names = T, sep = "\t",quote = F)

threshold <- 0.8 * nrow(data_fiber_diet)
filtered_data <- data_fiber_diet %>%
  dplyr::select(Row.names, 6:ncol(data_fiber_diet)) %>%
  select_if(~sum(. != 0) > threshold)
filtered_data <- bind_cols(data_fiber_diet %>% dplyr::select(1:5), filtered_data)

data_fiber_diet <- filtered_data[,-6]

p = 10

##############################In group
###Data preprocessing
#BC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("BI1","BI2","BI3","BI4","BI5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),]
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[,-1],2,function(col) f(col, t))
  s_f_deriv <- apply(otu_table[,-1],2,function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_bj <- as.matrix(F_d)
F_dt_bj <- as.matrix(F_dt)
log_F_dt_bj <- F_dt/F_d
otu_cor_bj <- t(otu_cor)
#SC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("SI1","SI2","SI3","SI4","SI5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),] 
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[,-1],2,function(col) f(col, t))
  s_f_deriv <- apply(otu_table[,-1],2,function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_sh <- as.matrix(F_d)
F_dt_sh <- as.matrix(F_dt)
log_F_dt_sh <- F_dt/F_d
otu_cor_sh <- t(otu_cor)
#HC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("HI1","HI2","HI4","HI5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),] 
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[,-1],2,function(col) f(col, t))
  s_f_deriv <- apply(otu_table[,-1],2,function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_hn <- as.matrix(F_d)
F_dt_hn <- as.matrix(F_dt)
log_F_dt_hn <- F_dt/F_d
otu_cor_hn <- t(otu_cor)
#GC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("GI1","GI2","GI3","GI4")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),] 
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[,-1],2,function(col) f(col, t))
  s_f_deriv <- apply(otu_table[,-1],2,function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_gd <- as.matrix(F_d)
F_dt_gd <- as.matrix(F_dt)
log_F_dt_gd <- F_dt/F_d
otu_cor_gd <- t(otu_cor)

#In group
F_In <-Reduce(rbind,list(F_d_bj,F_d_sh,F_d_hn,F_d_gd))
F_dt_In <-Reduce(rbind,list(F_dt_bj,F_dt_sh,F_dt_hn,F_dt_gd))
F_F_dt_In <-Reduce(rbind,list(log_F_dt_bj,log_F_dt_sh,log_F_dt_hn,log_F_dt_gd))

Case <-Reduce(cbind,list(otu_cor_bj,otu_cor_sh,otu_cor_hn,otu_cor_gd))
colnames(Case) <- paste0("Sample",seq(1,ncol(Case)))
Case <- Case[-1,]
Case <- adjustdata(Case)

###Parameter estimation
X <- as.matrix(F_In)
Y <- as.matrix(F_F_dt_In)
pe_r <- as.data.frame(matrix(nrow=(p+1)))
for (i in 1:p){
  set.seed(12)
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
write.table(A_ridge_write,file ="A_In.txt",row.names = F,col.names = T, sep = "\t",quote = F)
write.table(r_ridge,file ="r_In.txt",row.names = T,col.names = T, sep = "\t",quote = F)

###Driver prediction
In_Score <- Driver_Score(A_ridge,r_ridge)
#write.table(In_Score,file ="In_Score.txt",row.names = F,col.names = T, sep = "\t",quote = F)

In_driver <- In_Score[c(1:5),]
colnames(In_driver) <- c("Driver","Score")
In_driver$Group <- "In"
#write.table(In_driver,file ="In_driver.txt",row.names = F,col.names = T, sep = "\t",quote = F)


###################################Rs group
p = 10
###Data preprocessing
#BC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("BR1","BR2","BR3","BR4")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject== i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),] 
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[,-1],2,function(col) f(col, t))
  s_f_deriv <- apply(otu_table[,-1],2,function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_bj <- as.matrix(F_d)
F_dt_bj <- as.matrix(F_dt)
log_F_dt_bj <- F_dt/F_d
otu_cor_bj <- t(otu_cor)
#SC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("SR1","SR2","SR3","SR4","SR5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),] 
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[, -1], 2, function(col) f(col, t))
  s_f_deriv <- apply(otu_table[, -1], 2, function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_sh <- as.matrix(F_d)
F_dt_sh <- as.matrix(F_dt)
log_F_dt_sh <- F_dt/F_d
otu_cor_sh <- t(otu_cor)
#HC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("HR1","HR2","HR3","HR4","HR5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),] 
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[, -1], 2, function(col) f(col, t))
  s_f_deriv <- apply(otu_table[, -1], 2, function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_hn <- as.matrix(F_d)
F_dt_hn <- as.matrix(F_dt)
log_F_dt_hn <- F_dt/F_d
otu_cor_hn <- t(otu_cor)
#GC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("GR1","GR2","GR3","GR4","GR5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),] 
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[, -1], 2, function(col) f(col, t))
  s_f_deriv <- apply(otu_table[, -1], 2, function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_gd <- as.matrix(F_d)
F_dt_gd <- as.matrix(F_dt)
log_F_dt_gd <- F_dt/F_d
otu_cor_gd <- t(otu_cor)

#Rs group
F_Rs <-Reduce(rbind,list(F_d_bj,F_d_sh,F_d_hn,F_d_gd))
F_dt_Rs <-Reduce(rbind,list(F_dt_bj,F_dt_sh,F_dt_hn,F_dt_gd))
F_F_dt_Rs <-Reduce(rbind,list(log_F_dt_bj,log_F_dt_sh,log_F_dt_hn,log_F_dt_gd))

Case <-Reduce(cbind,list(otu_cor_bj,otu_cor_sh,otu_cor_hn,otu_cor_gd))
colnames(Case) <- paste0("Sample",seq(1,ncol(Case)))
Case <- Case[-1,]
Case <- adjustdata(Case)

###Parameter estimation
X <- as.matrix(F_Rs)
Y <- as.matrix(F_F_dt_Rs)
pe_r <- as.data.frame(matrix(nrow=(p+1)))
for (i in 1:p){
  set.seed(12)
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
write.table(A_ridge_write,file ="A_Rs.txt",row.names = F,col.names = T, sep = "\t",quote = F)
write.table(r_ridge,file ="r_Rs.txt",row.names = T,col.names = T, sep = "\t",quote = F)

###Driver prediction
Rs_Score <- Driver_Score(A_ridge,r_ridge)
#write.table(Rs_Score,file ="Rs_Score.txt",row.names = F,col.names = T, sep = "\t",quote = F)
Rs_driver <- Rs_Score[c(1:5),]
colnames(Rs_driver) <- c("Driver","Score")
Rs_driver$Group <- "Rs"
#write.table(Rs_driver,file ="Rs_driver.txt",row.names = F,col.names = T, sep = "\t",quote = F)


###################################Con group
p=10

###Data preprocessing
#BC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("BC1","BC2","BC3","BC4","BC5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject== i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),]
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[, -1], 2, function(col) f(col, t))
  s_f_deriv <- apply(otu_table[, -1], 2, function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_bj <- as.matrix(F_d)
F_dt_bj <- as.matrix(F_dt)
log_F_dt_bj <- F_dt/F_d
otu_cor_bj <- t(otu_cor)
#SC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("SC1","SC2","SC3","SC4","SC5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[, -1], 2, function(col) f(col, t))
  s_f_deriv <- apply(otu_table[, -1], 2, function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_sh <- as.matrix(F_d)
F_dt_sh <- as.matrix(F_dt)
log_F_dt_sh <- F_dt/F_d
otu_cor_sh <- t(otu_cor)
#HC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("HC1","HC2","HC3","HC4","HC5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),] 
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[, -1], 2, function(col) f(col, t))
  s_f_deriv <- apply(otu_table[, -1], 2, function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_hn <- as.matrix(F_d)
F_dt_hn <- as.matrix(F_dt)
log_F_dt_hn <- F_dt/F_d
otu_cor_hn <- t(otu_cor)
#GC
F_d <- data.frame()
F_dt <- data.frame()
C1 <- data.frame()
otu_cor <- data.frame()
for (i in c("GC1","GC2","GC3","GC4","GC5")){
  subject_n <- data_fiber_diet[data_fiber_diet$Subject==i,]
  subject_n <- subject_n[,-c(1:3,5)]
  subject_n_top <- subject_n[,c(1:(p+1))]
  otu_table <- subject_n_top[order(subject_n_top[,"Time"]),] 
  t <- otu_table[,1]
  n <- nrow(otu_table)
  s_f <- apply(otu_table[, -1], 2, function(col) f(col, t))
  s_f_deriv <- apply(otu_table[, -1], 2, function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
  otu_cor <- rbind(otu_cor,otu_table)
}
F_d_gd <- as.matrix(F_d)
F_dt_gd <- as.matrix(F_dt)
log_F_dt_gd <- F_dt/F_d
otu_cor_gd <- t(otu_cor)

#Con group
F_Con <-Reduce(rbind,list(F_d_bj,F_d_sh,F_d_hn,F_d_gd))
F_dt_Con <-Reduce(rbind,list(F_dt_bj,F_dt_sh,F_dt_hn,F_dt_gd))
F_F_dt_Con <-Reduce(rbind,list(log_F_dt_bj,log_F_dt_sh,log_F_dt_hn,log_F_dt_gd))

Con <-Reduce(cbind,list(otu_cor_bj,otu_cor_sh,otu_cor_hn,otu_cor_gd))
colnames(Con) <- paste0("Sample",seq(1,ncol(Con)))
Con <- Con[-1,]
Con <- adjustdata(Con)

###Parameter estimation
X <- as.matrix(F_Con)
Y <- as.matrix(F_F_dt_Con)
pe_r <- as.data.frame(matrix(nrow=(p+1)))
for (i in 1:p){
  set.seed(1)
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
write.table(A_ridge_write,file ="A_Con.txt",row.names = F,col.names = T, sep = "\t",quote = F)
write.table(r_ridge,file ="r_Con.txt",row.names = T,col.names = T, sep = "\t",quote = F)

###Driver prediction
Con_Score <- Driver_Score(A_ridge,r_ridge)
#write.table(Con_Score,file ="Con_Score.txt",row.names = F,col.names = T, sep = "\t",quote = F)

Con_driver <- Con_Score[c(1:5),]
colnames(Con_driver) <- c("Driver","Score")
Con_driver$Group <- "Con"
#write.table(Con_driver,file ="Con_driver.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###summary
summary <- rbind(Con_driver,Rs_driver,In_driver)
write.table(summary,file ="Fiber_driver_summary.txt",row.names = F,col.names = T, sep = "\t",quote = F)

