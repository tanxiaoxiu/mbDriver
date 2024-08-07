library(mbDriver)
library(MASS)
library(glmnet)
library(Rfast)
library(mgcv)
library(dplyr)
library(pheatmap)

setwd("~/mbDriver/example")

#Input data
abundance <- read.delim("absolute_abundance.txt",row.names = 1, sep = '\t', check.names = FALSE)
metadata <- read.delim("metadata.txt", row.names = 1, sep = '\t', check.names = FALSE)

#preprocessing
data_abundance <- rbind(abundance,Total=colSums(abundance))
data_abundance_sort <- data_abundance[,order(-as.numeric(data_abundance["Total",]))]
meta_abundance <- merge(metadata,data_abundance_sort,by = 'row.names', all = F)


p = 10
n = 5
group1 = "group1"
group2 = "group2"

#####################group1
###Data preprocessing
F_d <- data.frame()
F_dt <- data.frame()

for (i in 1:n) {
  subject_i <- meta_abundance[meta_abundance$Subject == i, ]
  subject_i <- subject_i[subject_i$Group == group1, ]
  subject_i <- subject_i[, -c(1, 2, 4)]
  subject_i_top <- subject_i[, c(1:(p + 1))]
  otu_table <- subject_i_top[order(subject_i_top$Time), ]
  t <- otu_table[, 1]
  s_f <- apply(otu_table[, -1], 2, function(col) f(col, t))
  s_f_deriv <- apply(otu_table[, -1], 2, function(col) f_deriv(col, t))
  F_d <- rbind(F_d, s_f)
  F_dt <- rbind(F_dt, s_f_deriv)
}
F_d <- as.matrix(F_d)
F_dt <- as.matrix(F_dt)
F_F_dt <- F_dt/F_d
filename_F <- paste0("F_",group1,".txt")
write.table(F_d,file = filename_F,row.names = F,col.names = T, sep = "\t",quote = F)
filename_F_F_dt <- paste0("F_F_dt_",group1,".txt")
write.table(F_F_dt,file = filename_F_F_dt,row.names = F,col.names = T, sep = "\t",quote = F)

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

#####################group2
###Data preprocessing
F_d <- data.frame()
F_dt <- data.frame()
for (i in 1:n){
  subject_i <- meta_abundance[meta_abundance$Subject==i,]
  subject_i <- subject_i[subject_i$Group == group2,]
  subject_i <- subject_i[,-c(1,2,4)]
  subject_i_top <- subject_i[,c(1:(p+1))]
  otu_table <- subject_i_top[order(subject_i_top[,"Time"]),]
  t <- otu_table[,1]
  s_f <- apply(otu_table[,-1],2,function(col) f(col, t))
  s_f_deriv <- apply(otu_table[,-1],2,function(col) f_deriv(col, t))
  F_d <- rbind(F_d,s_f)
  F_dt <- rbind(F_dt,s_f_deriv)
}
F_d <- as.matrix(F_d)
F_dt <- as.matrix(F_dt)
F_F_dt <- F_dt/F_d

filename_F <- paste0("F_",group2,".txt")
write.table(F_d,file = filename_F,row.names = F,col.names = T, sep = "\t",quote = F)
filename_F_F_dt <- paste0("F_F_dt_",group2,".txt")
write.table(F_F_dt,file = filename_F_F_dt,row.names = F,col.names = T, sep = "\t",quote = F)

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

###Summary
summary <- rbind(group1_driver,group2_driver)
write.table(summary,file ="driver_summary.txt",row.names = F,col.names = T, sep = "\t",quote = F)

###plot
library(ggplot2)
library(tidyr)
library(tibble)
library(gtable)

MAX_N_OTU = p
l_tb <- summary
l_tb %>% pivot_wider(names_from = "Group" ,values_from = "Score")

w_tb <- l_tb %>%
  pivot_wider(id_cols = "Driver",
              names_from = "Group",
              values_from = "Score",
              values_fill=0) %>%
  mutate(across(-Driver, function(x) log10(x + 1))) %>%
  column_to_rownames("Driver")

p1 <- pheatmap(w_tb, cluster_cols = F,angle_col=0,fontsize=14)
ggsave(filename="plot_driver.png",plot=p1,device="png",dpi=600,units="in",width=6,height=7)
