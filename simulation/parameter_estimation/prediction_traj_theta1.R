library(sampling)
library(tidyverse)
library(deSolve)
library(gtools)
library(MASS)
library(dplyr)
library(glmnet)
library(Rfast)
library(mgcv)

timepoint = "t8"  #"t8","t13","t18","t25"
path <- sprintf("~/mbDriver/simulation/parameter_estimation/trajectory/theta1/%s", timepoint)
setwd(path)

theta0 = 1 
p = 10 
n = 10 

### Generalized Lotka-Volterra model
GLV <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    x[x < 10^-5] <- 0 
    dxdt <- x * (r + A %*% x)
    list(dxdt)
  })
}

### general function to integrate GLV
integrate_GLV <- function(r, A, x0, maxtime = 50, steptime = 0.1){
  times <- seq(0, maxtime, by = steptime)
  parameters <- list(r = r, A = A)
  # solve numerically
  out <- ode(y = x0, times = times, 
             func = GLV, parms = parameters, 
             method = "ode45")
  return(out)
}

###calculate f
f <- function(X){ 
  Y<-c()
  X <- as.data.frame(X)
  fmla <- formula("X~s(t, m=2, k=5, bs = 'cr')")
  dat <- cbind(X,t)
  mod_nb <- gam(fmla, data = dat, family = nb(), method="REML")
  Y <- fitted.values(mod_nb)
  return(Y)
} 

###calculate f deriv
f_deriv <- function(X){
  Y<-c()
  X <- as.data.frame(X)
  fmla <- formula("X~s(t, m=2, k=5, bs = 'cr')")
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

### write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

#data
A0_tables <- list()
r0_tables <- list()
truth_tables <- list()

meta <- read.delim("metadata.txt", row.names = 1, sep = '\t', check.names = FALSE)
##
iterations = 100
for (iter in 1:iterations){
  set.seed(iter)
  ### generate A,r,x0
  r0 <- runif(p, min = 0, max = 0.2)
  r0 <- signif(r0,3)
  r0 <- as.matrix(r0)
  r0_write <- as.data.frame(r0)
  rownames(r0_write) <- paste('OTU',rownames(r0_write),sep='')
  colnames(r0_write) <- "r"
  A0 <- matrix(data=runif(p*p, min = -0.0005, max = 0.0005), nrow = p, ncol = p, byrow = FALSE, dimnames = NULL)
  matrx_len = dim(A0)[1] * dim(A0)[2]
  A0[sample(1:matrx_len,0.8*matrx_len)] <- 0
  diag(A0) <- rep(-0.001,times=p)
  A0_write <- as.data.frame(A0)
  rownames(A0_write) <- paste('OTU',rownames(A0_write),sep='')
  colnames(A0_write) <- rownames(A0_write)
  A0_tables[[iter]] <- A0_write
  r0_tables[[iter]] <- r0_write
  
  F_d <- as.data.frame(matrix(nrow=0,ncol=0))
  C1 <- as.data.frame(matrix(nrow=0,ncol=0))
  F_dt <- as.data.frame(matrix(nrow=0,ncol=0))
  count  <- as.data.frame(matrix(nrow=0,ncol=0))
  
  di_Y <- as.data.frame(matrix(nrow=0,ncol=p))
  colnames(di_Y) <- row.names(r0_write)
  di_C <- as.data.frame(matrix(nrow=0,ncol=p))
  colnames(di_C) <- row.names(r0_write)
  di_X <- as.data.frame(matrix(nrow=0,ncol=p))
  colnames(di_X) <- row.names(r0_write)
  
  #subject
  for (i in 1:n){
    set.seed(i)
    x0 <- runif(p, min = 100, max = 10000)
    x0 <- as.matrix(x0)
    otu_time <- integrate_GLV(r0, A0, x0)
    colnames(otu_time) <- paste("OTU",colnames(otu_time),sep = "")
    
    L <-nrow(otu_time) 
    E <- otu_time[,-1]
    
    otu_NB <- as.data.frame(matrix(nrow= L,ncol= p))
    
    for (j in 1:p){
      set.seed(j)
      otuj <-  as.data.frame(matrix(nrow=L,ncol=1))
      for (l in 1:L){
        otuj[l,] <-  rnegbin(1,mu=E[l,j],theta = theta0)
      }
      otu_NB[,j] <- otuj
    }
    
    otu_NB <- cbind(otu_time[,1],otu_NB)
    colnames(otu_NB) <- colnames(otu_time)
    
    #Sparse time
    t25 <-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25")
    t18 <-c("1","2","3","4","5","6","7","8","9","10","11","13","15","17","20","21","23","25")
    t13 <-c("1","2","3","4","6","7","11","15","17","19","21","23","25") 
    t8 <- c("1","2","3","7","11","15","20","25")
    if (timepoint == "t25") {
      mi <- t25
    } else if (timepoint == "t18") {
      mi <- t18
    } else if (timepoint == "t13") {
      mi <- t13
    } else if (timepoint == "t8") {
      mi <- t8
    } else {
      mi <- NULL
    }

    subject_sparse <- otu_NB[otu_NB$OTUtime %in% mi, ]
    subject_sparse_t <- subject_sparse[,-1]
    count <- rbind(count,subject_sparse_t)
    
    #spline
    subject_sparse_sp <- subject_sparse
    t <- subject_sparse_sp[,1]
    s_f <- apply(subject_sparse_sp[,-1],2,f)
    s_f_deriv <- apply(subject_sparse_sp[,-1],2,f_deriv)
    colnames(s_f) <- row.names(r0_write)
    colnames(s_f_deriv) <- row.names(r0_write)
    F_d <- rbind(F_d,s_f)
    F_dt <- rbind(F_dt,s_f_deriv)
  }
  
  F_d_meta <- cbind(meta,F_d)
  truth_tables[[iter]] <- F_d_meta
  
}

##save data
save(A0_tables, file = "A0_tables.Rdata")
save(r0_tables, file = "r0_tables.Rdata")
save(truth_tables, file = "truth_tables.Rdata")


### Generalized Lotka-Volterra model
GLV <- function(t, x, parameters){
  with(as.list(c(x, parameters)), {
    x[x < 10^-5] <- 0 
    dxdt <- x * (r + A %*% x)
    list(dxdt)
  })
}

### general function to integrate GLV
integrate_GLV <- function(r, A, x0, maxtime = 30, steptime = 0.1){
  times <- seq(0, maxtime, by = steptime)
  parameters <- list(r = r, A = A)
  # solve numerically
  out <- ode(y = x0, times = times, 
             func = GLV, parms = parameters, 
             method = "ode45")
  return(out)
}

### write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}


load("truth_tables.Rdata")

base_path1 <- "~/mbDriver/simulation/parameter_estimation/sparse"
base_path2 <- "~/mbDriver/simulation/parameter_estimation/other_methods"

load(paste0(base_path1, "/", timepoint, "/theta1/A0_tables.Rdata"))
load(paste0(base_path1, "/", timepoint, "/theta1/r0_tables.Rdata"))

# 读Ridge估计的参数
load(paste0(base_path1, "/", timepoint, "/theta1/sp_A_ridge_tables.Rdata"))
load(paste0(base_path1, "/", timepoint, "/theta1/sp_r_ridge_tables.Rdata"))

load(paste0(base_path2, "/", timepoint, "/theta1/BAL_list.Rdata"))
load(paste0(base_path2, "/", timepoint, "/theta1/BVS_list.Rdata"))
load(paste0(base_path2, "/", timepoint, "/theta1/MLRR_list.Rdata"))
load(paste0(base_path2, "/", timepoint, "/theta1/MLCRR_list.Rdata"))


#############K=1到10；
all_iter <- 100
n <- 10

r_rmse_list <- list()
corr_list <- list()

for (k in 1:n) {
  
  r_rmse <- data.frame(matrix(ncol = 5, nrow = all_iter))
  colnames(r_rmse) <- c("Ridge", "MLRR", "MLCRR", "BAL", "BVS")
  rownames(r_rmse) <- 1:all_iter
  
  corr <- data.frame(matrix(ncol = 5, nrow = all_iter))
  colnames(corr) <- c("Ridge", "MLRR", "MLCRR", "BAL", "BVS")
  rownames(corr) <- 1:all_iter
  
  if (timepoint == "t8") {
    invalid_iter <- c(2)
  } else if (timepoint == "t13") {
    invalid_iter <- c(27, 28, 30)
  } else if (timepoint == "t18") {
    invalid_iter <- c(18, 71)
  } else if (timepoint == "t25") {
    invalid_iter <- c(11, 43, 52, 60, 65, 77)
  } else {
    invalid_iter <- NULL  
  }
  
  for (iter in 1:all_iter) {
    if (is.null(BAL_list[[iter]]) || iter %in% invalid_iter) {
      next
    }
    
    truth_iter1 <- truth_tables[[iter]]
    truth_iter <- truth_iter1[truth_iter1$Subject == k, ]
    x0_truth <- as.matrix(as.numeric(truth_iter[1,-c(1,2)]))

    A_ridge <- as.matrix(sp_A_ridge_tables[[iter]])
    r_ridge <- as.matrix(sp_r_ridge_tables[[iter]])
    traj_pre_ridge <- integrate_GLV(r_ridge,A_ridge,x0_truth)
    traj_pre_ridge <- as.data.frame(traj_pre_ridge)
    traj_pre_ridge$time <- traj_pre_ridge$time + 1
    traj_ridge_time <- traj_pre_ridge[traj_pre_ridge$time %in% mi, ]
    
    A_MLRR <- as.matrix(MLRR_list[[iter]][,-1])
    r_MLRR <- as.matrix(MLRR_list[[iter]][,1])
    traj_pre_MLRR <- integrate_GLV(r_MLRR,A_MLRR,x0_truth)
    traj_pre_MLRR <- as.data.frame(traj_pre_MLRR)
    traj_pre_MLRR$time <- traj_pre_MLRR$time + 1
    traj_MLRR_time <- traj_pre_MLRR[traj_pre_MLRR$time %in% mi, ]
    
    A_MLCRR <- as.matrix(MLCRR_list[[iter]][,-1])
    r_MLCRR <- as.matrix(MLCRR_list[[iter]][,1])
    traj_pre_MLCRR <- integrate_GLV(r_MLCRR,A_MLCRR,x0_truth)
    traj_pre_MLCRR <- as.data.frame(traj_pre_MLCRR)
    traj_pre_MLCRR$time <- traj_pre_MLCRR$time + 1
    traj_MLCRR_time <- traj_pre_MLCRR[traj_pre_MLCRR$time %in% mi, ]
    
    A_BAL <- as.matrix(BAL_list[[iter]][,-1])
    r_BAL <- as.matrix(BAL_list[[iter]][,1])
    traj_pre_BAL <- integrate_GLV(r_BAL,A_BAL,x0_truth)
    traj_pre_BAL <- as.data.frame(traj_pre_BAL)
    traj_pre_BAL$time <- traj_pre_BAL$time + 1
    traj_BAL_time <- traj_pre_BAL[traj_pre_BAL$time %in% mi, ]
    
    A_BVS <- as.matrix(BVS_list[[iter]][,-1])
    r_BVS <- as.matrix(BVS_list[[iter]][,1])
    traj_pre_BVS <- integrate_GLV(r_BVS,A_BVS,x0_truth)
    traj_pre_BVS <- as.data.frame(traj_pre_BVS)
    traj_pre_BVS$time <- traj_pre_BVS$time + 1
    traj_BVS_time <- traj_pre_BVS[traj_pre_BVS$time %in% mi, ]
    
    truth_time <- as.data.frame(truth_iter)
    truth <- truth_time[,-c(1,2)]
    
    ridge_pre <- traj_ridge_time[,-1]
    MLRR_pre <- traj_MLRR_time[,-1]
    MLCRR_pre <- traj_MLCRR_time[,-1]
    BAL_pre <- traj_BAL_time[,-1]
    BVS_pre <- traj_BVS_time[,-1]
    
    #RMSE
    Y0 <- sqrt(mean(as.matrix(truth^2)))
    Ridge_rmse <- sqrt(mean(as.matrix((truth - ridge_pre)^2)))/Y0
    MLRR_rmse <- sqrt(mean(as.matrix((truth - MLRR_pre)^2)))/Y0
    MLCRR_rmse <- sqrt(mean(as.matrix((truth - MLCRR_pre)^2)))/Y0
    BAL_rmse <- sqrt(mean(as.matrix((truth - BAL_pre)^2)))/Y0
    BVS_rmse <- sqrt(mean(as.matrix((truth - BVS_pre)^2)))/Y0
    
    #Correlation
    Ridge_corr <- mean(diag(cor(truth, ridge_pre)))
    MLRR_corr <- mean(diag(cor(truth, MLRR_pre)))
    MLCRR_corr <- mean(diag(cor(truth, MLCRR_pre)))
    BAL_corr <- mean(diag(cor(truth, BAL_pre)))
    BVS_corr <- mean(diag(cor(truth, BVS_pre)))
    
    corr[iter, "Ridge"] <- Ridge_corr
    corr[iter, "MLRR"] <- MLRR_corr
    corr[iter, "MLCRR"] <- MLCRR_corr
    corr[iter, "BAL"] <- BAL_corr
    corr[iter, "BVS"] <- BVS_corr
    
    r_rmse[iter, "Ridge"] <- Ridge_rmse
    r_rmse[iter, "MLRR"] <- MLRR_rmse
    r_rmse[iter, "MLCRR"] <- MLCRR_rmse
    r_rmse[iter, "BAL"] <- BAL_rmse
    r_rmse[iter, "BVS"] <- BVS_rmse
  }
  
  r_rmse_list[[k]] <- r_rmse
  corr_list[[k]] <- corr
}


rmse_file_name <- paste0("r_rmse_", timepoint, ".Rdata")
corr_file_name <- paste0("corr_", timepoint, ".Rdata")
save(r_rmse_list, file = rmse_file_name)
save(corr_list, file = corr_file_name)
