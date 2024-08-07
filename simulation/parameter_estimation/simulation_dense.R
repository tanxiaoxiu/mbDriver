library(sampling)
library(tidyverse)
library(deSolve)
library(gtools)
library(MASS)
library(dplyr)
library(glmnet)
library(Rfast)
library(mgcv)

timepoint = "t8"   #"t8", "t13", "t18", "t25" 
theta0 = 1        # 1, 3, 5 
path <- paste0("~/mbDriver/simulation/parameter_estimation/dense/", timepoint, "/theta", theta0)
setwd(path)


p = 10 #p=10,15
n = 10 #n=10,15

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

##difference
#calculate deltaY
deltaY<- function(X1){ 
  Y<-c()
  for (j in 1:l){
    y <- (log(X1[j+1])-log(X1[j]))/(time[j+1]-time[j])
    Y <- c(Y,y)
  }
  return(Y)
} 


### write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}

#data
A0_tables <- list()
r0_tables <- list()
count_tables <- list()

#spline
sp_rmse_tables <- list()
sp_re_rmse_tables <- list()
sp_A_lse_tables <- list()
sp_r_lse_tables <- list()
sp_A_lasso_tables <- list()
sp_r_lasso_tables <- list()
sp_A_ridge_tables <- list()
sp_r_ridge_tables <- list()
sp_A_elastic_tables <- list()
sp_r_elastic_tables <- list()

#difference
di_rmse_tables <- list()
di_re_rmse_tables <- list()
di_A_lse_tables <- list()
di_r_lse_tables <- list()
di_A_lasso_tables <- list()
di_r_lasso_tables <- list()
di_A_ridge_tables <- list()
di_r_ridge_tables <- list()
di_A_elastic_tables <- list()
di_r_elastic_tables <- list()


##

invalid_iter = c(73,100)
iterations = 100 + length(invalid_iter)
for (iter in 1:iterations){
  if(iter %in% invalid_iter){
    next 
  }
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
  A0[sample(1:matrx_len,0.2*matrx_len)] <- 0
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
    
    L <-nrow(otu_time) #time point
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
    
    #difference
    subject_sparse_di <- subject_sparse
    time <- subject_sparse_di[,1]
    l <- nrow(subject_sparse_di) - 1
    subject_sparse_di[subject_sparse_di == 0] <- subject_sparse_di[subject_sparse_di == 0] + 1
    subject_sparse_di <- as.matrix(subject_sparse_di)
    
    Yi <- apply(subject_sparse_di[,-1],2,deltaY)
    
    Yi <- as.data.frame(Yi)
    Xi <- subject_sparse_di[,-1]
    Xi <- Xi[-nrow(Xi),]
    Xi <- as.data.frame(Xi)
    colnames(Yi) <- row.names(r0_write)
    colnames(Xi) <- row.names(r0_write)
    
    Ci <- cbind(rep(1,l),Xi)
    di_C <- dplyr::bind_rows(di_C,Ci)
    di_X <- dplyr::bind_rows(di_X,Xi)
    di_Y <- dplyr::bind_rows(di_Y,Yi)
    
  }

  N <- nrow(count)
  OTU_ID <- as.matrix(1:N)
  colnames(OTU_ID) <- "#OTU ID"
  counts <- cbind(OTU_ID,count)
  counts <- t(counts)
  count_tables[[iter]] <- counts
  
  ######################################spline
  F_d <- as.matrix(F_d)
  F_dt <- as.matrix(F_dt)
  log_F_dt <- F_dt/F_d
  
  C1 <- cbind(rep(1,nrow(F_d)),F_d)
  sp_X <- F_d
  sp_C <- as.matrix(C1)
  sp_Y <- as.matrix(log_F_dt)
  
  ##Parameter estimation: lse,lasso,ridge,elastic net
  #####LSE
  sp_pe_LSE <- solve( t(sp_C) %*% sp_C ) %*% t(sp_C) %*% sp_Y
  sp_pe_LSE <- t(sp_pe_LSE)
  sp_r_LSE <- as.matrix(sp_pe_LSE[,1])
  sp_A_LSE <- sp_pe_LSE[,2:(p+1)]
  sp_A_lse_tables[[iter]] <- sp_A_LSE
  sp_r_lse_tables[[iter]] <- sp_r_LSE
  
  ###lasso
  sp_pe_lasso <- as.data.frame(matrix(nrow=(p+1)))
  
  for (i in 1:p){
    set.seed(12345)
    Ys <- as.matrix(sp_Y[,i])
    pf <- c(rep(1,p))
    pf[i] <- 0 
    lasso_fit <- glmnet(sp_X, Ys,family = "gaussian",alpha=1)
    l_cv_fit <- cv.glmnet(sp_X, Ys,family = "gaussian",type.measure="mse")
    lasso_best <- glmnet(sp_X, Ys,family = "gaussian",lambda = l_cv_fit$lambda.min)
    coefficients <- coef(lasso_best,l_cv_fit$lambda.min)
    sp_pe_lasso[,i] <- as.matrix(coefficients)
  }
  
  sp_pe_lasso <- t(sp_pe_lasso)
  sp_r_lasso <- as.matrix(sp_pe_lasso[,1])
  sp_A_lasso <- sp_pe_lasso[,2:(p+1)]
  colnames(sp_A_lasso) <- colnames(sp_X)
  row.names(sp_A_lasso) <- colnames(sp_A_lasso)
  sp_A_lasso_tables[[iter]] <- sp_A_lasso
  sp_r_lasso_tables[[iter]] <- sp_r_lasso
  
  ###ridge
  sp_pe_r <- as.data.frame(matrix(nrow=(p+1)))
  for (i in 1:p){
    set.seed(12345)
    Yr <- as.matrix(sp_Y[,i])
    pf <- c(rep(1,p))
    pf[i] <- 0 
    ridge_fit <- glmnet(sp_X, Yr,family = "gaussian",alpha=0)
    r_cv_fit <- cv.glmnet(sp_X, Yr,family = "gaussian",type.measure="mse",alpha=0)
    ridge_best <- glmnet(sp_X, Yr,family = "gaussian",lambda = r_cv_fit$lambda.min,alpha=0)
    coefficients <- coef(ridge_best,r_cv_fit$lambda.min)
    sp_pe_r[,i] <- as.matrix(coefficients)
  }
  
  sp_pe_r <- t(sp_pe_r)
  sp_r_ridge <- as.matrix(sp_pe_r[,1])
  sp_A_ridge <- sp_pe_r[,2:(p+1)]
  colnames(sp_A_ridge) <- colnames(sp_X)
  row.names(sp_A_ridge) <- colnames(sp_A_ridge)
  sp_A_ridge_tables[[iter]] <- sp_A_ridge
  sp_r_ridge_tables[[iter]] <- sp_r_ridge
  
  ###elastic net
  sp_pe_e <- as.data.frame(matrix(nrow=(p+1)))
  
  for (i in 1:p){
    set.seed(12345)
    Ye <- as.matrix(sp_Y[,i])
    pf <- c(rep(1,p))
    pf[i] <- 0 
    elastic_fit <- glmnet(sp_X, Ye,family = "gaussian",alpha=0.5)
    e_cv_fit <- cv.glmnet(sp_X, Ye,family = "gaussian",type.measure="mse",alpha=0.5)
    elastic_best <- glmnet(sp_X, Ye,family = "gaussian",lambda = e_cv_fit$lambda.min,alpha=0.5)
    coefficients <- coef(elastic_best,e_cv_fit$lambda.min)
    sp_pe_e[,i] <- as.matrix(coefficients)
  }
  
  sp_pe_e <- t(sp_pe_e)
  sp_r_elastic <- as.matrix(sp_pe_e[,1])
  sp_A_elastic <- sp_pe_e[,2:(p+1)]
  colnames(sp_A_elastic) <- colnames(sp_X)
  row.names(sp_A_elastic) <- colnames(sp_A_elastic)
  sp_A_elastic_tables[[iter]] <- sp_A_elastic
  sp_r_elastic_tables[[iter]] <- sp_r_elastic
  
  ###RMSE
  sp_A_lasso <- as.matrix(sp_A_lasso)
  sp_A_ridge <- as.matrix(sp_A_ridge)
  sp_A_elastic <- as.matrix(sp_A_elastic)
  #A
  sp_rmse_A_LSE <- sqrt(sum((sp_A_LSE - A0)^2) / (ncol(A0)*nrow(A0)))
  sp_rmse_A_lasso <- sqrt(sum((sp_A_lasso - A0)^2) / (ncol(A0)*nrow(A0)))
  sp_rmse_A_ridge <- sqrt(sum((sp_A_ridge - A0)^2) / (ncol(A0)*nrow(A0)))
  sp_rmse_A_elastic <- sqrt(sum((sp_A_elastic - A0)^2) / (ncol(A0)*nrow(A0)))
  
  #r
  sp_rmse_r_LSE <- sqrt(sum((sp_r_LSE - r0)^2) / nrow(r0))
  sp_rmse_r_lasso <- sqrt(sum((sp_r_lasso - r0)^2) / nrow(r0))
  sp_rmse_r_ridge <- sqrt(sum((sp_r_ridge - r0)^2) / nrow(r0))
  sp_rmse_r_elastic <- sqrt(sum((sp_r_elastic - r0)^2) / nrow(r0))
  
  sp_RMSE <- matrix(c(sp_rmse_A_LSE, sp_rmse_A_lasso,sp_rmse_A_ridge,sp_rmse_A_elastic,
                      sp_rmse_r_LSE, sp_rmse_r_lasso, sp_rmse_r_ridge,sp_rmse_r_elastic),nrow=4, ncol=2, byrow=FALSE)
  row.names(sp_RMSE) <- c("LSE","lasso","ridge","elastic")
  colnames(sp_RMSE) <- c("RMSE_A","RMSE_r")
  
  sp_rmse_tables[[iter]] <- sp_RMSE
  
  ##relative RMSE
  #A
  rmse_A0 <- sqrt(sum((A0)^2) / (ncol(A0)*nrow(A0)))
  #r
  rmse_r0 <- sqrt(sum((r0)^2) / nrow(r0))
  
  sp_re_RMSE <- matrix(c(sp_rmse_A_LSE/rmse_A0, sp_rmse_A_lasso/rmse_A0,sp_rmse_A_ridge/rmse_A0,sp_rmse_A_elastic/rmse_A0,
                         sp_rmse_r_LSE/rmse_r0, sp_rmse_r_lasso/rmse_r0, sp_rmse_r_ridge/rmse_r0,sp_rmse_r_elastic/rmse_r0),nrow=4, ncol=2, byrow=FALSE)
  row.names(sp_re_RMSE) <- c("LSE","lasso","ridge","elastic")
  colnames(sp_re_RMSE) <- c("RMSE_A","RMSE_r")
  
  sp_re_rmse_tables[[iter]] <- sp_re_RMSE
  print(paste("spline:", iter))
  
  
  ######################################difference
  di_X <- as.matrix(di_X)
  di_C <- as.matrix(di_C)
  di_Y <- as.matrix(di_Y)
  ##Parameter estimation: lse,lasso,ridge,elastic net
  #####LSE
  di_pe_LSE <- solve( t(di_C) %*% di_C ) %*% t(di_C) %*% di_Y
  di_pe_LSE <- t(di_pe_LSE)
  di_r_LSE <- as.matrix(di_pe_LSE[,1])
  di_A_LSE <- di_pe_LSE[,2:(p+1)]
  di_A_lse_tables[[iter]] <- di_A_LSE
  di_r_lse_tables[[iter]] <- di_r_LSE
  
  
  ###lasso
  di_pe_lasso <- as.data.frame(matrix(nrow=(p+1)))
  for (i in 1:p){
    set.seed(12345)
    Ys <- as.matrix(di_Y[,i])
    pf <- c(rep(1,p))
    pf[i] <- 0 
    s_lasso_fit <- glmnet(di_X, Ys,family = "gaussian",alpha=1)
    s_cv_fit <- cv.glmnet(di_X, Ys,family = "gaussian",type.measure="mse")
    s_lasso_best <- glmnet(di_X, Ys,family = "gaussian",lambda = s_cv_fit$lambda.min)
    coefficients <- coef(s_lasso_best,s_cv_fit$lambda.min)
    di_pe_lasso[,i] <- as.matrix(coefficients)
  }
  
  di_pe_lasso <- t(di_pe_lasso)
  di_r_lasso <- as.matrix(di_pe_lasso[,1])
  di_A_lasso <- di_pe_lasso[,2:(p+1)]
  
  colnames(di_A_lasso) <- colnames(di_X)
  row.names(di_A_lasso) <- colnames(di_A_lasso)
  di_A_lasso_tables[[iter]] <- di_A_lasso
  di_r_lasso_tables[[iter]] <- di_r_lasso
  
  ###ridge
  di_pe_r <- as.data.frame(matrix(nrow=(p+1)))
  for (i in 1:p){
    set.seed(12345)
    Yr <- as.matrix(di_Y[,i])
    pf <- c(rep(1,p))
    pf[i] <- 0 
    ridge_fit <- glmnet(di_X, Yr,family = "gaussian",alpha=0)
    
    r_cv_fit <- cv.glmnet(di_X, Yr,family = "gaussian",type.measure="mse",alpha=0)
    
    ridge_best <- glmnet(di_X, Yr,family = "gaussian",lambda = r_cv_fit$lambda.min,alpha=0)
    coefficients <- coef(ridge_best,r_cv_fit$lambda.min)
    di_pe_r[,i] <- as.matrix(coefficients)
  }
  
  di_pe_r <- t(di_pe_r)
  di_r_ridge <- as.matrix(di_pe_r[,1])
  di_A_ridge <- di_pe_r[,2:(p+1)]
  colnames(di_A_ridge) <- colnames(di_X)
  row.names(di_A_ridge) <- colnames(di_A_ridge)
  di_A_ridge_tables[[iter]] <- di_A_ridge
  di_r_ridge_tables[[iter]] <- di_r_ridge
  
  ###elastic net
  di_pe_e <- as.data.frame(matrix(nrow=(p+1)))
  for (i in 1:p){
    set.seed(12345)
    Ye <- as.matrix(di_Y[,i])
    pf <- c(rep(1,p))
    pf[i] <- 0 
    elastic_fit <- glmnet(di_X, Ye,family = "gaussian",alpha=0.5)
    e_cv_fit <- cv.glmnet(di_X, Ye,family = "gaussian",type.measure="mse",alpha=0.5)
    elastic_best <- glmnet(di_X, Ye,family = "gaussian",lambda = e_cv_fit$lambda.min,alpha=0.5)
    coefficients <- coef(elastic_best,e_cv_fit$lambda.min)
    di_pe_e[,i] <- as.matrix(coefficients)
  }
  
  di_pe_e <- t(di_pe_e)
  di_r_elastic <- as.matrix(di_pe_e[,1])
  di_A_elastic <- di_pe_e[,2:(p+1)]
  colnames(di_A_elastic) <- colnames(di_X)
  row.names(di_A_elastic) <- colnames(di_A_elastic)
  di_A_elastic_tables[[iter]] <- di_A_elastic
  di_r_elastic_tables[[iter]] <- di_r_elastic
  
  ###RMSE
  di_A_lasso <- as.matrix(di_A_lasso)
  di_A_ridge <- as.matrix(di_A_ridge)
  di_A_elastic <- as.matrix(di_A_elastic)
  
  #A
  di_rmse_A_LSE <- sqrt(sum((di_A_LSE - A0)^2) / (ncol(A0)*nrow(A0)))
  di_rmse_A_lasso <- sqrt(sum((di_A_lasso - A0)^2) / (ncol(A0)*nrow(A0)))
  di_rmse_A_ridge <- sqrt(sum((di_A_ridge - A0)^2) / (ncol(A0)*nrow(A0)))
  di_rmse_A_elastic <- sqrt(sum((di_A_elastic - A0)^2) / (ncol(A0)*nrow(A0)))
  
  #r
  di_rmse_r_LSE <- sqrt(sum((di_r_LSE - r0)^2) / nrow(r0))
  di_rmse_r_lasso <- sqrt(sum((di_r_lasso - r0)^2) / nrow(r0))
  di_rmse_r_ridge <- sqrt(sum((di_r_ridge - r0)^2) / nrow(r0))
  di_rmse_r_elastic <- sqrt(sum((di_r_elastic - r0)^2) / nrow(r0))
  
  di_RMSE <- matrix(c(di_rmse_A_LSE, di_rmse_A_lasso,di_rmse_A_ridge,di_rmse_A_elastic,
                      di_rmse_r_LSE, di_rmse_r_lasso, di_rmse_r_ridge,di_rmse_r_elastic),nrow=4, ncol=2, byrow=FALSE)
  row.names(di_RMSE) <- c("LSE","lasso","ridge","elastic")
  colnames(di_RMSE) <- c("RMSE_A","RMSE_r")
  
  di_rmse_tables[[iter]] <- di_RMSE
  
  ##relative RMSE
  #A
  rmse_A0 <- sqrt(sum((A0)^2) / (ncol(A0)*nrow(A0)))
  #r
  rmse_r0 <- sqrt(sum((r0)^2) / nrow(r0))
  
  di_re_RMSE <- matrix(c(di_rmse_A_LSE/rmse_A0, di_rmse_A_lasso/rmse_A0,di_rmse_A_ridge/rmse_A0,di_rmse_A_elastic/rmse_A0,
                         di_rmse_r_LSE/rmse_r0, di_rmse_r_lasso/rmse_r0, di_rmse_r_ridge/rmse_r0,di_rmse_r_elastic/rmse_r0),nrow=4, ncol=2, byrow=FALSE)
  row.names(di_re_RMSE) <- c("LSE","lasso","ridge","elastic")
  colnames(di_re_RMSE) <- c("RMSE_A","RMSE_r")
  
  di_re_rmse_tables[[iter]] <- di_re_RMSE
  print(paste("difference:", iter))
}

##write spline
sp_rmse_tables <- lapply(sp_rmse_tables, function(matrix) {
  if (is.null(matrix)) {
    matrix <- matrix(0, nrow = 4, ncol = 2) 
  }
  matrix
})

sp_re_rmse_tables <- lapply(sp_re_rmse_tables, function(matrix) {
  if (is.null(matrix)) {
    matrix <- matrix(0, nrow = 4, ncol = 2)  
  }
  matrix
})

sp_mean_rmse_tables <- Reduce(`+`, sp_rmse_tables) / 100
sp_mean_re_rmse_tables <- Reduce(`+`, sp_re_rmse_tables) / 100
sp_mean_rmse_tables_write <- adjustdata(sp_mean_rmse_tables)
sp_mean_re_rmse_tables_write <- adjustdata(sp_mean_re_rmse_tables)

write.table(sp_mean_rmse_tables_write,file ="sp_mean_rmse_tables.txt",row.names = F,col.names = T, sep = "\t",quote = F)
write.table(sp_mean_re_rmse_tables_write,file ="sp_mean_re_rmse_tables.txt",row.names = F,col.names = T, sep = "\t",quote = F)

##write difference
di_rmse_tables <- lapply(di_rmse_tables, function(matrix) {
  if (is.null(matrix)) {
    matrix <- matrix(0, nrow = 4, ncol = 2)  
  }
  matrix
})

di_re_rmse_tables <- lapply(di_re_rmse_tables, function(matrix) {
  if (is.null(matrix)) {
    matrix <- matrix(0, nrow = 4, ncol = 2)  
  }
  matrix
})

di_mean_rmse_tables <- Reduce(`+`, di_rmse_tables) / 100
di_mean_re_rmse_tables <- Reduce(`+`, di_re_rmse_tables) / 100

di_mean_rmse_tables_write <- adjustdata(di_mean_rmse_tables)
di_mean_re_rmse_tables_write <- adjustdata(di_mean_re_rmse_tables)
write.table(di_mean_rmse_tables_write,file ="di_mean_rmse_tables.txt",row.names = F,col.names = T, sep = "\t",quote = F)
write.table(di_mean_re_rmse_tables_write,file ="di_mean_re_rmse_tables.txt",row.names = F,col.names = T, sep = "\t",quote = F)

##save data
save(A0_tables, file = "A0_tables.Rdata")
save(r0_tables, file = "r0_tables.Rdata")
save(count_tables, file = "count_tables.Rdata")

##save results of spline
save(sp_rmse_tables, file = "sp_rmse_tables.Rdata")
save(sp_re_rmse_tables, file = "sp_re_rmse_tables.Rdata")
save(sp_A_lse_tables, file = "sp_A_lse_tables.Rdata")
save(sp_r_lse_tables, file = "sp_r_lse_tables.Rdata")
save(sp_A_lasso_tables, file = "sp_A_lasso_tables.Rdata")
save(sp_r_lasso_tables, file = "sp_r_lasso_tables.Rdata")
save(sp_A_ridge_tables, file = "sp_A_ridge_tables.Rdata")
save(sp_r_ridge_tables, file = "sp_r_ridge_tables.Rdata")
save(sp_A_elastic_tables, file = "sp_A_elastic_tables.Rdata")
save(sp_r_elastic_tables, file = "sp_r_elastic_tables.Rdata")

##save results of difference
save(di_rmse_tables, file = "di_rmse_tables.Rdata")
save(di_re_rmse_tables, file = "di_re_rmse_tables.Rdata")
save(di_A_lse_tables, file = "di_A_lse_tables.Rdata")
save(di_r_lse_tables, file = "di_r_lse_tables.Rdata")
save(di_A_lasso_tables, file = "di_A_lasso_tables.Rdata")
save(di_r_lasso_tables, file = "di_r_lasso_tables.Rdata")
save(di_A_ridge_tables, file = "di_A_ridge_tables.Rdata")
save(di_r_ridge_tables, file = "di_r_ridge_tables.Rdata")
save(di_A_elastic_tables, file = "di_A_elastic_tables.Rdata")
save(di_r_elastic_tables, file = "di_r_elastic_tables.Rdata")


