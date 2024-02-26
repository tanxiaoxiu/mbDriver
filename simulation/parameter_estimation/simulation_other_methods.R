### write rownames of data
adjustdata <- function(data) {
  data<-cbind(rownames(data),data)
}


#####RMSE
timepoint <- "t8"  # "t8", "t13", "t18", "t25"
path <- paste0("~/mbDriver/simulation/parameter_estimation/other_methods/", timepoint, "/theta1")
setwd(path)

load("A0_tables.Rdata")
load("r0_tables.Rdata")
load("MLRR_list.Rdata")
load("MLCRR_list.Rdata")
load("BAL_list.Rdata")
load("BVS_list.Rdata")
load("time_list.Rdata")

rmse_tables <-  list()
re_rmse_tables <- list()

invalid_iter = c(3,16,22,23,24,32,40,46,47,54,63,64,66,68,70,72,75,81,82,83,86,88,92,96)
for (j in 1:100){
  if (j %in% invalid_iter) {
    next  
  }
  A0 <- A0_tables[[j]]
  A_mlrr <- MLRR_list[[j]][,-1]
  A_mlcrr <- MLCRR_list[[j]][,-1]
  A_bal <- BAL_list[[j]][,-1]
  A_bvs <- BVS_list[[j]][,-1]
  
  r0 <- r0_tables[[j]]
  r_mlrr <- MLRR_list[[j]][,1]
  r_mlcrr <- MLCRR_list[[j]][,1]
  r_bal <- BAL_list[[j]][,1]
  r_bvs <- BVS_list[[j]][,1]
  
  ##
  #A
  rmse_A_mlrr <- sqrt(sum((A_mlrr - A0)^2) / (ncol(A0)*nrow(A0)))
  rmse_A_mlcrr <- sqrt(sum((A_mlcrr - A0)^2) / (ncol(A0)*nrow(A0)))
  rmse_A_bal <- sqrt(sum((A_bal - A0)^2) / (ncol(A0)*nrow(A0)))
  rmse_A_bvs <- sqrt(sum((A_bvs - A0)^2) / (ncol(A0)*nrow(A0)))
  
  #r
  rmse_r_mlrr  <- sqrt(sum((r_mlrr - r0)^2) / nrow(r0))
  rmse_r_mlcrr <- sqrt(sum((r_mlcrr - r0)^2) / nrow(r0))
  rmse_r_bal <- sqrt(sum((r_bal - r0)^2) / nrow(r0))
  rmse_r_bvs <- sqrt(sum((r_bvs - r0)^2) / nrow(r0))
  
  RMSE <- matrix(c(rmse_A_mlrr, rmse_A_mlcrr,rmse_A_bal,rmse_A_bvs,
                   rmse_r_mlrr, rmse_r_mlcrr, rmse_r_bal,rmse_r_bvs),nrow=4, ncol=2, byrow=FALSE)
  row.names(RMSE) <- c("MLRR","MLCRR","BAL","BVS")
  colnames(RMSE) <- c("RMSE_A","RMSE_r")
  rmse_tables[[j]] <- RMSE
  
  ##相对RMSE
  #A
  rmse_A0 <- sqrt(sum((A0)^2) / (ncol(A0)*nrow(A0)))
  #r
  rmse_r0 <- sqrt(sum((r0)^2) / nrow(r0))
  
  re_RMSE <- matrix(c(rmse_A_mlrr/rmse_A0, rmse_A_mlcrr/rmse_A0, rmse_A_bal/rmse_A0, rmse_A_bvs/rmse_A0,
                      rmse_r_mlrr/rmse_r0, rmse_r_mlcrr/rmse_r0, rmse_r_bal/rmse_r0, rmse_r_bvs/rmse_r0),nrow=4, ncol=2, byrow=FALSE)
  row.names(re_RMSE) <- c("MLRR","MLCRR","BAL","BVS")
  colnames(re_RMSE) <- c("RMSE_A","RMSE_r")
  re_rmse_tables[[j]] <- re_RMSE
}


time_list <- lapply(time_list, function(matrix) {
  if (is.null(matrix)) {
    matrix <- matrix(0, nrow = 1, ncol = 4)  
  }
  matrix
})

rmse_tables <- lapply(rmse_tables, function(matrix) {
  if (is.null(matrix)) {
    matrix <- matrix(0, nrow = 4, ncol = 2)  
  }
  matrix
})

re_rmse_tables <- lapply(re_rmse_tables, function(matrix) {
  if (is.null(matrix)) {
    matrix <- matrix(0, nrow = 4, ncol = 2)  
  }
  matrix
})

length <- 100-length(c(3,16,22,23,24,32,40,46,47,54,63,64,66,68,70,72,75,81,82,83,86,88,92,96))
mean_time_list <- t(Reduce(`+`, time_list) / length)
mean_rmse_tables <- Reduce(`+`, rmse_tables) / length
mean_re_rmse_tables <- Reduce(`+`, re_rmse_tables) / length
colnames(mean_time_list) <- "time"
row.names(mean_time_list) <- c("MLRR","MLCRR","BAL","BVS")
mean_time_list_write <- adjustdata(mean_time_list)
mean_rmse_tables_write <- adjustdata(mean_rmse_tables)
mean_re_rmse_tables_write <- adjustdata(mean_re_rmse_tables)

save(rmse_tables, file = "rmse_tables.Rdata")
save(re_rmse_tables, file = "re_rmse_tables.Rdata")
write.table(mean_rmse_tables_write,file ="mdsine_mean_rmse_tables.txt",row.names = F,col.names = T, sep = "\t",quote = F)
write.table(mean_re_rmse_tables_write,file ="mdsine_mean_re_rmse_tables.txt",row.names = F,col.names = T, sep = "\t",quote = F)
write.table(mean_time_list_write,file ="mean_time_list.txt",row.names = F,col.names = T, sep = "\t",quote = F)
