#' Calculate f
#'
#' @param X Data frame of input variables
#' @return Fitted values
#' @export
f <- function(X, t){
  Y <- c()
  X <- as.data.frame(X)
  fmla <- as.formula("X ~ s(t, m=2, k=5, bs = 'cr')")
  dat <- cbind(X, t = t)
  mod_nb <- gam(fmla, data = dat, family = nb(), method="REML")
  Y <- fitted.values(mod_nb)
  return(Y)
}

#' Calculate f deriv
#'
#' @param X Data frame of input variables
#' @return Derivative of fitted values
#' @export
f_deriv <- function(X, t){
  Y <- c()
  X <- as.data.frame(X)
  fmla <- as.formula("X ~ s(t, m=2, k=5, bs = 'cr')")
  dat <- cbind(X, t = t)
  mod_nb <- gam(fmla, data = dat, family = nb(), method="REML")
  fitted.values <- fitted.values(mod_nb)
  newdat <- data.frame(t = t)
  X0 <- predict(mod_nb, newdata = newdat, type = "lpmatrix")
  eps <- 1e-7
  t <- t + eps
  newdat <- data.frame(t = t)
  X1 <- predict(mod_nb, newdata = newdat, type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  df <- Xp %*% coef(mod_nb)
  Y <- fitted.values * df
  return(Y)
}

#' Calculate D3 Score
#'
#' @param A Matrix A
#' @param r Vector r
#' @return Data frame with driver scores
#' @export
Driver_Score <- function(A, r) {
  x_star <- -solve(A) %*% r
  p <- length(r)
  D_square <- as.data.frame(matrix(nrow = p, ncol = 1))
  for (i in 1:p) {
    Ai <- A
    Ai[-i, i] <- 0
    ri <- r
    z_star <- -solve(Ai) %*% ri
    di <- x_star - z_star
    D_square[i, ] <- sum(as.numeric(di * di))
  }
  rownames(D_square) <- row.names(A)
  D_square <- cbind(rownames(D_square), D_square)
  D_square_sort <- D_square[order(-D_square[, 2]), ]
  return(D_square_sort)
}

#' Adjust Data
#'
#' @param data Data frame to adjust
#' @return Data frame with rownames as first column
#' @export
adjustdata <- function(data) {
  data <- cbind(rownames(data), data)
  return(data)
}

