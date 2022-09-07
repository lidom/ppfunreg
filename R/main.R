#' Estimation function
#'
#' @param Y Outcome
#' @param X Predictor 
#' @param grd Grid
#' @param rho Smoothingparameter
#' @export
ppfunreg <- function(Y, X, grd, rho = NULL, rho_rng = NULL){# Y=Y; X=X; grd=grid; rho = NULL
  ##
  Y_orig  <- Y
  X_orig  <- X
  n       <- base::ncol(Y)
  p       <- base::nrow(Y)
  a       <- base::min(grd)
  b       <- base::max(grd)
  ##
  if(is.null(rho_rng)) rho_rng <- c(1e-10, 1)
  ##  
  if(p > 101){
    base::warning("Just to let you know: large 'grd' vectors make the estimation procedure slow.")
  }
  ##
  ## centering the data:
  mean_Y  <- base::rowMeans(Y)
  mean_X  <- base::rowMeans(X)
  Y       <- base::apply(Y, 2, function(u) u - mean_Y)
  X       <- base::apply(X, 2, function(u) u - mean_X)
  ##
  ## covariances:
  covX_n  <- stats::cov(t(X))
  ##
  ## estimation of alpha star:
  alphaStar_hat   <- base::rowMeans(Y *X)/diag(covX_n)
  ## matplot(y=cbind(alphaStar, alphaStar_hat), x=grd, type="l")
  ##
  ## dimensions of delta(s,i,t): c(p,n,p)
  delta_ar <- base::sapply(1:p, FUN=function(t){
    X - base::matrix(covX_n[t,]/covX_n[t,t], nrow = p, ncol=n) *
      base::matrix(X[t,], nrow=p, ncol=n, byrow = TRUE)
  }, simplify = "array")
  
  # t <- 1
  # s <- 1
  # X[t,] %*% delta_ar[s,,t]
  # 
  # t        <- 3
  
  ## Need to start at t=3 (cubic splines)
  beta_hat_t <- base::sapply(3:p, FUN=function(t){
    tmp <- FunRegPoI::FunRegPoI(Y         = Y[t,],
                                X_mat     = delta_ar[1:t,,t],
                                grd       = grd[1:t], 
                                estimator = "CKS",
                                rho_rng   = rho_rng,
                                rho       = rho)
    return(c(tmp$coefficients$betaCurve,rep(NA,p-t)))})
  ##
  # clm1       <- beta_hat_t[,1]; clm1[2:3] <- NA
  # clm2       <- beta_hat_t[,1]; clm2[3]   <- NA
  ##
  beta_hat_t <- cbind(rep(NA,p),rep(NA,p),beta_hat_t)
  #
  
  # t <- 25
  # matplot(x=grd[1:t],
  #         y=cbind(beta_fun(t = grd[t], s = grd[1:t]),
  #                 na.omit(beta_t_mat[,t])), type="l")
  
  ## integral_0^t beta_hat(t,s)_* sigmaQuotient ds
  tmp <- sapply(3:p,
         FUN=function(t){
           result <- stats::na.omit(beta_hat_t[,t])  %*%
                           covX_n[t,1:t]/covX_n[t,t]  * diff(grd)[1]
           return(result)})
  
  alpha_hat    <- alphaStar_hat - c(0,NA,tmp)
  alpha_hat[2] <- mean(alpha_hat[c(1,3)])

  #   par(mfrow=c(2,1))
  # matplot(x=grid, y=cbind(alphaStar, alphaStar_hat), type="l")
  # matplot(x=grid, y=cbind(alpha_fun(t = grid), alpha_hat), type="l")
  # par(mfrow=c(1,1))
  
  
  ##
  result <- list("alpha_hat"     = alpha_hat, 
                 "alphaStar_hat" = alphaStar_hat, 
                 "beta_hat"      = beta_hat_t,
                 "grid"          = grd)
  class(result) <- "ppfunreg"
  return(result)
  
}



#' #' Visualizes estimation result from \link{ppfunreg}.
#' #'
#' #' @param x A 'confidence_band' object, the output from \link{confidence_band} funciton.
#' #' @param center Whether to include the functional estimate or not.
#' #' @param legendx position `x' of the legend. If NULL is passed, the legend will not be drawn (However, it may be added manually)
#' #' @param legendy position `y' of the legend.
#' #' @param ... Graphical parameters to be passed/overrided. If 'center' is TRUE, the first elements of 'col', 'lwd', 'lty' will be used for the estimate and the next ones will be used for the bands, but using the same values for one pair, i.e., lower and upper bounds.
#' #' @method plot confidence_band
#' #' @export
#' plot.ppfunreg <- function(x, center=TRUE, legendx="topleft", legendy=NULL, ...){
#'   par(mfrow=c(1,3))
#'   matplot(x = grid, y=cbind(alphaStar,           estim_results$alphaStar_hat), type="l",
#'           xlab = "", ylab = "", main = "alphaStar", col=c(1,2), lty=c(1,2), lwd=c(1,2))
#'   matplot(x = grid, y=cbind(alpha_fun(t = grid), estim_results$alpha_hat),     type="l",
#'           xlab = "", ylab = "", main = "alpha", col=c(1,2), lty=c(1,2), lwd=c(1,2))
#'   slct <- c(floor(p*1/4), floor(p*2/4), floor(p*3/4), p)
#'   matplot(x=grid,
#'           y=cbind(beta_fun(t = grid[slct[1]], s = grid),
#'                   beta_fun(t = grid[slct[2]], s = grid),
#'                   beta_fun(t = grid[slct[3]], s = grid),
#'                   beta_fun(t = grid[slct[4]], s = grid),
#'                   estim_results$beta_hat[,slct]),
#'           type="l", xlab = "", ylab = "", main = "beta", col=c(1,1,1,1,2,2,2,2), lty=c(1,1,1,1,2,2,2,2))
#'   par(mfrow=c(1,1))
  
  
  
  