#' Estimation function
#'
#' @param Y Outcome. A matrix with ncol(Y) = 'sample size' and nrow(Y) = 'number of grid 'grd' points' 
#' @param X Predictor. A matrix with ncol(X) = 'sample size' and nrow(X) = 'number of grid 'grd' points' 
#' @param grd Grid between 0 and 1.
#' @param rho Smoothing parameter. If left unspecified (rho = NULL), then rho is 
#' determined by Generalized Cross Validation (GCV).
#' @param rho_rng The range c(min(rho_rng), max(rho_rng)) is used for finding 
#' the GCV-optimal smoothing parameter rho, if rho = NULL. 
#' @export
ppfunreg <- function(Y, X, grd, rho = NULL, rho_rng = c(0, 100)){# Y=Y_mat; X=X; grd=grid; rho = NULL; rho_rng = c(0, 100)
  ##
  a       <- base::min(grd)
  b       <- base::max(grd)
  grd     <- (grd - a)/(b-a) # standardize grid to [0,1]
  ##
  Y_orig  <- Y
  X_orig  <- X
  n       <- base::ncol(Y)
  p       <- base::nrow(Y)
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
  
  ## Need to start at t=3 (cubic splines)
  beta_hat_t <- base::sapply(3:p, FUN=function(t){
    tmp <- beta_fun_estim(Y         = Y[t,],
                          X         = delta_ar[1:t,,t],
                          grd       = grd[1:t], 
                          rho       = rho,
                          rho_rng   = rho_rng)# t=50;Y=Y_mat[t,];X=delta_ar[1:t,,t];grd=grd[1:t];rho=rho;rho_rng=rho_rng
    ##
    return(c(tmp$beta_hat_fun(grd), rep(NA,p-t)))})
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

  # par(mfrow=c(2,1))
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

# when beta is estimated over [0,1], then we need to scale it back by (b-a) for [a,b]. 


beta_fun_estim <- function (Y, X, grd, rho, rho_rng) {
  ##
  a            <- base::min(grd)
  b            <- base::max(grd)
  grd          <- (grd - a)/(b-a) # standardize grid to [0,1]
  ##
  p            <- length(grd)
  n            <- length(as.vector(Y))
  ##
  polymat      <- as.matrix(cbind(1,grd)) # Polynomial t^0 and t^1 basis functions evaluated at 'grd'
  P_mat        <- polymat %*% base::solve(t(polymat) %*% polymat) %*% base::t(polymat)
  ##
  B_mat        <- splines2::naturalSpline(x              = grd, 
                                          knots          = grd[-c(1,length(grd))], 
                                          Boundary.knots = grd[c(1,length(grd))], 
                                          intercept      = TRUE, 
                                          derivs         = 0)
  grd4integr   <- seq(0, 1, len=1001)
  bd2_mat      <- splines2::naturalSpline(x              = grd4integr, 
                                          knots          = grd[-c(1,length(grd))], 
                                          Boundary.knots = grd[c(1,length(grd))], 
                                          intercept      = TRUE, 
                                          derivs         = 2)
  integr_bd2sq <- t(bd2_mat) %*% bd2_mat / length(grd4integr)
  A_m_star     <- B_mat %*% solve(t(B_mat) %*% B_mat) %*% integr_bd2sq %*% solve(t(B_mat) %*% B_mat) %*% t(B_mat)
  A_m          <- P_mat + p * A_m_star
  ##
  
  if (is.null(rho)) {
    optRhoViaGCV <- function(r){
      H_r <- (1/n) * t(X) %*% solve(tcrossprod(X)/(n * p) + r * A_m) %*% X
      gcv <- ( sum((Y - H_r %*% Y)^2)/n )/( (1 - sum(diag(H_r))/n )^2 )
      return(gcv)
    }
    # optRhoViaGCV <- Vectorize(optRhoViaGCV)
    # optRhoViaGCV(.5)
    # rr <- seq(.Machine$double.eps, 1e-5, len=100)
    # yy <- optRhoViaGCV(r=rr)
    # plot(y=yy, x=rr, type="b")
    # plot(y=yy[-c(1:20)], x=rr[-c(1:20)], type="b")
    optRho <- stats::optimize(f = optRhoViaGCV, interval = rho_rng)
    rho    <- optRho$minimum
  }
  ##
  alpha_hat    <- (1/n) * solve(tcrossprod(X)/(n * p) + rho * A_m) %*% X %*% Y
  
  #alpha_hat
  
  beta_hat_fun <- function(t){
    splines2::naturalSpline(x              = t, 
                            knots          = (b-a) * grd[-c(1,length(grd))] + a, 
                            Boundary.knots = (b-a) * grd[ c(1,length(grd))] + a, 
                            intercept      = TRUE, 
                            derivs         = 0) %*% solve(t(B_mat) %*% B_mat) %*% t(B_mat) %*% alpha_hat * (b-a)
  }
  ## cbind(beta_hat_fun( (b-a)*grd +a), alpha_hat * (b-a))
  results   <- list(
    "beta_hat_fun" = beta_hat_fun,
    "rho"          = rho)
  return(results)
}

# # t=p;Y=Y[t,];X=delta_ar[1:t,,t];grd=grd[1:t];rho=rho;rho_rng=rho_rng
# test <- beta_fun_estim(Y=Y, X=X, grd=grd, rho=1e-7, rho_rng=rho_rng)
# 
# plot(y=test$beta_hat_fun(grd), x=grd, type="l")
# 





# estBetaCraKneSa <- function (Y, X, A_m, X_B, rho, rho_rng) 
# {
#   n   <- length(Y)
#   p   <- nrow(X)
# 
#   if (is.null(rho)) {
#     optRhoViaGCV <- function(r){
#       H_r <- t(X) %*% solve(tcrossprod(X, X)/(n * p^2) + (r/p) * A_m) %*% X / (n*p)
#       ( sum((Y - H_r %*% Y)^2)/n )/((1 - sum(diag(H_r))/n)^2)
#     }
#     optRho <- stats::optim(stats::quantile(rho_rng, .001), optRhoViaGCV, method = "Brent", 
#                             lower = range(rho_rng)[1], upper = range(rho_rng)[2])
#     rho <- optRho$par
#     gcv <- optRho$value
#   }
#   else {
#     gcv <- NULL
#   }
#   betaEstimates          <- projEstimatorMatrix(X = X, X_B, Y, A_m = A_m, rho, n, p)
#   betaEstimates[["rho"]] <- rho
#   betaEstimates[["gcv"]] <- gcv
#   betaEstimates
# }
# 
# 
# projEstimatorMatrix <- function (X, X_B, Y, A_m, rho, n, p) 
# {
#   XtX_1Xt   <- 1/n * chol2inv(chol(1/(n * p) * tcrossprod(X) + rho * A_m)) %*% X
#   alpha_hat <- XtX_1Xt %*% Y
#   beta_hat  <- tcrossprod(X_B %*% backsolve(chol(crossprod(X_B, 
#                                                            X_B)), diag(1, ncol(X_B)))) %*% alpha_hat[1:p]
#   list(estBeta  = beta_hat, 
#        XtX1Xt   = XtX_1Xt, 
#        estAlpha = alpha_hat)
# }
# 
# calSecDerNatSpline <- function(grd) {
  
  ## m = 2 cubic smoothing splines (Crambes, Kneip, Sarda AOS 2009)
  
  ## --------------

  
  ## Check derivatives
  # B_test       <- splines2::naturalSpline(x = grd4integr, knots = grd[-c(1,length(grd))], Boundary.knots = grd[c(1,length(grd))], intercept = TRUE, derivs = 0)
  # B2_test      <- splines2::naturalSpline(x = grd4integr, knots = grd[-c(1,length(grd))], Boundary.knots = grd[c(1,length(grd))], intercept = TRUE, derivs = 1)
  # par(mfrow=c(1,2))
  # plot(x = grd4integr[1:50], y= B_test[,1][1:50], type="l")
  # plot(x = grd4integr[1:50], y=B2_test[,1][1:50], type="l")
  ## check interpol
  # b        <- splines2::naturalSpline(x = grd4integr, knots = grd[-c(1,length(grd))], Boundary.knots = grd[c(1,length(grd))], intercept = TRUE, derivs = 0)
  # yy       <- b %*% solve(t(B_mat) %*% B_mat) %*% t(B_mat)  %*% 1:101
  # plot(x = grd4integr, y=yy, type="l")
  # points(x = grd, y=1:101)
  ## --------------
  
  
  ## ORIGINAL CODE:
  # p <- length(grd)
  # if (any(!(range(grd) != c(0, 1)))) 
  #   grd <- seq(0, 1, len = length(grd))
  # order <- 4
  # grid.len.2d <- 10 * p
  # polynom <- paste(paste("grd^", c(0:1), sep = ""), collapse = " , ")
  # polynom <- paste("cbind(", polynom, ")", sep = "")
  # polymat <- eval(parse(text = polynom))
  # P_mat <- polymat %*% solve(t(polymat) %*% polymat) %*% t(polymat)
  # grd_2nd <- seq(0, 1, length = grid.len.2d)[1:(grid.len.2d - 
  #                                                 1)] + 1/(2 * grid.len.2d)
  # Aknots <- c(rep(range(grd), order), grd[2:(length(grd) - 
  #                                              1)])
  # B <- splineDesign(knots = Aknots, x = grd, ord = order, derivs = 0)
  # B_d2 <- splineDesign(knots = Aknots, x = grd_2nd, ord = order, 
  #                      derivs = rep(2, length(grd_2nd)))
  # C <- splineDesign(knots = Aknots, x = range(grd_2nd), ord = order, 
  #                   derivs = c(2, 2))
  # qr.c <- qr(t(C))
  # X_B <- as.matrix((t(qr.qty(qr.c, t(B))))[, -(1:2), drop = FALSE])
  # BB_d2 <- as.matrix((t(qr.qty(qr.c, t(B_d2))))[, -(1:2), drop = FALSE])
  # G_mat <- t(BB_d2) %*% BB_d2 * 1/length(grd_2nd)
  # btb <- chol2inv(chol(crossprod(X_B, X_B)))
  # A_star_m <- X_B %*% btb %*% G_mat %*% btb %*% t(X_B)
  # A_m <- P_mat + p * A_star_m
  # list(A_m = A_m, X_B = X_B)
#}



#' Estimation function for the classic function-on-function linear regression model 
#'
#' @param Y Outcome
#' @param X Predictor 
#' @param grd Grid
#' @param rho Smoothing parameter. If left unspecified (rho = NULL), then rho is 
#' determined by Generalized Cross Validation (GCV).
#' @param rho_rng The range c(min(rho_rng), max(rho_rng)) is used for finding 
#' the GCV-optimal smoothing parameter rho, if rho = NULL.
#' @export
ffreg <- function(Y, X, grd, rho = NULL, rho_rng = c(1e-10, 10)){# Y=Y; X=X; grd=grid; rho = NULL; rho_rng = NULL
  ##
  Y_orig  <- Y
  X_orig  <- X
  n       <- base::ncol(Y)
  p       <- base::nrow(Y)
  a       <- base::min(grd)
  b       <- base::max(grd)
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
  ## Need to start at t=3 (cubic splines)
  beta_hat_t <- base::sapply(3:p, FUN=function(t){
    tmp <- beta_fun_estim(Y         = Y[t,],
                          X         = X[1:t,],
                          grd       = grd[1:t], 
                          rho       = rho,
                          rho_rng   = rho_rng)
    return(c(tmp$estBeta,rep(NA,p-t)))})
  ##
  beta_hat_t <- cbind(rep(NA,p),rep(NA,p),beta_hat_t)
  ##
  result <- list("beta_hat"      = beta_hat_t,
                 "grid"          = grd,
                 "rho"           = rho)
  class(result) <- "ffreg"
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
  
  
  
  