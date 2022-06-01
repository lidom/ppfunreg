## library("usethis")
## create_package("ppfunreg")
## use_mit_license("Dominik Liebl")
## use_package("tidyfun", "ffscb")
## use_readme_md()
## use_git()




## devtools::install_github("christophrust/FunRegPoI")


## library("microbenchmark")
## library("fdapace")
library("ffscb")      # Data generating process:
library("FunRegPoI")

rm(list=ls())

set.seed(123)

# Generate a sample
p          <- 251
n          <- 600
grid       <- make_grid(p, rangevals=c(0,10))
## Mean function
mu         <- meanf_poly(grid,c(0,0)) 
names(mu)  <- grid
## Covariance function
param_X    <- c(1.1,1)
cov_X      <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=param_X)
## Generate (centered) functional random sample
X          <- make_sample(mu,cov_X,n)
X          <- apply(X, 2, function(u) u - rowMeans(X))
## Generate (centered) functional error term
param_eps  <- c(2,1)
cov_eps    <- make_cov_m(cov.f = covf.st.matern, grid=grid, cov.f.params=param_eps)
eps        <- make_sample(mu,cov_eps,n) * .0001
eps        <- apply(eps, 2, function(u) u - rowMeans(eps))
##
rowMeans(X); rowMeans(eps)
##

## Plot X and error functions
par(mfrow=c(1,2))
matplot(grid, X, type="l")
matplot(grid, eps, type="l")
par(mfrow=c(1,1))
dev.off()

## parameter function alpha
alpha_fun <- function(t) 1 + t * 5    # exp(-(1/(1-(t-0.5)^2))) 
## parameter function beta
beta_fun  <- function(t,s) 1 + t * 5 * s  # (1-t)^3 * (1-s)^2 + 1


par(mfrow=c(1,2))
plot(grid, alpha_fun(grid), type="l")
## plot of beta(t,s)
matplot(grid, cbind(beta_fun(grid, 0.5), beta_fun(grid, 0.25), beta_fun(grid, 0.75)), type="l")
par(mfrow=c(1,1))
dev.off()


## Generate dependent variable Y(t)
Y <- matrix(NA, p, n)
for(i in 1:n){ # i <- 1; t <- 10
  for(t in 1:p){
    Y[t,i] <- alpha_fun(grid[t]) * X[t,i] +
      sum(beta_fun(t = grid[t],
                   s = grid[1:t]) * X[1:t,i]) * diff(grid)[1] + eps[t,i]
  }
}

## Generate dependent variable Y(t)
## POSSIBLE ERRONEOUS:
## Y <- sapply(1:n, FUN= function(i){
##   int_beta_X <- c(0, sapply(1:(p-1), 
##                       FUN=function(t){
##                         sum(beta_fun(t = grid[t],
##                                      s = grid[1:t]) * X[1:t,i]) * diff(grid)[1]
##                       }))
##   alpha_fun(grid) * X[,i] + int_beta_X + eps[,i]
## })

matplot(grid, Y, type="l")

range(rowMeans(Y))

##
alphaTilde    <- vector("numeric", p)
alphaTilde[1] <- 0
for(t in 2:p){
  alphaTilde[t] <- sum(beta_fun(t = grid[t],
                                s = grid[1:t]) *
                       cov_X[t,1:t]/cov_X[t,t] ) * diff(grid)[1] # diff(seq(0,grid[t],len=t))[1]
}

## POSSIBLE ERRONEOUS:
## alphaTilde <- c(0, sapply(1:(p-1), 
##                      FUN=function(t){
##                        sum(beta_fun(t=grid[t],s=grid[1:t]) * 
##                              cov_X[t,1:t]/cov_X[t,t] ) * diff(grid)[1]
##                      }))

alphaStar <- alpha_fun(t = grid) + alphaTilde


matplot(cbind(alphaStar, alpha_fun(t = grid)), type="l")
dev.off()

## Y2  <- sapply(1:n, FUN=function(i){
##     int_beta_Xdelta <- c(0, sapply(1:(p-1), 
##                               FUN=function(t){
##       ## Trapezoidal rule
##       # sum(c(1,rep(2,times=(p-2)),1) * 
##       # beta_fun(t=grid[t], s=grid) * 
##       #   (X[,i] - cov_X[t,]/cov_X[t,t] * X[,i] ) )/ (p*2)
##                                 ##
##       sum(beta_fun(t=grid[t], s=grid[1:t]) *
##             (X[1:t,i] - (cov_X[t,1:t]/cov_X[t,t]) * X[t,i] ) ) * diff(grid)[1]
##     }))
##     ##
##   alphaStar * X[,i] + int_beta_Xdelta + eps[,i]
## })


## range(rowMeans(Y2))

## i <- 1
## i <- i +1
## matplot(cbind(Y[,i],Y2[,i]), type="l")


## ##############################
## Estimation
## ##############################

ppfunreg <- function(Y, X, grd){# Y=Y; X=X; grid=grid
  ##
  Y_orig  <- Y
  X_orig  <- X
  n       <- ncol(Y)
  p       <- nrow(Y)
  a       <- min(grid)
  b       <- max(grid)
  #p       <- 251
  #grid    <- seq(from = a, to = b, length.out = p)
  #t_seq   <- seq.int(from = 11, to = p, by=10)
  ##
  ## Interpolate data over grid of length p
  ## Y  <- apply(Y, 2, function(u){
  ##   stats::approx(x = seq(from = a, to = b, length.out = p_orig),
  ##                 y = u, xout = grid)$y})
  ## X  <- apply(X, 2, function(u){
  ##   stats::approx(x = seq(from = 0, to = 1, length.out = p_orig),
  ##                 y = u, xout = grid)$y}) 
  
  ## centering the data:
  mean_Y  <- rowMeans(Y)
  mean_X  <- rowMeans(X)
  Y       <- apply(Y, 2, function(u) u - mean_Y)
  X       <- apply(X, 2, function(u) u - mean_X)
  ##
  ## (cross-)covariances:
  covX_n  <- cov(t(X))
  covYX_n <- cov(t(Y),t(X))
  ##
  ## estimation of alpha star:
  alphaStar_hat        <- diag(covYX_n)/diag(covX_n)
  ## matplot(y=cbind(alphaStar, alphaStar_hat), x=grd, type="l")
  ##
  ## delta_ar <- array(0, dim = c(p,p,n))
  ## for(i in 1:n){
  ##   for(t in 1:p){
  ##     for(s in 1:p){
  ##       delta_ar[t,s,i] <- X[s,i] - (covX_n[t,s]/covX_n[t,t]) * X[t,i]
  ##     }
  ##   }
  ## }
  ## t <- 240
  ## s <- 20
  ## X[t,] %*% delta_ar[t,s,]

  ## dimensions of delta(s,i,t): c(p,n,p)
  delta_ar <- sapply(1:p, FUN=function(t){
    X - matrix(covX_n[t,]/covX_n[t,t], nrow = p, ncol=n) *
      matrix(X[t,], nrow=p, ncol=n, byrow = TRUE)
  }, simplify = "array")

  t <- 240
  s <- 20
  X[t,] %*% delta_ar[s,,t]

  t        <- 151
 
  funReg_t <- FunRegPoI(Y     = Y[t,],
                        # X_mat = delta_ar[t,1:t,],
                        X_mat = delta_ar[1:t,,t],
                        grd   = seq(0,1,len=t),
                        #grd   = grid[1:t],
                        estimator = "CKS", rho = 1e-7)
  ##
  matplot(x=grid[1:t],
          y=cbind(beta_fun(t = grid[t], s = grid[1:t]),
                  funReg_t$coefficients$betaCurve / (grid[t]) ), type="l")

  ## FunRegPoI:::CraKneSaEst
  ## FunRegPoI:::estBetaCraKneSa
  ## FunRegPoI:::projEstimatorMatrix
  
  par(mfrow=c(2,1))
  matplot(x=grid,
          y=cbind(alphaStar,
                  alphaStar_hat), type="l")
  par(mfrow=c(1,1))


  eigen_t <- lapply(t_seq, FUN = function(t){
    data_l   <- MakeFPCAInputs(IDs  = rep(1:n, each=t),
                               tVec = rep(grid[1:t],n),
                               yVec = delta_ar[1:t,,t])
    FPCA_res <- FPCA(Ly = data_l$Ly,
                     Lt = data_l$Lt,
                     optns = list("FVEthreshold"  = 0.999,
                                  "userMu"        = list("mu" = rep(0, times = t),
                                                         "t" = grid[1:t]),
                                  "dataType"      = 'Dense',
                                  "methodXi"      = 'IN',
                                  "methodSelectK" = 'FVE'))
    return(list("values" = FPCA_res$lambda,
                "efcts"  = FPCA_res$phi,
                "scores" = FPCA_res$xiEst))})
  ##
  beta_hat <- lapply(seq_len(length(t_seq)), # slct <- 1
                     FUN=function(slct){
                       t       <- t_seq[slct]
                       b_hat_t <- coef(
                         fastmatrix::ridge(Y[t,] ~ -1 + eigen_t[[slct]]$scores,
                                           method = "GCV"))
                       if(length(b_hat_t)==1){
                         beta_hat_t   <- eigen_t[[slct]]$efcts * b_hat_t
                       }else{ 
                         beta_hat_t   <- rowSums(
                           eigen_t[[slct]]$efcts * 
                           matrix(rep(b_hat_t, each = nrow(eigen_t[[slct]]$efcts)),
                                  nrow = nrow(eigen_t[[slct]]$efcts)))
                       }
                       return(list(
                         "beta_hat"= beta_hat_t,
                         "t"       = grid[t],
                         "s_seq"   = grid[1:t]
                         ))})
  ##
  alpha_hat  <- alphaStar_hat[t_seq] -
    sapply(seq_len(length(t_seq)),
           FUN=function(slct){
             t      <- t_seq[slct]
             result <- sum(beta_hat[[slct]]$beta_hat *
                           cov_X[t,1:t]/cov_X[t,t])  * diff(grid)[1]
             return(result)})
  alpha_hat <- c(0, alpha_hat)
  alpha_hat <- stats::approx(x = c(0,t_seq), y = alpha_hat, xout=grid)$y
  ##
  return(list("alphaStar_hat" = alphaStar_hat, 
              "beta_hat"      = beta_hat,
              "alpha_hat"     = alpha_hat))
}


estim_results <- ppfunreg(Y = Y, X = X)


par(mfrow=c(3,1))
matplot(cbind(alphaStar, estim_results$alphaStar_hat), type="l")
matplot(cbind(alpha_fun(t = grid), estim_results$alpha_hat), type="l")
matplot(cbind(beta_fun(s = estim_results$beta_hat[[1]]$t, t=estim_results$beta_hat[[1]]$s_seq),
              estim_results$beta_hat[[1]]$beta_hat), type="l")
par(mfrow=c(1,1))






## ppfunreg <- function(Y, X){# Y=Y; X=X
##   ##
##   Y_orig  <- Y
##   X_orig  <- X
##   n       <- ncol(Y)
##   p_orig  <- nrow(Y)
##   p       <- 251
##   grid    <- seq(from = 0, to = 1, length.out = p)
##   t_seq   <- seq.int(from = 11, to = p, by=10)
##   ##
##   ## Interpolate data over grid of length p
##   Y  <- apply(Y, 2, function(u){
##     stats::approx(x = seq(from = 0, to = 1, length.out = p_orig),
##                   y = u, xout = grid)$y})
##   X  <- apply(X, 2, function(u){
##     stats::approx(x = seq(from = 0, to = 1, length.out = p_orig),
##                   y = u, xout = grid)$y})  
##   ## centering the data:
##   mean_Y  <- rowMeans(Y)
##   mean_X  <- rowMeans(X)
##   Y       <- apply(Y, 2, function(u) u - mean_Y)
##   X       <- apply(X, 2, function(u) u - mean_X)
##   ##
##   ## (cross-)covariances:
##   covX_n  <- cov(t(X))
##   covYX_n <- cov(t(Y),t(X))
##   ##
##   ## estimation of alpha star:
##   alphaStar_hat        <- diag(covYX_n)/diag(covX_n)
##   ##
##   ## dimensions of delta(s,i,t): c(p,n,p)
##   delta_ar <- sapply(1:p, FUN=function(t){
##     X - matrix(covX_n[t,]/covX_n[t,t], nrow = p, ncol=n) *
##       matrix(X[t,], nrow=p, ncol=n, byrow = TRUE)
##   }, simplify = "array")
##   ## 
##   eigen_t <- lapply(t_seq, FUN = function(t){
##     data_l   <- MakeFPCAInputs(IDs  = rep(1:n, each=t),
##                                tVec = rep(grid[1:t],n),
##                                yVec = delta_ar[1:t,,t])
##     FPCA_res <- FPCA(Ly = data_l$Ly,
##                      Lt = data_l$Lt,
##                      optns = list("FVEthreshold"  = 0.999,
##                                   "userMu"        = list("mu" = rep(0, times = t),
##                                                          "t" = grid[1:t]),
##                                   "dataType"      = 'Dense',
##                                   "methodXi"      = 'IN',
##                                   "methodSelectK" = 'FVE'))
##     return(list("values" = FPCA_res$lambda,
##                 "efcts"  = FPCA_res$phi,
##                 "scores" = FPCA_res$xiEst))})
##   ##
##   beta_hat <- lapply(seq_len(length(t_seq)), # slct <- 1
##                      FUN=function(slct){
##                        t       <- t_seq[slct]
##                        b_hat_t <- coef(
##                          fastmatrix::ridge(Y[t,] ~ -1 + eigen_t[[slct]]$scores,
##                                            method = "GCV"))
##                        if(length(b_hat_t)==1){
##                          beta_hat_t   <- eigen_t[[slct]]$efcts * b_hat_t
##                        }else{ 
##                          beta_hat_t   <- rowSums(
##                            eigen_t[[slct]]$efcts * 
##                            matrix(rep(b_hat_t, each = nrow(eigen_t[[slct]]$efcts)),
##                                   nrow = nrow(eigen_t[[slct]]$efcts)))
##                        }
##                        return(list(
##                          "beta_hat"= beta_hat_t,
##                          "t"       = grid[t],
##                          "s_seq"   = grid[1:t]
##                          ))})
##   ##
##   alpha_hat  <- alphaStar_hat[t_seq] -
##     sapply(seq_len(length(t_seq)),
##            FUN=function(slct){
##              t      <- t_seq[slct]
##              result <- sum(beta_hat[[slct]]$beta_hat *
##                            cov_X[t,1:t]/cov_X[t,t])  * diff(grid)[1]
##              return(result)})
##   alpha_hat <- c(0, alpha_hat)
##   alpha_hat <- stats::approx(x = c(0,t_seq), y = alpha_hat, xout=grid)$y
##   ##
##   return(list("alphaStar_hat" = alphaStar_hat, 
##               "beta_hat"      = beta_hat,
##               "alpha_hat"     = alpha_hat))
## }


  ## estimating the covariance function
  ## Cov(delta(u,i,t), delta(v,i,t)) separately
  ## for each 0<t<=1 [therefore the start over 2:p]
  ## with 0 <= u,v <= t
  ## eigen_t <- lapply(2:p, FUN=function(t){# slct selects t in 2:p
  ##   ## Gamma_delta_t <- tcrossprod(delta_ar[1:t,,t], delta_ar[1:t,,t]) / n
  ##   Gamma_delta_t <- cov(t(delta_ar[1:t,,t]))
  ##   tmp           <- eigen(Gamma_delta_t, symmetric = TRUE) #/ sqrt(diff(grid)[1])
  ##   values        <- tmp$values[tmp$values > 0]
  ##   efcts         <- tmp$vectors[,tmp$values > 0]
  ##   scores        <- crossprod(efcts, delta_ar[1:t,,t])
  ##   return(list("values" = values,
  ##               "efcts"  = efcts,
  ##               "scores" = scores))
  ## })# tcrossprod(eigenfunctions_t[[2]], eigenfunctions_t[[2]]) * diff(grid)[1]


  ## Gamma_delta_t <- lapply(2:p, FUN=function(t){
  ##   tcrossprod(delta_ar[1:t,,t], delta_ar[1:t,,t]) / n
  ## })
  ## eigenfunctions_t <- lapply(1:(p-1), FUN=function(slct){# slct selects t in 2:p
  ##   eigen(Gamma_delta_t[[slct]], symmetric = TRUE)$vectors #/ sqrt(diff(grid)[1])
  ## })# tcrossprod(eigenfunctions_t[[2]], eigenfunctions_t[[2]]) * diff(grid)[1]
  ## ##
  ## eigenvalues_t <- lapply(1:(p-1), FUN=function(slct){# slct selects t in 2:p
  ##   eigen(Gamma_delta_t[[slct]], symmetric = TRUE, only.values = TRUE)$values #* diff(grid)[1]
  ## })
  ## scores_t <- lapply(1:(p-1), FUN=function(slct){# slct selects t in 2:p
  ##   t <- slct + 1
  ##   ## rth row = rth eigencomponent
  ##   ## ith col = ith observation
  ##   crossprod(eigenfunctions_t[[slct]], delta_ar[1:t,,t]) #* diff(grid)[1]
  ## })




  ## GCV
  #sapply(1:(p-1), FUN=function(slct){
  ## slct selects the list-elements corresponding to t in 2:p
  ## slct <- 1
  ## t            <- slct + 1
  ##   ##
  ## Lambda_t_inv <- diag(1/(eigenvalues_t[[slct]] + kappa))
  ## b_hat_t      <- Lambda_t_inv %*% scores_t[[slct]] %*% Y[t,] #/ n
  ## cbind(Y[t,], crossprod(scores_t[[slct]], b_hat_t))
  

  ## slct <- 2
   ## kappa <- 0
  ## trace_H <- sum(diag(
  ##   t(eigen_t[[slct]]$scores) %*% diag(1/(eigen_t[[slct]]$values + kappa)) %*% eigen_t[[slct]]$scores
  ## ))
  ## (1 - trace_H/n)^2

  ## Lambda_t_inv <- diag(1/(eigen_t[[slct]]$values + kappa))

  ## solve(cov(t(eigen_t[[slct]]$scores)))

  ## b_hat_t      <- Lambda_t_inv %*% eigen_t[[slct]]$scores %*% Y[t,] / n
  ## Y_t_hat      <- c(crossprod(b_hat_t, eigen_t[[slct]]$scores))



##
## Lambda_t_inv <- diag(1/(eigen_t[[slct]]$values + kappa))
## b_hat_t      <- Lambda_t_inv %*% eigen_t[[slct]]$scores %*% Y[t,] / n
## ## computing the estimate of beta(t,.)
## beta_hat_t   <- rowSums(
##   eigen_t[[slct]]$efcts * 
##   matrix(rep(b_hat_t, each = nrow(eigen_t[[slct]]$efcts)),
##          nrow = nrow(eigen_t[[slct]]$efcts))
## )
## if(n/t >= 10){
##   b_hat_t <- coef(stats::lm(Y[t,] ~ -1 + eigen_t[[slct]]$scores))
## }else{

my.fpca <- function(Ly, Lu, reconst_fcts = NULL, pev = 0.99, CEscores = FALSE, PACE = FALSE, PACE_E0 = FALSE, center = TRUE, maxbins = NULL){
  
  n      <- length(Ly)
  id_vec <- NULL
  for(i in 1:n){id_vec <- c(id_vec, rep(i, length(Ly[[i]])))}
  ##
  ydata  <-  data.frame(".id"    = id_vec, 
                        ".index" = unname(unlist(Lu)), 
                        ".value" = unname(unlist(Ly)))
  ##
  nbasis         <- 10
  maxbins        <- ifelse(is.null(maxbins), 1000, maxbins)
  useSymm        <- FALSE # if true, there will be no smoothing accross the diagonal
  makePD         <- FALSE # if true, the NP-estimate of cov is forced to be positive definite
  ##
  # Y (nobs x length(argvals)) with NA's if no data at a argvalue
  Y        <- irreg2mat(ydata, binning = TRUE, maxbins = maxbins) 
  argvals  <- as.numeric(colnames(Y))
  ##
  # functions to be reconstructed
  if(is.null(reconst_fcts)){reconst_fcts <- 1:n}
  Y.pred   <- Y[reconst_fcts,,drop=FALSE]
  # argvals of observed fragments to be reconstructed
  argvalsO <- vector("list", length = length(reconst_fcts))
  for(i in seq_len(length(reconst_fcts))){
    minmax        <- range(as.numeric(names(c(stats::na.omit(Y.pred[i,])))))
    argvalsO[[i]] <- argvals[argvals >= minmax[1] & argvals <= minmax[2]]
  }
  ##
  D      <- NCOL(Y)      # nobs per function
  I      <- NROW(Y)      # number of functions
  I.pred <- NROW(Y.pred) # number of functions to be reconstruced
  ##
  d.vec <- rep(argvals, each = I)
  id    <- rep(1:I, rep(D, I))
  
  ## MEAN ##############################################################
  if(center){
    ## mean
    gam0    <- mgcv::gam(as.vector(Y) ~ s(d.vec, k = nbasis))
    mu      <- mgcv::predict.gam(gam0, newdata = data.frame(d.vec = argvals))
    Y.tilde <- Y - matrix(mu, I, D, byrow = TRUE)
  }else{
    Y.tilde <- Y
    mu      <- rep(0, D)
  }
  ## plot(x=argvals,y=mu, type="l")
  ##
  
  ## COV ###############################################################
  ## 1. pointwise (at argvalues) sample covariance matrix (=naive cov-estimator)
  ## 2. smooth this matrix
  cov.sum = cov.count = cov.mean = matrix(0, D, D)
  for(i in 1:I){
    obs.points = which(!is.na(Y[i, ]))
    cov.count[obs.points, obs.points] <- cov.count[obs.points, obs.points] + 1
    cov.sum[  obs.points, obs.points] <- cov.sum[  obs.points, obs.points] + tcrossprod(Y.tilde[i, obs.points])
  }
  G.0       <- ifelse(cov.count == 0, NA, cov.sum/cov.count)
  diag.G0   <- diag(G.0)
  diag(G.0) <- NA
  if(!useSymm){
    row.vec <- rep(argvals, each = D)
    col.vec <- rep(argvals, D)
    npc.0   <- matrix(mgcv::predict.gam(mgcv::gam(as.vector(G.0) ~ te(row.vec, col.vec, k = nbasis),
                                                  weights = as.vector(cov.count)), 
                                        newdata = data.frame(row.vec = row.vec, col.vec = col.vec)), D, D)
    npc.0 = (npc.0 + t(npc.0))/2
  }else{
    use          <- upper.tri(G.0, diag = TRUE)
    use[2, 1]    <- use[ncol(G.0), ncol(G.0) - 1] <- TRUE
    usecov.count <- cov.count
    usecov.count[2, 1] <- usecov.count[ncol(G.0), ncol(G.0) - 1] <- 0
    usecov.count <- as.vector(usecov.count)[use]
    use          <- as.vector(use)
    vG.0         <- as.vector(G.0)[use]
    row.vec      <- rep(argvals, each = D)[use]
    col.vec      <- rep(argvals, times = D)[use]
    mCov         <- mgcv::gam(vG.0 ~ te(row.vec, col.vec, k = nbasis), weights = usecov.count)
    npc.0        <- matrix(NA, D, D)
    spred        <- rep(argvals, each = D)[upper.tri(npc.0, diag = TRUE)]
    tpred        <- rep(argvals, times = D)[upper.tri(npc.0, diag = TRUE)]
    # Estimated covariance function:
    smVCov       <- mgcv::predict.gam(mCov, newdata = data.frame(row.vec = spred, col.vec = tpred))
    npc.0[upper.tri(npc.0, diag = TRUE)] <- smVCov
    npc.0[lower.tri(npc.0)] <- t(npc.0)[lower.tri(npc.0)]
    # slct <- seq.int(1,length(argvals),len=25)
    # persp(z=npc.0[slct,slct],x=argvals[slct],y=argvals[slct])
  }
  if(makePD){
    npc.0 <- {
      tmp <- Matrix::nearPD(npc.0, corr = FALSE, keepDiag = FALSE, do2eigen = TRUE, trace = TRUE)
      as.matrix(tmp$mat)
    }
  }
  # Nnumerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch. 8)
  w          <- quadWeights(argvals, method = "trapezoidal")
  Wsqrt      <- diag(sqrt(w))
  Winvsqrt   <- diag(1/(sqrt(w)))
  V          <- Wsqrt %*% npc.0 %*% Wsqrt
  evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  evalues    <- replace(evalues, which(evalues <= 0), 0)
  npc        <- length(evalues[evalues>0])
  npc        <- ifelse(is.null(pev), npc, which(cumsum(evalues[evalues>0])/sum(evalues[evalues>0])>=pev)[1])
  efunctions <- matrix(Winvsqrt %*% eigen(V, symmetric = TRUE)$vectors[, seq(len = npc)], nrow = nrow(V), ncol = npc)
  evalues    <- eigen(V, symmetric = TRUE, only.values = TRUE)$values[1:npc]  # use correct matrix for eigenvalue problem
  # Estimated covariance function
  cov        <- efunctions %*% tcrossprod(diag(evalues, nrow = npc, ncol = npc), efunctions)
  # Numerical integration for estimation of sigma2
  T.len      <- argvals[D] - argvals[1]  # total interval length
  T1.min     <- min(which(argvals >= argvals[1] + 0.25 * T.len))  # left bound of narrower interval T1
  T1.max     <- max(which(argvals <= argvals[D] - 0.25 * T.len))  # right bound of narrower interval T1
  DIAG       <- (diag.G0 - diag(cov))[T1.min:T1.max]  # function values
  w2         <- quadWeights(argvals[T1.min:T1.max], method = "trapezoidal")
  sigma2     <- max(stats::weighted.mean(DIAG, w = w2, na.rm = TRUE), 0)
  ##
  
  ## PACE ################################################################################
  if(PACE | PACE_E0){
    scoresP           <- vector(mode = "list", length(reconst_fcts))
    # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch.8)
    wP                <- quadWeights(argvals, method = "trapezoidal")
    WsqrtP            <- diag(sqrt(wP))
    WinvsqrtP         <- diag(1/(sqrt(wP)))
    # Cov
    VP                <- WsqrtP %*% cov %*% WsqrtP
    evalP             <- eigen(VP, symmetric = TRUE, only.values = TRUE)$values
    evalP             <- replace(evalP, which(evalP <= 0), 0)
    npcP              <- length(evalP[evalP>0])  
    npcP              <- ifelse(is.null(pev), npcP, which(cumsum(evalP[evalP>0])/sum(evalP[evalP>0])>=pev)[1])
    efunctionsP       <- matrix(WinvsqrtP %*% eigen(VP, symmetric = TRUE)$vectors[, seq(len = npcP)], nrow = nrow(VP), ncol = npcP)
    evaluesP          <- evalP[1:npcP] 
    ##
    D.invP            <- diag(1/evaluesP, nrow = npcP, ncol = npcP)
    ##
    for(i in seq_len(length(reconst_fcts))){# i <- 1
      obs_locP        <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvals))
      Y.centP         <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, D))
      ##
      if(sigma2 == 0){sigma2 <- 1e-6}
      if(PACE_E0){sigma2  <- 0}
      if(length(obs_locP) < npcP){npcP <- length(obs_locP)}
      ZcurP           <- efunctionsP[obs_locP, 1:npcP, drop=FALSE]
      ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcP, 1:npcP]), silent = TRUE)
      # while(is.error(ZtZ_sD.invP)){# for preventing singularities
      #   sigma2          <- sigma2 + .Machine$double.eps
      #   ZtZ_sD.invP     <- try(solve(crossprod(ZcurP) + sigma2 * D.invP[1:npcP, 1:npcP]), silent = TRUE)
      # }
      if(is.error(ZtZ_sD.invP)){
        ZtZ_sD.invP     <- matrix(NA, npcP, npcP)
      }
      scoresP[[i]]    <- c(ZtZ_sD.invP %*% t(ZcurP) %*% c(stats::na.omit(Y.centP)))
    }
  } else {
    efunctionsP <- NA
    evaluesP    <- NA
    scoresP     <- NA  
  }
  ## End PACE ############################################################################
  
  
  ## computations for observed fragments
  muO          <- vector("list", length(reconst_fcts))
  scoresO      <- vector("list", length(reconst_fcts))
  CE_scoresO    <- vector("list", length(reconst_fcts))
  evaluesO     <- vector("list", length(reconst_fcts))
  efunctionsO  <- vector("list", length(reconst_fcts))
  efun_reconst <- vector("list", length(reconst_fcts))
  ##
  obs_argvalsO <- vector("list", length(reconst_fcts))
  ##
  for(i in seq_len(length(reconst_fcts))){# i <- 1
    # Numerical integration for calculation of eigenvalues (see Ramsay & Silverman, Ch.8)
    w                 <- quadWeights(argvalsO[[i]], method = "trapezoidal")
    Wsqrt             <- diag(sqrt(w))
    Winvsqrt          <- diag(1/(sqrt(w)))
    locO              <- match(argvalsO[[i]],argvals)
    # CovOO
    VO                <- Wsqrt %*% cov[locO,locO] %*% Wsqrt
    evalO             <- eigen(VO, symmetric = TRUE, only.values = TRUE)$values
    evalO             <- replace(evalO, which(evalO <= 0), 0)
    npcO              <- length(evalO[evalO>0])  
    npcO              <- ifelse(is.null(pev), npcO, which(cumsum(evalO[evalO>0])/sum(evalO[evalO>0])>=pev)[1])
    efunctionsO[[i]]  <- matrix(Winvsqrt %*% eigen(VO, symmetric = TRUE)$vectors[, seq(len = npcO)], nrow = nrow(VO), ncol = npcO)
    evaluesO[[i]]     <- evalO[1:npcO]  
    ##
    D.inv             <- diag(1/evaluesO[[i]], nrow = npcO, ncol = npcO)
    Z                 <- efunctionsO[[i]]
    Y.cent            <- c(Y.pred[i,,drop=FALSE] - matrix(mu, 1, D))
    obs_locO          <- match(names(c(stats::na.omit((Y.pred[i,])))), as.character(argvalsO[[i]]))
    obs_argvalsO[[i]] <- argvalsO[[i]][obs_locO]
    ## 
    if(CEscores){
      ## CEScores (i.e., PACE-Scores)
      if(sigma2           ==   0){sigma2 <- 1e-6}
      if(length(obs_locO) < npcO){npcO <- length(obs_locO)}
      ##
      Zcur           <- Z[obs_locO,1:npcO,drop=FALSE]
      ZtZ_sD.inv     <- solve(crossprod(Zcur) + sigma2 * D.inv[1:npcO,1:npcO])
      CE_scoresO[[i]] <- c(ZtZ_sD.inv %*% t(Zcur) %*% c(stats::na.omit(Y.cent)))
    } else {
      CE_scoresO[[i]] <- NA
    }
    ## Classical scores (via intergral approximation)
    scoresO[[i]] <- apply(X      = efunctionsO[[i]][obs_locO,,drop=FALSE], 
                          MARGIN = 2, 
                          FUN    = function(ef){pracma::trapz(y=ef*c(stats::na.omit(Y.cent)),x=obs_argvalsO[[i]])})
    ##
    muO[[i]]     <- mu[locO]
    ## ##################################################################
    ## Reconstructive eigenfunctions
    efun_reconst[[i]]  <- matrix(NA, nrow=length(argvals), ncol=npcO)
    ##
    for(k in seq_len(npcO)){
      efun_reconst[[i]][,k] <- apply(X   = cov[locO,,drop=FALSE], MARGIN = 2,
                                FUN = function(x){pracma::trapz(x=argvalsO[[i]], efunctionsO[[i]][,k] * x)})
      efun_reconst[[i]][,k] <- efun_reconst[[i]][,k] / evaluesO[[i]][k]
    }
    ## ##################################################################
  }

  ## Return results ##################################################
  ret.objects <- c("Y", "mu", "muO", "cov", "sigma2",
                   "argvals",    "argvalsO", "obs_argvalsO",
                   "CE_scoresO",  "scoresO", "scoresP",
                   "efunctions", "efunctionsO", "efun_reconst", "efunctionsP", 
                   "evalues",    "evaluesO", "evaluesP")
  ret         <- lapply(1:length(ret.objects), function(u) get(ret.objects[u]))
  names(ret)  <- ret.objects
  return(ret)
}

irreg2mat <- function(ydata, binning = FALSE, maxbins = 1000){
  ##
  ydata <- ydata[stats::complete.cases(ydata), ]
  nobs  <- length(unique(ydata$.id))
  newid <- as.numeric(as.factor(ydata$.id))
  bins  <- sort(unique(ydata$.index))
  if(binning && (length(bins) > maxbins)){
    binvalues <- seq((1 - 0.001 * sign(bins[1])) * bins[1], 
                     (1 + 0.001 * sign(bins[length(bins)])) * bins[length(bins)], 
                     l = maxbins + 1)
    bins      <- binvalues
    binvalues <- utils::head(stats::filter(binvalues, c(0.5, 0.5)), -1)
  }else{
    binvalues <- bins
    bins      <- c((1 - 0.001 * sign(bins[1])) * bins[1], bins[-length(bins)], 
                   (1 + 0.001 * sign(bins[length(bins)])) * bins[length(bins)])
    if(bins[1] == 0           ){bins[1]            <- -0.001}
    if(bins[length(bins)] == 0){bins[length(bins)] <-  0.001}
  }
  newindex         <- cut(ydata$.index, breaks = bins, include.lowest = TRUE)
  Y                <- matrix(NA, nrow = nobs, ncol = nlevels(newindex))
  colnames(Y)      <- binvalues
  attr(Y, "index") <- binvalues
  ## If there are more than one data-point within a bin, 
  ## then only one of these is used (the last one).
  Y[cbind(newid, as.numeric(newindex))] <- ydata$.value
  ##
  return(Y)
}

quadWeights <- function(argvals, method = "trapezoidal"){
  ret <- switch(method, 
                trapezoidal = {D <- length(argvals); 1/2 * c(argvals[2] - argvals[1], argvals[3:D] - argvals[1:(D - 2)], argvals[D] - argvals[D - 1])}, 
                midpoint    = c(0, diff(argvals)), 
                stop("function quadWeights: choose either trapezoidal or midpoint quadrature rule")
  )
  ##
  return(ret)
}
