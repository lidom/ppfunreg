#' Make grid
#'
#' @param p Length of the grid 
#' @param rangevals Endpoints of the grid 
#' @export
make_grid <- function(p=100,rangevals=c(0,1)){seq(rangevals[1],rangevals[2],len=p)}


#' Make discretized covariance function
#'
#' @param cov.f Covariance function
#' @param grid Evaluation-grid 
#' @param cov.f.params Parameters of the covariance function
#' @export
make_cov_m <- function(cov.f=covf.st.matern, grid, cov.f.params=NULL){  ### Make cov. matrix from cov. function.
  grid.size <- length(grid)
  cov.m     <- matrix(0,nrow=grid.size,ncol=grid.size)
  if (is.null(cov.f.params)) {
    for (i in (1:grid.size)){
      cov.m[i,]=sapply(grid, cov.f, x1=grid[i])
    }
  }
  else{
    for (i in (1:grid.size)){
      cov.m[i,]=sapply(grid, cov.f, x1=grid[i], params=cov.f.params)
    }
  }
  return(cov.m)
}

#' Make sample (for simulation)
#'
#' @param mean.v Mean-Vector (discretized mean function)
#' @param cov.m Covariance-Matrix (discretized covariance function)
#' @param N Number of functions
#' @param dist Distribution
#' @param ... further parameters for dist
#' @export
make_sample <- function(mean.v,cov.m,N,dist="rnorm",...){
  p <- length(mean.v)
  if (p != dim(cov.m)[1] | p != dim(cov.m)[2]) stop("Dimensions of mean vector and cov. matrix do not match")
  dist        <- get(dist, mode <- "function", envir <- parent.frame())
  Z           <- matrix(dist(N*p,...),nrow=p,ncol=N)
  eigen.cov.m <- eigen(cov.m);
  eigen.cov.m$values[eigen.cov.m$values<0] <- 0
  X           <- crossprod( t(eigen.cov.m$vectors), crossprod(diag(sqrt(eigen.cov.m$values)), Z)) + mean.v
  rownames(X) <- names(mean.v)
  colnames(X) <- paste0("x",c(1:N))
  return(X)
}


#' Matern Covariance Function
#'
#' @param x1 First argument of cov(x1,x2).
#' @param x2 Second argument of cov(x1,x2).
#' @param params Matern covariance function parameters: params=c(nu, sigma). 
#' @export
covf.st.matern <- function(x1, x2, params = c(1,1)){
  nu    <- params[1]
  sigma <- params[2]  
  l     <- 1
  d     <- sqrt(2*nu)*abs(x1-x2)/l
  if (d>0) {sigma^2 * 2^(1-nu) / gamma(nu) * d^nu * besselK(d,nu)} else {sigma^2}
}


#' Modified Matern Covariance Function (varying roughness parameter)
#'
#' @param x1 First argument of cov(x1,x2). Caution: It is assumed that 0<=x1<=1.
#' @param x2 Second argument of cov(x1,x2). Caution: It is assumed that 0<=x2<=1.
#' @param params Covariance function parameters: params=c(nu1, nu2, sigma). 
#' @export
covf.nonst.matern <- function(x1,x2,params=c(3/2, 1/2, 1)){
  nu    <- params[1] + sqrt(max(x1,x2)) * (params[2] - params[1])
  sigma <- params[3]  
  l     <- 1 
  d     <- sqrt(2*nu)*abs(x1-x2)/l
  if (d>0) {sigma^2 * 2^(1-nu) / gamma(nu) * d^nu * besselK(d,nu)} else {sigma^2}
}


#' Meanfunction (polynomial)
#'
#' @param x function argument
#' @param params Parameters: params=c(shift,scale). 
#' @examples 
#' curve(meanf_poly(x, c(0,2)), from=0, to=1, 
#' main="Meanfct Poly", ylab="",xlab="")
#' curve(meanf_poly(x, c(0,1)), from=0, to=1, 
#' lty=2, add=TRUE)
#' @export
meanf_poly <- function(x,params=c(0,1)){ 
  f <- params[1] + params[2]*(10*x^3 - 15*x^4 + 6*x^5)
  names(f) <- x
  return(f)
}



