#' Generic function for information matrix.
#'
#' Returns a information matrix for a given design
#' @param design Matrix. A design, each row is a design point with weight or
#'     repetition on the last entry.
#' @param t Numeric. Number of levels of treatments.
#' @param ... Other control parameter to be passed to methods
#' @return An information matrix.
infor_design <- function(design,t,...){
  UseMethod('infor_design',t)
}


#' @rdname infor_design
infor_design.default <- function(design,t){
  stop('Please specific your model')
}

#' @rdname infor_design
infor_design.dropout <- function(design, # matrix that stands for a design
                                   t, # number of treatment levels
                                   ... # dropout mechamism
                                   ){

  # information matrix at certain point at some period

  infor_point <- function(point, # design point, must be a matrix
                          t, # treatment number
                          p, # periods
                          l # longest time stay
  ){
    X <- matrix(0,nrow = p, ncol = p + t + t)
    X[,1:p] <- diag(rep(1, p))
    for (i in 1:p) {
      y <- point[i]
      X[i, p + y] <- 1
    }
    X[-1,p + t + 1:t] <- X[-p,p + 1:t]
    M <- diag(c(rep(1, l), rep(0, p - l)))
    out <- t(X) %*% (M - M %*% matrix(1/l, nrow = p, ncol = p) %*% M) %*% X
    return(out)
  }
  # expected information matrix at certain point

  infor_drop <- function(point, # design point
                         t, # number of treatment levels
                         ... # dropout mechanism
                         ){
    drop <- list(...)[[1]]
    p <- length(point)
    out <- matrix(0, ncol = 2*t + p,nrow = 2*t + p)
    for (i in 1:p) {
      tmp <- infor_point(point, t, p, i)
      out <- out + drop[i]*tmp
    }
    return(out)
  }

  # information matrix for design
  d <- dim(design)
  p <- d[2] - 1
  n <- d[1]
  infor1 <- matrix(0, ncol = 2*t + p, nrow = 2*t + p)
  weight <- design[, p + 1, drop = F]
  design <- design[, 1:p,drop = F]
  for (j in 1:n) {
    infor1 = infor1 + weight[j]*infor_drop(design[j,,drop = F], t,...)
  }
  return(infor1)
}


#' @rdname infor_design
infor_design.interference <- function(design, # matrix that stands for a design
                                      t, # number of treatment levels
                                      ... # sigma, covariance matrix
                                      ){

  # information matrix at certain point

  infor_point <- function(point,# design point, must be a matrix
                          t, # treatment number
                          sigma # var_cov matrix
  ){
    p <- length(point)
    if (missing(sigma)) {sigma <- diag(1,p)}
    X <- matrix(0,nrow = p,ncol = t + t + t)
    for (i in 1:p) {
      y <- point[i]
      X[i,y] <- 1 # direct treatment
    }
    X[2:p,t + 1:t] <- X[1:(p - 1),1:t] # left effects
    X[1:(p - 1),t + t + 1:t] <- X[2:p,1:t] # right effects
    sig_inv <- solve(sigma)
    out <- t(X) %*% (sig_inv - sig_inv %*% matrix(1/sum(sig_inv),ncol = p, nrow = p) %*% sig_inv) %*% X
    return(out)
  }

  # information matrix for a design

    d <- dim(design)
    p <- d[2] - 1
    n <- d[1]
    paras <- list(...)
    sigma <- paras[[1]]

  infor1 <- matrix(0,ncol = 3*t,nrow = 3*t)
  weight <- design[, p + 1, drop = F]
  design <- design[, 1:p, drop = F]
  for (j in 1:n) {
    infor1 = infor1 + weight[j]*infor_point(design[j,,drop = F],t,sigma)
  }
  return(infor1)
}


#' @rdname infor_design
infor_design.proportional <- function(design, # matrix that stands for a design
                                   t, # number of treatment levels
                                   ... # initial value of parameters and covariance
                                   ){
  
  # information matrix at certain design point

  infor_point <- function(point,# design point, must be a matrix
                          t,
                          sigma,
                          lambda,
                          tau
  ){
    p <- length(point)
    tau <- matrix(tau,ncol = 1)
    X <- matrix(0,nrow = p,ncol = p + t + 1)
    pr <- diag(rep(1,p)) # period effects::Pr
    Fd <- Rd <- matrix(0,ncol = t, nrow = p)
    for (i in 1:p) {
      y <- point[i]
      Fd[i,y] <- 1 # direct treatment :: Ld
    }
    Rd[2:p,] = Fd[1:(p - 1),]
    X[,1:p] <- pr
    X[,p + 1:t] <- Fd + lambda*Rd
    X[,p + t + 1] <- Rd %*% tau

    sig_inv <- solve(sigma)

    out <- t(X) %*% (sig_inv %*% diag(1,p) - matrix(1/sum(sig_inv),ncol = p,nrow = p) %*% sig_inv) %*% X
    return(out)
  }

  # information matrix for a design


    d <- dim(design)
    p <- d[2] - 1
    n <- d[1]
  paras <- list(...)
  sigma <- paras$sigma
  lambda <- paras$lambda
  tau <- paras$tau

  infor1 <- matrix(0,ncol = 1 + t + p,nrow = 1 + t + p)

  weight <- design[, p + 1, drop = F]
  design <- design[,1:p, drop = F]
  for (j in 1:n) {
    infor1 = infor1 + weight[j]*infor_point(design[j,,drop = F],t,sigma,lambda,tau)
  }
  return(infor1)
}


#' Generate contrast matrix
#'
#' A function return a matrix of contrast for given values of configuration.
#' @keywords internal
#' @param opt Numeric. \code{opt} = 0 means D-optimal, \code{opt} = 1 means A-optimal.
#' @param t Positive Integer. The number of levels of treatments.  of one of those
#'      classes 'dropout','proportional', or 'interference'.
#' @param p Postitive Integer. The number of periods, or number of plots in a block.
generate_contrast <- function(opt,t,p) {
  UseMethod('generate_contrast',t)
}


#' @rdname generate_contrast
generate_contrast.interference <- function(opt, t,p){
  if (opt == 0) {
    g_part <- matrix(0,nrow = t - 1,ncol = 3*t)
    g_part[,1:t] <- cbind(diag(1,t - 1),-1)
  }
  if (opt == 1) {
    g_part <- matrix(0,nrow = t,ncol = 3*t)
    g_part[,1:t] <- diag(1,t)
  }
  return(g_part)
}

#' @rdname generate_contrast
generate_contrast.dropout <- function(opt,t,p){
  if (opt == 0) {
    g_part <- matrix(0,nrow = t - 1 ,ncol = 2*t + p)
    g_part[,(t + 1):(2*t)] <- cbind(diag(1,t - 1),-1)
  }
  if (opt == 1) {
    g_part <- matrix(0,nrow = t,ncol = 2*t + p)
    g_part[,(t + 1):(2*t)] <- diag(1,t)
  }
  return(g_part)
}

#' @rdname generate_contrast
generate_contrast.proportional <- function(opt,t,p) {
  if (opt == 0) {
    g_part <- matrix(0,nrow = t - 1, ncol = 1 + t + p)
    g_part[,(t + 1):(2*t)] <- cbind(diag(1,t - 1),-1)
  }
  if (opt == 1) {
    g_part <- matrix(0,nrow = t,ncol = 1 + t + p )
    g_part[,(t + 1):(2*t)] <- diag(1,t)
  }
  return(g_part)
}












