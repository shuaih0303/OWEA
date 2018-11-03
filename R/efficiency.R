#' Efficiency generic function
#'
#' A generic function that returns the effciency for either exact designs to approximate
#' designs or exact design to a given design
#' @param exact_design A S3 object returned by \code{design} function.
#' @param ex Matrix. Design to be compared to. Default is \code{NULL}.
#' @export
#' @return Numeric. Relatve Efficiency.
#'     \item{If \code{ex} is given}{return relative effciency by \deqn{\Phi_{example}/\Phi_{exact\_design};}}
#'     \item{If \code{ex} is missing}{return relative efficiency by \deqn{\Phi_{approx\_design}/\Phi_{exact\_design}.}}
#' @seealso see examples in \code{\link{design}}.
eff <- function(exact_design, ex = NULL){
  UseMethod('eff', exact_design)
}


#' @rdname eff
#' @export
eff.default <- function(exact_design, ex = NULL){
  stop('No methods available for the design provided.')
}

#' @rdname eff
#' @export
eff.dropout <- function(exact_design, ex = NULL){
  list_design <- exact_design
  opt <- list_design$opt
  drop <- list_design$drop
  t <- structure(list_design$t, class = list_design$model)
  # require necessary quantities
  g_part <- generate_contrast.dropout(opt,list_design$t, list_design$p)

  cov_exact_design <- g_part %*%
    MASS::ginv(infor_design(list_design$exact_design,t,drop = drop)) %*% t(g_part)

  if (is.null(ex)) {
    cov_example <- g_part %*%
      MASS::ginv(list_design$n*infor_design(list_design$approx_design,t,drop = drop)) %*% t(g_part)
  } else{
    cov_example <- g_part %*%
      MASS::ginv(infor_design(ex,t,drop = drop)) %*% t(g_part)
  }
  phi_exact_design <- phi(cov_exact_design,opt)
  phi_example <- phi(cov_example,opt)
  efficiency <- phi_example/phi_exact_design

  # optimality
  if (opt == 0) {
    opt <- 'D-optimal'
  } else if (opt == 1) {
    opt <- 'A-optimal'
  }

  out <- list(Optimal_Criterion = opt, efficiency = efficiency)
  return(out)
}


#' @rdname eff
#' @export
eff.proportional <- function(exact_design, ex = NULL){
  list_design <- exact_design
  opt <- list_design$opt
  sigma <- list_design$sigma
  lambda <- list_design$lambda
  tau <- list_design$tau
  t <- structure(list_design$t, class = list_design$model)
  # require necessary quantities
  g_part <- generate_contrast.proportional(opt,list_design$t, list_design$p)

  cov_exact_design <- g_part %*%
    MASS::ginv(infor_design(list_design$exact_design,t,sigma = sigma,tau = tau,lambda = lambda)) %*% t(g_part)

  if ( is.null(ex)) {
    cov_example <- g_part %*%
      MASS::ginv(list_design$n*infor_design(list_design$approx_design,t,sigma = sigma,tau = tau,lambda = lambda)) %*% t(g_part)
  } else{
    cov_example <- g_part %*%
      MASS::ginv(infor_design(ex,t,sigma = sigma,tau = tau,lambda = lambda)) %*% t(g_part)
  }
  phi_exact_design <- phi(cov_exact_design,opt)
  phi_example <- phi(cov_example,opt)
  efficiency <- phi_example/phi_exact_design

  # optimality
  if (opt == 0) {
    opt <- 'D-optimal'
  } else if (opt == 1) {
    opt <- 'A-optimal'
  }

  out <- list(Optimal_Criterion = opt, efficiency = efficiency)
  return(out)
}

#' @rdname eff
#' @export
eff.interference <- function(exact_design, ex = NULL){
  list_design <- exact_design
  opt <- list_design$opt
  sigma <- list_design$sigma
  t <- structure(list_design$t, class = list_design$model)
  # require necessary quantities
  g_part <- generate_contrast.interference(opt,list_design$t, list_design$p)

  cov_exact_design <- g_part %*%
    MASS::ginv(infor_design(list_design$exact_design,t,sigma = sigma)) %*% t(g_part)

  if ( is.null(ex)) {
    cov_example <- g_part %*%
      MASS::ginv(list_design$n*infor_design(list_design$approx_design,t,sigma = sigma)) %*% t(g_part)
  } else{
    cov_example <- g_part %*%
      MASS::ginv(infor_design(ex,t,sigma = sigma)) %*% t(g_part)
  }
  phi_exact_design <- phi(cov_exact_design,opt)
  phi_example <- phi(cov_example,opt)
  efficiency <- phi_example/phi_exact_design

  # optimality
  if (opt == 0) {
    opt <- 'D-optimal'
  } else if (opt == 1) {
    opt <- 'A-optimal'
  }

  out <- list(Optimal_Criterion = opt, efficiency = efficiency)
  return(out)
}




#' Lower Bound Efficiency for Crossover-Dropout Model
#'
#' The function take S3 object of class 'dropout' as input and return its lower
#' bound of efficiency of exact design.
#' @param exact_design A object of class  returned by design function.
#' @return A list of relavent numerics.
#'     \item{optimal}{ Optimal Criterion}
#'     \item{lower.bound}{ Lower Bound of the exact design }
#'     \item{optimal.value}{ The value of objective function at optimal approxiamte design}
#' @export
#' @seealso see examples in \code{\link{design}}.

effLB <- function(exact_design) {
  if (class(exact_design) != 'dropout') {
    stop("This function is only for 'dropout' class!")
  }
  # a function to calculate informatrix at point without expectation.
  infor_point <- function(point, # design point, must be a matrix
                          t, # treatment number
                          p, # periods
                          l # longest time stay
  ){
    X <- matrix(0,nrow = p,ncol = p + t + t)
    X[,1:p] <- diag(rep(1,p))
    for (i in 1:p) {
      y <- point[i]
      X[i,p + y] <- 1
    }
    X[-1,p + t + 1:t] <- X[-p,p + 1:t]
    M <- diag(c(rep(1,l),rep(0,p - l)))
    out <- t(X) %*% (M - M %*% matrix(1/l,nrow = p, ncol = p) %*% M) %*% X
    return(out)
  }


  # unpack necessary parameters
  drop <- exact_design$drop
  p <- exact_design$p
  opt <- exact_design$opt
  t <- exact_design$t
  m <- min(which(drop != 0))
  design1 <- exact_design$exact_design
  design2 <- exact_design$approx_design # approximate design
  weight1 <- design1[, ncol(design1),drop = F]

  g_part <- generate_contrast.dropout(opt, t, p)


  # generate all realizations and probs#
  n1 <- nrow(design1)
  realize1 <- gtools::permutations(n = p - m + 1, r = n1, p:m, repeats.allowed = T)
  realize.prob1 <- apply(realize1, 1, function(x) {
    out <- drop[x[1]]
    for (i in 2:length(x)) {
      out <- out * drop[x[i]]
    }
    return(out)
  })



  # calculate the numerator and denomenator of efficiency formula
  phi1 <- phi2 <- 0
  covar2 <- g_part %*% MASS::ginv(exact_design$n*infor_design(design2, t, drop)) %*% t(g_part)
  if (opt == 0) {
    for (i in 1:NROW(realize1)) {
      infor1 <- 0
      for (j in 1:NROW(design1)) {
        infor1 <- infor1 + weight1[j] * infor_point(design1[j,- p - 1], t, p, realize1[i,
                                                                                 j])
      }
      covar1 <- g_part %*% MASS::ginv(infor1) %*% t(g_part)
      phi1 <- phi1 + realize.prob1[i] * det(covar1)
    }
    phi2 <- det(covar2)
    out <- list(efficiency.self = phi2/phi1, optimal.value = phi2)
  }


  if (opt == 1) {
    for (i in 1:NROW(realize1)) {
      infor1 <- 0
      for (j in 1:NROW(design1)) {
        infor1 <- infor1 + weight1[j] * infor_point(design1[j, - p - 1], t, p, realize1[i,j])
      }
      covar1 <- g_part %*% MASS::ginv(infor1) %*% t(g_part)
      phi1 <- phi1 + realize.prob1[i] * sum(diag(covar1))
    }
    phi2 <- sum(diag(covar2))

    # output two values, first is lowerbound, second is its optimal value of the
    # approximate design
    # optimality
    if (opt == 0) {
      opt <- 'D-optimal'
    } else if (opt == 1) {
      opt <- 'A-optimal'
    }

    out <- list(optimal = opt, lower.bound = phi2/phi1, optimal.value = phi2 )
  }
  #print(out)
  return(out)
}






