#' Function for objective function
#' @keywords internal
#' A function returns value of objective function
#' @param covar Matrix, a covariance matrix for parameters at certain design
#' @param opt Integer. \code{opt=0} means D-optimal, \code{opt=1} means A-optimal.
phi <- function(covar, opt)
{
  if (opt == 0)
  {
    out <- det(covar)
  } else
  {
    out <- sum(diag(covar))^opt
  }
  return(out)
}


#' Function for General Equivalence Theorem
#'
#' A function that examines general equivalence theorem and return the maximum
#' of drectional derivative as long as its associated design point.
#' @keywords internal
#' @param opt Integer. \code{opt=0} means D-optimal, \code{opt=1} means A-optimal.
#' @param space Matrix. Design space, each row is a design point.
#' @param opt_infor Matrix. Information matrix of current design.
#' @param t Integer. Number of levels of treatment.
#' @param g_part Matrix. Contrast matrix. Its column numbers must match the row
#'      number of information matrix.
#' @param infor_design Function. A function for calculating information matrix.
#' @param ... other necessary control parameters.
#' @return A vector stands for a design point with its last entry being the maximum
#'     of directional derivative.
verify_equiv <- function(opt, space, opt_infor, t, g_part, infor_design, ...)
{
  temp <- -1e+05
  d <- dim(space)
  nx <- d[1]
  np <- d[2]
  result <- matrix(-10000, 1, np + 1)
  coeff <- 1
  
  
  temp_inv <- g_part %*% MASS::ginv(opt_infor)
  if (opt == 0)
  {
    part <- t(temp_inv) %*% solve(temp_inv %*% t(g_part)) %*% temp_inv
    diff1 <- sum(diag(opt_infor %*% part))
  }
  if (opt == 1)
  {
    var_cov <- temp_inv %*% t(g_part)
    part <- t(temp_inv) %*% temp_inv
    diff1 <- sum(diag(opt_infor %*% part))
    coeff <- (1/nrow(g_part))^(1/opt) * sum(diag(var_cov))^(1/opt - 1)
  }
  for (i in 1:nx)
  {
    x <- cbind(space[i, , drop = F], 1)
    inforx <- infor_design(x, t, ...)
    diff <- (sum(diag(inforx %*% part)) - diff1) * coeff
    if (diff - temp > 0)
    {
      temp <- diff
      result[1, 1:np] <- x[, 1:np, drop = F]
      result[1, np + 1] <- temp
    }
  }
  return(result)
}

#' Implementation of Newton's method, Part 1.
#'
#' A function that optimizes weights for a given design using newton's method
#' @keywords internal
#' @inheritParams verify_equiv
#' @param w0 Numeric Vector. Initial weights.
#' @param x Matrix, the design at current stage.
#' @param ... other necessary control parameters
#' @return A numeric vector of optimized weights.
weight1 <- function(opt, x, t, g_part, w0, infor_design, ...)
{
  n <- nrow(x)
  xn <- cbind(x[n, , drop = F], 1)
  last_infor <- infor_design(xn, t, ...)
  p <- dim(last_infor)[1]
  dim_g <- nrow(g_part)
  dinfor_m <- infor_m <- NULL
  for (i in 1:(n - 1))
  {
    xi <- NULL
    xi <- cbind(x[i, , drop = F], 1)
    temp_infor <- infor_design(xi, t, ...)
    dinfor_m <- cbind(dinfor_m, temp_infor - last_infor)
    infor_m <- cbind(infor_m, temp_infor)
  }
  infor_m <- cbind(infor_m, last_infor)
  
  weight <- w0[-n]
  diff <- 1
  repli <- 1
  indic <- 1
  delta <- 1
  
  while ((diff > 1e-14) & (indic > 0.5) & (repli < 40))
  {
    
    last_infor <- infor_m[, p * (n - 1) + 1:p]
    infor <- (1 - sum(weight)) * last_infor
    for (i in 1:(n - 1))
    {
      infor <- infor + weight[i] * infor_m[, p * (i - 1) + 1:p]
    }
    
    d1w <- matrix(0, nrow = n - 1, 1)
    d2w <- diag(rep(1, n - 1))
    temp_inv <- MASS::ginv(infor)
    var_cov <- g_part %*% temp_inv %*% t(g_part)
    if (opt == 0)
    {
      inv_var <- solve(var_cov)
      tr_inv_var1 <- dim(var_cov)[1]^(1/p - 1)
      tr_inv_var2 <- dim(var_cov)[1]^(1/p - 2)
    }
    if (opt == 1)
    {
      inv_var <- diag(rep(1, dim(var_cov)[1]))
      tr_inv_var1 <- sum(diag(var_cov))^(1/p - 1)
      tr_inv_var2 <- sum(diag(var_cov))^(1/p - 2)
    }
    part1 <- g_part %*% temp_inv %*% dinfor_m
    part2 <- -g_part %*% temp_inv %*% dinfor_m[, 1:p] %*% temp_inv %*% t(g_part)
    
    if (n > 2)
    {
      for (i in 2:(n - 1))
      {
        part2 <- cbind(part2, -g_part %*% temp_inv %*% dinfor_m[, p * (i -
                                                                         1) + 1:p] %*% temp_inv %*% t(g_part))
      }
    }
    
    for (i in 1:(n - 1))
    {
      d1w[i, 1] <- sum(diag(inv_var %*% part2[, dim_g * (i - 1) + 1:dim_g]))
      for (j in 1:i)
      {
        temp1 <- inv_var %*% part1[, p * (i - 1) + 1:p] %*% temp_inv
        temp2 <- part1[, p * (j - 1) + 1:p]
        if (opt == 0)
        {
          d2w[i, j] <- sum(diag(2 * temp1 %*% t(temp2) - inv_var %*% part2[,
                                                                           dim_g * (i - 1) + 1:dim_g] %*% inv_var %*% part2[, dim_g * (j -
                                                                                                                                         1) + 1:dim_g]))
        }
        if (opt == 1)
        {
          d2w[i, j] <- 2 * sum(diag(temp1 %*% t(temp2)))
        }
        if (d2w[i, j] > 1e+16)
        {
          d2w[i, j] <- 1e+16
        }
        if (d2w[i, j] < -1e+16)
        {
          d2w[i, j] <- -1e+16
        }
        d2w[j, i] <- d2w[i, j]
      }
    }
    # print(list(d1w=d1w,d2w=d2w))
    new_weight = weight - delta * MASS::ginv(d2w) %*% d1w
    
    if ((min(new_weight) < 0) | (sum(new_weight) > 1))
    {
      if (delta > 1e-06)
      {
        delta <- delta/2
      } else
      {
        indic <- 0
      }
    } else
    {
      diff = sum((new_weight - weight)^2)
      repli <- repli + 1
      weight <- new_weight
    }
  }
  temp <- 1 - sum(new_weight)
  weight <- c(new_weight, temp)
  return(weight)
}

#' Implementation of Newton's method, part 2.
#'
#' A function that removing boundary points after newton's method.
#' @keywords internal
#' @inheritParams weight1
#' @param x Matrix, a design returned by \code{weight1}.
#' @return Matrix, an optimal design
weight2 <- function(opt, x, t, g_part, w0, infor_design, ...)
{
  new_x <- x
  n1 <- nrow(new_x)
  if (n1 > 1)
  {
    opt_weight <- matrix(0, nrow(x), 1)
    weight <- weight1(opt, new_x, t, g_part, w0, infor_design, ...)
    while ((nrow(new_x) > 1) & (min(weight) < 1e-05))
    {
      if (min(weight) < 1e-05)
      {
        d_point <- which.min(weight)
        new_x <- new_x[-d_point, , drop = F]
        new_w0 <- weight[-d_point]
      }
      if (nrow(new_x) > 1)
      {
        weight <- weight1(opt, new_x, t, g_part, new_w0, infor_design, ...)
      }
    }
  } else
  {
    weight = 1
  }
  n1 <- nrow(new_x)
  if (n1 == 1)
  {
    weight <- 1
  }
  opt_weight <- weight
  design <- cbind(new_x, opt_weight)
  return(design)
}

#' Design Generator for Three Models
#'
#' Construct optimal approximate designs as well as efficient exact designs for
#' crossover model with subject dropout, crossover model with proportional residual
#' effect, and interference model.
#'
#' @param model an model indicator, must be one of 'dropout', 'proportional', or 'interference'.
#' @param n Positive Integer, total number of observations needed.
#' @param opt Integer. optimal criterion indicator, opt = 0 means D-opt, opt = 1 means A-opt
#' @param t Positive interger,number or levels of treatment, the default coding is integer from 1 to t
#' @param p Numeric, number of periods for crossover model or number of blocks for intereference model
#' @param ... other necessary control parameters required by specific model
#'  For crossover with dropout, \code{drop}, a numeric vector of dropout mechanism
#'  For crossover proportional, \code{lambda},value of proportion cofficient in proportional model
#'  and \code{sigma}, assumed covariance matrix.
#'  For interference model, \code{sigma}, assumed covariance matrix.
#' @param max_iter a positive integer. Controls maximum iteration time of exchange. Default is 40.
#' @return A S3 object of one of classes 'dropout', 'proportional' or 'interference'.
#'     \item{model}{the model name}
#'     \item{n}{total number of observations of exact design}
#'     \item{opt}{optimal criterion}
#'     \item{t}{number of levels of treaments}
#'     \item{p}{number of periods or plots in a block}
#'     \item{...}{other inputs}
#'     \item{initial_design}{a randomly chosen design as a starting point for
#'         newton's method}
#'     \item{exact_design}{an exact design rounded from approximate design}
#'     \item{approx_design}{optimal approximate design}
#'     \item{verify_equivalence}{result of general equivalence theorem, the last
#'         entry is the value of directional derivative}
#'     \item{time}{computing time for approximate design}
#' @examples
#' # NOTE: max_iter is usually set to 40. 
#' # Here max_iter = 5 is for demenstration only.
#' # crossover dropout model
#' ## D-optimal
#' 
#' example1 <- design('dropout',10,0,3,3,drop=c(0,0,0.5), max_iter = 5)
#' summary(example1)
#' eff(example1) # efficiency from rounding
#' effLB(example1) # obtain lower bound of efficiency
#' 
#' ## A-optimal
#' design('dropout',10,1,3,3,drop=c(0,0,0.5), max_iter = 5)
#' 
#' 
#' # proportional model
#' ## D-optimal
#' design('proportional',10,0,3,3, sigma = diag(1,3),tau = matrix(sqrt(1+3),
#'     nrow=3, ncol=1),lambda = 0.2, max_iter = 5)
#' 
#' ## A-optimal
#' design('proportional',10,1,3,3, sigma = diag(1,3), tau = matrix(sqrt(1+3),
#'     nrow=3, ncol=1),lambda = 0.2, max_iter = 5)
#'
#'
#' # interference model
#' ## D-optimal
#' design('interference',10,0,3,3, sigma = diag(1,3), max_iter = 5)
#' 
#' ## A-optimal
#' design('interference',10,1,3,3, sigma = diag(1,3), max_iter = 5)
#' 
#' @seealso \code{\link{eff}}, \code{\link{effLB}}, \code{\link{summary}}
#' @export
design <- function(model = c("dropout", "proportional", "interference"), n, opt, t,
                   p, ..., max_iter = 40)
{
  t <- as.numeric(t)
  t <- structure(t,class = model)
  n <- as.numeric(n)
  opt <- as.numeric(opt)
  p <- as.numeric(p)
  
  model <- match.arg(model)
  g_part <- generate_contrast(opt, t, p)
  space <- gtools::permutations(t, p, repeats.allowed = TRUE)
  space <- space[apply(space,1,function(x){length(unique(x))>=2}),]
  index_initial <- sample(1:NROW(space), t + 1)  # initial selection
  x0 <- space[index_initial, , drop = F]  # initial design
  
  w0 <- matrix(1/(t + 1), nrow = t + 1, ncol = 1)
  
  temp <- weight2(opt, x0, t, g_part, w0, infor_design, ...)
  opt_infor <- infor_design(temp, t, ...)
  temp1 <- verify_equiv(opt, space, opt_infor, t, g_part, infor_design, ...)
  iter <- 1
  # cat('at iteration ', iter, ' equiv is','\n') print(temp1)
  tt <- system.time({
    while ((temp1[length(temp1)] > 1e-06) & (iter < max_iter + 1))
    {
      newx <- rbind(temp[, -ncol(temp), drop = F], temp1[, -ncol(temp), drop = F])
      nn <- nrow(newx)
      w0 <- 1/nrow(newx) * matrix(1, nrow = nrow(newx), ncol = 1)
      temp <- weight2(opt, newx, t, g_part, w0, infor_design, ...)
      # cat('current design is\n') print(temp)
      opt_infor <- infor_design(temp, t, ...)
      temp1 <- verify_equiv(opt, space, opt_infor, t, g_part, infor_design,
                            ...)
      iter <- iter + 1
      # cat('at iteration ', iter, ' equiv is','\n') print(temp1)
      #if (iter == max_iter)
      #{
      #  print("Maximum Iter Reached")
      #}
    }
  })
  #cat("computing time is", tt, "\n")
  # print(temp1) print(temp)
  
  # simple rounding
  opt_design_ <- temp
  weight_opt <- temp[, ncol(temp)]
  weight_opt1 <- round(weight_opt * n)
  weight_opt2 <- weight_opt * n - weight_opt1
  weight_opt3 <- matrix(0, nrow = NROW(weight_opt1), ncol = NCOL(weight_opt1))
  diffn <- n - sum(weight_opt1)
  if (abs(diffn) > 0)
  {
    for (i in 1:abs(diffn))
    {
      if (diffn > 0)
      {
        # need extra points
        tmp.ind <- which.max(weight_opt2)
        weight_opt3[tmp.ind] <- 1
        weight_opt2[tmp.ind] <- 0
      }
      if (diffn < 0)
      {
        # need remove points
        tmp.ind <- which.min(weight_opt2)
        weight_opt3[tmp.ind] <- -1
        weight_opt2[tmp.ind] <- 0
      }
    }
  }
  
  # summarying exact design and weights #
  weight_exact <- weight_opt1 + weight_opt3
  tmp.zero <- (weight_exact == 0)
  design_exact <- opt_design_[!tmp.zero, -ncol(opt_design_), drop = F]
  weight_exact <- weight_exact[!tmp.zero]
  
  
  
  
  # comprehensive rounding #
  for(i in 1:NROW(space)){
    phi.old <- phi(g_part%*%MASS::ginv(infor_design(cbind(design_exact,weight_exact),t, ...))%*%t(g_part),opt)
    for(j in 1:NROW(design_exact)){
      design_tmp <- rbind(design_exact,space[i,])
      weight_tmp <- weight_exact
      weight_tmp[j] <- weight_tmp[j] - 1
      weight_tmp <- c(weight_tmp,1)
      tmp.zero <- weight_tmp==0
      weight_tmp <- weight_tmp[!tmp.zero]
      design_tmp <- design_tmp[!tmp.zero,]
      phi.new <- phi(g_part%*%MASS::ginv(infor_design(cbind(design_tmp,weight_tmp),
                                                      t,...))%*%t(g_part),opt)
      if(phi.new - phi.old < 0){
        tmp.zero <- weight_tmp==0
        weight_exact <- weight_tmp[!tmp.zero]
        design_exact <- design_tmp[!tmp.zero,,drop = F]
      }               
    }
  } 
  
  
  exact_design <- cbind(design_exact, weight_exact)
  
  out <- list(model = model, n = n, opt = opt, t = t, p = p,..., initial_design = x0,
              exact_design = exact_design,
              approx_design = temp, verify_equivalence = temp1,
              time = tt)
  out <- structure(out, class = model)
  return(out)
}



