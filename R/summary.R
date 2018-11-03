#' Summary method for S3 object
#'
#' Return summary info for S3 object return by \code{design} function.
#' @param object A S3 object of class 'dropout', 'proportional', or 'interference'.
#' @param ... other control parameters, but usually not necessary.
#' @return A list of key info.
#'     \item{exact_design}{exact design and its repetitions}
#'     \item{approximate_design}{approximate design and its weights}
#'     \item{computing_time}{computing time for approximate design}
#' @rdname summary
#' @seealso see examples in \code{\link{design}}.
#' @export
summary.dropout <- function(object, ...) {
  ans <- list()
  exact_design <- object$exact_design
  approx_design <- object$approx_design
  colnames(exact_design)[NCOL(exact_design)] <- 'Repetitions'
  colnames(approx_design)[NCOL(exact_design)] <- 'Weights'
  Time <- object$time
  model <- object$model
  if (object$opt == 0) {criterion <- 'D-optimal'
  } else if (object$opt == 1) {
    criterion <- 'A-optimal'
  } else {criterion = 'Unknown Optimal Criterion'}

  out.title <- paste(criterion, 'designs for', model, 'model with', object$t, 'treatments',
                 object$p, 'periods', 'dropout mechanism',toString(object$drop),':',
                 sep = ' ')
  ans <- list(exact_design = exact_design, approximate_design = approx_design,
              computing_time = Time )
  cat(out.title,'\n')
  return(ans)
}

#' @rdname summary
#' @export
summary.proportional <- function(object,...) {
  ans <- list()
  exact_design <- object$exact_design
  approx_design <- object$approx_design
  colnames(exact_design)[NCOL(exact_design)] <- 'Repetitions'
  colnames(approx_design)[NCOL(exact_design)] <- 'Weights'
  Time <- object$time
  model <- object$model
  if (object$opt == 0) {criterion <- 'D-optimal'
  } else if (object$opt == 1) {
    criterion <- 'A-optimal'
  } else {criterion = 'Unknown Optimal Criterion'}
  out.title1 <- paste(criterion, 'designs for', model, 'model with',object$t,'treatments',
                 object$p, 'periods', sep = ' ')
  out.title2 <- paste('proportional parameter', object$lambda, 'and',
                      'initial treatment effects', toString(object$tau), sep = ' ')
  out.title3 <- 'and assumed variance covariance matrix:'
  ans <- list(exact_design = exact_design, approximate_design = approx_design,
              computing_time = Time )
  cat(out.title1,'\n')
  cat(out.title2,'\n')
  cat(out.title3,'\n')
  print(object$sigma)
  cat('\n')
  return(ans)
}

#' @rdname summary
#' @export
summary.interference <- function(object,...) {
  ans <- list()
  exact_design <- object$exact_design
  approx_design <- object$approx_design
  colnames(exact_design)[NCOL(exact_design)] <- 'Repetitions'
  colnames(approx_design)[NCOL(exact_design)] <- 'Weights'
  Time <- object$time
  model <- object$model
  if (object$opt == 0) {criterion <- 'D-optimal'
  } else if (object$opt == 1) {
    criterion <- 'A-optimal'
  } else {criterion = 'Unknown Optimal Criterion'}
  out.title1 <- paste(criterion, 'designs for', model, 'model with', object$t, 'treatments',
                  object$p, 'blocks')
  out.title2 <- 'and assumed variance covariance matrix:'
  ans <- list(exact_design = exact_design, approximate_design = approx_design,
              computing_time = Time )
  #cat(out.title1,'\n')
  #cat(out.title2,'\n')
  #print(object$sigma)
  cat('\n')
  return(ans)
}

