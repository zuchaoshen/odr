#' Optimal sample allocation calculation for two-level MRTs detecting main effects
#'
#' @description The optimal design of two-level
#'     multisite randomized trials (MRTs) detecting main effects is to calculate
#'     the sample allocation that minimizes the variance of a
#'     treatment effect under a fixed budget, which is approximately the optimal
#'     sample allocation that maximizes statistical power under a fixed budget.
#'     The optimal design parameters include
#'     the level-one sample size per site (\code{n})
#'     and the proportion of level-one unit to be assigned to treatment (\code{p}).
#'     This function solves the optimal \code{n} and/or \code{p}
#'     with and without a constraint.
#'
#' @inheritParams od.4
#' @inheritParams power.2m
#' @inheritParams power.3m
#'
#' @param m Total budget, default is the total costs of sampling 60
#'     sites.
#' @param plot.by Plot the variance by \code{n} and/or \code{p};
#'     default value is plot.by = list(n = "n", p = "p").
#' @param plab The plot label for \code{p},
#'     default value is "Proportion Level-1 Units in Treatment: p".
#' @param verbose Logical; print the values of \code{n}
#'    and \code{p} if TRUE, otherwise not; default value is TRUE.
#' @return
#'     Unconstrained or constrained optimal sample allocation
#'     (\code{n} and \code{p}).
#'     The function also returns the variance of the treatment effect,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.2m
#' @references
#'    Shen, Z., & Kelcey, B. (in press). Optimal sample
#'    allocation in multisite randomized trials.
#'    The Journal of Experimental Education.
#'    <https://doi.org/10.1080/00220973.2020.1830361>
#'
#' @examples
#' # Unconstrained optimal design #---------
#'   myod1 <- od.2m(icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
#'               c1 = 1, c2 = 10, c1t = 10,
#'               varlim = c(0, 0.005))
#'   myod1$out # n = 20, p =0.37
#' # Plots by p
#'   myod1 <- od.2m(icc = 0.2, omega = 0.02,
#'               r12 = 0.5, r22m = 0.5,
#'               c1 = 1, c2 = 10, c1t = 10,
#'               varlim = c(0, 0.005), plot.by = list(p = 'p'))
#'
#' # Constrained optimal design with p = 0.5 #---------
#'   myod2 <- od.2m(icc = 0.2, omega = 0.02,
#'               r12 = 0.5, r22m = 0.5,
#'               c1 = 1, c2 = 10, c1t = 10,
#'               varlim = c(0, 0.005), p = 0.5)
#'   myod2$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod2)
#'   myre$re # RE = 0.86
#'
#' # Constrained optimal design with n = 5 #---------
#'   myod3 <- od.2m(icc = 0.2, omega = 0.02,
#'               r12 = 0.5, r22m = 0.5, c1 = 1, c2 = 10,
#'               c1t = 10, varlim = c(0, 0.005), n = 5)
#'   myod3$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod3)
#'   myre$re # RE = 0.79
#'
#' # Constrained n and p, no calculation performed #---------
#'   myod4 <- od.2m(icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
#'               c1 = 1, c2 = 10, c1t = 10,
#'               varlim = c(0, 0.005), p = 0.5, n = 10)
#'   myod4$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod4)
#'   myre$re # RE = 0.84
#'
od.2m <- function(n = NULL, p = NULL, icc = NULL,
                 r12 = NULL, r22m = NULL,
                 c1 = NULL, c2 = NULL,
                 c1t = NULL, omega = NULL,
                 m = NULL, plots = TRUE, plot.by = NULL,
                 nlim = NULL, plim = NULL, varlim = NULL,
                 nlab = NULL, plab = NULL, varlab = NULL,
                 vartitle = NULL, verbose = TRUE, iter = 100, tol = 1e-10) {
  funName <- "od.2m"
  designType <- "two-level MRTs"
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (sum(sapply(list(icc, r12, r22m,
                      c1, c2, c1t, omega),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc', 'r12', 'r22m',
         'c1', 'c2', 'c1t', and 'omega' must be specified")
  if (sum(sapply(list(icc), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1)
    stop("'icc' must be numeric in [0, 1]")
    if (sum(sapply(list(r12, r22m), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1)
    stop("'r12' and 'r22m' must be numeric in [0, 1]")
  if (sum(sapply(list(c1, c2, c1t), function(x) {
    NumberCheck(x)})) >= 1)
    stop("'c1', 'c2', and 'c1t' must be numeric")
  if (!is.null(plot.by) && !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n'))")
  if (!is.numeric(iter) || iter < 2)
    stop("'iter' must be numeric with iter >= 2")
  par <- list(icc = icc,
              r12 = r12, r22m = r22m,
              c1 = c1, c2 = c2, c1t = c1t, omega = omega,
              n = n, p = p, iter = iter)
  if (is.null(n)) {
    n.expr <- quote({
      sqrt(((1 - icc) * (1 - r12)) /
        (p * (1 - p) * omega * (1 - r22m)) *
        c2 / ((1 - p) * c1 + p * c1t))
    })
  } else {
    n.expr <- ({n})
  }
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  nlim <- limFun(x = nlim, y = c(2, 50))
  plim <- limFun(x = plim, y = c(0, 1))
  varlim <- limFun(x = varlim, y = c(0, 0.05))
  if (is.null(p)) {
    p.expr <- quote({
     (c1t - c1) * omega * (1 - r22m) * n^2 * p^2 * (1 - p)^2 +
        (1 - icc) * (1 - r12) * n * p^2 -
        (1 - 2 * p) * (1 - icc) * (1 - r12) * (n * c1 + c2)
    })
  }
  if (!is.null(n)) {
    if (!is.numeric(n) || n <= 0)
      stop("constrained 'n' must be numeric with n > 0")
  } else {
    n <- sample(2:50, 1)
  }
  if (!is.null(p)) {
    if (!is.numeric(p) || any(p <= 0 | p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p.constr <- p
  } else {
    p.constr <- NULL
    p <- stats::runif(1, min = 0, max = 1)
  }
  nn <- pp <- NULL
  for (i in 1:iter) {
    if (is.null(p.constr)) {
      pp[i] <- stats::uniroot(function(p)
        eval(p.expr), plim)$root
      p <- pp[i]
    } else {
      pp[i] <- p
    }
    nn[i] <- eval(n.expr); n <- nn[i]
  }
  if (!is.null(par$n) && !is.null(par$p)) {
    cat("===============================\n",
        "Both n and p are constrained, there is no calculation from other parameters",
        ".\n===============================\n", sep = "")
  }
  if (verbose) {
    if (!is.null(par$n)) {
      cat("The constrained level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    } else {
      cat("The optimal level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    }
    if (!is.null(par$p)) {
      cat("The constrained proportion of level-1 units in treatment (p) is ", p, ".\n", "\n", sep = "")
    } else {
      cat("The optimal proportion of level-1 units in treatment (p) is ", p, ".\n", "\n" ,sep = "")
    }
  }
  if (nn[iter] - nn[iter-1] <= tol && pp[iter] - pp[iter-1] <= tol) {
    p <- pp[iter]
    nn <- pp <- NULL
  } else {
    cat("===============================\n",
        "The solutions are not converged to specified tolerance,
        please specify a large numer of 'iter' to replace the default value of 100",
        ".\n===============================\n", sep = "")
  }
  m <- ifelse(!is.null(m), m, 60 * (p * c1t * n + (1 - p) * c1 * n + c2 ))
  var.expr <- quote({
    (omega * (1 - r22m) * n * p * (1 - p) + (1 - icc) * (1 - r12)) /
      (p * (1 - p) * n * (m / (p * c1t * n + (1 - p) * c1 * n + c2)))
  })
  Var <- eval(var.expr)
  par <- c(par, list(m = m))
  out <- list(n = n, p = p, var = Var)
  od.out <- list(funName = funName, designType = designType,
                 par = par, out = out)
  labFun <- function(x, y) {
    if (!is.null(x) && length(x) == 1 && is.character(x)) {x} else {y}
  }
  nlab <- labFun(x = nlab, y = "Level-One Sample Size: n")
  plab <- labFun(x = plab, y = "Proportion: p")
  varlab <- labFun(x = varlab, y = "Variance")
  vartitle <- labFun(x = vartitle, y = "")
  plotbyFun <- function(x, y) {
    if (!is.null(x) && is.list(x)) {x} else {y}
  }
  plot.by <- plotbyFun(x = plot.by, y = list(n = "n", p = "p"))
  nrange <- seq(nlim[1], nlim[2], by = 1)
  prange <- seq(plim[1] + 0.05, plim[2] - 0.05, by = 0.01)
  if (length(plot.by) == 2) figure <- par(mfrow = c (1, 2))
  if (length(plot.by) == 1) figure <- par(mfrow = c (1, 1))
  if (plots == TRUE) {
    if (!is.null(plot.by$n)) {
      plot.y <- NULL
      for (n in nrange)
        plot.y <- c(plot.y, eval(var.expr))
      graphics::plot(nrange, plot.y,
           type = "l", lty = 1,
           xlim = nlim, ylim = varlim,
           xlab = nlab, ylab = varlab,
           main = vartitle, col = "black")
      n <- out$n
      graphics::abline(v = n, lty = 2, col = "Blue")
    }
    if (!is.null(plot.by$p)) {
      plot.y <- NULL
      for (p in prange)
        plot.y <- c(plot.y, eval(var.expr))
      graphics::plot(prange, plot.y,
           type = "l", lty = 1,
           xlim = plim, ylim = varlim,
           xlab = plab, ylab = varlab,
           main = vartitle, col = "black")
      p <- out$p
      graphics::abline(v = p, lty = 2, col = "Blue")
    }
  }
  par(figure)
  return(od.out)
}
