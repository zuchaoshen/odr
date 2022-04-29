#' Optimal sample allocation calculation for three-level MRTs detecting main effects
#'
#' @description The optimal design of three-level
#'     multisite randomized trials (MRTs) is to calculate
#'     the optimal sample allocation that minimizes the variance of
#'     treatment effect under fixed budget, which is approximately the optimal
#'     sample allocation that maximizes statistical power under a fixed budget.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n}),
#'     the level-2 sample size per level-3 unit (\code{J}),
#'     and the proportion of level-2 unit to be assigned to treatment (\code{p}).
#'     This function solves the optimal \code{n}, \code{J} and/or \code{p}
#'     with and without constraints.
#'
#' @inheritParams od.4
#' @inheritParams power.3m
#' @inheritParams power.4m
#' @param m Total budget, default is the total costs of sampling 60
#'     level-3 units.
#' @param plot.by Plot the variance by \code{n}, \code{J} and/or \code{p};
#'     default value is plot.by = list(n = "n", J = "J", p = "p").
#' @param plab The plot label for \code{p},
#'     default value is "Proportion Level-2 Units in Treatment: p".
#' @param verbose Logical; print the values of \code{n}, \code{J},
#'    and \code{p} if TRUE, otherwise not; default value is TRUE.
#' @return
#'     Unconstrained or constrained optimal sample allocation
#'     (\code{n}, \code{J}, and \code{p}).
#'     The function also returns the variance of the treatment effect,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.3m
#'
#' @references
#'   Shen, Z., & Kelcey, B. (in press). Optimal sampling ratios in three-level
#'   multisite experiments. Journal of Research on Educational Effectiveness.
#'
#' @examples
#' # Unconstrained optimal design #---------
#'   myod1 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005))
#'   myod1$out # output
#' # Plots by p and J
#'   myod1 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005), plot.by = list(p = 'p', J = 'J'))
#'
#' # Constrained optimal design with p = 0.5 #---------
#'   myod2 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005), p = 0.5)
#'   myod2$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod2)
#'   myre$re # RE = 0.81
#'
#' # Constrained optimal design with n = 5 #---------
#'   myod3 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005), n = 5)
#'   myod3$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod3)
#'   myre$re # RE = 0.89
#'
#' # Constrained n, J and p, no calculation performed #---------
#'   myod4 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005), p = 0.5, n = 15, J = 20)
#'   myod4$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod4)
#'   myre$re # RE = 0.75
#'
od.3m <- function(n = NULL, J = NULL, p = NULL, icc2 = NULL, icc3 = NULL,
                 r12 = NULL, r22 = NULL, r32m = NULL,
                 c1 = NULL, c2 = NULL, c3 = NULL,
                 c1t = NULL, c2t = NULL, omega = NULL,
                 m = NULL, plots = TRUE, plot.by = NULL,
                 nlim = NULL, Jlim = NULL, plim = NULL, varlim = NULL,
                 nlab = NULL, Jlab = NULL, plab = NULL, varlab = NULL,
                 vartitle = NULL,verbose = TRUE, iter = 100, tol = 1e-10) {
  funName <- "od.3m"
  designType <- "three-level MRTs"
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (sum(sapply(list(icc2, icc3, r12, r22, r32m,
                      c1, c2, c3, c1t, c2t, omega),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc2', 'icc3', 'r12', 'r22', 'r32m',
         'c1', 'c2', 'c3', 'c1t', 'c2t', and 'omega' must be specified")
    if (sum(sapply(list(icc2, icc3, r12, r22, r32m, omega), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1)
    stop("'icc2', 'icc3', 'r12', 'r22', 'r32m', and 'omega' must be numeric in [0, 1]")
  if (sum(sapply(list(c1, c2, c3, c1t, c2t), function(x) {
    NumberCheck(x) || x < 0})) >= 1)
    stop("'c1', 'c2', 'c3', 'c1t', and 'c2t' must be numeric in [0, inf)")
  if (!is.null(plot.by) && !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n', J = 'J'))")
  if (!is.numeric(iter) || iter < 2)
    stop("'iter' must be numeric with iter >= 2")
  par <- list(icc2 = icc2, icc3 = icc3, r12 = r12, r22 = r22, r32m = r32m,
              c1 = c1, c2 = c2, c3 = c3,
              c1t =c1t, c2t = c2t,  omega = omega,
              n = n, J = J, p = p, iter = iter)
  if (is.null(n)) {
    n.expr <- quote({
      sqrt(((1 - icc2 - icc3) * (1 - r12)) /
        (p * (1 - p) *J * omega * (1 - r32m) + icc2 * (1 - r22)) *
        ((1 - p) * J * c2 + p * J * c2t + c3) /
       ((1 - p) * c1 * J + p * c1t * J))
    })
  } else {
    n.expr <- ({n})
  }
  if (is.null(J)) {
    J.expr <- quote({
      sqrt((n * icc2 * (1 - r22) + (1 - icc2 - icc3) * (1 - r12)) /
             (p * (1 - p) *n * omega * (1 - r32m)) *
             c3 / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t)))
    })
  } else {
    J.expr <- ({J})
  }

  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  nlim <- limFun(x = nlim, y = c(2, 50))
  Jlim <- limFun(x = Jlim, y = c(2, 50))
  plim <- limFun(x = plim, y = c(0, 1))
  varlim <- limFun(x = varlim, y = c(0, 0.05))
  if (is.null(p)) {
    p.expr <- quote({
     (n * J * omega * (1 - r32m) * p * (1 - p) + n * icc2 * (1 - r22) +
        (1 - icc2 - icc3) * (1 - r12)) * (J * (c1t * n + c2t) -
        J * (c1 * n + c2)) * p * (1 - p) -
        (1 - 2 * p) * ((1 - p) * J * (c1 * n + c2) + p * J * (c1t * n + c2t) + c3) *
        (n * icc2 * (1 - r22) + (1 - icc2 - icc3) *(1 - r12))
    })
  }
  if (!is.null(n)) {
    if (!is.numeric(n) || n <= 0)
      stop("constrained 'n' must be numeric with n > 0")
  } else {
    n <- sample(2:50, 1)
  }
  if (!is.null(J)) {
    if (!is.numeric(J) || J <= 0)
      stop("constrained 'J' must be nu meric with J > 0")
  } else {
    J <- sample(2:50, 1)
  }
  if (!is.null(p)) {
    if (!is.numeric(p) || any(p <=0 | p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p.constr <- p
  } else {
    p.constr <- NULL
    p <- stats::runif(1, min = 0, max = 1)
  }
  nn <- JJ <- pp <- NULL
  for (i in 1:iter) {
    if (is.null(p.constr)) {
      pp[i] <- stats::uniroot(function(p)
        eval(p.expr), plim)$root
      p <- pp[i]
    } else {
      pp[i] <- p
    }
    n <- eval(n.expr); nn[i] <- n
    J <- eval(J.expr); JJ[i] <- J
  }
  if (!is.null(par$n) && !is.null(par$J) && !is.null(par$p)) {
    cat("===============================\n",
        "All of n, J and p are constrained, there is no calculation from other parameters",
        ".\n===============================\n", sep = "")
  }
  if (verbose) {
    if (!is.null(par$n)) {
      cat("The constrained level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    } else {
      cat("The optimal level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    }
    if (!is.null(par$J)) {
      cat("The constrained level-2 sample size per level-3 unit (J) is ", J, ".\n", sep = "")
    } else {
      cat("The optimal level-2 sample size per level-3 unit (J) is ", J, ".\n", sep = "")
    }
    if (!is.null(par$p)) {
      cat("The constrained proportion of level-2 units in treatment (p) is ", p, ".\n", "\n", sep = "")
    } else {
      cat("The optimal proportion of level-2 units in treatment (p) is ", p, ".\n", "\n" ,sep = "")
    }
  }
  if (nn[iter] - nn[iter-1] <= tol && JJ[iter] - JJ[iter-1] <= tol &&
      pp[iter] - pp[iter-1] <= tol) {
    p <- pp[iter]
    nn <- JJ <- pp <- NULL
  } else {
    cat("===============================\n",
        "The solutions are not converged to specified tolerance,
        please specify a large numer of 'iter' to replace the default value of 100",
        ".\n===============================\n", sep = "")
  }
  m <- ifelse(!is.null(m), m, 60 * (p * (c1t * n * J + c2t * J) +
                                      (1 - p) * (c1 * n * J + c2 * J) + c3))
  var.expr <- quote({
    K <- m / ((1 - p) * (c1 * n * J + c2 * J ) +
                    p * (c1t * n * J + c2t * J) + c3)
    (omega * (1 - r32m) * n * J * p * (1 - p) + icc2 * (1 - r22) * n +
        (1 - icc2 - icc3) * (1 - r12)) / (p * (1 - p) * n * J * K)
  })
  Var <- eval(var.expr)
  par <- c(par, list(m = m))
  out <- list(n = n, J = J, p = p, var = Var)
  od.out <- list(funName = funName, designType = designType,
                 par = par, out = out)
  labFun <- function(x, y) {
    if (!is.null(x) && length(x) == 1 && is.character(x)) {x} else {y}
  }
  nlab <- labFun(x = nlab, y = "Level-1 Sample Size: n")
  Jlab <- labFun(x = Jlab, y = "Level-2 Sample Size: J")
  plab <- labFun(x = plab, y = "Proportion Level-2 Units in Treatment: p")
  varlab <- labFun(x = varlab, y = "Variance")
  vartitle <- labFun(x = vartitle, y = "")
  plotbyFun <- function(x, y) {
    if (!is.null(x) && is.list(x)) {x} else {y}
  }
  plot.by <- plotbyFun(x = plot.by, y = list(n = "n", J = "J", p = "p"))
  nrange <- seq(nlim[1], nlim[2], by = 1)
  Jrange <- seq(Jlim[1], Jlim[2], by = 1)
  prange <- seq(plim[1] + 0.05, plim[2] - 0.05, by = 0.01)
  if (length(plot.by) == 3) figure <- par(mfrow = c (1, 3))
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
    if (!is.null(plot.by$J)) {
      plot.y <- NULL
      for (J in Jrange)
        plot.y <- c(plot.y, eval(var.expr))
      graphics::plot(Jrange, plot.y,
           type = "l", lty = 1,
           xlim = Jlim, ylim = varlim,
           xlab = Jlab, ylab = varlab,
           main = vartitle, col = "black")
      J <- out$J
      graphics::abline(v = J, lty = 2, col = "Blue")
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
