#' Optimal sample allocation calculation for four-level CRTs
#'
#' @description The optimal design of four-level
#'     cluster randomized trials (CRTs) is to choose
#'     the sample allocation that minimizes the variance of
#'     treatment effect under fixed budget and cost structure.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n}),
#'     the level-2 sample size per level-3 unit (\code{J}),
#'     the level-3 sample size per level-4 unit (\code{K}),
#'     and the proportion of level-4 clusters/groups to be assigned to treatment (\code{p}).
#'     This function solves the optimal \code{n}, \code{J}, \code{K} and/or \code{p}
#'     with and without constraints.
#'
#' @inheritParams power.4
#' @param m total budget, default value is the total costs of sampling 60
#'     level-4 units across treatment conditions.
#' @param plots logical, provide variance plots if TRUE, otherwise not; default value is TRUE.
#' @param plot.by specify variance plot by \code{n}, \code{J}, \code{K} and/or \code{p};
#'     default value is plot.by = list(n = "n", J = "J", K = 'K', p = "p").
#' @param nlim the plot range for n, default value is c(2, 50).
#' @param Jlim the plot range for J, default value is c(2, 50).
#' @param Klim the plot range for K, default value is c(2, 50).
#' @param plim the plot range for p, default value is c(0, 1).
#' @param varlim the plot range for variance, default value is c(0, 0.05).
#' @param nlab the plot label for \code{n},
#'     default value is "Level-1 Sample Size: n".
#' @param Jlab the plot label for \code{J},
#'     default value is "Level-2 Sample Size: J".
#' @param Klab the plot label for \code{K},
#'     default value is "Level-3 Sample Size: K".
#' @param plab the plot label for \code{p},
#'     default value is "Proportion Level-4 Units in Treatment: p".
#' @param varlab the plot label for variance,
#'     default value is "Variance".
#' @param vartitle the title of variance plot, default value is NULL.
#' @param verbose logical; print the values of \code{n}, \code{J},
#'    \code{K}, and \code{p} if TRUE,
#'    otherwise not; default value is TRUE.
#' @param iter number of iterations; default value is 100.
#' @param tol tolerance for convergence; default value is 1e-10.
#' @return
#'     unconstrained or constrained optimal sample allocation
#'     (\code{n}, \code{J}, \code{K}, and \code{p}).
#'     The function also returns the variance of treatment effect,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.4
#'
#' @references
#'   Shen, Z. (2019). Optimal sample allocation in multilevel experiments
#'   (Doctoral dissertation). University of Cincinnati, Cincinnati, OH.
#'
#' @examples
#' # unconstrained optimal design #---------
#'   myod1 <- od.4(icc2 = 0.2, icc3 = 0.1, icc4 = 0.05,
#'               r12 = 0.5, r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'               c1 = 1, c2 = 5, c3 = 25, c4 = 125,
#'               c1t = 1, c2t = 50, c3t = 250, c4t = 2500,
#'               varlim = c(0, 0.01))
#'   myod1$out # output
#' # plots by p and K
#'   myod1 <- od.4(icc2 = 0.2, icc3 = 0.1, icc4 = 0.05,
#'               r12 = 0.5, r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'               c1 = 1, c2 = 5, c3 = 25, c4 = 125,
#'               c1t = 1, c2t = 50, c3t = 250, c4t = 2500,
#'               varlim = c(0, 0.01), plot.by = list(p = 'p', K = 'K'))
#'
#' # constrained optimal design with p = 0.5 #---------
#'   myod2 <- od.4(icc2 = 0.2, icc3 = 0.1, icc4 = 0.05, p = 0.5,
#'               r12 = 0.5, r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'               c1 = 1, c2 = 5, c3 = 25, c4 = 125,
#'               c1t = 1, c2t = 50, c3t = 250, c4t = 2500,
#'               varlim = c(0, 0.01))
#'   myod2$out
#' # relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod2)
#'   myre$re # RE = 0.78
#'
#' # constrained optimal design with K = 20 #---------
#'   myod3 <- od.4(icc2 = 0.2, icc3 = 0.1, icc4 = 0.05,  K = 20,
#'               r12 = 0.5, r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'               c1 = 1, c2 = 5, c3 = 25, c4 = 125,
#'               c1t = 1, c2t = 50, c3t = 250, c4t = 2500,
#'               varlim = c(0, 0.01))
#'   myod3$out
#' # relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod3)
#'   myre$re # RE = 0.67
#'
#' # constrained n, J, K and p, no calculation performed #---------
#'   myod4 <- od.4(icc2 = 0.2, icc3 = 0.1, icc4 = 0.05,
#'               r12 = 0.5, n = 10, J = 10, K = 20, p = 0.5,
#'               r22 = 0.5, r32 = 0.5, r42 = 0.5,
#'               c1 = 1, c2 = 5, c3 = 25, c4 = 125,
#'               c1t = 1, c2t = 50, c3t = 250, c4t = 2500,
#'               varlim = c(0, 0.01))
#'   myod4$out
#' # relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod4)
#'   myre$re # RE = 0.27
#'
od.4 <- function(n = NULL, J = NULL, K = NULL, p = NULL, icc2 = NULL, icc3 = NULL, icc4 = NULL,
                 r12 = NULL, r22 = NULL, r32 = NULL, r42 = NULL,
                 c1 = NULL, c2 = NULL, c3 = NULL, c4 = NULL,
                 c1t = NULL, c2t = NULL, c3t = NULL, c4t = NULL,
                 m = NULL, plots = TRUE, plot.by = NULL,
                 nlim = NULL, Jlim = NULL, Klim = NULL, plim = NULL, varlim = NULL,
                 nlab = NULL, Jlab = NULL, Klab = NULL, plab = NULL, varlab = NULL,
                 vartitle = NULL,verbose = TRUE, iter = 100, tol = 1e-10) {
  funName <- "od.4"
  designType <- "four-level CRTs"
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (sum(sapply(list(icc2, icc3, icc4, r12, r22, r32, r42,
                      c1, c2, c3, c4, c1t, c2t, c3t, c4t),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc2', 'icc3', 'icc4', 'r12', 'r22', 'r32', 'r42',
         'c1', 'c2', 'c3', 'c4', 'c1t', 'c2t', 'c3t', and 'c4t' must be specified")
  if (sum(sapply(list(icc2, icc3, icc4), function(x) {
    NumberCheck(x) || any(0 >= x | x >= 1)
  })) >= 1)
    stop("'icc2', 'icc3', and 'icc4' must be numeric in (0, 1)")
    if (sum(sapply(list(r12, r22, r32, r42), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1)
    stop("'r12', 'r22', 'r32', and 'r42' must be numeric in [0, 1)")
  if (sum(sapply(list(c1, c2, c3, c4, c1t, c2t, c3t, c4t), function(x) {
    NumberCheck(x) || x < 0})) >= 1)
    stop("'c1', 'c2', 'c3', 'c4', 'c1t', 'c2t', 'c3t', and 'c4t' must be numeric in [0, inf)")
  if (!is.null(plot.by) && !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n'))")
  if (!is.numeric(iter) || iter < 2)
      stop("specified 'iter' must be numeric with iter >= 2")
  par <- list(icc2 = icc2, icc3 = icc3, icc4 = icc4,
              r12 = r12, r22 = r22, r32 = r32, r42 = r42,
              c1 = c1, c2 = c2, c3 = c3, c4 = c4,
              c1t =c1t, c2t = c2t, c3t = c3t, c4t = c4t,
              n = n, J = J, K = K, p = p, iter = iter)
  if (is.null(n)) {
    n.expr <- quote({
      sqrt(((1 - icc2 - icc3 - icc4) * (1 - r12)) /
        (icc4 * (1 - r42) * J * K +
               icc3 * (1 - r32) * J + icc2 * (1 - r22)) *
        ((1 - p) * (c4 + c3 * K + c2 * J * K) +
               p * (c4t + c3t * K + c2t * J * K)) /
        ((1 - p) * c1 * J * K + p * c1t * J * K))
    })
  } else {
    n.expr <- ({n})
  }
  if (is.null(J)) {
    J.expr <- quote({
      sqrt((n * icc2 * (1 - r22) + (1 - icc2 - icc3 - icc4) * (1 - r12)) /
        (n * K * icc4 * (1 - r42) + n * icc3 * (1 - r32)) *
        ((1 - p) * (c4 + c3 * K) + p * (c4t + c3t * K)) /
        ((1 - p) * (c2 * K + c1 * n * K) + p * (c2t * K + c1t * n * K)))
    })
  } else {
    J.expr <- ({J})
  }
  if (is.null(K)) {
    K.expr <- quote({
      sqrt((n * J * icc3 * (1 - r32)  + n * icc2 * (1 - r22) +
             (1 - icc2 - icc3 - icc4) * (1 - r12)) /
        (n * J * icc4 * (1 - r42)) *
        ((1 - p) * c4 + p * c4t) /
        ((1 - p) * (c3 + c2 * J + c1 * n * J) + p * (c3t + c2t * J + c1t * n * J)))
    })
  } else {
    K.expr <- ({K})
  }
  if (is.null(p)) {
    p.expr <- quote({
      sqrt((c4 + c3 * K + c2 * J * K + c1 * n * J * K) /
             (c4t + c3t * K + c2t * J * K + c1t * n * J * K)) /
        (1 + sqrt((c4 + c3 * K + c2 * J * K + c1 * n * J * K) /
                          (c4t + c3t * K + c2t * J * K + c1t * n * J * K)))
    })
  } else {
    p.expr <- ({p})
  }
  if (!is.null(n)) {
    if (!is.numeric(n) || n <= 0)
      stop("constrained 'n' must be numeric n > 0")
  } else {
    n <- sample(2:50, 1)
  }
  if (!is.null(J)) {
    if (!is.numeric(J) || J <= 0)
      stop("constrained 'J' must be numeric J > 0")
  } else {
    J <- sample(2:50, 1)
  }
  if (!is.null(K)) {
    if (!is.numeric(K) || J <= 0)
      stop("constrained 'K' must be numeric with K > 0")
  } else {
    K <- sample(2:50, 1)
  }
  if (!is.null(p)) {
    if (!is.numeric(p) || any(p <=0 | p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
  } else {
    p <- sqrt((c4 + c3 + c2 + c1 ) / (c4t + c3t + c2t + c1t)) /
      (1 + sqrt((c4 + c3 + c2 + c1 ) / (c4t + c3t + c2t + c1t)))
  }
  nn <- JJ <- pp <- KK <- NULL
  for (i in 1:iter) {
    n <- eval(n.expr); nn[i] <- n
    J <- eval(J.expr); JJ[i] <- J
    K <- eval(K.expr); KK[i] <- K
    p <- eval(p.expr); pp[i] <- p
  }
  if (!is.null(par$n) && !is.null(par$J) && !is.null(par$K) && !is.null(par$p)) {
    cat("===============================\n",
        "All of n, J, K, and p are constrained, there is no calculation from other parameters",
        ".\n===============================\n", sep = "")
  }
  if (verbose == TRUE) {
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
    if (!is.null(par$K)) {
      cat("The constrained level-3 sample size per level-4 unit (K) is ", K, ".\n", sep = "")
    } else {
      cat("The optimal level-3 sample size per level-4 unit (K) is ", K, ".\n", sep = "")
    }
    if (!is.null(par$p)) {
      cat("The constrained proportion of level-4 units in treatment (p) is ", p, ".\n", "\n", sep = "")
    } else {
      cat("The optimal proportion of level-4 units in treatment (p) is ", p, ".\n", "\n" ,sep = "")
    }
  }
  if (nn[iter] - nn[iter-1] <= tol && JJ[iter] - JJ[iter-1] <= tol &&
      KK[iter] - KK[iter-1] <= tol && pp[iter] - pp[iter-1] <= tol) {
    nn <- JJ <- pp <- KK <- NULL
  } else {
    cat("===============================\n",
        "The solutions are not converged to specified tolerance,
        please specify a large numer of 'iter' to replace the default value of 100",
        ".\n===============================\n", sep = "")
  }
  m <- ifelse(!is.null(m), m, 60 * (p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t) +
                                      (1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4)))
  var.expr <- quote({
    L <- m / ((1 - p) * (c1 * n * J * K + c2 * J * K + c3 * K + c4) +
                    p * (c1t * n * J * K + c2t * J * K + c3t * K + c4t))
    (icc4 * (1 - r42) + icc3 * (1 - r32) / K + icc2 * (1 - r22) / (J * K) +
        (1 - icc2 - icc3- icc4) * (1 - r12) / (n * J * K)) / (p * (1 - p) * L)
  })
  Var <- eval(var.expr)
  par <- c(par, list(m = m))
  out <- list(n = n, J = J, K = K, p = p, var = Var)
  od.out <- list(funName = funName, designType = designType,
                 par = par, out = out)
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  nlim <- limFun(x = nlim, y = c(2, 50))
  Jlim <- limFun(x = Jlim, y = c(2, 50))
  Klim <- limFun(x = Klim, y = c(2, 50))
  plim <- limFun(x = plim, y = c(0, 1))
  varlim <- limFun(x = varlim, y = c(0, 0.05))
  labFun <- function(x, y) {
    if (!is.null(x) && length(x) == 1 && is.character(x)) {x} else {y}
  }
  nlab <- labFun(x = nlab, y = "Level-1 Sample Size: n")
  Jlab <- labFun(x = Jlab, y = "Level-2 Sample Size: J")
  Klab <- labFun(x = Klab, y = "Level-3 Sample Size: K")
  plab <- labFun(x = plab, y = "Proportion Level-4 Units in Treatment: p")
  varlab <- labFun(x = varlab, y = "Variance")
  vartitle <- labFun(x = vartitle, y = "")
  plotbyFun <- function(x, y) {
    if (!is.null(x) && is.list(x)) {x} else {y}
  }
  plot.by <- plotbyFun(x = plot.by, y = list(n = "n", J = "J", K = 'K', p = "p"))
  nrange <- seq(nlim[1], nlim[2], by = 1)
  Jrange <- seq(Jlim[1], Jlim[2], by = 1)
  Krange <- seq(Klim[1], Klim[2], by = 1)
  prange <- seq(plim[1] + 0.05, plim[2] - 0.05, by = 0.01)
  if (length(plot.by) == 4) figure <- par(mfrow = c (2, 2))
  if (length(plot.by) == 3) figure <- par(mfrow = c (1, 3))
  if (length(plot.by) == 2) figure <- par(mfrow = c (1, 2))
  if (length(plot.by) == 1) figure <- par(mfrow = c (1, 1))
  if (plots) {
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
    if (!is.null(plot.by$K)) {
      plot.y <- NULL
      for (K in Krange)
        plot.y <- c(plot.y, eval(var.expr))
      graphics::plot(Krange, plot.y,
           type = "l", lty = 1,
           xlim = Klim, ylim = varlim,
           xlab = Klab, ylab = varlab,
           main = vartitle, col = "black")
      K <- out$K
      graphics::abline(v = K, lty = 2, col = "Blue")
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
