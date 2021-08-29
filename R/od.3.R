#' Optimal sample allocation calculation for three-level CRTs
#'
#' @description The optimal design of three-level
#'     cluster randomized trials (CRTs) is to choose
#'     the sample allocation that minimizes the variance of
#'     treatment effect under fixed budget and cost structure.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n}),
#'     the level-2 sample size per level-3 unit (\code{J}),
#'     and the proportion of level-3 clusters/groups to be assigned to treatment (\code{p}).
#'     This function solves the optimal \code{n}, \code{J} and/or \code{p}
#'     with and without constraints.
#'
#' @inheritParams power.3
#' @inheritParams od.4
#' @param m total budget, default is the total costs of sampling 60
#'     level-3 units across treatment conditions.
#' @param plot.by specify variance plot by \code{n}, \code{J}
#'     and/or \code{p}; default is plot.by = list(n = "n", J = "J", p = "p").
#' @param plab the plot label for p, default is "Proportion Level-3
#'     Units in Treatment: p"
#' @param verbose logical; print the values of \code{n}, \code{J}, and \code{p} if TRUE,
#'    otherwise not; default is TRUE.
#'
#' @return
#'     unconstrained or constrained optimal sample allocation
#'     (\code{n}, \code{J}, and \code{p}).
#'     The function also returns the variance of treatment effect,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.3
#'
#' @references
#'   Shen, Z., & Kelcey, B. (2020).
#'   Optimal sample allocation under unequal costs in cluster-randomized trials.
#'   Journal of Educational and Behavioral Statistics, 45(4): 446–474.
#'
#'   Shen, Z. (2019). Optimal sample allocation in multilevel experiments
#'   (Doctoral dissertation). University of Cincinnati, Cincinnati, OH.
#'
#' @examples
#' # unconstrained optimal design #---------
#'   myod1 <- od.3(icc2 = 0.2, icc3 = 0.1, r12 = 0.5, r22 = 0.5, r32 = 0.5,
#'              c1 = 1, c2 = 5, c3 = 25, c1t = 1, c2t = 50, c3t = 250,
#'              varlim = c(0.005, 0.025))
#'   myod1$out # output
#' # plots by p and J
#'   myod1 <- od.3(icc2 = 0.2, icc3 = 0.1, r12 = 0.5, r22 = 0.5, r32 = 0.5,
#'              c1 = 1, c2 = 5, c3 = 25, c1t = 1, c2t = 50, c3t = 250,
#'              varlim = c(0.005, 0.025), plot.by = list(p = 'p', J = 'J'))
#'
#' # constrained optimal design with J = 20 #---------
#'   myod2 <- od.3(icc2 = 0.2, icc3 = 0.1, r12 = 0.5, r22 = 0.5, r32 = 0.5, J = 20,
#'              c1 = 1, c2 = 5, c3 = 25, c1t = 1, c2t = 50, c3t = 250,
#'              varlim = c(0, 0.025))
#'   myod2$out
#' # relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod2)
#'   myre$re # RE = 0.53
#'
#' # constrained optimal design with p = 0.5 #---------
#'   myod3 <- od.3(icc2 = 0.2, icc3 = 0.1, r12 = 0.5, r22 = 0.5, r32 = 0.5, p = 0.5,
#'              c1 = 1, c2 = 5, c3 = 25, c1t = 1, c2t = 50, c3t = 250,
#'              varlim = c(0.005, 0.025))
#'   myod3$out
#' # relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod3)
#'   myre$re # RE = 0.84
#'
#' # constrained n, J and p, no calculation performed #---------
#'   myod4 <- od.3(icc2 = 0.2, icc3 = 0.1, r12 = 0.5, r22 = 0.5, r32 = 0.5, n = 10, J = 10, p = 0.5,
#'              c1 = 1, c2 = 5, c3 = 25, c1t = 1, c2t = 50, c3t = 250,
#'              varlim = c(0, 0.025))
#'   myod4$out
#' # relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod4)
#'   myre$re # RE = 0.61
#'
od.3 <- function(n = NULL, J = NULL, p = NULL, icc2 = NULL, icc3 = NULL, r12 = NULL, r22 = NULL,
                 r32 = NULL, c1 = NULL, c2 = NULL, c3 = NULL, c1t = NULL, c2t = NULL, c3t = NULL,
                 m = NULL, plots = TRUE, plot.by = NULL,
                 nlim = NULL, Jlim = NULL, plim = NULL, varlim = NULL,
                 nlab = NULL, Jlab = NULL, plab = NULL, varlab = NULL,
                 vartitle = NULL,verbose = TRUE, iter = 100, tol = 1e-10) {
  funName <- "od.3"
  designType <- "three-level CRTs"
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (sum(sapply(list(icc2, icc3, r12, r22, r32, c1, c2, c3, c1t, c2t, c3t),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc2', 'icc3', 'r12', 'r22', 'r32', 'c1', 'c2', 'c3',
         'c1t', 'c2t', 'c3t' must be specified")
  if (sum(sapply(list(icc2, icc3), function(x) {
    NumberCheck(x) || any(0 >= x | x >= 1)
  })) >= 1)
    stop("'icc2', 'icc3' must be numeric in (0, 1)")
    if (sum(sapply(list(r12, r22, r32), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1)
    stop("'r12', 'r22', 'r32' must be numeric in [0, 1]")
  if (sum(sapply(list(c1, c2, c3, c1t, c2t, c3t), function(x) {
    NumberCheck(x) || x < 0})) >= 1)
    stop("'c1', 'c2', 'c3', 'c1t', 'c2t', 'c3t' must be numeric in [0, inf)")
  if (!is.null(plot.by) && !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n'))")
    if (!is.numeric(iter) || iter < 2)
      stop("specified 'iter' must be numeric with iter >= 2")
  par <- list(icc2 = icc2, icc3 = icc3, r12 = r12, r22 = r22, r32 = r32,
              c1 = c1, c2 = c2, c3 = c3, c1t =c1t, c2t = c2t, c3t = c3t,
              n = n, J = J, p = p, iter = iter)
  if (is.null(n)) {
    n.expr <- quote({
      sqrt(((1 - icc2 - icc3) * (1 - r12)) /
        (icc3 * (1 - r32) * J + icc2 * (1 - r22)) *
        ((1 - p) * (c3 + c2 * J) + p * (c3t + c2t * J)) /
        ((1 - p) * c1 * J + p * c1t * J))
    })
  } else {
    n.expr <- quote(n)
  }
  if (is.null(J)) {
    J.expr <- quote({
      sqrt((n * icc2 * (1 - r22) + (1 - icc2 - icc3) * (1 - r12)) /
        (n * icc3 * (1 - r32)) *
        ((1 - p) * c3 + p * c3t) /
        ((1 - p) * (c2 + c1 * n) + p * (c2t + c1t * n)))
    })
  } else {
    J.expr <- quote(J)
  }
  if (is.null(p)) {
    p.expr <- quote({
      sqrt((c3 + c2 * J + c1 * n * J) / (c3t + c2t * J + c1t * n * J)) /
        (1 + sqrt((c3 + c2 * J + c1 * n * J) / (c3t + c2t * J + c1t * n * J)))
    })
  } else {
    p.expr <- quote(p)
  }
  if (!is.null(n)) {
    if (!is.numeric(n) || n <= 0)
      stop("constrained 'n' must be numeric with n > 0")
  } else {
    n <- sample(2:50, 1)
  }
  if (!is.null(J)) {
    if (!is.numeric(J) || J <= 0)
      stop("constrained 'J' must be numeric with J > 0")
  } else {
    J <- sample(2:50, 1)
  }
  if (!is.null(p)) {
    if (!is.numeric(p) || any(p <=0 | p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
  } else {
    p <- sqrt((c3 + c2 + c1 ) / (c3t + c2t + c1t)) / (1 + sqrt((c3 + c2 + c1 ) / (c3t + c2t + c1t)))
  }
  nn <- JJ <- pp <- NULL
  for (i in 1:iter) {
    n <- eval(n.expr); nn[i] <- n
    J <- eval(J.expr); JJ[i] <- J
    p <- eval(p.expr); pp[i] <- p
  }
  if (!is.null(par$n) && !is.null(par$J) && !is.null(par$p)) {
    cat("===============================\n",
        "All of n, J and p are constrained, there is no calculation from other parameters",
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
    if (!is.null(par$p)) {
      cat("The constrained proportion of level-3 units in treatment (p) is ", p, ".\n", "\n", sep = "")
    } else {
      cat("The optimal proportion of level-3 units in treatment (p) is ", p, ".\n", "\n" ,sep = "")
    }
  }
  if (nn[iter] - nn[iter-1] <= tol && JJ[iter] - JJ[iter-1] <= tol && pp[iter] - pp[iter-1] <= tol) {
    nn <- JJ <- pp <- NULL
  } else {
    cat("===============================\n",
        "The solutions are not converged to specified tolerance,
        please specify a large numer of 'iter' to replace the default value of 100",
        ".\n===============================\n", sep = "")
  }
  m <- ifelse(!is.null(m), m, 60 * (p * (c1t * n * J + c2t * J + c3t) + (1 - p) * (c1 * n * J + c2 * J + c3)))
  var.expr <- quote({
    K <- m / ((1 - p) * (c1 * n * J + c2 * J + c3)
                   + p * (c1t * n * J + c2t * J + c3t))
    (icc3 * (1 - r32) + icc2 * (1 - r22) / J + (1 - icc2 - icc3) * (1 - r12) / (n * J))/ (p * (1 - p) * K)
  })
  Var <- eval(var.expr)
  par <- c(par, list(m = m))
  out <- list(n = n, J = J, p = p, var = Var)
  od.out <- list(funName = funName, designType = designType,
                 par = par, out = out)
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  nlim <- limFun(x = nlim, y = c(2, 50))
  Jlim <- limFun(x = Jlim, y = c(2, 50))
  plim <- limFun(x = plim, y = c(0, 1))
  varlim <- limFun(x = varlim, y = c(0, 0.05))
  labFun <- function(x, y) {
    if (!is.null(x) && length(x) == 1 && is.character(x)) {x} else {y}
  }
  nlab <- labFun(x = nlab, y = "Level-1 Sample Size: n")
  Jlab <- labFun(x = Jlab, y = "Level-2 Sample Size: J")
  plab <- labFun(x = plab, y = "Proportion Level-3 Units in Treatment: p")
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
