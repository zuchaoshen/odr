#' Optimal sample allocation calculation for two-level CRTs
#'
#' @description The optimal design of two-level
#'     cluster randomized trials (CRTs) is to choose
#'     the sample allocation that minimizes the variance of
#'     treatment effect under fixed budget and cost structure.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n})
#'     and the proportion of level-2 clusters/groups to be assigned to treatment (\code{p}).
#'     This function solves the optimal \code{n} and/or \code{p}
#'     with and without constraints.
#'
#' @inheritParams power.2
#' @inheritParams od.4
#' @param m total budget, default value is the total costs of sampling 60
#'     level-2 units across treatment conditions.
#' @param plot.by specify variance plot by \code{n} and/or \code{p}; default value is
#'     plot.by = list(n = "n", p = "p").
#' @param plab the plot label for p, default value is "Proportion Level-2 Units in Treatment: p"
#' @param verbose logical; print the values of \code{n} and \code{p} if TRUE,
#'    otherwise not; default value is TRUE.
#'
#' @return
#'     unconstrained or constrained optimal sample allocation (\code{n} and \code{p}).
#'     The function also returns the variance of treatment effect,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.2
#'
#' @references
#'   Shen, Z., & Kelcey, B. (2018, April). Optimal design of cluster
#'   randomized trials under condition- and unit-specific cost structures. Roundtable
#'   discussion presented at American Educational Research Association (AERA)
#'   annual conference, New York City, NY;
#'
#'   Shen, Z., & Kelcey, B. (2020). Optimal sample allocation under unequal
#'   costs in cluster-randomized trials. Journal of Educational
#'   and Behavioral Statistics, 45(4): 446–474.
#'
#'   Shen, Z. (2019). Optimal sample allocation in multilevel experiments
#'   (Doctoral dissertation). University of Cincinnati, Cincinnati, OH.
#'
#' @examples
#' # unconstrained optimal design #---------
#'   myod1 <- od.2(icc = 0.2, r12 = 0.5, r22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'               varlim = c(0.01, 0.02))
#'   myod1$out # output
#' # plot by p
#'   myod1 <- od.2(icc = 0.2, r12 = 0.5, r22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'               varlim = c(0.01, 0.02), plot.by = list(p = 'p'))
#'
#' # constrained optimal design with n = 20 #---------
#'   myod2 <- od.2(icc = 0.2, r12 = 0.5, r22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'               n = 20, varlim = c(0.005, 0.025))
#'   myod2$out
#' # relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod2)
#'   myre$re # RE = 0.88
#'
#' # constrained optimal design with p = 0.5 #---------
#'   myod3 <- od.2(icc = 0.2, r12 = 0.5, r22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'              p = 0.5, varlim = c(0.005, 0.025))
#'   myod3$out
#' # relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod3)
#'   myre$re # RE = 0.90
#'
#' # constrained n and p, no calculation performed #---------
#'   myod4 <- od.2(icc = 0.2, r12 = 0.5, r22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'               n = 20, p = 0.5, varlim = c(0.005, 0.025))
#'   myod4$out
#' # relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod4)
#'   myre$re # RE = 0.83
#'
od.2 <- function(n = NULL, p = NULL, icc = NULL, r12 = NULL, r22 = NULL,
                 c1 = NULL, c2 = NULL, c1t = NULL, c2t = NULL, m = NULL,
                 plots = TRUE, plot.by = NULL,
                 nlim = NULL, plim = NULL, varlim = NULL,
                 nlab = NULL, plab = NULL, varlab = NULL,
                 vartitle = NULL,verbose = TRUE) {
  funName <- "od.2"
  designType <- "two-level CRTs"
  if (sum(sapply(list(icc, r12, r22, c1, c2, c1t, c2t),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc', 'r12', 'r22', 'c1', 'c2',
         'c1t', 'c2t' must be specified")
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (NumberCheck(icc) || any(0 >= icc | icc >= 1))
    stop("'icc' must be numeric in (0, 1)")
  if (sum(sapply(list(r12, r22), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1)
    stop("'r12', 'r22' must be numeric in [0, 1)")
  if (sum(sapply(list(c1, c2, c1t, c2t), function(x) {
    NumberCheck(x) || x < 0})) >= 1)
    stop("'c1', 'c2', 'c1t', 'c2t' must be numeric in [0, inf)")
  if (!is.null(plot.by) && !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n'))")
  if (c1 == 0 && c1t == 0 && is.null(n) && is.null(p))
    stop("when c1 and c1t are both zero, one of n or p must be constrained,
         please specify a value for n or p")
  if (c2 == 0 && c2t == 0 && is.null(n) && is.null(p))
    stop("when c2 and c2t are both zero, one of n or p must be constrained,
         please specify a value for n or p")
  derivative <- quote({
    n <- sqrt((1 - icc) * (1 - r12) / (icc * (1 - r22) )) *
      sqrt(((1 - p) * c2 + p * c2t)/(( 1 - p) * c1 + p * c1t))
    p - sqrt((c1 * n + c2) / (c1t * n + c2t)) /
      (1 + sqrt((c1 * n + c2) / (c1t * n + c2t)))
  })
  par <- list(icc = icc, r12 = r12, r22 = r22, c1 = c1, c2 = c2,
              c1t =c1t, c2t = c2t, n = n, p = p)
  if (is.null(p) && is.null(n)) {
    p <- stats::uniroot(function(p) eval(derivative), interval = c(0, 1))$root
    n <- sqrt((1 - icc) * (1 - r12) / (icc * (1 - r22) )) *
         sqrt(((1 - p) * c2 + p * c2t)/(( 1 - p) * c1 + p * c1t))
  } else if (!is.null(n) && is.null(p)) {
    if (!is.numeric(n) || n <= 0)
      stop("constrained 'n' must be numeric with n > 0")
    p <- sqrt((c1 * n + c2) / (c1t * n + c2t)) /
      (1 + sqrt((c1 * n + c2)/(c1t * n + c2t)))
  } else if (!is.null(p) && is.null(n)) {
    if (!is.numeric(p) || any(p <=0 | p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    n <- sqrt(( 1 - icc) * (1 - r12) / (icc * (1 - r22))) *
         sqrt(((1 - p) * c2 + p * c2t) / ((1 - p) * c1 + p * c1t))
  } else if (!is.null(p) && !is.null(n)) {
    if (!is.numeric(n) || n <= 0)
      stop("constrained 'n' must be numeric with n > 0")
    if (!is.numeric(p) || any(p <=0 | p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
      cat("===============================\n",
        "Both p and n are constrained, there is no calculation from other parameters",
        ".\n===============================\n", sep = "")
    }
  if (verbose == TRUE) {
    if (!is.null(par$n)) {
      cat("The constrained level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    } else {
      cat("The optimal level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    }
    if (!is.null(par$p)) {
      cat("The constrained proportion of level-2 units in treatment (p) is ", p, ".\n", "\n", sep = "")
    } else {
      cat("The optimal proportion of level-2 units in treatment (p) is ", p, ".\n", "\n" ,sep = "")
    }
  }
  m <- ifelse(!is.null(m), m, 60 * (p * (c1t * n + c2t) + (1 - p) * (c1 * n + c2)))
  var.expr <- quote({
    J <- m / ((1 - p) * (c1 * n + c2)
              + p * (c1t * n + c2t))
    (icc * (1 - r22) + (1 - icc) * (1 - r12) / n )/ (p * (1 - p) * J)
  })
  Var <- eval(var.expr)
  par <- c(par, list(m = m))
  out <- list(n = n, p = p, var = Var)
  od.out <- list(funName = funName, designType = designType,
                 par = par, out = out)
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  nlim <- limFun(x = nlim, y = c(2, 50))
  plim <- limFun(x = plim, y = c(0, 1))
  varlim <- limFun(x = varlim, y = c(0, 0.05))
  labFun <- function(x, y) {
    if (!is.null(x) && length(x) == 1 && is.character(x)) {x} else {y}
  }
  nlab <- labFun(x = nlab, y = "Level-1 Sample Size: n")
  plab <- labFun(x = plab, y = "Proportion Level-2 Units in Treatment: p")
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
