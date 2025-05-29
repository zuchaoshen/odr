#' Budget and/or sample size, power, MDES calculation for
#' three-level MRTs detecting main effects
#'
#' @description This function can calculate required budget for desired power,
#'     power or minimum detectable effect size (MDES) under fixed budget
#'     for three-level multisite randomized trials (MRTs).
#'     It also can perform conventional power analyses
#'     (e.g., required sample size, power, and MDES calculation).
#'
#' @inheritParams power.4
#' @inheritParams power.3
#' @inheritParams power.4m
#'
#' @param expr Returned objects from function \code{\link{od.3m}}; default is NULL;
#'     if \code{expr} is specified, parameter values of \code{icc2}, \code{icc3},
#'     \code{r12}, \code{r22}, \code{r32m},
#'     \code{c1}, \code{c2}, \code{c3},
#'     \code{c1t}, \code{c2t}, \code{p}, \code{n}, and \code{J}
#'     used or solved in function \code{\link{od.3m}} will
#'     be passed to current function;
#'     only the values of \code{p}, \code{n}, and/or \code{J} that specified or solved in
#'     function \code{\link{od.3m}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint Specify the constrained values of \code{p}, \code{n}, and/or \code{J},
#'     in list format to overwrite
#'     those from \code{expr}; default value is NULL.
#' @param r32m The proportion of variance of site-specific treatment effect explained by covariates.
#' @param c3 The cost of sampling one level-3 unit (site).
#' @param p The proportion of level-2 units to be assigned to treatment.
#' @param q The number of covariates at level 3.
#' @param mlim The range for searching the root of budget (\code{m}) numerically,
#'     default is the costs sampling \code{Klim} level-3 units
#'     or c(4 * Kcost, 1e+10 * Kcost) with Kcost =
#'     ((1 - p) * (c1 * n * J + c2 * J) +
#'     p * (c1t * n * J + c2t * J) + c3.
#'
#' @return Required budget (and/or required level-3 sample size), statistical power, or MDES
#'     depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.3m
#'
#' @references
#'   Shen, Z., & Kelcey, B. (2022). Optimal sampling ratios in three-level
#'   multisite experiments. Journal of Research on Educational Effectiveness,
#'   15(1), 130-150.
#'   <https://doi.org/10.1080/19345747.2021.1953200>
#'
#' @examples
#' # Unconstrained optimal design #---------
#'   myod1 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005))
#'   myod1$out # n = 13.1, J = 15.3, p = 0.23
#'
#' # ------- Power analyses by default considering costs and budget -------
#' # Required budget and sample size
#'   mym.1 <- power.3m(expr = myod1, d = 0.2, q = 1, power = 0.8)
#'   mym.1$out  # m = 15491, K = 13.6
#'   # mym.1$par  # parameters and their values used for the function
#' # Or, equivalently, specify every argument in the function
#'   mym.1 <- power.3m(d = 0.2, power = 0.8, q = 1,
#'                  icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'                   r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'                   c1 = 1, c2 = 5,
#'                   c1t = 1, c2t = 200, c3 = 200,
#'                   n = 13, J = 15, p = 0.23)
#' # Required budget and sample size with constrained p
#'   mym.2 <- power.3m(expr = myod1, d = 0.2, q = 1, power = 0.8,
#'                  constraint = list(p = 0.5))
#'   mym.2$out  # m = 21072, K = 10.9
#' # Required budget and sample size with constrained p and n
#'   mym.3 <- power.3m(expr = myod1, d = 0.2, q = 1, power = 0.8,
#'                  constraint = list(p = 0.5, n = 20))
#'   mym.3$out  # m = 21252, K = 10.4
#'
#' # Power calculation
#'   mypower <- power.3m(expr = myod1, q = 1, d = 0.2, m = 15491)
#'   mypower$out  # power = 0.80
#' # Power calculation under constrained p (p = 0.5)
#'   mypower.1 <- power.3m(expr = myod1, q = 1, d = 0.2, m = 15491,
#'                  constraint = list(p = 0.5))
#'   mypower.1$out  # power = 0.62
#'
#' # MDES calculation
#'   mymdes <- power.3m(expr = myod1, q = 1, power = 0.80, m = 15491)
#'   mymdes$out  # d = 0.20
#'
#'
#' # ------- Conventional power analyses with cost.model = FALSE-------
#' # Required sample size
#'   myK <- power.3m(cost.model = FALSE, expr = myod1, d = 0.2, q = 1, power = 0.8)
#'   myK$out  # K = 13.6
#'   # myK$par  # parameters and their values used for the function
#' # Or, equivalently, specify every argument in the function
#'   myK <- power.3m(cost.model = FALSE, d = 0.2, power = 0.8, q = 1,
#'                   icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'                   r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'                   c1 = 1, c2 = 5,
#'                   c1t = 1, c2t = 200, c3 = 200,
#'                   n = 13, J = 15, p = 0.23)
#'
#' # Power calculation
#'   mypower1 <- power.3m(cost.model = FALSE, expr = myod1, K = 13.6, d = 0.2, q = 1)
#'   mypower1$out  # power = 0.80
#'
#' # MDES calculation
#'   mymdes1 <- power.3m(cost.model = FALSE, expr = myod1, K = 13.6, power = 0.8, q = 1)
#'   mymdes1$out  # d = 0.20
#'
power.3m <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                    sig.level = 0.05, two.tailed = TRUE,
                    d = NULL, power = NULL, m = NULL,
                    n = NULL, J = NULL, K = NULL, p = NULL,
                    icc2 = NULL, icc3 = NULL,
                    r12 = NULL, r22 = NULL, r32m = NULL, q = NULL,
                    c1 = NULL, c2 = NULL, c3 = NULL,
                    c1t = NULL, c2t = NULL, omega = NULL,
                    dlim = NULL, powerlim = NULL, Klim = NULL, mlim = NULL,
                    rounded = TRUE) {
  funName <- "power.3m"
  designType <- "three-level MRTs"
  if (cost.model == TRUE) {
    if (sum(sapply(list(m, d, power), is.null)) != 1)
      stop("exactly one of 'm', 'd', and 'power' must be NULL
           when cost.model is TRUE")
    if (!is.null(K))
      stop("'K' must be NULL when cost.model is TRUE")
  } else {
    if (sum(sapply(list(K, d, power), is.null)) != 1)
      stop("exactly one of 'K', 'd', and 'power' must be NULL
           when cost.model is FALSE")
    if (!is.null(m))
      stop("'m' must be NULL when cost.model is FALSE")
  }
  if (!is.null(expr)) {
    if (expr$funName != "od.3m") {
      stop("'expr' can only be NULL or
           the return from the function of 'od.3m'")
    } else {
      if (sum(sapply(list(icc2, icc3, r12, r22, r32m, c1, c2, c3,
                          c1t, c2t, omega, n, J, p), function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'icc2', 'icc3', 'r12', 'r22', 'r32m',
             'c1', 'c2', 'c3', 'c1t', 'c2t', 'omega', 'n', 'J', and 'p'
             have been specified in expr of 'od.3m'")
      icc2 <- expr$par$icc2
      icc3 <- expr$par$icc3
      r12 <- expr$par$r12
      r22 <- expr$par$r22
      r32m <- expr$par$r32m
      c1 <- expr$par$c1
      c2 <- expr$par$c2
      c3 <- expr$par$c3
      c1t <- expr$par$c1t
      c2t <- expr$par$c2t
      omega <- expr$par$omega
      if (rounded == TRUE) {
        n <- round(expr$out$n, 0)
        J <- round(expr$out$J, 0)
        p <- round(expr$out$p, 2)
      } else {
        n <- expr$out$n
        J <- expr$out$J
        p <- expr$out$p
      }
    }
  } else {
    if (!is.null(constraint))
      stop("'constraint' must be NULL when 'expr' is NULL")
  }
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (!is.null(constraint) && !is.list(constraint))
    stop("'constraint' must be in list format
         (e.g., constraint = list(p = 0.5, n = 20))")
  if (length(constraint) > 3)
    stop("'constraint' must be limited to 'n', 'J', and/or 'p'")
  if (!is.null(constraint$n)) {
    if(NumberCheck(constraint$n) || constraint$n <= 0)
      stop("constrained 'n' must be numeric with n > 0")
    n <- constraint$n
  }
  if (!is.null(constraint$J)) {
    if(NumberCheck(constraint$J) || constraint$J <= 0)
      stop("constrained 'J' must be numeric with J > 0")
    J <- constraint$J
  }
  if (!is.null(constraint$p)) {
    if(NumberCheck(constraint$p) ||
       any (0 >= constraint$p | constraint$p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p <- constraint$p
  }
  if (sum(sapply(list(p, power, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'p', 'power', and 'sig.level'
                 must be numeric in (0, 1)")
  if (sum(sapply(list(icc2, icc3, r12, r22, r32m), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1) stop("'icc2', 'icc3', 'r12', 'r22', and 'r32m' must be numeric in [0, 1]")
  if (cost.model == TRUE){
    if (sum(sapply(list(c1, c2, c3, c1t, c2t), function(x) {
      NumberCheck(x) || x < 0})) >= 1)
      stop("'c1', 'c2', 'c3', 'c1t', 'c2t' must be numeric in [0, Inf)")
    if (NumberCheck(m))
      stop("'m' must be numeric in [0, Inf)")
  }
  if (NumberCheck(q) | q < 0)
    stop("'q' must be numeric with q >= 0")
  if (NumberCheck(n) || n <= 0)
    stop("'n' must be numeric with n > 0")
  if (NumberCheck(J) || J <= 0)
    stop("'J' must be numeric with J > 0")
  if (NumberCheck(d) || any(0 > d | d > 5))
    stop("'d' must be numeric in [0, 5],
         please transfer negative effect size to positive one if needed")
  par <- list(cost.model = cost.model,
              sig.level = sig.level,
              two.tailed = two.tailed,
              d = d, icc2 = icc2, icc3 = icc3,
              r12 = r12, r22 = r22, r32m = r32m,
              c1 = c1, c2 = c2, c3 = c3,
              c1t = c1t, c2t = c2t, omega = omega,
              n = n, J = J, p = p,
              q = q, m = m, power = power)
  tside <- ifelse(two.tailed == TRUE, 2, 1)
  if (cost.model == TRUE) {
    Kcost <- ((1 - p) * (c1 * n * J + c2 * J) +
                p * (c1t * n * J + c2t * J) + c3)
    if (two.tailed == TRUE) {
      pwr.expr <- quote({
        K <- m / ((1 - p) * (c1 * n * J + c2 * J) +
                    p * (c1t * n * J + c2t * J) + c3);
        lambda <- d * sqrt((p * (1 - p) * n * J *  m / ((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3)) /
          (p * (1 - p) * n * J * omega * (1 - r32m) +
                 n * icc2 * (1 - r22)  +
                 (1 - icc2 - icc3) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = (K - q - 1)),
               df = (K - q - 1) , lambda) +
          pt(qt(sig.level / tside, df = (K - q - 1) ),
             df = (K - q - 1) , lambda)
      })
    } else {
      pwr.expr <- quote({
        K <- m / ((1 - p) * (c1 * n * J + c2 * J) +
                    p * (c1t * n * J + c2t * J) + c3);
        lambda <- d * sqrt((p * (1 - p) * n * J *  m / ((1 - p) * (c1 * n * J + c2 * J) +
                                p * (c1t * n * J + c2t * J) + c3)) /
                             (p * (1 - p) * n * J * omega * (1 - r32m) +
                                n * icc2 * (1 - r22)  +
                                (1 - icc2 - icc3) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = (K - q - 1)),
               df = (K - q - 1), lambda)
      })
    }
  } else {
    if (two.tailed == TRUE) {
      pwr.expr <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J * K) /
                             (p * (1 - p) * n * J * omega * (1 - r32m) +
                                n * icc2 * (1 - r22)  +
                                (1 - icc2 - icc3) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = (K - q - 1) ),
               df = (K - q - 1) , lambda) +
          pt(qt(sig.level / tside, df = (K - q - 1)),
             df = (K - q - 1), lambda)
      })
    } else {
      pwr.expr <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J * K) /
          (p * (1 - p) * n * J * omega * (1 - r32m) +
                 n * icc2 * (1 - r22)  +
                 (1 - icc2 - icc3) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = (K - q - 1)),
               df = (K - q - 1), lambda)
      })
    }
  }
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  Klim <- limFun(x = Klim, y = c(4, 1e+10))
  powerlim <- limFun(x = powerlim, y = c(1e-10, 1 - 1e-10))
  dlim <- limFun(x = dlim, y = c(0, 5))
  if(cost.model == TRUE) {
    mlim <- limFun(x = mlim, y = c(Klim[1] * Kcost, Klim[2] * Kcost))
    if (is.null(power)) {
      out <- list(power = eval(pwr.expr))
    } else if (is.null(m)) {
      out <- list(m = stats::uniroot(function(m)
        eval(pwr.expr) - power, mlim)$root)
      out <- c(out, list(K = out$m / ((1 - p) * (c1 * n * J + c2 * J) +
                                        p * (c1t * n * J + c2t * J) + c3)))
    } else if (is.null(d)) {
      out <- list(d = stats::uniroot(function(d)
        eval(pwr.expr) - power, dlim)$root)
    }
  } else {
    if (is.null(power)) {
      out <- list(power = eval(pwr.expr))
    } else if (is.null(K)) {
      out <- list(K = stats::uniroot(function(K)
        eval(pwr.expr) - power, Klim)$root)
    } else if (is.null(d)) {
      out <- list(d = stats::uniroot(function(d)
        eval(pwr.expr) - power, dlim)$root)
    }
  }
  power.out <- list(funName = funName,
                    designType = designType,
                    par = par, out = out)
  return(power.out)
}
