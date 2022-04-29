#' Budget and/or sample size, power, MDES calculation for two-level
#' MRTs detecting main effects
#'
#' @description This function can calculate required budget for desired power,
#'     power or minimum detectable effect size (MDES) under fixed budget
#'     for two-level multisite randomized trials (MRTs).
#'     It also can perform conventional power analyses
#'     (e.g., required sample size, power, and MDES calculation).
#'
#' @inheritParams power.4
#' @inheritParams power.2
#' @inheritParams power.3m
#'
#' @param expr Returned objects from function \code{\link{od.2m}}; default is NULL;
#'     if \code{expr} is specified, parameter values of \code{icc},
#'     \code{r12}, \code{r22m},
#'     \code{c1}, \code{c2},
#'     \code{c1t}, \code{p}, and \code{n}
#'     used or solved in function \code{\link{od.2m}} will
#'     be passed to current function;
#'     only the values of \code{p} and \code{n} that specified or solved in
#'     function \code{\link{od.2m}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint Specify the constrained values of \code{p} and/or \code{n}
#'     in list format to overwrite those from \code{expr}; default value is NULL.
#' @param r22m The proportion of variance of site-specific treatment effect explained by covariates.
#' @param c2 The cost of sampling one level-2 unit.
#' @param p The proportion of level-1 units to be assigned to treatment.
#' @param q The number of covariates at level 2.
#' @param mlim The range for searching the root of budget (\code{m}) numerically,
#'     default is the costs sampling \code{Jlim} level-2 units
#'     or c(4 * Jcost, 1e+10 * Jcost) with Jcost =
#'     (1 - p) * c1 * n + p * c1t * n + c2.
#'
#' @return Required budget (and/or required level-2 sample size), statistical power, or MDES
#'     depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.2m
#'
#' @references
#'    Shen, Z., & Kelcey, B. (in press). Optimal sample
#'    allocation in multisite randomized trials. The Journal of Experimental Education.
#'    <https://doi.org/10.1080/00220973.2020.1830361>
#'
#' @examples
#' # Unconstrained optimal design #---------
#'   myod1 <- od.2m(icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
#'               c1 = 1, c2 = 10, c1t = 10,
#'               varlim = c(0, 0.005))
#'   myod1$out # n = 19.8, p = 0.37
#'
#' # ------- Power analyses by default considering costs and budget -------
#' # Required budget and sample size
#'   mym.1 <- power.2m(expr = myod1, d = 0.2, q = 1, power = 0.8)
#'   mym.1$out  # m = 2019, J = 20.9
#'   # mym.1$par  # parameters and their values used for the function
#' # Or, equivalently, specify every argument in the function
#'   mym.1 <- power.2m(d = 0.2, power = 0.8, q = 1,
#'                  icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
#'                  c1 = 1, c2 = 10, c1t = 10,
#'                  n = 20, p = 0.37)
#' # Required budget and sample size with constrained p
#'   mym.2 <- power.2m(expr = myod1, d = 0.2, q = 1, power = 0.8,
#'                  constraint = list(p = 0.5))
#'   mym.2$out  # m = 2373, J = 19.8
#' # Required budget and sample size with constrained p and n
#'   mym.3 <- power.2m(expr = myod1, d = 0.2, q = 1, power = 0.8,
#'                  constraint = list(p = 0.5, n = 5))
#'   mym.3$out  # m = 2502, J = 66.7
#'
#' # Power calculation
#'   mypower <- power.2m(expr = myod1, q = 1, d = 0.2, m = 2019)
#'   mypower$out  # power = 0.80
#' # Power calculation under constrained p (p = 0.5)
#'   mypower.1 <- power.2m(expr = myod1, q = 1, d = 0.2, m = 2019,
#'                  constraint = list(p = 0.5))
#'   mypower.1$out  # power = 0.72
#'
#' # MDES calculation
#'   mymdes <- power.2m(expr = myod1, q = 1, power = 0.80, m = 2019)
#'   mymdes$out  # d = 0.20
#'
#'
#' # ------- Conventional power analyses with cost.model = FALSE-------
#' # Required sample size
#'   myJ <- power.2m(cost.model = FALSE, expr = myod1, d = 0.2, q = 1, power = 0.8)
#'   myJ$out  # J = 6.3
#'   # myL$par  # parameters and their values used for the function
#' # Or, equivalently, specify every argument in the function
#'   myJ <- power.2m(cost.model = FALSE, d = 0.2, power = 0.8, q = 1,
#'                  icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
#'                  c1 = 1, c2 = 10, c1t = 10,
#'                  n = 20, p = 0.37)
#'
#' # Power calculation
#'   mypower1 <- power.2m(cost.model = FALSE, expr = myod1, J = 6.3, d = 0.2, q = 1)
#'   mypower1$out  # power = 0.80
#'
#' # MDES calculation
#'   mymdes1 <- power.2m(cost.model = FALSE, expr = myod1, J = 6.3, power = 0.8, q = 1)
#'   mymdes1$out  # d = 0.20
#'
power.2m <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                    sig.level = 0.05, two.tailed = TRUE,
                    d = NULL, power = NULL, m = NULL,
                    n = NULL, J = NULL, p = NULL,
                    icc = NULL, r12 = NULL, r22m = NULL, q  = NULL,
                    c1 = NULL, c2 = NULL, c1t = NULL,  omega = NULL,
                    dlim = NULL, powerlim = NULL, Jlim = NULL, mlim = NULL,
                    rounded = TRUE) {
  funName <- "power.2m"
  designType <- "two-level MRTs"
  if (cost.model == TRUE) {
    if (sum(sapply(list(m, d, power), is.null)) != 1)
      stop("exactly one of 'm', 'd', and 'power' must be NULL
           when cost.model is TRUE")
    if (!is.null(J))
      stop("'J' must be NULL when cost.model is TRUE")
  } else {
    if (sum(sapply(list(J, d, power), is.null)) != 1)
      stop("exactly one of 'J', 'd', and 'power' must be NULL
           when cost.model is FALSE")
    if (!is.null(m))
      stop("'m' must be NULL when cost.model is FALSE")
  }
  if (!is.null(expr)) {
    if (expr$funName != "od.2m") {
      stop("'expr' can only be NULL or
           the return from the function of 'od.2m'")
    } else {
      if (sum(sapply(list(icc, r12, r22m, c1, c2,
                          c1t, omega, n, p), function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'icc', 'r12', 'r22m',
             'c1', 'c2', 'c1t', 'omega', 'n', and 'p'
             have been specified in expr of 'od.2m'")
      icc <- expr$par$icc
      r12 <- expr$par$r12
      r22m <- expr$par$r22m
      c1 <- expr$par$c1
      c2 <- expr$par$c2
      c1t <- expr$par$c1t
      omega <- expr$par$omega
      if (rounded == TRUE) {
        n <- round(expr$out$n, 0)
        p <- round(expr$out$p, 2)
      } else {
        n <- expr$out$n
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
  if (length(constraint) > 2)
    stop("'constraint' must be limited to 'n' and/or 'p'")
  if (!is.null(constraint$n)) {
    if(NumberCheck(constraint$n) || constraint$n <= 0)
      stop("constrained 'n' must be numeric with n > 0")
    n <- constraint$n
  }
  if (!is.null(constraint$p)) {
    if(NumberCheck(constraint$p) ||
       any (0 >= constraint$p | constraint$p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p <- constraint$p
  }
  if (sum(sapply(list(icc, p, power, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'icc', 'p', 'power', and 'sig.level'
                 must be numeric in (0, 1)")
  if (sum(sapply(list(r12, r22m), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'r12', and 'r22m' must be numeric in [0, 1)")
  if (cost.model == TRUE){
    if (sum(sapply(list(c1, c2, c1t), function(x) {
      NumberCheck(x) || x < 0})) >= 1)
      stop("'c1', 'c2', 'c1t' must be numeric in [0, Inf)")
    if (NumberCheck(m))
      stop("'m' must be numeric in [0, Inf)")
  }
  if (NumberCheck(q ) | q  < 0)
    stop("'q' must be numeric with q  >= 0")
  if (NumberCheck(n) || n <= 0)
    stop("'n' must be numeric with n > 0")
  if (NumberCheck(d) || any(0 > d | d > 5))
    stop("'d' must be numeric in [0, 5],
         please transfer negative effect size to positive one if needed")
  par <- list(cost.model = cost.model, sig.level = sig.level,
              two.tailed = two.tailed, d = d, icc = icc,
              r12 = r12, r22m = r22m, c1 = c1, c2 = c2,
              c1t = c1t, omega = omega, n = n, p = p,
              q = q, m = m, power = power)
  tside <- ifelse(two.tailed == TRUE, 2, 1)
  if (cost.model == TRUE) {
    Jcost <- ((1 - p) * c1 * n +  p * c1t * n + c2)
    if (two.tailed == TRUE) {
      pwr.expr <- quote({
        J <- m / ((1 - p) * c1 * n +  p * c1t * n + c2);
        lambda <- d * sqrt((p * (1 - p) * n * m / ((1 - p) * c1 * n +  p * c1t * n + c2)) /
          (p * (1 - p) * n * omega * (1 - r22m) + (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = J - q - 1),
               df = J - q - 1, lambda) +
          pt(qt(sig.level / tside, df = J - q - 1),
             df = J - q - 1, lambda)
      })
    } else {
      pwr.expr <- quote({
        J <- m / ((1 - p) * c1 * n +  p * c1t * n + c2);
        lambda <- d * sqrt((p * (1 - p) * n * m / ((1 - p) * c1 * n +  p * c1t * n + c2)) /
                             (p * (1 - p) * n * omega * (1 - r22m) + (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, J - q - 1),
               df = J - q - 1, lambda)
      })
    }
  } else {
    if (two.tailed == TRUE) {
      pwr.expr <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J) /
                             (p * (1 - p) * n * omega * (1 - r22m) + (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = J - q - 1),
               df = J - q - 1, lambda) +
          pt(qt(sig.level / tside, df = J - q - 1),
             df = J - q - 1, lambda)
      })
    } else {
      pwr.expr <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J) /
                             (p * (1 - p) * n * omega * (1 - r22m) + (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, J - q - 1),
               df = J - q - 1, lambda)
      })
    }
  }
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  Jlim <- limFun(x = Jlim, y = c(4, 1e+10))
  powerlim <- limFun(x = powerlim, y = c(1e-10, 1 - 1e-10))
  dlim <- limFun(x = dlim, y = c(0, 5))
  if(cost.model == TRUE) {
    mlim <- limFun(x = mlim, y = c(Jlim[1] * Jcost, Jlim[2] * Jcost))
    if (is.null(power)) {
      out <- list(power = eval(pwr.expr))
    } else if (is.null(m)) {
      out <- list(m = stats::uniroot(function(m)
        eval(pwr.expr) - power, mlim)$root)
      out <- c(out, list(J = out$m / ((1 - p) * c1 * n + p * c1t * n + c2)))
    } else if (is.null(d)) {
      out <- list(d = stats::uniroot(function(d)
        eval(pwr.expr) - power, dlim)$root)
    }
  } else {
    if (is.null(power)) {
      out <- list(power = eval(pwr.expr))
    } else if (is.null(J)) {
      out <- list(J = stats::uniroot(function(J)
        eval(pwr.expr) - power, Jlim)$root)
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
