#' Budget and/or sample size, power, MDES calculation for individual randomized controlled trials
#'
#' @description This function can calculate required budget for desired power,
#'     power or minimum detectable effect size (MDES) under fixed budget
#'     for individual randomized controlled trials (RCTs).
#'     It also can perform conventional power analyses
#'     (e.g., required sample size, power, and MDES calculation).
#'
#' @inheritParams power.4
#' @param expr returned object from function \code{\link{od.1}}; default value is NULL;
#'     if \code{expr} is specified, parameter values of \code{r12},
#'     \code{c1}, \code{c1t}, and \code{p}
#'     used or solved in function \code{\link{od.1}} will
#'     be passed to the current function;
#'     only the value of \code{p} that specified or solved in
#'     function \code{\link{od.1}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint specify the constrained value of
#'     \code{p} in list format to overwrite that
#'     from \code{expr}; default value is NULL.
#' @param r12 the proportion of outcome variance explained by covariates.
#' @param c1 the cost of sampling one unit in control condition.
#' @param c1t the cost of sampling one unit in treatment condition.
#' @param n the total sample size.
#' @param p the proportion of individuals to be assigned to treatment.
#' @param q the number of covariates.
#' @param nlim the range for searching the root of sample size (\code{n}) numerically,
#'     default value is c(4, 10e10)
#' @param mlim the range for searching the root of budget (\code{m}) numerically,
#'     default value is the costs sampling \code{nlim} units across treatment conditions
#'     or c(4 * ncost, 10e10 * ncost) with ncost = ((1 - p) * c1 + p * c1t)
#' @param rounded logical; round \code{p} that is from functions \code{od.1}
#'     to two decimal places if TRUE,
#'     otherwise no rounding; default value is TRUE.
#'
#' @return Required budget (or required sample size), statistical power, or MDES
#'     depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.1
#'
#' @references
#'   Shen, Z. (2019). Optimal sample allocation in multilevel experiments
#'   (Doctoral dissertation). University of Cincinnati, Cincinnati, OH.
#'
#' @examples
#' # unconstrained optimal design
#'   myod1 <- od.1(r12 = 0.5, c1 = 1, c1t = 5, varlim = c(0, 0.2))
#'   myod1$out   # p = 0.31
#'
#' # ------- power analyses by default considering costs and budget -------
#' # required budget and sample size
#'   mym.1 <- power.1(expr = myod1, d = 0.2, q = 1, power = 0.8)
#'   mym.1$out  # m = 1032 n = 461
#'   # mym.1$par  # parameters and their values used for the function
#' # or equivalently, specify every argument in the function
#'   mym.1 <- power.1(d = 0.2, power = 0.8, c1 = 1, c1t = 5,
#'                   r12 = 0.5, p = 0.31, q = 1)
#' # required budget and sample size with constrained p
#'   mym.2 <- power.1(expr = myod1, d = 0.2, q = 1, power = 0.8,
#'                constraint = list(p = 0.5))
#'   mym.2$out  # m = 1183, n = 394
#'
#' # Power calculation
#'   mypower <- power.1(expr = myod1, q = 1, d = 0.2, m = 1032)
#'   mypower$out  # power = 0.80
#' # Power calculation under constrained p (p = 0.5)
#'   mypower.1 <- power.1(expr = myod1, q = 1, d = 0.2, m = 1032,
#'                constraint = list(p = 0.5))
#'   mypower.1$out  # power = 0.74
#'
#' # MDES calculation
#'   mymdes <- power.1(expr = myod1, q = 1, power = 0.80, m = 1032)
#'   mymdes$out  # d = 0.20
#'
#'
#' # ------- conventional power analyses with cost.model = FALSE-------
#' # Required sample size n
#'   myn <- power.1(cost.model = FALSE, expr = myod1, d = 0.2, q = 1, power = 0.8)
#'   myn$out  # n = 461
#'   # myn$par  # parameters and their values used for the function
#' # or equivalently, specify every argument in the function
#'   myn <- power.1(cost.model = FALSE, d = 0.2, power = 0.8,
#'                   r12 = 0.5, p = 0.31, q = 1)
#'
#' # Power calculation
#'   mypower1 <- power.1(cost.model = FALSE, expr = myod1, n = 461, d = 0.2, q = 1)
#'   mypower1$out  # power = 0.80
#'
#' # MDES calculation
#'   mymdes1 <- power.1(cost.model = FALSE, expr = myod1, n = 461, power = 0.8, q = 1)
#'   mymdes1$out  # d = 0.20
#'
power.1 <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                    sig.level = 0.05, two.tailed = TRUE,
                    d = NULL, power = NULL, m = NULL,
                    n = NULL, p = NULL,
                    r12 = NULL, q = NULL,
                    c1 = NULL, c1t = NULL,
                    dlim = NULL, powerlim = NULL, nlim = NULL, mlim = NULL,
                    rounded = TRUE) {
  funName <- "power.1"
  designType <- "individual RCTs"
  if (cost.model == TRUE) {
    if (sum(sapply(list(m, d, power), is.null)) != 1)
      stop("exactly one of 'm', 'd', and 'power' must be NULL
           when cost.model is TRUE")
    if (!is.null(n))
      stop("'n' must be NULL when cost.model is TRUE")
  } else {
    if (sum(sapply(list(n, d, power), is.null)) != 1)
      stop("exactly one of 'n', 'd', and 'power' must be NULL
           when cost.model is FALSE")
    if (!is.null(m))
      stop("'m' must be NULL when cost.model is FALSE")
  }
  if (!is.null(expr)) {
    if (expr$funName != "od.1") {
      stop("'expr' can only be NULL or
            the return from the function of 'od.1'")
    } else {
      if (sum(sapply(list(r12, c1, c1t, p),
                     function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'r12', 'c1', 'c1t', 'p'
             have been specified in expr of 'od.1'")
      r12 <- expr$par$r12
      c1 <- expr$par$c1
      c1t <- expr$par$c1t
      if (rounded == TRUE) {
        p <- round(expr$out$p, 2)
      } else {
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
         (e.g., constraint = list(p = 0.5))")
  if (length(constraint) > 1)
    stop("'constraint' must be limited to 'p'")
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
  if (sum(sapply(list(r12), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'r12' must be numeric in [0, 1)")
  if (cost.model == TRUE){
   if (sum(sapply(list(c1, c1t), function(x) {
    NumberCheck(x) || x < 0})) >= 1)
    stop("'c1', 'c1t' must be numeric in [0, Inf)")
   if (NumberCheck(m))
    stop("'m' must be numeric in [0, Inf)")
   }
  if (NumberCheck(q) | q < 0)
   stop("'q' must be numeric in [0, 10e3]")
  if (NumberCheck(d) || any(0 > d | d > 5))
    stop("'d' must be numeric in [0, 5],
         please transfer negative effect size to positive one if needed")
  if (r12 > 0 && q == 0)
    stop("'q' must be q >= 1 when r12 != 0")
  par <- list(cost.model = cost.model,
              sig.level = sig.level,
              two.tailed = two.tailed,
              d = d, r12 = r12,
              c1 = c1, c1t = c1t,
              n = n, p = p,
              q = q, m = m, power = power)
  tside <- ifelse(two.tailed == TRUE, 2, 1)
  if (cost.model == TRUE) {
    if (two.tailed == TRUE) {
      pwr.expr <- quote({
        n <- m / ((1 - p) * c1 + p * c1t);
        lambda <- d * sqrt(p * (1 - p) * n) /
          sqrt(1 - r12);
        1 - pt(qt(1 - sig.level / tside, df = n - q - 2) ,
               df = n - q - 2, lambda) +
          pt(qt(sig.level / tside, df = n - q - 2),
             df = n - q - 2, lambda)
      })
    } else {
      pwr.expr <- quote({
        n <- m / ((1 - p) * c1 + p * c1t);
        lambda <- d * sqrt(p * (1 - p) * n) /
          sqrt(1 - r12);
        1 - pt(qt(1 - sig.level / tside, df = n - q - 2),
               df = n - q - 2, lambda)
      })
    }
  } else {
    if (two.tailed == TRUE) {
      pwr.expr <- quote({
        lambda <- d * sqrt(p * (1 - p) * n) /
          sqrt(1 - r12);
        1 - pt(qt(1 - sig.level / tside, df = n - q - 2),
               df = n - q - 2, lambda) +
          pt(qt(sig.level / tside, df = n - q - 2),
             df = n - q - 2, lambda)
      })
    } else {
      pwr.expr <- quote({
        lambda <- d * sqrt(p * (1 - p) * n) /
          sqrt(1 - r12);
        1 - pt(qt(1 - sig.level / tside, df = n - q - 2),
               df = n - q - 2, lambda)
      })
    }
  }
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  nlim <- limFun(x = nlim, y = c(4, 10e10))
  powerlim <- limFun(x = powerlim, y = c(1e-10, 1 - 1e-10))
  dlim <- limFun(x = dlim, y = c(0, 5))
  if(cost.model == TRUE){
    if (is.null(power)) {
      out <- list(power = eval(pwr.expr))
    } else if (is.null(m)) {
      ncost <- ((1 - p) * c1 + p * c1t)
      mlim <- limFun(x = mlim, y = c(nlim[1] * ncost, nlim[2] * ncost))
      out <- list(m = stats::uniroot(function(m)
        eval(pwr.expr) - power, mlim)$root)
      out <- c(out, list(n = out$m / (((1 - p) * c1
                                       + p * c1t))))
    } else if (is.null(d)) {
      out <- list(d = stats::uniroot(function(d)
        eval(pwr.expr) - power, dlim)$root)
    }
  } else {
    if (is.null(power)) {
      out <- list(power = eval(pwr.expr))
    } else if (is.null(n)) {
      out <- list(n = stats::uniroot(function(n)
        eval(pwr.expr) - power, nlim)$root)
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

