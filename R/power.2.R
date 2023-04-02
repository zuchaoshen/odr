#' Budget and/or sample size, power, MDES calculation for
#' two-level CRTs detecting main effects
#'
#' @description This function can calculate required budget for desired power,
#'     power or minimum detectable effect size (MDES) under fixed budget
#'     for two-level cluster randomized trials (CRTs).
#'     It also can perform conventional power analyses
#'     (e.g., required sample size, power, and MDES calculation).
#'
#' @inheritParams power.4
#' @param expr Returned object from function \code{\link{od.2}}; default is NULL;
#'     if \code{expr} is specified, parameter values of \code{icc},
#'     \code{r12}, \code{r22},
#'     \code{c1}, \code{c2}, \code{c1t}, \code{c2t}, \code{n}, and \code{p}
#'     used or solved in function \code{\link{od.2}} will
#'     be passed to the current function;
#'     only the values of \code{n} and \code{p} that specified or solved in
#'     function \code{\link{od.2}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint Specify the constrained values of \code{n}
#'     and/or \code{p} in list format to overwrite those
#'     from \code{expr}; default is NULL.
#' @param icc The unconditional intraclass correlation coefficient (ICC) in population or in
#'     each treatment condition.
#' @param J The total level-2 sample size.
#' @param p The proportion of level-2 clusters/units to be assigned to treatment.
#' @param q The number of level-2 covariates.
#' @param Jlim The range for searching the root of level-2 sample size (\code{J}) numerically,
#'     default is c(4, 10e10).
#' @param mlim The range for searching the root of budget (\code{m}) numerically,
#'     default is the costs sampling \code{Jlim} level-2 units across treatment conditions
#'     or c(4 * Jcost, 10e10 * Jcost), with Jcost = ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t)).
#' @param rounded Logical; round \code{n} and \code{p} that are from functions \code{od.2}
#'     to integer and two decimal places, respectively if TRUE,
#'     otherwise no rounding; default value is TRUE.
#'
#' @return Required budget (and/or required level-2 sample size), statistical power, or MDES
#'     depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.2
#'
#' @references
#'   Shen, Z., & Kelcey, B. (2020). Optimal sample allocation under unequal
#'   costs in cluster-randomized trials. Journal of Educational
#'   and Behavioral Statistics, 45(4): 446–474.
#'   <https://doi.org/10.3102/1076998620912418>
#' @examples
#' # Unconstrained optimal design
#'   myod1 <- od.2(icc = 0.2, r12 = 0.5, r22 = 0.5, c1 = 1, c2 = 5, c1t = 1, c2t = 50)
#'   myod1$out   # n = 8.9, p = 0.33
#'
#' # ------- Power analyses by default considering costs and budget -------
#' # Required budget and sample size
#'   mym.1 <- power.2(expr = myod1, d = 0.2, q = 1, power = 0.8)
#'   mym.1$out  # m = 3755, J = 130.2
#'   #mym.1$par  # parameters and their values used for the function
#' # Or, equivalently, specify every argument in the function
#'   mym.1 <- power.2(d = 0.2, power = 0.8, icc = 0.2,
#'                  c1 = 1, c2 = 5, c1t = 1, c2t = 50,
#'                   r12 = 0.5, r22 = 0.5, n = 9, p = 0.33, q = 1)
#' # Required budget and sample size with constrained p
#'   mym.2 <- power.2(expr = myod1, d = 0.2, q = 1, power = 0.8,
#'                constraint = list(p = 0.5))
#'   mym.2$out  # m = 4210, J = 115.3
#' # Required budget and sample size with constrained p and n
#'   mym.3 <- power.2(expr = myod1, d = 0.2, q = 1, power = 0.8,
#'                constraint = list(p = 0.5, n = 20))
#'   mym.3$out  # m = 4568, J = 96.2
#'
#' # Power calculation
#'   mypower <- power.2(expr = myod1, q = 1, d = 0.2, m = 3755)
#'   mypower$out  # power = 0.80
#' # Power calculation under constrained p (p = 0.5)
#'   mypower.1 <- power.2(expr = myod1, q = 1, d = 0.2, m = 3755,
#'                constraint = list(p = 0.5))
#'   mypower.1$out  # power = 0.75
#'
#' # MDES calculation
#'   mymdes <- power.2(expr = myod1, q = 1, power = 0.80, m = 3755)
#'   mymdes$out  # d = 0.20
#'
#'
#' # ------- Conventional power analyses with cost.model = FALSE-------
#' # Required J
#'   myJ <- power.2(cost.model = FALSE, expr = myod1, d = 0.2, q = 1, power = 0.8)
#'   myJ$out  # J = 130.2
#'   #myJ$par  # parameters and their values used for the function
#' # Or, equivalently, specify every argument in the function
#'   myJ <- power.2(cost.model = FALSE, d = 0.2, power = 0.8, icc = 0.2,
#'                   r12 = 0.5, r22 = 0.5, n = 9, p = 0.33, q = 1)
#'
#' # Power calculation
#'   mypower1 <- power.2(cost.model = FALSE, expr = myod1, J = 130, d = 0.2, q = 1)
#'   mypower1$out  # power = 0.80
#'
#' # MDES calculation
#'   mymdes1 <- power.2(cost.model = FALSE, expr = myod1, J = 130, power = 0.8, q = 1)
#'   mymdes1$out  # d = 0.20
#'
power.2 <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                    sig.level = 0.05, two.tailed = TRUE,
                    d = NULL, power = NULL, m = NULL,
                    n = NULL, J = NULL, p = NULL,
                    icc = NULL, r12 = NULL, r22 = NULL, q = NULL,
                    c1 = NULL, c2 = NULL, c1t = NULL, c2t = NULL,
                    dlim = NULL, powerlim = NULL, Jlim = NULL, mlim = NULL,
                    rounded = TRUE) {
  funName <- "power.2"
  designType <- "two-level CRTs"
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
    if (expr$funName != "od.2") {
      stop("'expr' can only be NULL or
            the return from the function of 'od.2'")
    } else {
      if (sum(sapply(list(icc, r12, r22, c1, c2, c1t, c2t, n, p),
                     function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'icc', 'r12', 'r22', 'c1', 'c2',
             'c1t', 'c2t', 'n', 'p'
             have been specified in expr of 'od.2'")
      icc <- expr$par$icc
      r12 <- expr$par$r12
      r22 <- expr$par$r22
      c1 <- expr$par$c1
      c2 <- expr$par$c2
      c1t <- expr$par$c1t
      c2t <- expr$par$c2t
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
         (e.g., constraint = list(p = 0.5))")
  if (length(constraint) > 2)
    stop("'constraint' must be limited to 'n' and/or 'p'")
  if (!is.null(constraint$n)) {
    if(NumberCheck(constraint$n) || constraint$n <= 0)
      stop("constrained 'n' must be numeric in (0, 10e10)")
    n <- constraint$n
  }
  if (!is.null(constraint$p)) {
    if(NumberCheck(constraint$p) ||
       any (0 >= constraint$p | constraint$p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p <- constraint$p
  }
  if (sum(sapply(list(icc, p, power, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1) stop("'icc', 'p', 'power', and 'sig.level'
                 must be numeric in [0, 1]")
  if (sum(sapply(list(r12, r22), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1) stop("'r12', 'r22' must be numeric in [0, 1]")
  if (cost.model == TRUE){
   if (sum(sapply(list(c1, c2, c1t, c2t), function(x) {
    NumberCheck(x) })) >= 1)
    stop("'c1', 'c2', 'c1t', 'c2t' must be numeric")
   if (NumberCheck(m))
    stop("'m' must be numeric")
   }
  if (NumberCheck(q) | q < 0)
   stop("'q' must be numeric in [0, 10e3]")
  if (NumberCheck(n) || n <= 0)
    stop("'n' must be numeric in (0, 10e10)")
  if (NumberCheck(d) || any(0 > d | d > 5))
    stop("'d' must be numeric in [0, 5],
         please transfer negative effect size to positive one if needed")
  par <- list(cost.model = cost.model,
              sig.level = sig.level,
              two.tailed = two.tailed,
              d = d, icc = icc, r12 = r12, r22 = r22,
              c1 = c1, c2 = c2, c1t = c1t, c2t = c2t,
              n = n, J = J, p = p,
              q = q, m = m, power = power)
  tside <- ifelse(two.tailed == TRUE, 2, 1)
  if (cost.model == TRUE) {
    if (two.tailed == TRUE) {
      pwr.expr <- quote({
        J <- m / ((1 - p) * (c1 * n + c2)
                  + p * (c1t * n + c2t));
        lambda <- d * sqrt(p * (1 - p) * J) /
          sqrt(icc * (1 - r22) + (1 - icc) * (1 - r12) / n);
        1 - pt(qt(1 - sig.level / tside, df = J - q - 2) ,
               df = J - q - 2, lambda) +
          pt(qt(sig.level / tside, df = J - q - 2),
             df = J - q - 2, lambda)
      })
    } else {
      pwr.expr <- quote({
        J <- m / ((1 - p) * (c1 * n + c2)
                  + p * (c1t * n + c2t));
        lambda <- d * sqrt(p * (1 - p) * J) /
          sqrt(icc * (1 - r22) + (1 - icc) * (1 - r12) / n);
        1 - pt(qt(1 - sig.level / tside, df = J - q - 2),
               df = J - q - 2, lambda)
      })
    }
  } else {
    if (two.tailed == TRUE) {
      pwr.expr <- quote({
        lambda <- d * sqrt(p * (1 - p) * J) /
          sqrt(icc * (1 - r22) + (1 - icc) * (1 - r12) / n);
        1 - pt(qt(1 - sig.level / tside, df = J - q - 2),
               df = J - q - 2, lambda) +
          pt(qt(sig.level / tside, df = J - q - 2),
             df = J - q - 2, lambda)
      })
    } else {
      pwr.expr <- quote({
        lambda <- d * sqrt(p * (1 - p) * J) /
          sqrt(icc * (1 - r22) + (1 - icc) * (1 - r12) / n);
        1 - pt(qt(1 - sig.level / tside, df = J - q - 2),
               df = J - q - 2, lambda)
      })
    }
  }
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  Jlim <- limFun(x = Jlim, y = c(4, 10e10))
  powerlim <- limFun(x = powerlim, y = c(1e-10, 1 - 1e-10))
  dlim <- limFun(x = dlim, y = c(0, 5))
  if(cost.model == TRUE){
    if (is.null(power)) {
      out <- list(power = eval(pwr.expr))
    } else if (is.null(m)) {
      Jcost <- ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))
      mlim <- limFun(x = mlim, y = c(Jlim[1] * Jcost, Jlim[2] * Jcost))
      out <- list(m = stats::uniroot(function(m)
        eval(pwr.expr) - power, mlim)$root)
      out <- c(out, list(J = out$m / (((1 - p) * (c1 * n + c2)
                                       + p * (c1t * n + c2t)))))
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

