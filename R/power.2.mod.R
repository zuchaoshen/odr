#' Statistical power, sample size (and/or budget), minimum detectable moderator
#' effect size calculation for two-level cluster-randomized trials (CRTs)
#' detecting moderation effects
#'
#' @description This function can calculate power, required sample size/budget
#'     for desired power, or minimum detectable moderation effect size (MDMES)
#'     under a fixed budget in two-level CRTs.
#'     It also can perform conventional power analyses
#'     (e.g., required sample size, power, and MDMES calculation).
#'
#' @inheritParams od.2.mod
#' @inheritParams power.2
#' @param power Statistical power.mod for a moderation effect.
#' @param expr Returned objects from function \code{\link{od.2.mod}}; default is NULL;
#'     if \code{expr} is specified, parameter values of \code{icc},
#'     \code{r12}, \code{r22}, \code{r12m}, \code{r22m},
#'     \code{c1}, \code{c2},
#'     \code{c1t}, \code{c2t}, \code{p}, and \code{n}
#'     used or solved in function \code{\link{od.2.mod}} will
#'     be passed to the current function;
#'     only the values of \code{p} and \code{n} that specified or solved in
#'     function \code{\link{od.2.mod}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint Specify the constrained values of \code{p} and/or \code{n}
#'     in list format to overwrite those from \code{expr}; default value is NULL.
#' @param gammalim The range for numerically solving the root of standardized
#'     moderation effect (gamma). Default is c(0, 5).
#' @param mlim The range for numerically solving the root of budget (\code{m}).
#'     The default is NULL, which mlim = Jlim times the costs for each site and
#'     its members.
#' @param Jlim The range for numerically solving the root of
#'     the sample size requirement(\code{J}).
#' @return Required budget (and/or required level-2 sample size), statistical power,
#'     or MDMES depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#' @export power.2.mod
#' @examples
#' myod <- od.2.mod(icc = .2, r12 = .5, r22 = .5,
#'                  c1 = 10, c1t = 100, c2 = 50, c2t = 500,
#'                  gamma = 0.2, d = 0.2)
#' mypower <- power.2.mod(expr = myod, m=myod$out$m, gamma = 0.2); mypower$out
#' mym <- power.2.mod(expr = myod, power.mod = .80, gamma = 0.2); mym$out
#'
power.2.mod <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                          sig.level = 0.05, two.tailed = TRUE,
                          gamma = NULL, power.mod = NULL, m = NULL,
                          n = NULL, J = NULL, p = NULL,
                          icc = NULL, r12 = NULL, r22 =NULL,
                          r12m = NULL, r22m = NULL, q.mod = 1,
                          c1 = NULL, c2 = NULL, c1t = NULL, c2t = NULL,
                          gammalim = c(0, 5), powerlim = c(1e-10, 1 - 1e-10),
                          Jlim = c(5.5, 1e+10),
                          binary = TRUE,
                          mlim = NULL,
                          rounded = TRUE, Q = 0.5) {
  funName <- "power.2.mod"
  designType <- "2-2-1 moderation in 2-level CRTs"
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (!is.null(expr)) {
    if (expr$funName %in% c("od.2", "od.2.mod", "od.2.only.mod")) {
      if (sum(sapply(list(icc, r12, r22, r12m, r22m, c1, c2,
                          c1t, c2t), function(x) {!is.null(x)})) >= 1){
        stop("parameters of 'icc', 'r12', 'r22', 'r12m', 'r22m',
             'c1', 'c2', 'c1t', and 'c2t'
             have been specified in expr of 'od.2.mod'
             or a similar od function")
      }
      icc <- expr$par$icc
      r12 <- expr$par$r12
      r22 <- expr$par$r22
      r12m <- expr$par$r12m
      r22m <- expr$par$r22m
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
      if(expr$funName %in% c("od.2", "od.2.only.mod")){
      }else{q.mod <- expr$par$q.mod}
      binary <- expr$par$binary
      if (is.null(n)){
        if (rounded) {n <- round(expr$out$n, 0)} else {n <- expr$out$n}
      }
      if(is.null(p)){
        if (rounded) {p <- round(expr$out$p, 2)} else { p <- expr$out$p}
      }
    } else{
      stop("'expr' can only be NULL or
           the return from the function of 'od.2' or a similar od function
             ('od.2.mod', 'od.2.only.mod')")
    }
  }
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
  if (isTRUE(cost.model)) {
    if (sum(sapply(list(m, gamma, power.mod), is.null)) != 1)
      stop("exactly one of 'm', 'gamma', and 'power.mod' must be NULL
           when cost.model is TRUE")
    if (!is.null(J)){stop("'J' must be NULL when cost.model is TRUE")}
    if (sum(sapply(list(c1, c2, c1t, c2t), function(x) {
      NumberCheck(x)})) >= 1)
      stop("'c1', 'c2', 'c1t', 'c2t' must be numeric")
    if (NumberCheck(m))
      stop("'m' must be numeric")
  }
  if (sum(sapply(list(icc, r12, r22,
                      c1, c1t, c2, c2t),
                 function(x) is.null(x))) >= 1){
    stop("All of 'icc', 'r12', 'r22',
         'c1', 'c1t', 'c2', and 'c2t' must be specified")}
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
   if (sum(sapply(list(icc, r12, r22, p, power.mod, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1) stop("'icc', 'r12', 'r22', 'p', 'power.mod', and 'sig.level'
                 must be numeric in [0, 1]")
  if (NumberCheck(q.mod ) | q.mod  < 0)
    stop("'q.mod' must be numeric with q  >= 0")
  if (NumberCheck(n) || n <= 0)
    stop("'n' must be numeric with n > 0")
  if (NumberCheck(gamma) || any(0 > gamma | gamma > 5))
    stop("'gamma' must be numeric in [0, 5],
         please transfer negative effect size to positive one if needed")
  if(is.null(r12m)) {r12m = r12}
  if(is.null(r22m)) {r22m = r22}
  par <- list(cost.model = cost.model, sig.level = sig.level,
              two.tailed = two.tailed, gamma = gamma, icc = icc,
              Q = Q, J = J, binary = binary,
              r12 = r12, r22 = r22, c1 = c1, c2 = c2,
              c1t = c1t, c2t = c2t, n = n, p = p,
              q.mod = q.mod, m = m, power.mod = power.mod)

  tside <- ifelse(two.tailed == TRUE, 2, 1)
  if(binary){var.mod <- Q*(1-Q)} else {var.mod <- 1}
    if(cost.model){
      if (two.tailed) {
        pwr.mod <- quote({
          J <- m / ((1 - p) * (c1 * n + c2) +  p *(c1t * n + c2t));
        lambda <- gamma/sqrt((icc*(1-r22m)+(1-icc)*(1-r12m)/n)/
                               (p*(1-p)*var.mod*J));
        1 - pt(qt(1 - sig.level/tside, df = J - q.mod - 4),
               df = J - q.mod - 4, lambda)
        })
      } else {
        pwr.mod <- quote({
          J <- m / ((1 - p) * (c1 * n + c2) +  p *(c1t * n + c2t));
          lambda <- gamma/sqrt((icc*(1-r22m)+(1-icc)*(1-r12m)/n)/
                                 (p*(1-p)*var.mod*J));
          1 - pt(qt(1 - sig.level/tside, df = J - q.mod - 4),
                 df = J - q.mod - 4, lambda)
        })}
    } else {
      if (two.tailed) {
        pwr.mod <- quote({
          lambda <- gamma/sqrt((icc*(1-r22m)+(1-icc)*(1-r12m)/n)/
                                 (p*(1-p)*var.mod*J));
          1 - pt(qt(1 - sig.level/tside, df = J - q.mod - 4),
                 df = J - q.mod - 4, lambda)
        })
      } else {
        pwr.mod <- quote({
          lambda <- gamma/sqrt((icc*(1-r22m)+(1-icc)*(1-r12m)/n)/
                                 (p*(1-p)*var.mod*J));
          1 - pt(qt(1 - sig.level/tside, df = J - q.mod - 4),
                 df = J - q.mod - 4, lambda)
        })}
    }

  if(cost.model) {
    Jcost <- (1 - p) *(c1 * n + c2) + p * (c1t * n + c2t)
    if(is.null(mlim)){  mlim <-  c(Jlim[1] * Jcost, Jlim[2] * Jcost)}

    if (is.null(power.mod)) {
      out <- list(power.mod = eval(pwr.mod))
    } else if (is.null(m)) {
      out <- list(m = stats::uniroot(function(m)
        eval(pwr.mod) - power.mod, mlim)$root)
      out <- c(out, list(J = out$m / ((1 - p) * (c1 * n + c2) +
                                        p * (c1t * n + c2t))))
    } else if (is.null(gamma)) {
      out <- list(gamma = stats::uniroot(function(gamma)
        eval(pwr.mod) - power.mod, gammalim)$root)
    }
  } else {
    if (is.null(power.mod)) {
      out <- list(power.mod = eval(pwr.mod))
    } else if (is.null(J)) {
      out <- list(J = stats::uniroot(function(J)
        eval(pwr.mod) - power.mod, Jlim)$root)
    } else if (is.null(gamma)) {
      out <- list(gamma = stats::uniroot(function(gamma)
        eval(pwr.mod) - power.mod, gammalim)$root)
    }
  }
  return(list(funName = funName,
              designType = designType,
              par = par, out = out))
}
