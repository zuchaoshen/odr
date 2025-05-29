#' Budget and/or sample size, power, MDES calculation for two-level
#' MRTs detecting moderation effects with moderators at level 1
#'
#' @description This function can calculate required budget for desired power,
#'     power or minimum detectable effect size (MDES) under fixed budget
#'     for two-level multisite randomized trials (MRTs)
#'     detecting moderation effects with moderators at level 1.
#'     It also can perform conventional power analyses
#'     (e.g., required sample size, power, and MDES calculation).
#'
#' @inheritParams od.2m.111m
#' @inheritParams power.2m
#' @param binary Logical; binary moderator if TURE and continuous moderator if
#'     FALSE. Default is TRUE.
#' @param omega The treatment-by-site variance of the outcome.
#' @param power Statistical power.
#' @param q The number of additional covariates at level 1 beyond the treatment
#'     indicator, covariate, and the interaction between the moderator and
#'     the treatment. Default is 1.
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
#' @param c2 The cost of sampling one level-2 unit (site).
#' @param p The proportion of level-1 units to be assigned to treatment.
#' @param J The number of sites.
#' @param gammalim The range for searching the root of standardized
#'     moderation effect (gamma). Default is c(0, 5).
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
#' @export power.2m.111m
#'
#' @examples
#' myod <- od.2m.111m(icc = .2, r12 = .5, r22m = .5,
#'                    c1 = 10, c1t = 100, c2 = 50, omega = .01, gamma = 0.1)
#' mypowwer <- power.2m.111m(expr = myod, gamma = .1, power = .8)


power.2m.111m <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                          sig.level = 0.05, two.tailed = TRUE,
                          omega = NULL,
                          gamma = NULL, power = NULL, m = NULL,
                          n = NULL, J = NULL, p = NULL,
                          icc = NULL, r12 = NULL, q  = 1,
                          c1 = NULL, c2 = NULL, c1t = NULL,
                          gammalim = c(0, 5), powerlim = c(1e-10, 1 - 1e-10),
                          Jlim = c(6, 1e+10),
                          mlim = NULL,
                          rounded = TRUE,  binary = TRUE, Q = 0.5,
                          random.slope = TRUE) {
  funName <- "power.2m.111m"
  designType <- "two-level MRTs with individaul-level moderators"
  if (!is.null(expr)) {
    if (expr$funName %in% c("od.2m", "od.2m.111m")) {
      if (sum(sapply(list(icc, r12, c1, c2,
                          c1t, omega, n, p), function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'icc', 'r12',
             'c1', 'c2', 'c1t', 'omega', 'n', and 'p'
             have been specified in expr of 'od.2m.111m'
             or a similar od function")
      icc <- expr$par$icc
      r12 <- expr$par$r12
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
    } else  {
      stop("'expr' can only be NULL or
           the return from the function of 'od.2m' or a similar od function
             ('od.2m.111m')")
    }
  } else {
    if (!is.null(constraint))
      stop("'constraint' must be NULL when 'expr' is NULL")
  }
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (sum(sapply(list(icc, r12,
                      c1, c2, c1t, omega),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc', 'r12',
         'c1', 'c2', 'c1t', and 'omega' must be specified")
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
  if (sum(sapply(list(icc, r12, p, power, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1) stop("'icc', 'r12', 'p', 'power', and 'sig.level'
                 must be numeric in [0, 1]")

  if (cost.model == TRUE){
    if (sum(sapply(list(c1, c2, c1t), function(x) {
      NumberCheck(x)})) >= 1)
      stop("'c1', 'c2', 'c1t' must be numeric")
    if (NumberCheck(m))
      stop("'m' must be numeric")
  }
  if (NumberCheck(q ) | q  < 0)
    stop("'q' must be numeric with q  >= 0")
  if (NumberCheck(n) || n <= 0)
    stop("'n' must be numeric with n > 0")
  if (NumberCheck(gamma) || any(0 > gamma | gamma > 5))
    stop("'gamma' must be numeric in [0, 5],
         please transfer negative effect size to positive one if needed")


  par <- list(cost.model = cost.model, sig.level = sig.level,
              two.tailed = two.tailed, gamma = gamma, icc = icc,
              Q = Q, J = J,
              r12 = r12, c1 = c1, c2 = c2,
              c1t = c1t, omega = omega, n = n, p = p,
              q = q, m = m, power = power)

  tside <- ifelse(two.tailed == TRUE, 2, 1)
  if(binary){
    var.mod <- Q*(1-Q)
  } else {
    var.mod <- 1
  }

    if(random.slope){#power for moderation
    if(cost.model){
      if (two.tailed) {
        pwr2.expr <- quote({
          J <- m / ((1 - p) * c1 * n +  p * c1t * n + c2);
          lambda <- gamma * sqrt((p * (1 - p) * n * J*var.mod) /
                                   (p * (1 - p) * n * omega * var.mod +
                                      (1 - icc) * (1 - r12)));
          1 - pt(qt(1 - sig.level / tside, df = J - 1),
                 df = J - 1, lambda) +
            pt(qt(sig.level / tside, df = J - 1),
               df = J - 1, lambda)
        })
      } else {
        pwr2.expr <- quote({
          J <- m / ((1 - p) * c1 * n +  p * c1t * n + c2);
          lambda <- gamma * sqrt((p * (1 - p) * n * J*var.mod) /
                                   (p * (1 - p) * n * omega * var.mod +
                                      (1 - icc) * (1 - r12)));
          1 - pt(qt(1 - sig.level / tside, J - 1),
                 df = J - 1, lambda)
        })}
    } else {
      if (two.tailed) {
        pwr2.expr <- quote({
          lambda <- gamma * sqrt((p * (1 - p) * n * J*var.mod) /
                                   (p * (1 - p) * n * omega * var.mod +
                                      (1 - icc) * (1 - r12)));
          1 - pt(qt(1 - sig.level / tside, df = J - 1),
                 df = J - 1, lambda) +
            pt(qt(sig.level / tside, df = J - 1),
               df = J - 1, lambda)
        })
      } else {
        pwr2.expr <- quote({
          lambda <- gamma * sqrt((p * (1 - p) * n * J*var.mod) /
                                   (p * (1 - p) * n * omega * var.mod +
                                      (1 - icc) * (1 - r12)));
          1 - pt(qt(1 - sig.level / tside, J - 1),
                 df = J - 1, lambda)
        })}
    }

  } else {
    if (cost.model){
    if (two.tailed) {
      pwr2.expr <- quote({
        J <- m / ((1 - p) * c1 * n +  p * c1t * n + c2);
        lambda <- gamma * sqrt((p * (1 - p) * n * J*var.mod) /
                                 ((1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = J*(n-1) - q - 3),
               df = J*(n-1) - q - 3, lambda) +
          pt(qt(sig.level / tside, df = J*(n-1) - q - 3),
             df = J*(n-1) - q - 3, lambda)
      })
    } else {
      pwr2.expr <- quote({
        J <- m / ((1 - p) * c1 * n +  p * c1t * n + c2);
        lambda <- gamma * sqrt((p * (1 - p) * n * J*var.mod) /
                                 ((1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, J*(n-1) - q - 3),
               df = J*(n-1) - q - 3, lambda)
      })}
    } else {
    if (two.tailed) {
      pwr2.expr <- quote({
        lambda <- gamma * sqrt((p * (1 - p) * n * J*var.mod) /
                                 ((1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = J*(n-1) - q - 3),
               df = J*(n-1) - q - 3, lambda) +
          pt(qt(sig.level / tside, df = J*(n-1) - q - 3),
             df = J*(n-1) - q - 3, lambda)
      })
    } else {
      pwr2.expr <- quote({
        lambda <- gamma * sqrt((p * (1 - p) * n * J*var.mod) /
                                 ((1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, J*(n-1) - q - 3),
               df = J*(n-1) - q - 3, lambda)
      })}
  }}



  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  gammalim <- limFun(x = gammalim, y = c(0, 5))
  if(cost.model == TRUE) {
    Jcost <- ((1 - p) * c1 * n +  p * c1t * n + c2)
    mlim <- limFun(x = mlim, y = c(Jlim[1] * Jcost, Jlim[2] * Jcost))
    if (is.null(power)) {
      out <- list(power = eval(pwr2.expr))
    } else if (is.null(m)) {
      out <- list(m = stats::uniroot(function(m)
        eval(pwr2.expr) - power, mlim)$root)
      out <- c(out, list(J = out$m / ((1 - p) * c1 * n + p * c1t * n + c2)))
    } else if (is.null(gamma)) {
      out <- list(gamma = stats::uniroot(function(gamma)
        eval(pwr2.expr) - power, gammalim)$root)
    }
  } else {
    if (is.null(power)) {
      out <- list(power = eval(pwr2.expr))
    } else if (is.null(J)) {
      out <- list(J = stats::uniroot(function(J)
        eval(pwr2.expr) - power, Jlim)$root)
    } else if (is.null(gamma)) {
      out <- list(gamma = stats::uniroot(function(gamma)
        eval(pwr2.expr) - power, gammalim)$root)
    }
  }
  power.out <- list(funName = funName,
                    designType = designType,
                    par = par, out = out)
  return(power.out)
}

