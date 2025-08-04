#' Budget and/or sample size, power, MDES calculation for single-level
#'     randomized controlled trials (RCTs) investigating
#'     moderation effects  (1-1-1m)
#' @description This function can calculate required budget for desired power,
#'     the minimum detectable effect size, and
#'     statistical power under a fixed budget in
#'     randomized controlled trials (RCTs) probing moderation effects.
#'     It also can perform conventional power analyses
#'     (e.g., required sample size calculation, minimum detectable effect size
#'     calculation, and power calculation).
#' @param expr Returned object from function \code{\link{od.1.111m}}; default value is NULL;
#'     if \code{expr} is specified, parameter values of \code{a}, \code{b},
#'     \code{c}, \code{ct}, and \code{p}
#'     used or solved in function \code{\link{od.1.111m}} will
#'     be passed to the current function;
#'     only the values of \code{p} that specified or solved in
#'     function \code{\link{od.1.111m}} can be overwritten
#'     if \code{constraint} is specified.
#' @param cost.model Logical; power analyses accommodating costs and budget
#'     (e.g., required budget for a desired power, power under fixed budget)
#'     if TRUE. Otherwise, conventional power analyses are performed
#'     (e.g., required sample size and power calculation); default value is TRUE.
#' @param binary Logical. The moderator is binary if TRUE and continuous if
#'     FALSE. Default is TRUE.
#' @param gamma Moderated treatment effect.
#' @param Q The proportion of individuals in one group the binary moderator.
#'     Default value is 0.5, which requires the minimum number of individuals
#'     to achieve a targeted power. Change it as necessary.
#' @param n Total number of individuals.
#' @param m Total budget.
#' @param p The proportion of individuals assigned to the experimental group.
#' @param sig.level Significance level, default value is .05.
#' @param two.tailed Logical; two-tailed tests if TRUE,
#'     otherwise one-tailed tests; default value is TRUE.
#' @param r.yx Within-treatment correlation between the outcome (y) and
#'     the covariate (x) for continuous moderators. Within-treatment
#'     within-moderator correlation between the outcome (y) and
#'     the covariate (x) for binary moderators.
#' @param r.mx Within-treatment correlation between the moderator (m) and
#'     the covariate (x), if specified, for continuous moderators.
#' @param r.ym Within-treatment correlation between the outcome (y) and
#'     the moderator (m), if specified, for continuous moderators.
#' @param c1 The cost of sampling one unit in control condition.
#' @param c1t The cost of sampling one unit in treatment condition.
#' @param constraint If specified, the constrained value of
#'     \code{p} in a list format (e.g., constraint = list(p = 0.5))
#'     will overwrite that
#'     from \code{expr}; default value is NULL.
#' @param q.mod The number of covariates in the moderation model (besides the
#'     treatment, moderator, and their interaction term). The default value is 1.
#' @param gammalim  The range for identifying the root of moderation
#'     effect size (\code{gamma}) numerically,
#'     default value is c(0.005, 5).
#' @param power Statistical power.
#' @param powerlim  The range for identifying the root of power
#'     (\code{power}) numerically,
#'     default value is c(0.0001, 0.9999).
#' @param mlim The range for identifying the root of budget (\code{m}) numerically,
#'     default value is the costs sampling \code{nlim} units.
#' @param nlim The range for identifying the root of sample size (\code{n})
#'     numerically. Default is c(20, 1e7).
#' @return Required budget (\code{m}) or required sample size (\code{n}),
#'     statistical power(\code{power}),
#'     minimum detectable moderation effect size (\code{gamma}),
#'     depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.1.111m
#'
#' @examples
#' # Optimal design and power analyses accommodating costs and budget
#' myod <- od.1.111m(d =.1, gamma = .2, r12 = .50,
#'                  c1 = 10, c1t = 100)
#' myod
#' N <- power.1.111m(expr = myod, power = .8)
#' N$out
#'
#' # Conventional power analyses
#' # Required sample size for a binary moderator
#' N <- power.1.111m(cost.model = FALSE, gamma = .2, power = .8, p =.5)
#' N
#'
#' # Required sample size for a continuous moderator
#' N <- power.1.111m(cost.model = FALSE,
#'                   gamma = .2, power = .8, p =.5, binary = FALSE)
#' N
#'
power.1.111m <- function(cost.model = TRUE, expr = NULL,
                        constraint = NULL,
                        sig.level = 0.05, two.tailed = TRUE,
                        gamma = NULL,
                        binary = TRUE,
                        power = NULL, m = NULL,
                        n = NULL, p = NULL, Q = 0.5,
                        c1 = NULL, c1t = NULL,
                        r.yx = 0, r.mx = 0, r.ym = 0,
                        q.mod = 1, gammalim = c(0.005, 5),
                        powerlim = c(0.0001, 0.9999), nlim = c(20, 1e7), mlim = NULL
                        ) {
  funName <- "power.1.111m"
  designType <- "1-1-1 moderation in single-level RCTs"
  if (!is.null(expr)) {
    if (expr$funName != "od.1.111m") {
      stop("'expr' can only be NULL or
            the return from the function of 'od.1.111m'")
    } else {
      if (sum(sapply(list(gamma, c1, c1t, p),
                     function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'gamma', 'c1', 'c1t',
          and 'p'
             have been specified in expr of 'od.1.111m'")
      gamma <- expr$par$gamma
      Q <- expr$par$Q
      c1 <- expr$par$c1
      c1t <- expr$par$c1t
      r.yx <- expr$par$r.yx
      r.mx <- expr$par$r.mx
      r.ym <- expr$par$r.ym
      q.mod <- expr$par$q.mod
      two.tailed <- expr$par$two.tailed
      sig.level <- expr$par$sig.level
      binary <- expr$par$binary
      p <- expr$out$p
    }
  } else {
    if (!is.null(constraint))
      stop("'constraint' must be NULL when 'expr' is NULL")
  }
  if (cost.model) {
    if (sum(sapply(list(m, gamma, power), is.null)) != 1)
      stop("exactly one of 'm', 'gamma', and 'power' must be NULL
           when cost.model is TRUE")
    if (!is.null(n))
      stop("'n' must be NULL when cost.model is TRUE")
  } else {
    if (sum(sapply(list(n, gamma, power), is.null)) != 1)
      stop("exactly one of 'n', 'gamma', and 'power' must be NULL
           when cost.model is FALSE")
    if (!is.null(m))
      stop("'m' must be NULL when cost.model is FALSE")
  }

  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (!is.null(constraint) && !is.list(constraint))
    stop("'constraint' must be in list format
         (e.g., constraint = list(p = 0.5))")
  if (length(constraint) > 1)
    stop("'constraint' must be limited to 'p' ")
  if (!is.null(constraint$p)) {
    if(NumberCheck(constraint$p) ||
       any (0 >= constraint$p | constraint$p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p <- constraint$p
  }
  if (sum(sapply(list(p, power, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'p', 'power', and 'sig.level' must be numeric in (0, 1]")

  if (cost.model){
    if (sum(sapply(list(c1, c1t), function(x) {
      NumberCheck(x) || x < 0})) >= 1)
      stop("'c1', 'c1t' must be numeric")
    if (NumberCheck(m))
      stop("'m' must be numeric in [0, Inf)")
  }
  if (NumberCheck(gamma) || any(-5 > gamma | gamma > 5))
    stop("'gamma' must be numeric in [-5, 5]")
  par <- list(cost.model = cost.model,
              expr = expr,
              sig.level = sig.level,
              two.tailed = two.tailed,
              binary = binary,
              gamma = gamma,
              r.yx = r.yx, r.mx = r.mx, r.ym = r.ym,
              c1 = c1, c1t = c1t,
              n = n, p = p, Q = Q, funName = funName,
              q.mod = q.mod, m = m, power = power)
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  mlim <- limFun(x = mlim, y = c(6 * ((1 - p) * c1 + p * c1t),
                                 1e6 * ((1 - p) * c1 + p * c1t)))
  if (cost.model){
    if (binary){
      if (two.tailed) {
        pwr <- quote({
          n <- m/(c1*(1-p)+c1t*p);
          1 - pt(qt(1-sig.level/2, df = n-q.mod-4),
                 df = n-q.mod-4, gamma/sqrt((1-r.yx^2)/(n*(p*(1-p)*Q*(1-Q))))) +
            pt(qt(sig.level/2, df = n-q.mod-4),
               df = n-q.mod-4, gamma/sqrt((1-r.yx^2)/(n*(p*(1-p)*Q*(1-Q)))))
        })
      } else {
        pwr <- quote({
          n <- m/(c1*(1-p)+c1t*p);
          1 - pt(qt(1-sig.level, df = n-q.mod-4), df = n-q.mod-4,
                 gamma/sqrt((1-r.yx^2)/(n*(p*(1-p)*Q*(1-Q)))))
        })
      }} else {
        if (two.tailed) {
            pwr <- quote({
              n <- m/(c1*(1-p)+c1t*p);
              lambda <-gamma/sqrt((1 - r.yx^2 - r.ym^2 - r.mx^2 +
                                     2*r.yx*r.ym*r.mx)/(n*(p*(1-p))*(1 - r.mx^2)));
              1 - pt(qt(1 - sig.level/2, df = n-q.mod-4),
                     df = n-q.mod-4, lambda) +
                pt(qt(sig.level/2, df = n-q.mod-4),
                   df = n-q.mod-4, lambda)
            })
         }else{
            pwr <- quote({
              n <- m/(c1*(1-p)+c1t*p);
              lambda <-gamma/sqrt((1 - r.yx^2 - r.ym^2 - r.mx^2 +
                                     2*r.yx*r.ym*r.mx)/(n*(p*(1-p))*(1 - r.mx^2)));
              1 - pt(qt(1 - sig.level, df = n-q.mod-4),
                     df = n-q.mod-4, lambda)
            })
        }
      }
   # cost model power formulas
  } else {
    if (binary){
      if (two.tailed) {
      pwr <- quote({
        1 - pt(qt(1-sig.level/2, df = n-q.mod-4),
                df = n-q.mod-4, gamma/sqrt((1-r.yx^2)/(n*(p*(1-p)*Q*(1-Q))))) +
           pt(qt(sig.level/2, df = n-q.mod-4),
              df = n-q.mod-4, gamma/sqrt((1-r.yx^2)/(n*(p*(1-p)*Q*(1-Q)))))
      })
    } else {
      pwr <- quote({
        1 - pt(qt(1-sig.level, df = n-q.mod-4), df = n-q.mod-4,
               gamma/sqrt((1-r.yx^2)/(n*(p*(1-p)*Q*(1-Q)))))
      })
      }
      } else {
    if (two.tailed) {
        pwr <- quote({
          lambda <-gamma/sqrt((1 - r.yx^2 - r.ym^2 - r.mx^2 +
                                 2*r.yx*r.ym*r.mx)/(n*(p*(1-p))*(1 - r.mx^2)));
          1 - pt(qt(1 - sig.level/2, df = n-q.mod-4),
                 df = n-q.mod-4, lambda) +
            pt(qt(sig.level/2, df = n-q.mod-4),
               df = n-q.mod-4, lambda)
        })
    }else{
        pwr <- quote({
          lambda <-gamma/sqrt((1 - r.yx^2 - r.ym^2 - r.mx^2 +
                                 2*r.yx*r.ym*r.mx)/(n*(p*(1-p))*(1 - r.mx^2)));
            1 - pt(qt(1 - sig.level, df = n-q.mod-4),
                   df = n-q.mod-4, lambda)
        })
      }
     }
   }
      if (is.null(par$power)) {
      out <- list(power = eval(pwr))
    } else if (is.null(par$gamma)) {
      out <- list(gamma = stats::uniroot(
        function(gamma) eval(pwr) - power, gammalim)$root)
    }  else if (is.null(par$m) & is.null(par$n)) {
      if (cost.model) {
        out <- list(m = stats::uniroot(
          function(m) eval(pwr) - power, mlim)$root)
        out <- c(out, list(n = out$m/((1 - p) *c1  + p * c1t)))
      } else {
        out <- list(n = stats::uniroot(function(n) eval(pwr) -
                                               power, nlim)$root)
      }
    }
   power.out <- list(funName = funName, designType = designType,
                     pwr = pwr,
                      par = par, out = out)
    return(power.out)
  }





