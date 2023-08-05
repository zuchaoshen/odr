#' Budget and/or sample size, power, MDES calculation for MRTs investigating
#'     mediation effects with individual-level mediators
#' @description This function can calculate required budget for desired power and
#'     power under a fixed budget
#'     for multisite-randomized trials (MRTs) with individual mediators
#'     probing mediation effects.
#'     It also can perform conventional power analyses
#'     (e.g., required sample size and power calculation).
#' @inheritParams od.1.111
#' @inheritParams power.4
#' @inheritParams od.4
#' @param expr Returned object from function \code{\link{od.1.111}}; default value is NULL;
#'     if \code{expr} is specified, parameter values of \code{a}, \code{b},
#'     \code{c}, \code{ct}, and \code{p}
#'     used or solved in function \code{\link{od.1.111}} will
#'     be passed to the current function;
#'     only the values of \code{p} that specified or solved in
#'     function \code{\link{od.1.111}} can be overwritten
#'     if \code{constraint} is specified.
#' @param cost.model Logical; power analyses accommodating costs and budget
#'     (e.g., required budget for a desired power, power under fixed budget)
#'     if TRUE. Otherwise, conventional power analyses are performed
#'     (e.g., required sample size and power calculation); default value is TRUE.
#' @param constraint If specified, the constrained value of
#'     \code{p} in a list format (e.g., constraint = list(p = 0.5))
#'     will overwrite that
#'     from \code{expr}; default value is NULL.
#' @param mlim The range for identifying the root of budget (\code{m}) numerically,
#'     default value is the costs sampling \code{nlim} units.
#' @param alim The range for identifying the root of a path
#'     effect (\code{a}) numerically. Default value is c(0, 4).
#' @param blim The range for identifying the root of b path within-treatment
#'     correlation between the mediator and outcome (\code{b}) numerically.
#'     Default value is (.01, .99), if (\code{b}) is negative, please re-code
#'     the outcome or mediator to make it positive.
#' @return Required budget (or required sample size), statistical power,
#'     (\code{a}) , or (\code{b})  depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.1.111
#'
#' @examples
#' # Optimal design and power analyses accommodating costs and budget
#' myod <- od.1.111(a = .3, b = .5, c = 10, ct = 100)
#' # myod
#' mypower <- power.1.111(expr = myod, power = .8)
#' #mypower
#'
#' # Conventional power analyses
#' mypower <- power.1.111(cost.model = FALSE, a = .3, b = .5,
#'                        power = .8, p =.5)
#' #mypower
#' mypower <- power.1.111(cost.model = FALSE, n = 350, b = .5,
#'                        power = .8, p =.5)
#' #mypower
#'

power.1.111 <- function(cost.model = TRUE, expr = NULL,
                        constraint = NULL,
                        sig.level = 0.05, two.tailed = TRUE,
                        a = NULL, b = NULL,
                        power = NULL, m = NULL, test = "joint",
                        n = NULL, p = NULL,
                        c = NULL, ct = NULL,
                        r.yx = 0, r.mx = 0, r.mw = 0,
                        q.a = 0, q.b = 0, max.iter = 300,
                        alim = c(0, 4), blim = c(.01, .99),
                        powerlim = NULL, nlim = c(6, 1e7), mlim = NULL
                        ) {
  funName <- "power.1.111"
  designType <- "1-1-1 mediation in single-level RCTs"
  if (!is.null(expr)) {
    if (expr$funName != "od.1.111") {
      stop("'expr' can only be NULL or
            the return from the function of 'od.1.111'")
    } else {
      if (sum(sapply(list(c, ct, n, p),
                     function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'a', 'b', 'c', 'ct',
            'n', and 'p'
             have been specified in expr of 'od.1.111'")
      if (is.null(a)){a <- expr$par$a}
      if (is.null(b)){b <- expr$par$b}
      c <- expr$par$c
      ct <- expr$par$ct
      r.yx <- expr$par$r.yx
      r.mx <- expr$par$r.mx
      r.mw <- expr$par$r.mw
      q.a <- expr$par$q.a
      q.b <- expr$par$q.b

      if(is.null(test)){
        test <- expr$par$test
      } else {
        if (test != expr$par$test)
        {cat('Tests are different in the power and the od functions ',
             '(', test, ' and ', expr$par$test,
             '). \n', 'The power analysis is for the ',
             test, ' test.', sep = "")}
      }
      p <- expr$out$p
    }
  } else {
    if (!is.null(constraint))
      stop("'constraint' must be NULL when 'expr' is NULL")
  }
  if (cost.model) {
    if (sum(sapply(list(m, a, b, power), is.null)) != 1)
      stop("exactly one of 'm', 'a', 'b' and 'power' must be NULL
           when cost.model is TRUE")
    if (!is.null(n))
      stop("'n' must be NULL when cost.model is TRUE")
  } else {
    if (sum(sapply(list(n, a, b, power), is.null)) != 1)
      stop("exactly one of 'n', 'a', 'b', and 'power' must be NULL
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

  if (sum(sapply(list(b), function(x) {
    NumberCheck(x) || any(-1 > x | x > 1)
  })) >= 1) stop("'b' must be numeric in [-1, 1]")
  if (cost.model == TRUE){
    if (sum(sapply(list(c, ct), function(x) {
      NumberCheck(x) || x < 0})) >= 1)
      stop("'c', 'ct' must be numeric")
    if (NumberCheck(m))
      stop("'m' must be numeric in [0, Inf)")
  }
  if (NumberCheck(a) || any(-5 > a | a > 5))
    stop("'a' must be numeric in [-5, 5]")
  par <- list(cost.model = cost.model,
              sig.level = sig.level,
              two.tailed = two.tailed,
              a = a, b = b,
              r.yx = r.yx, r.mx = r.mx, r.mw = r.mw,
              c = c, ct = ct,
              n = n, p = p, funName = funName,
              q.a = q.a, q.b = q.b, m = m, power = power)
  tside <- ifelse(two.tailed == TRUE, 2, 1)
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  nlim <- limFun(x = nlim, y = c(5, 1e6))
  powerlim <- limFun(x = powerlim, y = c(5e-10, 1 - 5e-10))
  mlim <- limFun(x = mlim, y = c(5 * ((1 - p) * c + p * ct),
                                 1e6 * ((1 - p) * c + p * ct)))


    if (cost.model){
      if (test == "joint" | test == "Joint" | test == "JOINT") {
        if (two.tailed) {
          pwr <- quote({
            B <- (b-r.yx*r.mx)/(1-r.mx^2);
            n <- m/(p*ct + (1-p)*c);
            se.a <- sqrt((1-r.mw^2)/(p*(1-p)*n));
            se.B <- sqrt((1-b^2-r.mx^2-r.yx^2+2*b*r.mx*r.yx)/
                           (n*(1-r.mx^2)*(1-r.mx^2)));
            (1 - pt(qt(1 - sig.level/tside, df = n-q.a-2),
                    df = n-q.a-2, a/se.a) +
               pt(qt(sig.level/tside, df = n-q.a-2),
                  df = n-q.a-2, a/se.a)) *
              (1 - pt(qt(1 - sig.level/tside, df = n-q.b-3),
                      df = n-q.b-3, B/se.B) +
                 pt(qt(sig.level/tside, df = n-q.b-3),
                    df = n-q.b-3, B/se.B))
          })
        } else {
          pwr <- quote({
            B <- (b-r.yx*r.mx)/(1-r.mx^2);
            n <- m/(p*ct + (1-p)*c);
            se.a <- sqrt((1-r.mw^2)/(p*(1-p)*n));
            se.B <- sqrt((1-b^2-r.mx^2-r.yx^2+2*b*r.mx*r.yx)/
                           (n*(1-r.mx^2)*(1-r.mx^2)));
            (1 - pt(qt(1 - sig.level/tside, df = n-q.a-2),
                    df = n-q.a-2, a/se.a)) *
              (1 - pt(qt(1 - sig.level/tside, df = n-q.b-3),
                      df = n-q.b-3, B/se.B))
          })
        }
      } else if (test == "sobel" | test == "Sobel" | test == "SOBEL"){

        if (two.tailed) {
          pwr <- quote({
            B <- (b-r.yx*r.mx)/(1-r.mx^2);
            n <- m/(p*ct + (1-p)*c);
            se.a <- sqrt((1-r.mw^2)/(p*(1-p)*n));
            se.B <- sqrt((1-b^2-r.mx^2-r.yx^2+2*b*r.mx*r.yx)/
                           (n*(1-r.mx^2)*(1-r.mx^2)));
            z.sobel <- a*B/sqrt(a^2*se.B^2+B^2*se.a^2);
            1-pnorm(qnorm(1-sig.level/tside)-z.sobel)+
              pnorm(qnorm(sig.level/tside)-z.sobel)
          })
        } else {
          pwr <- quote({
            B <- (b-r.yx*r.mx)/(1-r.mx^2);
            n <- m/(p*ct + (1-p)*c);
            se.a <- sqrt((1-r.mw^2)/(p*(1-p)*n));
            se.B <- sqrt((1-b^2-r.mx^2-r.yx^2+2*b*r.mx*r.yx)/
                           (n*(1-r.mx^2)*(1-r.mx^2)));
            z.sobel <- a*B/sqrt(a^2*se.B^2+B^2*se.a^2);
            1-pnorm(qnorm(1-sig.level/tside)-z.sobel)
          })
        }
      }
    } else{
      if (test == "joint" | test == "Joint" | test == "JOINT") {
        if (two.tailed) {
          pwr <- quote({
            B <- (b-r.yx*r.mx)/(1-r.mx^2);
            se.a <- sqrt((1-r.mw^2)/(p*(1-p)*n));
            se.B <- sqrt((1-b^2-r.mx^2-r.yx^2+2*b*r.mx*r.yx)/
                           (n*(1-r.mx^2)*(1-r.mx^2)));
            (1 - pt(qt(1 - sig.level/tside, df = n-q.a-2),
                    df = n-q.a-2, a/se.a) +
               pt(qt(sig.level/tside, df = n-q.a-2),
                  df = n-q.a-2, a/se.a)) *
              (1 - pt(qt(1 - sig.level/tside, df = n-q.b-3),
                      df = n-q.b-3, B/se.B) +
                 pt(qt(sig.level/tside, df = n-q.b-3),
                    df = n-q.b-3, B/se.B))
          })
        } else {
          pwr <- quote({
            B <- (b-r.yx*r.mx)/(1-r.mx^2);
            se.a <- sqrt((1-r.mw^2)/(p*(1-p)*n));
            se.B <- sqrt((1-b^2-r.mx^2-r.yx^2+2*b*r.mx*r.yx)/
                           (n*(1-r.mx^2)*(1-r.mx^2)));
            (1 - pt(qt(1 - sig.level/tside, df = n-q.a-2),
                    df = n-q.a-2, a/se.a)) *
              (1 - pt(qt(1 - sig.level/tside, df = n-q.b-3),
                      df = n-q.b-3, B/se.B))
          })
        }
      } else if (test == "sobel" | test == "Sobel" | test == "SOBEL"){

        if (two.tailed) {
          pwr <- quote({
            B <- (b-r.yx*r.mx)/(1-r.mx^2);
            se.a <- sqrt((1-r.mw^2)/(p*(1-p)*n));
            se.B <- sqrt((1-b^2-r.mx^2-r.yx^2+2*b*r.mx*r.yx)/
                           (n*(1-r.mx^2)*(1-r.mx^2)));
            z.sobel <- a*B/sqrt(a^2*se.B^2+B^2*se.a^2);
            1-pnorm(qnorm(1-sig.level/tside)-z.sobel)+
              pnorm(qnorm(sig.level/tside)-z.sobel)
          })
        } else {
            pwr <- quote({
            B <- (b-r.yx*r.mx)/(1-r.mx^2);
            se.a <- sqrt((1-r.mw^2)/(p*(1-p)*n));
            se.B <- sqrt((1-b^2-r.mx^2-r.yx^2+2*b*r.mx*r.yx)/
                           (n*(1-r.mx^2)*(1-r.mx^2)));
            z.sobel <- a*B/sqrt(a^2*se.B^2+B^2*se.a^2);
            1-pnorm(qnorm(1-sig.level/tside)-z.sobel)
          })
        }
      }
    }

    if (is.null(power)) {
      out <- list(power = eval(pwr))
    } else if (is.null(a)) {
      out <- list(a = stats::uniroot(
        function(a) eval(pwr) - power, alim)$root)
    }  else if (is.null(b)) {
      out <- list(b = stats::uniroot(
        function(b) eval(pwr) - power, blim)$root)
    } else if (is.null(m) & is.null(n)) {
      if (cost.model) {
        out <- list(m = stats::uniroot(
          function(m) eval(pwr) - power, mlim)$root)
        out <- c(out, list(n = out$m/((1 - p) *c  + p * ct)))
      } else {
        out <- list(n = stats::uniroot(function(n) eval(pwr) -
                                               power, nlim)$root)
      }
    }
   power.out <- list(funName = funName, designType = designType,
                      par = par, test=test, out = out)
    return(power.out)
  }





