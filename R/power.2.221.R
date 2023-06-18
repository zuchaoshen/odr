#' Budget and/or sample size, power calculation for CRTs probing
#'     mediation effects with cluster-level mediators
#' @description This function can calculate required budget for desired power and
#'     power under a fixed budget
#'     for experimental studies with group mediators probing mediation effects.
#'     It also can perform conventional power analyses
#'     (e.g., required sample size and power calculation).
#' @inheritParams power.2
#' @inheritParams od.2.221
#' @param expr returned object from function \code{\link{od.2.221}};
#'     default value is NULL;
#'     if \code{expr} is specified, parameter values of \code{a}, \code{b},
#'     \code{c1}, \code{c1t}, and \code{p}
#'     used or solved in function \code{\link{od.2.221}} will
#'     be passed to the current function;
#'     only the values of \code{p} and \code{n} that specified or solved in
#'     function \code{\link{od.2.221}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint specify the constrained value of
#'     \code{p} and/or \code{n} in a list format to overwrite that/those
#'     from \code{expr}; default value is NULL.
#' @param mlim the range for searching the root of budget (\code{m}) numerically,
#'     default value is the costs sampling \code{nlim}
#'     units across treatment conditions
#'     or c(4 * ncost, 10e10 * ncost) with ncost = ((1 - p) * c1 + p * c1t)
#'
#' @return Required budget (or required sample size), statistical power, or MDES
#'     depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.2.221
#'


power.2.221 <- function(cost.model = TRUE, expr = NULL, constraint = NULL,
                        sig.level = 0.05, two.tailed = TRUE,
                        a = NULL, b = NULL,
                        test = "joint",
                        n = NULL, p = NULL,
                        power= NULL, J = NULL, m = NULL,
                        c1 = NULL, c1t = NULL,
                        c2 = NULL, c2t = NULL,
                        r2m = r2m, r.yx = 0, r.mw = 0, r.yw = 0,
                        icc = NULL,  q = 0,
                        q.a = 0, q.b = 0,
                        powerlim = NULL, Jlim = NULL, mlim = NULL) {
  funName <- "power.2.221"
  designType <- "2-2-1 mediation in two-level CRTs"
  if (cost.model == TRUE) {
    if (sum(sapply(list(m, power), is.null)) != 1)
      stop("exactly one of 'm' and 'power' must be NULL
           when cost.model is TRUE")
    if (!is.null(J))
      stop("'J' must be NULL when cost.model is TRUE")
  } else {
    if (sum(sapply(list(J, power), is.null)) != 1)
      stop("exactly one of 'J' and 'power' must be NULL
           when cost.model is FALSE")
    if (!is.null(m))
      stop("'m' must be NULL when cost.model is FALSE")
  }
  if (!is.null(expr)) {
    if (expr$funName != "od.2.221") {
      stop("'expr' can only be NULL or
            the return from the function of 'od.2.221'")
    } else {
      if (sum(sapply(list(a, b, c1, c1t, c2, c2t, n, p),
                     function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'a', 'b', 'c1', 'c1t', 'c2', 'c2t', 'n', 'p'
             have been specified in expr of 'od.2.221'")
      a <- expr$par$a
      b <- expr$par$b
      c1 <- expr$par$c1
      c1t <- expr$par$c1t
      c2 <- expr$par$c2
      c2t <- expr$par$c2t
      icc <- expr$par$icc
      r2m <- expr$par$r2m
      r.yx <- expr$par$r.yx
      r.yw <- expr$par$r.yw
      r.mw <- expr$par$r.mw
      q.a <- expr$par$q.a
      q.b <- expr$par$q.b
      if(!is.null(test)){
        if (test != expr$par$test)
        {cat('Tests are different in power and the od function ',
             '(', test, ' and ', expr$par$test,
            '). \n', 'The power analysis is for the ', test,
            ' test.', sep = "")}
      } else {
        test <- expr$par$test
      }
      n <- expr$out$n
      p <- expr$out$p
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
    stop("'constraint' must be limited to 'p' and 'n'")
  if (!is.null(constraint$p)) {
    if(NumberCheck(constraint$p) ||
       any (0 >= constraint$p | constraint$p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p <- constraint$p
  }
  if (!is.null(constraint$n)) {
    if(NumberCheck(constraint$n) ||
       (0 >= constraint$n))
      stop("constrained 'n' must be numeric in (0.05, 1e10)")
    n <- constraint$n
  }
  if (sum(sapply(list(p, power, sig.level), function(x) {
    NumberCheck(x) || any(0 > x | x >= 1)
  })) >= 1) stop("'p', 'power', and 'sig.level' must be numeric in (0, 1]")

  if (sum(sapply(list(b), function(x) {
    NumberCheck(x) || any(-1 >= x | x >= 1)
  })) >= 1) stop("'b' must be numeric in (-1, 1)")
  if (cost.model == TRUE){
    if (sum(sapply(list(c1, c1t, c2, c2t), function(x) {
      NumberCheck(x) || x < 0})) >= 1)
      stop("'c1', 'c1t', 'c2', 'c2t' must be numeric in [0, Inf)")
    if (NumberCheck(m))
      stop("'m' must be numeric in [0, Inf)")
  }
  if (NumberCheck(a) || any(-5 > a | a > 5))
    stop("'a' must be numeric in [-5, 5]")
  par <- list(cost.model = cost.model,
              sig.level = sig.level,
              two.tailed = two.tailed,
              a = a, b = b,
              r2m = r2m, r.yx = r.yx, r.yw = r.yw, r.mw = r.mw,
              c1 = c1, c1t = c1t, c2 = c2, c2t = c2t,
              n = n, p = p, J = J, funName = funName,
              q.a = q.a, q.b = q.b, m = m, power = power)
  tside <- ifelse(two.tailed == TRUE, 2, 1)
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  Jlim <- limFun(x = Jlim, y = c(6, 1e6))
  powerlim <- limFun(x = powerlim, y = c(5e-10, 1 - 5e-10))
  mlim <- limFun(x = mlim, y = c(4 * ((1 - p) * (c1 * n + c2) +
                                        p * (c1t * n + c2t)),
                                 1e10 * ((1 - p) * (c1 * n + c2) +
                                           p * (c1t * n + c2t))))

  if (test == "joint" | test == "Joint" | test == "JOINT"){
    if (cost.model){
      # when cost.model is true
      if (two.tailed) {
        pwr <- quote({
          J <- m / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t));
          lambda.a <- a/sqrt((1-r2m)/(p*(1-p)*J));
          lambda.B <- ((b-r.yw*r.mw)/(1-r.mw^2))/
            sqrt((icc*(1-b^2-r.mw^2-r.yw^2+2*b*r.mw*r.yw)/
                    (1-r.mw^2)+(1-icc)*(1-r.yx^2))/
                   (J*(1-r.mw^2)));
          (1 - pt(qt(1-sig.level/tside, df = J-q.a-2),
                  df = J-q.a-2, lambda.a) +
             pt(qt(sig.level/tside, df = J-q.a-2),
                df = J-q.a-2, lambda.a)) *
            (1 - pt(qt(1-sig.level/tside, df = J-q.b-3),
                    df = J-q.b-3, lambda.B) +
               pt(qt(sig.level/tside, df = J-q.b-3),
                  df = J-q.b-3, lambda.B))
        })
      } else {
        pwr <- quote({
          J <- m / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t));
          lambda.a <- a/sqrt((1-r2m)/(p*(1-p)*J));
          lambda.B <- ((b-r.yw*r.mw)/(1-r.mw^2))/
            sqrt((icc*(1-b^2-r.mw^2-r.yw^2+2*b*r.mw*r.yw)/
                    (1-r.mw^2)+(1-icc)*(1-r.yx^2))/
                   (J*(1-r.mw^2)));
          (1 - pt(qt(1-sig.level/tside, df = J-q.a-2),
                  df = J-q.a-2, lambda.a)) *
            (1 - pt(qt(1-sig.level/tside, df = J-q.b-3),
                    df = J-q.b-3, lambda.B))
        })
      }
    } else {
      # when cost.model is not true
      if (two.tailed) {
        pwr <- quote({
          lambda.a <- a/sqrt((1-r2m)/(p*(1-p)*J));
          lambda.B <- ((b-r.yw*r.mw)/(1-r.mw^2))/
            sqrt((icc*(1-b^2-r.mw^2-r.yw^2+2*b*r.mw*r.yw)/
                    (1-r.mw^2)+(1-icc)*(1-r.yx^2))/
                   (J*(1-r.mw^2)));
            (1 - pt(qt(1-sig.level/tside, df = J-q.a-2),
                  df = J-q.a-2, lambda.a) +
             pt(qt(sig.level/tside, df = J-q.a-2),
                df = J-q.a-2, lambda.a)) *
            (1 - pt(qt(1-sig.level/tside, df = J-q.b-3),
                    df = J-q.b-3, lambda.B) +
               pt(qt(sig.level/tside, df = J-q.b-3),
                  df = J-q.b-3, lambda.B))
        })
      } else {
        pwr <- quote({
          lambda.a <- a/sqrt((1-r2m)/(p*(1-p)*J));
          lambda.B <- ((b-r.yw*r.mw)/(1-r.mw^2))/
            sqrt((icc*(1-b^2-r.mw^2-r.yw^2+2*b*r.mw*r.yw)/
                    (1-r.mw^2)+(1-icc)*(1-r.yx^2))/
                   (J*(1-r.mw^2)));
          (1 - pt(qt(1-sig.level/tside, df = J-q.a-2),
                  df = J-q.a-2, lambda.a)) *
            (1 - pt(qt(1-sig.level/tside, df = J-q.b-3),
                    df = J-q.b-3, lambda.B))
        })
      }
    }
    }

    if (is.null(par$power)) {
      out = list(power = eval(pwr))
    } else if (is.null(par$m) & is.null(par$J)) {
      if(cost.model){
        out = list(
        m = stats::uniroot(function(m)
          eval(pwr) - power, mlim)$root)
        out <- list(m = out$m, J = out$m / (((1 - p) * (c1 * n + c2)
                                  + p * (c1t * n + c2t))))
      } else {
        out = list(J = stats::uniroot(function(J)
          eval(pwr) - power, Jlim)$root)
      }
    }
      results <- list(funName = funName,
               designType = designType,
               par = par, test =test, out = out)
    return(results)

}



