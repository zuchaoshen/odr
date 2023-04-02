#' Budget and/or sample size, power, MDES calculation for MRTs investigating
#'     mediation effects with individual-level mediators
#' @description This function can calculate required budget for desired power and
#'     power under a fixed budget
#'     for multisite-randomized trials (MRTs) with individual mediators
#'     probing mediation effects.
#'     It also can perform conventional power analyses
#'     (e.g., required sample size and power calculation).
#' @inheritParams power.2m
#' @inheritParams od.2m.111
#' @param expr returned object from function \code{\link{od.2m.111}}; default value is NULL;
#'     if \code{expr} is specified, parameter values of \code{a}, \code{b},
#'     \code{c1}, \code{c1t}, and \code{p}
#'     used or solved in function \code{\link{od.2m.111}} will
#'     be passed to the current function;
#'     only the values of \code{p} and \code{n} that specified or solved in
#'     function \code{\link{od.2m.111}} can be overwritten
#'     if \code{constraint} is specified.
#' @param constraint specify the constrained value of
#'     \code{p} and/or \code{n} in a list format to overwrite that/those
#'     from \code{expr}; default value is NULL.
#' @param mlim the range for searching the root of budget (\code{m}) numerically,
#'     default value is the costs sampling \code{nlim} units.
#'
#' @return Required budget (or required sample size), statistical power, or MDES
#'     depending on the specification of parameters.
#'     The function also returns the function name, design type,
#'     and parameters used in the calculation.
#'
#' @export power.2m.111
#'

power.2m.111 <- function(cost.model = TRUE, expr = NULL, 
                         constraint = NULL,
                        sig.level = 0.05, two.tailed = TRUE,
                        a = NULL, b = NULL,
                        power = NULL, m = NULL, test = NULL,
                        n = NULL, p = NULL,
                        c1 = NULL, c1t = NULL, 
                        c2 = NULL, 
                        r12 = 0, r22m = 0, r12m = 0,
                        icc.m = NULL, omega = NULL,
                        icc = NULL, J = NULL, q = 0,
                        q.a = 0, q.b = 0, max.iter = 300,
                        powerlim = NULL, Jlim = NULL, mlim = NULL,
                        rounded = TRUE) {
  funName <- "power.2m.111"
  designType <- "1-1-1 mediation in 2-level MRTs"
  if (cost.model) {
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
    if (expr$funName != "od.2m.111") {
      stop("'expr' can only be NULL or
            the return from the function of 'od.2.221'")
    } else {
      if (sum(sapply(list(a, b, c1, c1t, c2,icc, icc.m, omega, n, p),
                     function(x) {!is.null(x)})) >= 1)
        stop("parameters of 'a', 'b', 'c1', 'c1t', 'c2', 'icc',
            'icc.m', 'omega', 'n', and 'p'
             have been specified in expr of 'od.2m.111'")
      a <- expr$par$a
      b <- expr$par$b
      c1 <- expr$par$c1
      c1t <- expr$par$c1t
      c2 <- expr$par$c2
      icc <- expr$par$icc
      icc.m <-expr$par$icc.m
      r12 <- expr$par$r12
      r22 <- expr$par$r22
      r12m <- expr$par$r12m
      r22m <- expr$par$r22m
      omega <- expr$par$omega
      q.a <- expr$par$q.a
      q.b <- expr$par$q.b
      
      if(is.null(test)){
        test <- expr$par$test
      } else {
        if (test != expr$par$test)
        {cat('Tests are different in the power and the od functions ', '(', test, ' and ', expr$par$test,
             '). \n', 'The power analysis is for the ', test, ' test.', sep = "")}
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
    NumberCheck(x) || any(-1 > x | x > 1)
  })) >= 1) stop("'b' must be numeric in [-1, 1]")
  if (cost.model == TRUE){
    if (sum(sapply(list(c1, c1t, c2), function(x) {
      NumberCheck(x) || x < 0})) >= 1)
      stop("'c1', 'c1t', 'c2' must be numeric")
    if (NumberCheck(m))
      stop("'m' must be numeric in [0, Inf)")
  }
  if (NumberCheck(a) || any(-5 > a | a > 5))
    stop("'a' must be numeric in [-5, 5]")
  par <- list(cost.model = cost.model,
              sig.level = sig.level,
              two.tailed = two.tailed,
              a = a, b = b, icc = icc, icc.m = icc.m,
              omega = omega,
              r12 = r12, r12m = r12m, r22m = r22m,
              c1 = c1, c1t = c1t, c2 = c2,
              n = n, p = p, J = J, funName = funName,
              q = q, m = m, power = power)
  tside <- ifelse(two.tailed == TRUE, 2, 1)
  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  Jlim <- limFun(x = Jlim, y = c(4, 1e6))
  powerlim <- limFun(x = powerlim, y = c(5e-10, 1 - 5e-10))
  mlim <- limFun(x = mlim, y = c(4 * ((1 - p) * (c1 * n) + p * c1t * n + c2),
                                 1e10 * ((1 - p) * (c1 * n) + p * c1t * n + c2)))
  if (test == "sobel" | test == "Sobel" | test == "SOBEL"){
    if (cost.model){
      if (two.tailed){
        pwr.sobel <- quote({
          var.ab <- 
            a^2*((1-icc)*(1-r12)/
                   ((m / ((1 - p) * n * c1+ p * n * c1t + c2))*n*(1-icc.m)*(1-r12m)))+
            b^2*((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                   (p*(1-p)*n*(m / ((1 - p) * n * c1+ p * n * c1t + c2))));
          1 - pnorm(qnorm(1 - sig.level/tside),a*b/sqrt(var.ab))+
         pnorm(qnorm(sig.level/tside),a*b/sqrt(var.ab))
        })
      }else{
        pwr.sobel <- quote({
          var.ab <- 
            a^2*((1-icc)*(1-r12)/
                   ((m / ((1 - p) * n * c1+ p * n * c1t + c2))*n*(1-icc.m)*(1-r12m)))+
            b^2*((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                   (p*(1-p)*n*(m / ((1 - p) * n * c1+ p * n * c1t + c2))));
          1 - pnorm(qnorm(1 - sig.level/tside),a*b/sqrt(var.ab))
        })
      }
    } else{
      if (two.tailed){
        pwr.sobel <- quote({
          var.ab <- 
            a^2*((1-icc)*(1-r12)/
                   (J*n*(1-icc.m)*(1-r12m)))+
            b^2*((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                   (p*(1-p)*n*J));
          1 - pnorm(qnorm(1 - sig.level/tside),a*b/sqrt(var.ab))+
          pnorm(qnorm(sig.level/tside),a*b/sqrt(var.ab))
        })
      }else{
        pwr.sobel <- quote({
          var.ab <- 
            a^2*((1-icc)*(1-r12)/
                   (J*n*(1-icc.m)*(1-r12m)))+
            b^2*((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                   (p*(1-p)*n*J));
          1 - pnorm(qnorm(1 - sig.level/tside),
                    a*b/sqrt(var.ab))
        })
      }
    }
    if (is.null(power)) {
      sobel.out <- list(power = eval(pwr.sobel))
    } else if (is.null(m) & is.null(J)) {
      if (cost.model) {
        sobel.out <- list(m = stats::uniroot(
          function(m) eval(pwr.sobel) - power, mlim)$root)
        sobel.out <- c(sobel.out, 
                       list(J = sobel.out$m/
                              ((1 - p) *c1 * n  + p * c1t * n+ c2)))
      } else {
        sobel.out <- list(J = stats::uniroot(function(J) eval(pwr.sobel) - 
                                               power, Jlim)$root)
      }
    }
    power.out <- list(funName = funName, designType = designType, 
                      par = par, sobel.out = sobel.out)
    return(power.out)
  }
  if (test == "joint" | test == "Joint" | test == "JOINT"){
    if (cost.model){
      # when cost.model is true for the Joint test
      if (two.tailed) {
        pwr.joint <- quote({
          J <- m / ((1 - p) * n * c1+ p * n * c1t + c2);
          lambda.a <- a/sqrt((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                               (p*(1-p)*n*J));
          lambda.b <- b/sqrt((1-icc)*(1-r12)/
                               (J*n*(1-icc.m)*(1-r12m)));
          (1 - pt(qt(1 - sig.level/tside, df = J-q.a-1),
                  df = J-q.a-1, lambda.a) +
             pt(qt(sig.level/tside, df = J-q.a-1),
                df = J-q.a-1, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J*(n-1) -q.b),
                    df = J*(n-1) -q.b, lambda.b) +
               pt(qt(sig.level/tside, df = J*(n-1) -q.b),
                  df = J*(n-1) -q.b, lambda.b))
        })
      } else {
        pwr.joint <- quote({
          J <- m / ((1 - p) * n * c1+ p * n * c1t + c2);
          lambda.a <- a/sqrt((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                               (p*(1-p)*n*J));
          lambda.b <- b/sqrt((1-icc)*(1-r12)/
                               (J*n*(1-icc.m)*(1-r12m)));
          (1 - pt(qt(1 - sig.level/tside, df = J-q.a-1),
                  df = J-q.a-1, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J*(n-1) - q.b),
                    df = J*(n-1) - q.b, lambda.b))
        })
      }
    } else {
      # when cost.model is not true for the Joint test
      if (two.tailed) {
        pwr.joint <- quote({
          lambda.a <- a/sqrt((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                               (p*(1-p)*n*J));
          lambda.b <- b/sqrt((1-icc)*(1-r12)/
                               (J*n*(1-icc.m)*(1-r12m)));
          (1 - pt(qt(1 - sig.level/tside, df = J-q.a-1),
                  df = J-q.a-1, lambda.a) +
             pt(qt(sig.level/tside, df = J-q.a-1),
                df = J-q.a-1, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J*(n-1) -q.b),
                    df = J*(n-1) -q.b, lambda.b) +
               pt(qt(sig.level/tside, df = J*(n-1) -q.b),
                  df = J*(n-1) -q.b, lambda.b))
        })
      } else {
        pwr.joint <- quote({
          lambda.a <- a/sqrt((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                               (p*(1-p)*n*J));
          lambda.b <- b/sqrt((1-icc)*(1-r12)/
                               (J*n*(1-icc.m)*(1-r12m)));
          (1 - pt(qt(1 - sig.level/tside, df = J-q.a-1),
                  df = J-q.a-1, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J*(n-1) -q.b),
                    df = J*(n-1) -q.b, lambda.b))
        })
      }
    }

    if (is.null(power)) {
      joint.out <- list(power = eval(pwr.joint))
    } else if (is.null(m) & is.null(J)) {
      if(cost.model){
        joint.out <- list(m = stats::uniroot(function(m)
          eval(pwr.joint) - power, mlim)$root);
        joint.out <- c(joint.out, list(J = joint.out$m / ((1 - p) *n* c1
                                                + p *n* c1t  + c2)))
      } else {
        joint.out <- list(J = stats::uniroot(function(J)
          eval(pwr.joint) - power, Jlim)$root);
      }

    }
    power.out <- list(funName = funName,
                      designType = designType,
                      par = par, joint.out = joint.out)
    return(power.out)
  }

}





