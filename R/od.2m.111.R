#' Optimal sample allocation calculation for two-level multisite-randomized
#'     trials investigating mediation effects with individual-level mediators
#'     (1-1-1)
#'
#' @description The optimal design of two-level
#'     multisite-randomized trials (MRTs) probing mediation effects with
#'     individual-level mediators, for the
#'     Sobel test, is to calculate
#'     the optimal sample allocation that minimizes the variance of
#'     a mediation effect under a fixed budget. For the joint significance test, it is to identify
#'     the optimal sample allocation that requires the minimum budget
#'     to achieve certain power level.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n})
#'     and the proportion of level-1 individuals/units to be assigned to treatment (\code{p}).
#'     This function solves the optimal \code{n} and/or \code{p}
#'     with and without a constraint.
#'
#' @inheritParams power.2m.111
#' @inheritParams od.4m
#' @param a The treatment effect on the mediator.
#' @param b The within treatment correlation between the outcome and the mediator.
#' @param q.b The number of covariates in the outcome model (except the treatment indicator
#'     and the mediator).
#' @param q.a The number of covariates at the individual level of the mediator model
#'     (except the treatment indicator).
#' @param omega The treatment-by-site variance of the outcome.
#' @param test The type of test will be used to detect mediation effects. Default is
#'     the joint significance test (i.e., test = "joint").
#'     Another choice is the Sobel test by specifying the argument as test = "sobel".
#' @param power Statistical power specified, default is .80.
#' @param r12m The proportion of within treatment mediator variance at the level one
#'     explained by covariates.
#' @param r22m The proportion of treatment-by-site variance 
#'     explained by covariates.
#' @param r12 The proportion of within treatment individual-level outcome variance
#'     explained by covariates.
#' @param icc.m The intraclass correlation coefficient for the mediator.
#' @param e Maximum error value used when solution quality used as
#'     the stopping criterion, default is 1e-10.
#' @param max.value Maximal value of optimization when used as
#'     the stopping criterion. Default is infinite.
#' @param d.p The initial sampling domains for p. Default is c(0.10, 0.50).
#' @param d.n The initial sampling domain for n. Default is c(4, 500).
#' @param max.iter Maximal number of function evaluations when used as
#'     the stopping criterion. Default is 200.
#' @param n.of.archive  Size of the solution archive, default is 100.
#' @param q Locality of the search (0,1), default is 0.0001.
#' @param xi  Convergence pressure (0, Inf), suggested: (0, 1), default is 0.5.
#' @param verbose Print out evaluation process if TRUE, default is TRUE.
#' @param Jlim The range for J to search for a numerical solution. Default is c(3, 10e4).
#' @param n.of.ants Number of ants used in each iteration after the initialization
#'     of power analysis for calculating required budget, default value is 10.
#' @param iter number of iteration used for solving roots in the Sobel test.
#' @param tol convergence tolerance.
#'
#' @return
#'     Unconstrained or constrained optimal sample allocation (\code{n} and \code{p}).
#'     The function also returns statistical power,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.2m.111
#'

od.2m.111 <- function(a = NULL, b = NULL, icc.m = NULL, icc = NULL,
                     c1 = NULL, c1t = NULL, c2 = NULL, m = NULL,
                     r12m = 0, r22m = 0, r12 = 0,
                     omega = 0.01,
                     q.a = 0, q.b = 3,
                     test = "joint",
                     n = NULL, p = NULL, iter = 100,
                     tol = 1e-11, power = 0.80,
                     d.p = c(0.1, 0.5), d.n = c(5, 50),
                     sig.level = 0.05, two.tailed = TRUE,
                     plots = TRUE, 
                     nlim = c(4, 100), plim = c(.01, .99),
                     varlim = c(0, 0.001),
                     nlab = NULL, plab = NULL, varlab = NULL,
                     vartitle = NULL,
                     Jlim = c(3, 10e4), verbose = TRUE,
                     max.value = Inf, max.iter = 300,  e = 1e-10,
                     n.of.ants = 10, n.of.archive = 50, q = 0.0001,
                     xi = 0.5,
                     plot.by = list(n = "n", p = "p")
) {
  funName <- "od.2m.111"
  designType <- "1-1-1 mediation in 2-level MRTs"
  par <- list(a = a, b = b, icc = icc, icc.m = icc.m,
              r12 = r12, r12m = r12m, r22m = r22m,
              omega = omega,
              c1 = c1, c2 = c2, c1t =c1t,
              n = n, p = p, m = m,
              q.a = q.a, q.b = q.b,
              sig.level = sig.level, two.tailed = two.tailed,
              test = test,
              max.iter = max.iter,
              n.of.ants = n.of.ants, n.of.archive = n.of.archive,
              q = q,
              xi = xi
  )
  if (sum(sapply(list(icc, icc.m, r12m, r22m, r12, c1, c2, c1t),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc', 'icc.m', 'r12m', 'r22m', 'r12', 'c1', 'c2',
         'c1t' must be specified")
  NumberCheck <- function(x) {!is.null(x) & !is.numeric(x)}
  if (NumberCheck(icc) | any(0 > icc | icc > 1))
    stop("'icc' must be numeric in [0, 1]")
  if (sum(sapply(list(r12m, r22m, r12), function(x) {
    NumberCheck(x) | any(0 > x | x > 1)
  })) >= 1)
    stop("'r12m', 'r22m', 'r12' must be numeric in [0, 1]")
  if (sum(sapply(list(c1, c2, c1t), function(x) {
    NumberCheck(x) | x < 0})) >= 1)
    stop("'c1', 'c2', 'c1t' must be numeric")
  if (!is.null(plot.by) & !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n'))")
  if (c1 == 0 & c1t == 0 & is.null(par$n) & is.null(par$p))
    stop("when c1 and c1t are both zero, n and/or p must be constrained,
         please specify a value for n or p")
  if (c2 == 0 & is.null(par$n) & is.null(par$p))
    stop("when c2 is zero, n and/or p must be constrained,
         please specify a value for n or p")

  labFun <- function(x, y) {
    if (!is.null(x) & length(x) == 1 & is.character(x)) {x} else {y}
  }

  tside <- ifelse(two.tailed == TRUE, 2, 1)

# Sobel Test----  
  if (test == "sobel" | test == "Sobel" | test == "SOBEL") {
    
    var.expr <- quote({
      J <- m / (c1 * (1 - p) * n + c1t * p * n + c2);
      a^2 * ((1 - icc) *(1 - r12)/(J * n * (1 - icc.m) * (1 - r12m))) +
        b^2 * (omega * (1 - r22m) + (1 - icc.m) * (1 - r12m) / (n * p * (1 - p)))/J
    })
      
    if (is.null(par$p)){
      p.expr <- quote({
        (a^2 * (1 - icc) * (1 - r12) * ( c1t * n - c1 * n) *
          (p - p^2) +
          b^2 * (omega * (1 - r22m) * n * p * (1 - p) + 
                   (1 - icc.m) * (1 - r12m)) *
          (c1t * n - c1 * n) * (1 - icc.m) * (1 - r12m)) * (p * (1 - p)) - 
          (1 - 2 * p) * (b^2 * (1 - icc.m) * (1 - r12m) * (c1 * (1 - p) * n + 
           c1t * p * n + c2) * (1 - icc.m) * (1 - r22m))
      })
    } else {
      if (!is.numeric(par$p) || any(par$p <= 0 | par$p >= 1))
        stop("constrained 'p' must be numeric in (0, 1)")
      p.expr <- par$p
    }
    
    if (is.null(par$n)){
      n.expr <- quote({
        sqrt((a^2 * (1 - icc) * (1 - r12) * c2 * p * (1 - p) + 
                b^2 * (1 - icc.m) * (1 - r12m) * c2 * (1 - icc.m) * (1 - r12m))/
          (b^2 * omega* (1 - r22m) * p * (1 - p) * (c1 * (1 - p) +
            c1t * p) * 
             (1 - icc.m) * (1 - r12m))
            )
      })
      n <- sample(2:50, 1)
    } else {
      if (!is.numeric(par$n) || par$n <= 0)
        stop("constrained 'n' must be numeric with n > 0")
      n.expr <- par$n
    }
    
    nn <- pp <- NULL
    for (i in 1:iter) {
      if (is.null(par$p)) {
        pp[i] <- p <- stats::uniroot(function(p)
          eval(p.expr), plim)$root
      } else {
        pp[i] <- par$p
      }
      if (is.null(par$n)){
        nn[i] <- n <- eval(n.expr)
      } else{
        nn[i] <- par$n
      }
      
    }
    if (!is.null(par$n) && !is.null(par$p)) {
      cat("===============================\n",
          "Both n and p are constrained, there is no calculation from other parameters",
          ".\n===============================\n", sep = "")
    }
    if (verbose) {
      if (!is.null(par$n)) {
        cat("The constrained level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
      } else {
        cat("The optimal level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
      }
      if (!is.null(par$p)) {
        cat("The constrained proportion of level-1 units in treatment (p) is ", p, ".\n", "\n", sep = "")
      } else {
        cat("The optimal proportion of level-1 units in treatment (p) is ", p, ".\n", "\n" ,sep = "")
      }
    }
    if (nn[iter] - nn[iter-1] <= tol && pp[iter] - pp[iter-1] <= tol) {
      p <- pp[iter]
      nn <- pp <- NULL
    } else {
      cat("===============================\n",
          "The solutions are not converged to specified tolerance,
        please specify a large numer of 'iter' to replace the default value of 100",
          ".\n===============================\n", sep = "")
    }
    m <- ifelse(is.null(par$m), 50*(c1t*p*n+(1-p)*n*c1+c2),par$m)
    Var <- eval(var.expr)
    par <- c(par, list(m = m))
    out <- list(n = n, p = p, var = Var)
    od.out <- list(funName = funName, designType = designType,
                   par = par, out = out)
    labFun <- function(x, y) {
      if (!is.null(x) && length(x) == 1 && is.character(x)) {x} else {y}
    }
    nlab <- labFun(x = nlab, y = "Level-One Sample Size: n")
    plab <- labFun(x = plab, y = "Proportion: p")
    varlab <- labFun(x = varlab, y = "Variance")
    vartitle <- labFun(x = vartitle, y = "")
    plotbyFun <- function(x, y) {
      if (!is.null(x) && is.list(x)) {x} else {y}
    }
    plot.by <- plotbyFun(x = plot.by, y = list(n = "n", p = "p"))
    nrange <- seq(nlim[1], nlim[2], by = 1)
    prange <- seq(plim[1] + 0.05, plim[2] - 0.05, by = 0.01)
    if (length(plot.by) == 2) figure <- par(mfrow = c (1, 2))
    if (length(plot.by) == 1) figure <- par(mfrow = c (1, 1))
    if (plots == TRUE) {
     if (!is.null(plot.by$n)) {
        plot.y <- NULL
        for (n in nrange)
          plot.y <- c(plot.y, eval(var.expr))
        graphics::plot(nrange, plot.y,
                       type = "l", lty = 1,
                       xlim = nlim, ylim = varlim,
                       xlab = nlab, ylab = varlab,
                       main = vartitle, col = "black")
        n <- out$n
        graphics::abline(v = n, lty = 2, col = "Blue")
      }
      if (!is.null(plot.by$p)) {
        plot.y <- NULL
        for (p in prange)
          plot.y <- c(plot.y, eval(var.expr))
        graphics::plot(prange, plot.y,
                       type = "l", lty = 1,
                       xlim = plim, ylim = varlim,
                       xlab = plab, ylab = varlab,
                       main = vartitle, col = "black")
        p <- out$p
        graphics::abline(v = p, lty = 2, col = "Blue")
      }
    }
    par(figure)
    return(od.out) 
  }
  
  
  
# Joint test----
  if (test == "joint" | test == "Joint" | test == "JOINT") {
    if (is.null(par$p) & is.null(par$n)) {
      n.of.opt.pars <- 2
      if (verbose) {cat('The ACO algorithm started initilization..', ".\n", sep = "")}
      e.abs <- e # absolute error
      e.rel <- e # relative error
      # initiate parameters
      eval <- 0
      last.impr <- max.iter
      design.pars <- data.frame()
      outcome <- vector()
      max.X <- rep(NA, n.of.opt.pars)
      max.y <- -Inf
      p.X <- vector()
      pp <- data.frame(v = numeric(), sd = numeric(), gr = numeric());
      outcome <- NULL

      #Power under joint significance test for two-tailed test
      if (two.tailed) {
        pwr <- quote({
          lambda.a <- a/sqrt((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                               (p*(1-p)*n*J));
          lambda.b <- b/sqrt((1-icc)*(1-r12)/
                               (J*n*(1-icc.m)*(1-r12m)));
          (1 - pt(qt(1 - sig.level/tside, df = J-q.a-1),
                  df = J-q.a-1, lambda.a) +
             pt(qt(sig.level/tside, df = J-q.a-1),
                df = J-q.a-1, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J*(n-1)-q.b),
                    df = J*(n-1)-q.b, lambda.b) +
               pt(qt(sig.level/tside, df = J*(n-1) -q.b),
                  df = J*(n-1) -q.b, lambda.b))
        })
      } else {
        pwr <- quote({
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

      n.of.initial <- round(sqrt(n.of.archive), 0)
      n.initial <- seq(from = d.n[1], to = d.n[2], length = n.of.initial)
      p.initial <- seq(from = d.p[1], to = d.p[2], length = n.of.initial)
      n.of.archive <- n.of.initial^2
      nl <- matrix(NA, n.of.archive, n.of.archive-1)
      X <- NULL
      p.X <- NULL
      y <- NULL
      budget <- NULL
      for (n in n.initial){
        for (p in p.initial){
          X <- rbind(X, c(p, n))
          p.X <- rbind(p.X, c(p, n))
          J <- stats::uniroot(function(J) eval(pwr) - power, c(Jlim[1], Jlim[2]))$root
          m <- (p*n*c1t + (1-p)*n*c1 + c2) * J
          y <- c(y, 1/m)
          budget <- c(budget, m)

        }
      }
      pp <- rbind(pp, data.frame(v = y, sd = 0, gr = 0, m = budget))
      pp$gr <- rank(-pp$v, ties.method = "random")

      for (i in 1:n.of.archive){
        nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]
      }
      #  colnames(p.X) <- c("p", "n")
      #  colnames(X) <- c("p", "n")
      #  p.X <- as.data.frame(p.X)
      n.iter <- n.of.archive

      if (verbose)
      {cat('The ACO algorithm finished initilization of ', n.of.archive, ' analyses',".\n", sep = "")}

      while (TRUE) { # the algorithm will stop if one of the stop criteria is met
        dist.mean <- p.X
        # the algorithm will stop if it converges (no more change in sensitivity parameter values)
        if (sum(apply(dist.mean, 2, stats::sd)) == 0) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X[1], n = max.X[2])))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants, nl, q, n.of.archive, xi)
        # the algorithm will stop if it converges (no more available random samples)
        if (length(o.X) == 0) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X[1], n = max.X[2])))
        }
        #X <- o.X
        X <- NULL
        for (i in 1:n.of.ants){ # exclude unreasonable values
          if (sum((0.01 < o.X[i, 1] & o.X[i, 1] < 0.99),
                  (1 < o.X[i, 2]  && o.X[i, 2] < 10000)) == n.of.opt.pars) {
            X <- rbind(X, o.X[i,])
          }
        }

        if(length(X)>0) {

          p.X <- rbind(p.X, X)
          dim(X) <- c(length(X)/n.of.opt.pars, n.of.opt.pars)

          for (j in 1:dim(X)[1]) { # redo power analysis with n.of.ants times for those reasonable

            n.iter <- n.iter + 1
            p <- X[j, 1]
            n <- X[j, 2]
            if (verbose) {cat('Number of tried evaluations is ', n.iter, ".\n", sep = "")}
            J <- stats::uniroot(function(J) eval(pwr) - power, c(4, 10e4))$root
            m <- (p*n*c1t + (1-p)*n*c1 + c2) * J
            y <- c(y, 1/m)
            pp <- rbind(pp, data.frame(v = 1/m, sd = 0, gr = 0, m = m))
          }
        }

        # recalculate the rank
        pp$gr <- rank(-pp$v, ties.method = "random") # calculate the rank of the solutions
        idx.final <- pp$gr <= n.of.archive
        pp <- pp[idx.final,]
        p.X <- p.X[idx.final,]
        y <- y[idx.final]
        dim(p.X) <- c(length(p.X)/n.of.opt.pars, n.of.opt.pars)
        for (i in 1:n.of.archive) {nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]}

        # check if the required accuracy have been obtained
        if (max(y, na.rm = TRUE) > max.y) {
          max.y <- max(y, na.rm = TRUE)
          max.X <- p.X[which.max(y), ]
          # colnames(max.X) <- c(phan.names, "eval")
          last.impr <- eval}

        if ((abs(max.y - max.value) < abs(e.rel * max.value + e.abs)) |
            (max.y > max.value)) {
          return(list(archive = pp, archive.design.pars = p.X, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X[1], n = max.X[2])))
        }

        # check if the maximum allowed number of objective function
        # evaluations has not been exceeded
        if (n.iter >= max.iter) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X[1], n = max.X[2])))
        }

      }
    } else if (!is.null(par$n) & is.null(par$p)){
      n.of.opt.pars <- 1
      if (verbose) {cat('The ACO algorithm started initilization..', ".\n", sep = "")}
      e.abs <- e # absolute error
      e.rel <- e # relative error
      # initiate parameters
      # eval <- 0
      # iter <- 0
      last.impr <- max.iter
      design.pars <- data.frame()
      outcome <- vector()
      max.X <- rep(NA, n.of.opt.pars)
      max.y <- -Inf
      p.X <- vector()
      pp <- data.frame(v = numeric(), sd = numeric(), gr = numeric());
      outcome <- NULL

      #Power under joint significance test for two-tailed test
      if (two.tailed) {
        pwr <- quote({
          lambda.a <- a/sqrt((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                               (p*(1-p)*n*J));
          lambda.b <- b/sqrt((1-icc)*(1-r12)/
                               (J*n*(1-icc.m)*(1-r12m)));
          (1 - pt(qt(1 - sig.level/tside, df = J*(n-1) -q.a),
                  df = J*(n-1) -q.a, lambda.a) +
             pt(qt(sig.level/tside, df = J*(n-1) -q.a),
                df = J*(n-1) -q.a, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J*(n-1) -q.b),
                    df = J*(n-1) -q.b, lambda.b) +
               pt(qt(sig.level/tside, df = J*(n-1) -q.b),
                  df = J*(n-1) -q.b, lambda.b))
        })
      } else {
        pwr <- quote({
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

      n.of.initial <- round(n.of.archive, 0)
      p.initial <- seq(from = d.p[1], to = d.p[2], length = n.of.initial)
      n.of.archive <- n.of.initial
      nl <- matrix(NA, n.of.archive, n.of.archive-1)
      X <- NULL
      p.X <- NULL
      y <- NULL
      budget <- NULL
      for (p in p.initial){
        X <- rbind(X, p)
        p.X <- rbind(p.X, p)
        J <- stats::uniroot(function(J) eval(pwr) - power, c(4, 10e4))$root
        m <- (p*n*c1t + (1-p)*n*c1 + c2) * J
        y <- c(y, 1/m)
        budget <- c(budget, m)
      }

      pp <- rbind(pp, data.frame(v = y, sd = 0, gr = 0, m = budget))
      pp$gr <- rank(-pp$v, ties.method = "random")

      for (i in 1:n.of.archive){
        nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]
      }
      #  colnames(p.X) <- c("p", "n")
      #  colnames(X) <- c("p", "n")
      #  p.X <- as.data.frame(p.X)
      n.iter <- n.of.archive

      if (verbose)
      {cat('The ACO algorithm finished initilization of ', n.of.archive, ' analyses',".\n", sep = "")}

      while (TRUE) { # the algorithm will stop if one of the stop criteria is met
        dist.mean <- p.X
        # the algorithm will stop if it converges (no more change in sensitivity parameter values)
        if (sum(apply(dist.mean, 2, stats::sd)) == 0) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X, n = par$n)))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants, nl, q, n.of.archive, xi)
        # the algorithm will stop if it converges (no more available random samples)
        if (length(o.X) == 0) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X, n = par$n)))
        }
        #X <- o.X
        X <- NULL
        for (i in 1:n.of.ants){ # exclude unreasonable values
          if (sum((0.001 < o.X[i, 1] & o.X[i, 1] < 0.999)) == n.of.opt.pars) {
            X <- rbind(X, o.X[i,])
          }
        }

        if(length(X)>0) {

          p.X <- rbind(p.X, X)
          dim(X) <- c(length(X)/n.of.opt.pars, n.of.opt.pars)

          for (j in 1:dim(X)[1]) { # redo power analysis with n.of.ants times for those reasonable

            n.iter <- n.iter + 1
            p <- X[j, 1]
            if (verbose) {cat('Number of tried evaluations is ', n.iter, ".\n", sep = "")}
            J <- stats::uniroot(function(J) eval(pwr) - power, c(4, 10e4))$root
            m <- (p*n*c1t + (1-p)*n*c1 + c2) * J
            y <- c(y, 1/m)
            pp <- rbind(pp, data.frame(v = 1/m, sd = 0, gr = 0, m = m))
          }
        }

        # recalculate the rank
        pp$gr <- rank(-pp$v, ties.method = "random") # calculate the rank of the solutions
        idx.final <- pp$gr <= n.of.archive
        pp <- pp[idx.final,]
        p.X <- p.X[idx.final,]
        y <- y[idx.final]
        dim(p.X) <- c(length(p.X)/n.of.opt.pars, n.of.opt.pars)
        for (i in 1:n.of.archive) {nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]}

        # check if the required accuracy have been obtained
        if (max(y, na.rm = TRUE) > max.y) {
          max.y <- max(y, na.rm = TRUE)
          max.X <- p.X[which.max(y), ]
          # colnames(max.X) <- c(phan.names, "eval")
          last.impr <- eval}

        if ((abs(max.y - max.value) < abs(e.rel * max.value + e.abs)) |
            (max.y > max.value)) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X, n = par$n)))
        }

        # check if the maximum allowed number of objective function
        # evaluations has not been exceeded
        if (n.iter >= max.iter) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X, n = par$n)))
        }

      }
    } else if (is.null(par$n) & !is.null(par$p)){
      n.of.opt.pars <- 1
      if (verbose) {cat('The ACO algorithm started initilization..', ".\n", sep = "")}
      e.abs <- e # absolute error
      e.rel <- e # relative error
      # initiate parameters
      # eval <- 0
      # iter <- 0
      last.impr <- max.iter
      design.pars <- data.frame()
      outcome <- vector()
      max.X <- rep(NA, n.of.opt.pars)
      max.y <- -Inf
      p.X <- vector()
      pp <- data.frame(v = numeric(), sd = numeric(), gr = numeric());
      outcome <- NULL

      #Power under joint significance test for two-tailed test
      if (two.tailed) {
        pwr <- quote({
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
                  df = J*(n-1)-q.b, lambda.b))
        })
      } else {
        pwr <- quote({
          lambda.a <- a/sqrt((omega*(1-r22m)*p*(1-p)*n + (1-icc.m)*(1-r12m))/
                               (p*(1-p)*n*J));
          lambda.b <- b/sqrt((1-icc)*(1-r12)/
                               (J*n*(1-icc.m)*(1-r12m)));
          (1 - pt(qt(1 - sig.level/tside, df = J-q.a-1),
                  df = J-q.a-1, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J*(n-1)-q.b),
                    df = J*(n-1)-q.b, lambda.b))
        })
      }

      n.of.initial <- round(n.of.archive, 0)
      n.initial <- seq(from = d.n[1], to = d.n[2], length = n.of.initial)
      n.of.archive <- n.of.initial
      nl <- matrix(NA, n.of.archive, n.of.archive-1)
      X <- NULL
      p.X <- NULL
      y <- NULL
      budget <- NULL

      for (n in n.initial){
        X <- rbind(X, n)
        p.X <- rbind(p.X, n)
        J <- stats::uniroot(function(J) eval(pwr) - power, c(4, 10e4))$root
        m <- (p*n*c1t + (1-p)*n*c1 + c2) * J
        y <- c(y, 1/m)
        budget <- c(budget, m)
      }

      pp <- rbind(pp, data.frame(v = y, sd = 0, gr = 0, m = budget))
      pp$gr <- rank(-pp$v, ties.method = "random")

      for (i in 1:n.of.archive){
        nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]
      }
      n.iter <- n.of.archive

      if (verbose)
      {cat('The ACO algorithm finished initilization of ', n.of.archive, ' analyses',".\n", sep = "")}

      while (TRUE) { # the algorithm will stop if one of the stop criteria is met
        dist.mean <- p.X
        # the algorithm will stop if it converges (no more change in sensitivity parameter values)
        if (sum(apply(dist.mean, 2, stats::sd)) == 0) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = par$p, n = max.X)))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants, nl, q, n.of.archive, xi)
        # the algorithm will stop if it converges (no more available random samples)
        if (length(o.X) == 0) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = par$p, n = max.X)))
        }
        #X <- o.X
        X <- NULL
        for (i in 1:n.of.ants){ # exclude unreasonable values
          if (sum((0.3 < o.X[i, 1]  && o.X[i, 1] < 10000)) == n.of.opt.pars) {

            X <- rbind(X, o.X[i,])
          }
        }

        if(length(X)>0) {

          p.X <- rbind(p.X, X)
          dim(X) <- c(length(X)/n.of.opt.pars, n.of.opt.pars)

          for (j in 1:dim(X)[1]) { # redo power analysis with n.of.ants times for those reasonable

            n.iter <- n.iter + 1
            #p <- X[j, 1]
            n <- X[j, 1]
            if (verbose) {cat('Number of tried evaluations is ', n.iter, ".\n", sep = "")}
            J <- stats::uniroot(function(J) eval(pwr) - power, c(5, 10e4))$root
            m <- (p*n*c1t + (1-p)*n*c1 + c2) * J
            y <- c(y, 1/m)
            pp <- rbind(pp, data.frame(v = 1/m, sd = 0, gr = 0, m = m))
          }
        }

        # recalculate the rank
        pp$gr <- rank(-pp$v, ties.method = "random") # calculate the rank of the solutions
        idx.final <- pp$gr <= n.of.archive
        pp <- pp[idx.final,]
        p.X <- p.X[idx.final,]
        y <- y[idx.final]
        dim(p.X) <- c(length(p.X)/n.of.opt.pars, n.of.opt.pars)
        for (i in 1:n.of.archive) {nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]}

        # check if the required accuracy have been obtained
        if (max(y, na.rm = TRUE) > max.y) {
          max.y <- max(y, na.rm = TRUE)
          max.X <- p.X[which.max(y), ]
          # colnames(max.X) <- c(phan.names, "eval")
          last.impr <- eval}

        if ((abs(max.y - max.value) < abs(e.rel * max.value + e.abs)) |
            (max.y > max.value)) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = par$p, n = max.X)))
        }

        # check if the maximum allowed number of objective function
        # evaluations has not been exceeded
        if (n.iter >= max.iter) {
          return(list(archive = pp, archive.design.pars = p.X, n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = par$p, n = max.X)))
        }

      }
    } else if (!is.null(par$n) & !is.null(par$p)) {
      cat("===============================\n",
          "There is no calculation performed
        because both p and n are contrained",
          ".\n===============================\n", sep = "")
      return(list(par = par, funName = funName,
                  designType = designType,
                  out = c(p = par$p, n = par$n)))
    }
  }
}

