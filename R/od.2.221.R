#' Optimal sample allocation calculation for two-level CRTs probing
#'     mediation effects with cluster-level mediators
#'
#' @description The optimal design of two-level
#'     cluster randomized trials (CRTs) probing mediation effects with
#'     cluster-level mediators is to to identify
#'     the optimal sample allocation that requires the minimum budget
#'     to achieve certain power level.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n})
#'     and the proportion of level-2 clusters/groups to be assigned to
#'     treatment (\code{p}).
#'     This function solves the optimal \code{n} and/or \code{p}
#'     with and without a constraint.
#'
#' @inheritParams power.2.221
#' @inheritParams od.2
#' @inheritParams od.1.111
#' @param a The treatment effect on the mediator.
#' @param b The within treatment correlation between the outcome and
#'     the mediator at the cluster level.
#' @param q.a The number of covariates in the mediator model
#'     (except the treatment indicator).
#' @param q.b The  number of covariates in the outcome model
#'     at the cluster level
#'     (except the treatment indicator and the mediator).
#' @param r.yx The correlation between the outcome and the covariate at the
#'     individual level.
#' @param r.yw The correlation between the outcome and the covariate at the
#'     cluster level.
#' @param r.mw The correlation between the mediator and the covariate at the
#'     cluster level.
#' @param r2m The proportion of mediator variance explained by covariates in
#'     the mediator model.
#' @param test The type of test will be used to detect mediation effects.
#'     Default is the joint significance test (i.e., test = "joint").
#'     The other choice is the Sobel test
#'     by specifying the argument as test = "sobel".
#' @param two.tailed Two tailed test, the default value is TRUE.
#' @param sig.level Significance level or type I error rate, default value is 0.05.
#' @param power Statistical power specified. The default value is .80.
#' @param e Maximum error value used when solution quality used as
#'     the stopping criterion, default is 1e-10.
#' @param max.value Maximal value of optimization when used as
#'     the stopping criterion. Default is -Inf.
#' @param d.p The initial sampling domains for p. Default is c(0.1, 0.5).
#' @param d.n The initial sampling domain for n. Default is c(2, 100).
#' @param max.iter Maximal number of function evaluations when used as
#'     the stopping criterion.
#' @param n.of.archive  Size of the solution archive, default is 100.
#' @param q Locality of the search (0,1), default is 0.0001.
#' @param xi  Convergence pressure (0, Inf), suggested: (0, 1), default is 0.5.
#' @param verbose Print out evaluation process if TRUE, default is TRUE.
#' @param Jlim The range for J to solve for a numerical solution.
#'     Default is c(max(q.a, q.b)+4, 1e6).
#' @param nrange The range of the individual-level sample size per cluster
#'     that used to exclude unreasonable values. Default value is c(1.5, 10000).
#' @param n.of.ants Number of ants used in each iteration after the initialization
#'     of power analysis for calculating required budget, default value is 10.
#' @param tol convergence tolerance.
#'
#'
#' @return
#'     Unconstrained or constrained optimal sample allocation
#'     (\code{n} and \code{p}).
#'     The function also returns the variance of a mediation effect
#'     or statistical power,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.2.221
#'


od.2.221 <- function(a = NULL, b = NULL, n = NULL,
                     p = NULL, icc = NULL,
                     c1 = NULL, c1t = NULL, c2 = NULL, c2t = NULL,
                     m = NULL,
                     r2m = 0, r.yx = 0, r.mw = 0, r.yw = 0,
                     q.a = 0, q.b = 0,
                     test = "joint",
                     tol = 1e-11, power = 0.80,
                     d.p = c(0.1, 0.5), d.n = c(2, 100),
                     sig.level = 0.05, two.tailed = TRUE,
                     Jlim = c(max(q.a, q.b)+4, 1e6),
                     verbose = TRUE, nrange = c(1.5, 10000),
                     max.value = Inf, max.iter = 300,  e = 1e-10,
                     n.of.ants = 10, n.of.archive = 50, q = 0.0001,
                     xi = 0.5
                 ) {
  funName <- "od.2.221"
  designType <- "2-2-1 mediation in 2-level CRTs"
  par <- list(a = a, b = b, n = n, p = p, icc = icc,
              r2m = r2m, r.yx = r.yx, r.yw = r.yw, r.mw = r.mw,
              c1 = c1, c2 = c2, c1t =c1t, c2t = c2t,
              m = m, q.a = q.a, q.b = q.b,
              sig.level = sig.level, two.tailed = two.tailed,
              test = test,
              max.iter = max.iter,
              n.of.ants = n.of.ants, n.of.archive = n.of.archive,
              q = q,
              xi = xi
              )
  if (sum(sapply(list(a, b, icc, r.yx, r.yw, r.mw, c1, c2, c1t, c2t),
                 function(x) is.null(x))) >= 1)
    stop("All of 'a', 'b', 'icc', 'r.yx', 'r.yw', 'r.mw', 'c1', 'c2',
         'c1t', 'c2t' must be specified")
  NumberCheck <- function(x) {!is.null(x) & !is.numeric(x)}
  if (NumberCheck(icc) | any(0 > icc | icc > 1))
    stop("'icc' must be numeric in [0, 1]")
  if (sum(sapply(list(r.yx, r.yw, r.mw), function(x) {
    NumberCheck(x) | any(-1 > x | x > 1)
  })) >= 1)
    stop("'r.yx', 'r.yw', 'r.mw' must be numeric in [-1, 1]")
  if (sum(sapply(list(c1, c2, c1t, c2t), function(x) {
    NumberCheck(x) | x < 0})) >= 1)
    stop("'c1', 'c2', 'c1t', 'c2t' must be numeric in [0, inf)")
  if (c1 == 0 & c1t == 0 & is.null(n) & is.null(p))
    stop("when c1 and c1t are both zero, one of n or p must be constrained,
         please specify a value for n or p")
  if (c2 == 0 & c2t == 0 & is.null(n) & is.null(p))
    stop("when c2 and c2t are both zero, one of n or p must be constrained,
         please specify a value for n or p")
  labFun <- function(x, y) {
    if (!is.null(x) & length(x) == 1 & is.character(x)) {x} else {y}
  }
  plotbyFun <- function(x, y) {
    if (!is.null(x) & is.list(x)) {x} else {y}
  }
  tside <- ifelse(two.tailed == TRUE, 2, 1)
  if (test == "joint" | test == "Joint" | test == "JOINT") {
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

    if (is.null(par$p) & is.null(par$n)) {
      n.of.opt.pars <- 2
      if (verbose) {cat('The ACO algorithm started initilization..',
                        ".\n", sep = "")}
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
          J <- stats::uniroot(function(J) eval(pwr) - power,
                              c(max(q.a, q.b) + 4, 1e6))$root
          m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
          y <- c(y, 1/m)
          budget <- c(budget, m)

        }
      }
      pp <- rbind(pp, data.frame(v = y, sd = 0, gr = 0, m = budget))
      pp$gr <- rank(-pp$v, ties.method = "random")

      for (i in 1:n.of.archive){
        nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]
      }
      n.iter <- n.of.archive

      if (verbose)
      {cat('The ACO algorithm finished initilization of ', n.of.archive,
           ' analyses',".\n", sep = "")}

      while (TRUE) { # the algorithm will stop if one of the criteria is met
        dist.mean <- p.X
        # the algorithm will stop if it converges
        if (sum(apply(dist.mean, 2, stats::sd)) == 0) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X[1], n = max.X[2])))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants,
                               nl, q, n.of.archive, xi)

        if (length(o.X) == 0) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X[1], n = max.X[2])))
        }
        X <- NULL
        for (i in 1:n.of.ants){ # exclude unreasonable values
          if (sum((0.001 < o.X[i, 1] & o.X[i, 1] < 0.999),
                  (nrange[1] < o.X[i, 2]  && o.X[i, 2] < nrange[2])) == n.of.opt.pars) {
            X <- rbind(X, o.X[i,])
          }
        }
        if(length(X)>0) {
          p.X <- rbind(p.X, X)
          dim(X) <- c(length(X)/n.of.opt.pars, n.of.opt.pars)

          for (j in 1:dim(X)[1]) {
            # redo power analysis with n.of.ants times for those reasonable
            n.iter <- n.iter + 1
            p <- X[j, 1]
            n <- X[j, 2]
            if (verbose) {cat('Number of tried evaluations is ', n.iter,
                              ".\n", sep = "")}
            J <- stats::uniroot(function(J) eval(pwr) - power, Jlim)$root
            m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
            y <- c(y, 1/m)
            pp <- rbind(pp, data.frame(v = 1/m, sd = 0, gr = 0, m = m))
          }
        }

        # recalculate the rank
        pp$gr <- rank(-pp$v, ties.method = "random")
        idx.final <- pp$gr <= n.of.archive
        pp <- pp[idx.final,]
        p.X <- p.X[idx.final,]
        y <- y[idx.final]
        dim(p.X) <- c(length(p.X)/n.of.opt.pars, n.of.opt.pars)
        for (i in 1:n.of.archive)
          {nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]}

        # check if the required accuracy have been obtained
        if (max(y, na.rm = TRUE) > max.y) {
          max.y <- max(y, na.rm = TRUE)
          max.X <- p.X[which.max(y), ]
          last.impr <- eval}

        if ((abs(max.y - max.value) < abs(e.rel * max.value + e.abs)) |
            (max.y > max.value)) {
          return(list(archive = pp, archive.design.pars = p.X,
                      archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X[1], n = max.X[2])))
        }

        # check if the maximum allowed number of objective function
        # evaluations has not been exceeded
        if (n.iter >= max.iter) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X[1], n = max.X[2])))
        }

      }
    } else if (!is.null(par$n) & is.null(par$p)){
      n.of.opt.pars <- 1
      if (verbose) {cat('The ACO algorithm started initilization..',
                        ".\n", sep = "")}
      e.abs <- e # absolute error
      e.rel <- e # relative error
      last.impr <- max.iter
      design.pars <- data.frame()
      outcome <- vector()
      max.X <- rep(NA, n.of.opt.pars)
      max.y <- -Inf
      p.X <- vector()
      pp <- data.frame(v = numeric(), sd = numeric(), gr = numeric());
      outcome <- NULL

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
        J <- stats::uniroot(function(J) eval(pwr) - power, Jlim)$root
        m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
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
      {cat('The ACO algorithm finished initilization of ', n.of.archive,
           ' analyses',".\n", sep = "")}

      while (TRUE) { # the algorithm will stop if one of the criteria is met
        dist.mean <- p.X
        # the algorithm will stop if it converges
        if (sum(apply(dist.mean, 2, stats::sd)) == 0) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X, n = par$n)))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank,
                               n.of.ants, nl, q, n.of.archive, xi)

        if (length(o.X) == 0) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X, n = par$n)))
        }
        X <- NULL
        for (i in 1:n.of.ants){ # exclude unreasonable values
          if (sum((0.001 < o.X[i, 1] & o.X[i, 1] < 0.999)) == n.of.opt.pars) {
            X <- rbind(X, o.X[i,])
          }
        }

        if(length(X)>0) {
          p.X <- rbind(p.X, X)
          dim(X) <- c(length(X)/n.of.opt.pars, n.of.opt.pars)

          for (j in 1:dim(X)[1]) {
            # redo power analysis with n.of.ants times for those reasonable

            n.iter <- n.iter + 1
            p <- X[j, 1]
            if (verbose) {cat('Number of tried evaluations is ', n.iter,
                              ".\n", sep = "")}
            J <- stats::uniroot(function(J) eval(pwr) - power, Jlim)$root
            m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
            y <- c(y, 1/m)
            pp <- rbind(pp, data.frame(v = 1/m, sd = 0, gr = 0, m = m))
          }
        }

        # recalculate the rank
        pp$gr <- rank(-pp$v, ties.method = "random")
        idx.final <- pp$gr <= n.of.archive
        pp <- pp[idx.final,]
        p.X <- p.X[idx.final,]
        y <- y[idx.final]
        dim(p.X) <- c(length(p.X)/n.of.opt.pars, n.of.opt.pars)
        for (i in 1:n.of.archive)
          {nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]}

        # check if the required accuracy have been obtained
        if (max(y, na.rm = TRUE) > max.y) {
          max.y <- max(y, na.rm = TRUE)
          max.X <- p.X[which.max(y), ]
          last.impr <- eval}

        if ((abs(max.y - max.value) < abs(e.rel * max.value + e.abs)) |
            (max.y > max.value)) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X, n = par$n)))
        }

        # check if the maximum allowed number of objective function
        # evaluations has not been exceeded
        if (n.iter >= max.iter) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X, n = par$n)))
        }

      }
    } else if (is.null(par$n) & !is.null(par$p)){
      n.of.opt.pars <- 1
      if (verbose) {cat('The ACO algorithm started initilization..',
                        ".\n", sep = "")}
      e.abs <- e # absolute error
      e.rel <- e # relative error
      last.impr <- max.iter
      design.pars <- data.frame()
      outcome <- vector()
      max.X <- rep(NA, n.of.opt.pars)
      max.y <- -Inf
      p.X <- vector()
      pp <- data.frame(v = numeric(), sd = numeric(), gr = numeric());
      outcome <- NULL

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
        J <- stats::uniroot(function(J) eval(pwr) - power, Jlim)$root
        m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
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
      {cat('The ACO algorithm finished initilization of ',
           n.of.archive, ' analyses',".\n", sep = "")}

      while (TRUE) { # the algorithm will stop if one of the criteria is met
        dist.mean <- p.X
        if (sum(apply(dist.mean, 2, stats::sd)) == 0) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = par$p, n = max.X)))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank,
                               n.of.ants, nl, q, n.of.archive, xi)
        if (length(o.X) == 0) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = par$p, n = max.X)))
        }
        X <- NULL
        for (i in 1:n.of.ants){ # exclude unreasonable values
          if (sum((nrange[1] < o.X[i, 1]  && o.X[i, 1] < nrange[2])) == n.of.opt.pars) {
            X <- rbind(X, o.X[i,])
          }
        }

        if(length(X)>0) {

          p.X <- rbind(p.X, X)
          dim(X) <- c(length(X)/n.of.opt.pars, n.of.opt.pars)

          for (j in 1:dim(X)[1]) {
            # redo power analysis with n.of.ants times for those reasonable

            n.iter <- n.iter + 1
            n <- X[j, 1]
            if (verbose) {cat('Number of tried evaluations is ',
                              n.iter, ".\n", sep = "")}
            J <- stats::uniroot(function(J) eval(pwr) - power, Jlim)$root
            m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
            y <- c(y, 1/m)
            pp <- rbind(pp, data.frame(v = 1/m, sd = 0, gr = 0, m = m))
          }
        }

        # recalculate the rank
        pp$gr <- rank(-pp$v, ties.method = "random")
        idx.final <- pp$gr <= n.of.archive
        pp <- pp[idx.final,]
        p.X <- p.X[idx.final,]
        y <- y[idx.final]
        dim(p.X) <- c(length(p.X)/n.of.opt.pars, n.of.opt.pars)
        for (i in 1:n.of.archive)
          {nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]}

        # check if the required accuracy have been obtained
        if (max(y, na.rm = TRUE) > max.y) {
          max.y <- max(y, na.rm = TRUE)
          max.X <- p.X[which.max(y), ]
          last.impr <- eval}

        if ((abs(max.y - max.value) < abs(e.rel * max.value + e.abs)) |
            (max.y > max.value)) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = par$p, n = max.X)))
        }

        # check if the maximum allowed number of objective function
        # evaluations has not been exceeded
        if (n.iter >= max.iter) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
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

