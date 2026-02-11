#' Optimal sample allocation identification for two-level multisite randomized
#'     trials (MRTs) investigating main and moderation effects
#'
#' @description The optimal design of two-level
#'     MRTs probing main and moderation effects
#'     identify the optimal sample allocations.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n})
#'     and the proportion of level-1 individuals/units assigned to
#'     the experimental group (\code{p}).
#'     This function solves the optimal \code{n} and/or \code{p}
#'     with and without a constraint.
#' @inheritParams od.2m
#' @param omega The treatment-by-site variance of the outcome.
#' @param power Statistical power specified for the main effect. The default is .80.
#' @param power.mod Statistical power for the moderation
#'     effect. The default is .80.
#' @param gamma The standardized moderated treatment effect.
#' @param Q The proportion of units in one group for the binary moderator.
#'     Default is 0.5.
#' @param binary Logical; The moderator is binary if TRUE, and continuous if
#'     FALSE. The default is TRUE.
#' @param mod.level The level of the moderator is at. The moderator is at level 1
#'     if mod.level is 1, and at level 2 if mod.level is 2. The default is
#'     mod.level = 1.
#' @param q The number of covariates in the model detecting the main/average
#'     treatment effect. The default is 1.
#' @param q.mod The number of predictors at the moderator level
#'     in the moderation model.
#' @param q.aco The locality of the search (0, 1)
#' @export od.2m.mod
#' @return
#'     Unconstrained or constrained optimal sample allocation (\code{n} and \code{p}).
#'     The function also returns statistical power formulas,
#'     function name, design type,
#'     and parameters used in the calculation.
#' @examples
#' myod <- od.2m.mod(icc = .2, r12 = .5, r22m = .5,
#'                   c1 = 10, c1t = 100, c2 = 50, omega = .01,
#'                   gamma = 0.1)
#' myod$out

od.2m.mod <- function(n = NULL, p = NULL, icc = NULL,
                       r12 = NULL, r22m = NULL,
                       c1 = NULL, c2 = NULL,
                       c1t = NULL, omega = NULL,
                       m = NULL, plots = TRUE, plot.by = list(n = "n", p = "p"),
                       verbose = TRUE, iter = 100,
                       tol = 1e-10, q = 1, q.mod = 1, d = 0.1, gamma = 0.1,
                       power = .80, power.mod = .80, mod.level = 2,
                       d.p = c(0.1, 0.5), d.n = c(2, 1000),
                       sig.level = 0.05, two.tailed = TRUE,
                       Jlim = c(2.5, 1e+10), binary = TRUE,
                       nrange = c(2, 10000),
                       Q = 0.5,
                       max.value = Inf, max.iter = 300,  e = 1e-10,
                       n.of.ants = 10, n.of.archive = 50, q.aco = 0.0001,
                       xi = 0.5) {
  funName <- "od.2m.mod"
  designType <- "two-level MRTs with moderators"
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (sum(sapply(list(icc, r12, r22m,
                      c1, c2, c1t, omega),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc', 'r12', 'r22m',
         'c1', 'c2', 'c1t', and 'omega' must be specified")
  if (sum(sapply(list(icc), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1)
    stop("'icc' must be numeric in [0, 1]")
  if (sum(sapply(list(r12, r22m), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1)
    stop("'r12' and 'r22m' must be numeric in [0, 1]")
  if (sum(sapply(list(c1, c2, c1t), function(x) {
    NumberCheck(x)})) >= 1)
    stop("'c1', 'c2', and 'c1t' must be numeric")
  if (!is.numeric(iter) || iter < 2)
    stop("'iter' must be numeric with iter >= 2")
  par <- list(icc = icc,
              r12 = r12, r22m = r22m,
              c1 = c1, c2 = c2, c1t = c1t, omega = omega,
              m = m, Q = Q,
              n = n, p = p, iter = iter, gamma = gamma, binary = binary,
              power = power, power.mod = power.mod, mod.level = mod.level,
              d = d, q = q, q.mod = q.mod,
              sig.level = sig.level, two.tailed = two.tailed,
              max.iter = max.iter,
              n.of.ants = n.of.ants, n.of.archive = n.of.archive,
              q.aco = q.aco,
              xi = xi
  )

  tside <- ifelse(two.tailed == TRUE, 2, 1)
    if (two.tailed) {#power for main
      pwr.main <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J) /
                             (p * (1 - p) * n * omega * (1 - r22m) +
                                (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = J - q - 1),
               df = J - q - 1, lambda) +
          pt(qt(sig.level / tside, df = J - q - 1),
             df = J - q - 1, lambda)
      })
     } else {
      pwr.main <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J) /
                             (p * (1 - p) * n * omega * (1 - r22m) +
                                (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, J - q - 1),
               df = J - q - 1, lambda)
      })}
  if(binary){var.mod <- Q*(1-Q)} else {var.mod <- 1}

  if(mod.level==1){#moderator at level 1
    if (two.tailed) {
      pwr.mod <- quote({
        lambda <- gamma*sqrt((p * (1 - p) * n * J*var.mod) /
                             ((1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = J*(n-1)-q.mod),
               df = J*(n-1)-q.mod, lambda) +
          pt(qt(sig.level / tside, df = J*(n-1)-q.mod),
             df = J*(n-1)-q.mod, lambda)
      })
    } else {
      pwr.mod <- quote({
        lambda <- gamma*sqrt((p * (1 - p) * n * J*var.mod) /
                               ((1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, J*(n-1)-q.mod),
               df = J*(n-1)-q.mod, lambda)
      })}
  } else if (mod.level==2) {#moderator at level 2
    if (two.tailed) {
      pwr.mod <- quote({
        lambda <- gamma*sqrt((p * (1 - p) * n * J*var.mod) /
                               (omega*(1-r22m)*p * (1 - p) * n +
                                  (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = J-q.mod-1),
               df = J-q.mod-1, lambda) +
          pt(qt(sig.level / tside, df = J-q.mod-1),
             df = J-q.mod-1, lambda)
      })
    } else {
      pwr.mod <- quote({
        lambda <- gamma*sqrt((p * (1 - p) * n * J*var.mod) /
                               (omega*(1-r22m)*p * (1 - p) * n +
                                  (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, J-q.mod-1),
               df = J-q.mod-1, lambda)
      })}
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
          J <- stats::uniroot(function(J)
            eval(pwr.main) - power, Jlim)$root
          m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
          J.mod <- stats::uniroot(function(J)
            eval(pwr.mod) - power.mod, Jlim)$root
          m.mod <- J.mod *((1 - p) * c1 * n + p * c1t * n + c2)
          m <- max(m, m.mod)
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
        if (sum(apply(dist.mean, 2, stats::sd)) < e) {
          m = 1/max.y; p = max.X[1]; n = max.X[2];
          J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- c("p", "n");
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants,
                               nl, q.aco, n.of.archive, xi)
        if (length(o.X) == 0) {
          m = 1/max.y; p = max.X[1]; n = max.X[2];
          J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- c("p", "n");
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
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
            if (verbose) {
              if(n.iter==n.of.archive+1){
                cat('The number of iterations is ', n.iter, sep = "")
              }else if((n.iter>n.of.archive+1)&(n.iter < max.iter)){
                if (n.iter %in% c(seq(100, max.iter, by=100), max.iter)){
                  cat(n.iter,
                      " ", sep = "")
                }else{
                  cat(".", sep = "")
                }
              }
            }
            J <- stats::uniroot(function(J)
              eval(pwr.main) - power, Jlim)$root
            m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
            J.mod <- stats::uniroot(function(J)
              eval(pwr.mod) - power.mod, Jlim)$root
            m.mod <- J.mod *((1 - p) * c1 * n + p * c1t * n + c2)
            m <- max(m, m.mod)
            y <- c(y, 1/m)
            budget <- c(budget, m)
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
          m = 1/max.y; p = max.X[1]; n = max.X[2];
          J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- c("p", "n");
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
        }

        # check if the maximum allowed number of objective function
        # evaluations has not been exceeded
        if (n.iter >= max.iter) {
          m = 1/max.y; p = max.X[1]; n = max.X[2];
          J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- c("p", "n");
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
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
        J <- stats::uniroot(function(J)
          eval(pwr.main) - power, Jlim)$root
        m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
        J.mod <- stats::uniroot(function(J)
          eval(pwr.mod) - power.mod, Jlim)$root
        m.mod <- J.mod *((1 - p) * c1 * n + p * c1t * n + c2)
        m <- max(m, m.mod)
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
        if (sum(apply(dist.mean, 2, stats::sd)) < e) {
          m = 1/max.y; p = max.X;
          n = par$n; J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- "p";
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank,
                               n.of.ants, nl, q.aco, n.of.archive, xi)

        if (length(o.X) == 0) {
          m = 1/max.y; p = max.X;
          n = par$n; J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- "p";
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
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
            if (verbose) {
              if(n.iter==n.of.archive+1){
                cat('The number of iterations is ', n.iter, sep = "")
              }else if((n.iter>n.of.archive+1)&(n.iter < max.iter)){
                if (n.iter %in% c(seq(100, max.iter, by=100), max.iter)){
                  cat(n.iter,
                      " ", sep = "")
                }else{
                  cat(".", sep = "")
                }
              }
            }
            J <- stats::uniroot(function(J)
              eval(pwr.main) - power, Jlim)$root
            m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
            J.mod <- stats::uniroot(function(J)
              eval(pwr.mod) - power.mod, Jlim)$root
            m.mod <- J.mod *((1 - p) * c1 * n + p * c1t * n + c2)
            m <- max(m, m.mod)
            y <- c(y, 1/m)
            budget <- c(budget, m)
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
          m = 1/max.y; p = max.X;
          n = par$n; J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- "p";
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
        }

        # check if the maximum allowed number of objective function
        # evaluations has not been exceeded
        if (n.iter >= max.iter) {
          m = 1/max.y; p = max.X;
          n = par$n; J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- "p";
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
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
        J <- stats::uniroot(function(J)
          eval(pwr.main) - power, Jlim)$root
        m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
        J.mod <- stats::uniroot(function(J)
          eval(pwr.mod) - power.mod, Jlim)$root
        m.mod <- J.mod *((1 - p) * c1 * n + p * c1t * n + c2)
        m <- max(m, m.mod)
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
        if (sum(apply(dist.mean, 2, stats::sd)) < e) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = par$p,
                                 n = max.X,
                                 J = 1/max.y/((1 - p) * c1 * n + p * c1t * n + c2))))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank,
                               n.of.ants, nl, q.aco, n.of.archive, xi)
        if (length(o.X) == 0) {
          m = 1/max.y; p = par$p; n = max.X;
          J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- "n";
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
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
            if (verbose) {
              if(n.iter==n.of.archive+1){
                cat('The number of iterations is ', n.iter, sep = "")
              }else if((n.iter>n.of.archive+1)&(n.iter < max.iter)){
                if (n.iter %in% c(seq(100, max.iter, by=100), max.iter)){
                  cat(n.iter,
                      " ", sep = "")
                }else{
                  cat(".", sep = "")
                }
              }
            }
            J <- stats::uniroot(function(J)
              eval(pwr.main) - power, Jlim)$root
            m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
            J.mod <- stats::uniroot(function(J)
              eval(pwr.mod) - power.mod, Jlim)$root
            m.mod <- J.mod *((1 - p) * c1 * n + p * c1t * n + c2)
            m <- max(m, m.mod)
            y <- c(y, 1/m)
            budget <- c(budget, m)
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
          m = 1/max.y; p = par$p; n = max.X;
          J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- "n";
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
        }

        # check if the maximum allowed number of objective function
        # evaluations has not been exceeded
        if (n.iter >= max.iter) {
          m = 1/max.y; p = par$p; n = max.X;
          J = m/((1 - p) * c1 * n + p * c1t * n + c2);
          colnames(p.X) <- "n";
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = m, p = p, n = n,
                                 J = J)))
        }
      }
    } else if (!is.null(par$n) & !is.null(par$p)) {
      cat("===============================\n",
          "There is no optimization performed
        because both p and n are contrained",
          ".\n===============================\n", sep = "")
      J <- stats::uniroot(function(J)
        eval(pwr.main) - power, Jlim)$root
      m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
      J.mod <- stats::uniroot(function(J)
        eval(pwr.mod) - power.mod, Jlim)$root
      m.mod <- J.mod *((1 - p) * c1 * n + p * c1t * n + c2)
      m <- max(m, m.mod)
      J <- max(J, J.mod)
      return(list(par = par, funName = funName,
                  designType = designType,
                  out = list(m = m,
                             p = par$p,
                             n = par$n, J = J)))
    }
}

