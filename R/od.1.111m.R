#' Jointly optimal sample allocation identification for single-level
#'     randomized controlled trials (RCTs) investigating
#'     main and moderation effects (1-1-1m)
#'
#' @description The optimal design of single-level RCTs
#'     probing main and moderation effects is to identify the jointly optimal sample
#'     allocation that use the minimum budget to achieve targeted
#'     statistical power for both effects. The optimal design parameter
#'     is the proportion of
#'     individuals/units assigned to the experimental condition.
#'     This function uses the ant colony optimization algorithm
#'     to identify the optimal \code{p}.
#'
#' @inheritParams power.1.111m
#' @param power.main Statistical power specified for the main effect.
#'     The default value is .80.
#' @param power.mod Statistical power specified for the moderation effect.
#'     The default value is .80.
#' @param d.p The initial sampling domain for p. Default is c(0.1, 0.5).
#' @param e Maximum error value used when solution quality used as
#'     the stopping criterion. The default value is 1e-10.
#' @param max.value Maximal value of optimization when used as
#'     the stopping criterion. Default is infinite.
#' @param d The standardized main effect size.
#' @param r12 The proportion of within-treatment outcome variance explained
#'     by covariates in the model that estimated the main effect.
#' @param q.main The number of covariates in the model estimating the main
#'     effect (besides the treatment, moderator). The default value is 1.
#' @param d.p The initial sampling domains for p. Default is c(0.10, 0.50).
#' @param max.iter Maximal number of function evaluations when used as
#'     the stopping criterion. The default value is 300.
#' @param n.of.archive  Size of the solution archive, default is 20.
#' @param q Locality of the search (0,1). The default value is 0.0001.
#' @param xi  Convergence pressure (0, Inf), suggested: (0, 1).
#'     The default value is 0.5.
#' @param verbose Print out evaluation process if TRUE. The default value is TRUE.
#' @param n.of.ants Number of ants used in each iteration after
#'    the initialization stage. The default value is 10.
#'
#' @return
#'     Unconstrained or constrained optimal sample allocation \code{p}).
#'     The function also returns statistical power for
#'     main and moderation effects,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.1.111m
#' @examples
#' myod <- od.1.111m(d =.1, gamma = .2, r12 = .50,
#'                  c1 = 10, c1t = 100)
#' myod


od.1.111m <- function(d = NULL, gamma = NULL, n = NULL, Q = .50,
                      p = NULL, binary = TRUE,
                      c1 = NULL, c1t = NULL,
                      r12 = NULL, r.yx = 0, r.mx = 0, r.ym = 0,
                      m = NULL,
                      q.main = 1, q.mod = 1,
                      power.mod = 0.80, power.main = 0.80,
                      d.p = c(0.1, 0.5),
                      sig.level = 0.05, two.tailed = TRUE,
                      verbose = TRUE, nlim = c(20, 1e7),
                      max.value = Inf, max.iter = 300,  e = 1e-10,
                      n.of.ants = 10, n.of.archive = 50, q = 0.0001,
                      xi = 0.5) {
  funName <- "od.1.111m"
  designType <- "1-1-1 moderation in single-level RCTs"
  par <- list(d = d, gamma = gamma, n = n, p = p, Q = Q,
              r12 = r12, r.yx = r.yx, r.mx = r.mx, r.ym = r.ym,
              c1 = c1, c1t = c1t,
              m = m, q.mod = q.mod, q.main = q.main,
              sig.level = sig.level, two.tailed = two.tailed,
              binary = binary,
              max.iter = max.iter,
              n.of.ants = n.of.ants, n.of.archive = n.of.archive,
              q = q,
              xi = xi)
  if (sum(sapply(list(d, gamma, c1, c1t),
                 function(x) is.null(x))) >= 1)
    stop("All of 'd', 'gamma', 'c1',
         'c1t' must be specified")
  NumberCheck <- function(x) {!is.null(x) & !is.numeric(x)}
  if (sum(sapply(list(c1, c1t), function(x) {
    NumberCheck(x) | x < 0})) >= 1)
    stop("'c1', 'c1t' must be numeric in [0, inf)")
  if (c1 == 0 & c1t == 0 & is.null(n) & is.null(p))
    stop("when c1 and c1t are both zero, p must be constrained,
         please specify a value for p")
  labFun <- function(x, y) {
    if (!is.null(x) & length(x) == 1 & is.character(x)) {x} else {y}
  }
  plotbyFun <- function(x, y) {
    if (!is.null(x) & is.list(x)) {x} else {y}
  }
  tside <- ifelse(two.tailed == TRUE, 2, 1)

  if (binary){
    if (two.tailed) {
      pwr.mod <- quote({
        1 - pt(qt(1-sig.level/2, df = n-q.mod-4),
               df = n-q.mod-4, gamma/sqrt((1-r.yx^2)/(n*(p*(1-p)*Q*(1-Q))))) +
          pt(qt(sig.level/2, df = n-q.mod-4),
             df = n-q.mod-4, gamma/sqrt((1-r.yx^2)/(n*(p*(1-p)*Q*(1-Q)))))
      })
    } else {
      pwr.mod <- quote({
        1 - pt(qt(1-sig.level, df = n-q.mod-4), df = n-q.mod-4,
               gamma/sqrt((1-r.yx^2)/(n*(p*(1-p)*Q*(1-Q)))))
      })
    }
  } else {
    if (two.tailed) {
      pwr.mod <- quote({
        lambda <-gamma/sqrt((1 - r.yx^2 - r.ym^2 - r.mx^2 +
                               2*r.yx*r.ym*r.mx)/(n*(p*(1-p))*(1 - r.mx^2)));
        1 - pt(qt(1 - sig.level/2, df = n-q.mod-4),
               df = n-q.mod-4, lambda) +
          pt(qt(sig.level/2, df = n-q.mod-4),
             df = n-q.mod-4, lambda)
      })
    }else{
      pwr.mod <- quote({
        lambda <-gamma/sqrt((1 - r.yx^2 - r.ym^2 - r.mx^2 +
                               2*r.yx*r.ym*r.mx)/(n*(p*(1-p))*(1 - r.mx^2)));
        1 - pt(qt(1 - sig.level, df = n-q.mod-4),
               df = n-q.mod-4, lambda)
      })
    }
  }

  if (two.tailed) {
    pwr.main <- quote({
      lambda <- d * sqrt(p * (1 - p) * n) /
        sqrt(1 - r12);
      1 - pt(qt(1 - sig.level / tside, df = n - q.main - 2),
             df = n - q.main - 2, lambda) +
        pt(qt(sig.level / tside, df = n - q.main - 2),
           df = n - q.main - 2, lambda)
    })
  } else {
    pwr.main <- quote({
      lambda <- d * sqrt(p * (1 - p) * n) /
        sqrt(1 - r12);
      1 - pt(qt(1 - sig.level / tside, df = n - q.main - 2),
             df = n - q.main - 2, lambda)
    })
  }

  par <- c(par, pwr.main = pwr.main, pwr.mod = pwr.mod)
  if(!is.null(par$p)){d.p[1] = par$p; prange[1] = par$p}

 if (is.null(par$p)){
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
      n.mod <- stats::uniroot(function(n) eval(pwr.mod) -
                                power.mod, nlim)$root
      n.main <- stats::uniroot(function(n) eval(pwr.main) -
                                 power.main, nlim)$root
      n <- max(n.mod, n.main)
      m <- p * c1t * n + (1 - p) * c1*n
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
                    out = list(m = 1/max.y, p = max.X)))
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
                    out = list(m = 1/max.y, p = max.X)))
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
          n.mod <- stats::uniroot(function(n) eval(pwr.mod) -
                                    power.mod, nlim)$root
          n.main <- stats::uniroot(function(n) eval(pwr.main) -
                                     power.main, nlim)$root
          n <- max(n.mod, n.main)
          m <- p * c1t * n + (1 - p) * c1*n
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
                    out = list(m = 1/max.y, p = max.X)))
      }

      # check if the maximum allowed number of objective function
      # evaluations has not been exceeded
      if (n.iter >= max.iter) {
        return(list(archive = pp, archive.design.pars = p.X,
                    n.iter = n.iter, par = par, funName = funName,
                    designType = designType,
                    out = list(m = 1/max.y, p = max.X)))
      }

    }
  } else if (!is.null(par$p)) {
    cat("===============================\n",
        "There is no calculation performed
        because p is contrained",
        ".\n===============================\n", sep = "")
    return(list(par = par, funName = funName,
                designType = designType,
                out = c(p = par$p)))
  }

}
