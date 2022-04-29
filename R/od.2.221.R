#' Optimal sample allocation calculation for two-level CRTs probing
#'     mediation effects with cluster-level mediators
#'
#' @description The optimal design of two-level
#'     cluster randomized trials (CRTs) probing mediation effects with
#'     cluster-level mediators, for the
#'     Sobel test, is to calculate
#'     the optimal sample allocation that minimizes the variance of
#'     a mediation effect under a fixed budget. For the joint significance test, it is to identify
#'     the optimal sample allocation that requires the minimum budget
#'     to achieve certain power level.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n})
#'     and the proportion of level-2 clusters/groups to be assigned to treatment (\code{p}).
#'     This function solves the optimal \code{n} and/or \code{p}
#'     with and without a constraint.
#'
#' @inheritParams power.2.221
#' @inheritParams od.2
#' @param a The treatment effect on the mediator.
#' @param b The within treatment correlation between the outcome and the mediator.
#' @param q.b The number of covaraites in the mediator model (except the treatment indicator).
#' @param q.a The  number of covaraites in the outcome model (except the treatment indicator and the mediator).
#' @param test The type of test will be used to detect mediation effects. Default is
#'     the joint significance test (i.e., test = "joint").
#'     Other choices are the Sobel test and Monte Carlo confidence interval
#'     test by specifying the argument as test = "sobel" or test = "mcci".
#' @param power.joint Statistical power specified for the joint significance test,
#'     default is .80.
#' @param r2m The proportion of within treatment mediator variance explained by covariates.
#' @param r12 The proportion of within treatment individual-level outcome variance
#'     explained by covariates.
#' @param r22 The proportion of within treatment group-level outcome variance
#'     explained by covariates and the mediator.
#' @param e Maximum error value used when solution quality used as
#'     the stopping criterion, default is 1e-10.
#' @param max.value Maximal value of optimization when used as
#'     the stopping criterion. Default is infinite.
#' @param d.p The initial sampling domains for p. Default is c(0.03, 0.97).
#' @param d.n The initial sampling domain for n. Default is c(0.5, 500).
#' @param max.iter Maximal number of function evaluations when used as
#'     the stopping criterion.
#' @param n.of.archive  Size of the solution archive, default is 100.
#' @param n.rep Number of replications in MCCI power calculation.
#' @param q Locality of the search (0,1), default is 0.0001.
#' @param xi  Convergence pressure (0, Inf), suggested: (0, 1), default is 0.5.
#' @param verbose Print out evaluation process if TRUE, default is TRUE.
#' @param Jlim The range for J to search for a numerical solution. Default is c(4, 10e4).
#' @param n.of.ants Number of ants used in each iteration after the initialization
#'     of power analysis for calculating required budget, default value is 10.
#' @param iter number of iteration used for solving roots in Sobel test.
#' @param tol convergence tolerance.
#'
#'
#' @return
#'     Unconstrained or constrained optimal sample allocation (\code{n} and \code{p}).
#'     The function also returns the variance of a mediation effect or statistical power,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.2.221
#'


od.2.221 <- function(a = NULL, b = NULL, icc = NULL,
                 c1 = NULL, c1t = NULL, c2 = NULL, c2t = NULL, m = NULL,
                 r2m = 0, r12 = 0, r22 = 0,
                 q.a = 0, q.b = 0,
                 test = "joint",
                 n = NULL, p = NULL, iter = 100,
                 tol = 1e-11, power.joint = 0.80,
                 d.p = c(0.1, 0.5), d.n = c(1, 50),
                 sig.level = 0.05, two.tailed = TRUE,
                 plots = TRUE, plot.by = NULL,
                 n.rep = 1000,
                 nlim = c(0, 50), plim = c(0, 1), varlim = c(0, 0.05),
                 Jlim = c(4, 10e4),
                 nlab = NULL, plab = NULL, varlab = NULL,
                 vartitle = NULL,verbose = TRUE,
                 max.value = Inf, max.iter = 150,  e = 1e-10,
                 n.of.ants = 10, n.of.archive = 50, q = 0.0001,
                 xi = 0.5
                 ) {
  funName <- "od.2.221"
  designType <- "2-2-1 mediation in 2-level CRTs"
  par <- list(a = a, b = b, icc = icc,
              r12 = r12, r22 = r22, r2m = r2m,
              c1 = c1, c2 = c2, c1t =c1t, c2t = c2t,
              n = n, p = p, m = m,
              q.a = q.a, q.b = q.b,
              sig.level = sig.level, two.tailed = two.tailed,
              test = test,
              max.iter = max.iter,
              n.of.ants = n.of.ants, n.of.archive = n.of.archive,
              q = q,
              xi = xi
              )
  if (sum(sapply(list(icc, r2m, r12, r22, c1, c2, c1t, c2t),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc', 'r2m', 'r12', 'r22', 'c1', 'c2',
         'c1t', 'c2t' must be specified")
  NumberCheck <- function(x) {!is.null(x) & !is.numeric(x)}
  if (NumberCheck(icc) | any(0 > icc | icc > 1))
    stop("'icc' must be numeric in [0, 1]")
  if (sum(sapply(list(r2m, r12, r22), function(x) {
    NumberCheck(x) | any(0 > x | x > 1)
  })) >= 1)
    stop("'r2m', 'r12', 'r22' must be numeric in [0, 1]")
  if (sum(sapply(list(c1, c2, c1t, c2t), function(x) {
    NumberCheck(x) | x < 0})) >= 1)
    stop("'c1', 'c2', 'c1t', 'c2t' must be numeric in [0, inf)")
  if (!is.null(plot.by) & !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n'))")
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

  if (test == "sobel" | test == "Sobel" | test == "SOBEL") {
  out <- NULL
  if (is.null(par$p)) {
  p.expr <- quote({
    a^2 * (icc * (1 - r22) + (1 - icc) *
             (1 - r12)/n) *
      ((c1t * n + c2t) - (c1 * n + c2)) *
      p^2 * (1 - p)^2 +
      b^2 * (1-r2m)^2 *((c1t * n + c2t) - (c1 * n + c2)) *
      p * (1 - p) - (1 - 2*p) * b^2 * (1 - r2m)^2 *(
        (1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))
  })}

  if (is.null(par$n)) {
    n.expr <- quote({
      n - sqrt((a^2 * (1 - icc) * (1 - r12) * ((1 - p) * (c1 * n + c2) +
          p * (c1t *n + c2t)) * p * (1 - p))/
            (a^2 * p * (1 - p) * (icc * (1 - r22) + (1 - icc) * (1 - r12)/n) *
                ((1 - p) * c1 + p * c1t) +
                b^2 * (1 - r2m)^2 *((1 - p) * c1 + p * c1t)))
            })
  }

  if (!is.null(par$p) & !is.null(par$n)) {
    cat("===============================\n",
        "There is no calculation performed
        because both p and n are contrained",
        ".\n===============================\n", sep = "")
  } else if (is.null(par$p) & !is.null(par$n)) {
    out$p <- p <- stats::uniroot(function(p)
      eval(p.expr), c(0.05, 0.95))$root
  } else if (!is.null(par$p) & is.null(par$n)) {
    out$n <- n <- stats::uniroot(function(n)
      eval(n.expr), c(0.3, 10000))$root
  } else if (is.null(par$p) & is.null(par$n)) {
   nn <- pp <- NULL; n <- sample(1:100, 1)
   for (i in 1:iter) {
     out$p <- p <- stats::uniroot(function(p)
       eval(p.expr), c(0.05, 0.95))$root; pp[i] <- p
      out$n <- n <- stats::uniroot(function(n)
         eval(n.expr), c(0.3, 10000))$root; nn[i] <- n
   }
     if (abs(nn[iter] - nn[iter-1]) <= tol & abs(pp[iter] - pp[iter-1]) <= tol) {
       nn <- pp <- NULL
     } else {
       cat("===============================\n",
           "The solutions are not converged to the specified tolerance,
        please specify a large numer of 'iter' ",
           ".\n===============================\n", sep = "")
     }
    }

   m <- ifelse(!is.null(m), m, 60 * (p * (c2t + n * c1t) + (1 - p) *
                                         (c2 + n * c1)))
   var.expr <- quote({
       a^2 * (icc - b^2 + (1 - icc)/n)/(m / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))) +
         b^2 / (p * (1 - p) * (m / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t))))
     })

     var.ab <- eval(var.expr)
     par <- c(par, list(m = m))
     out <- list(n = n, p = p, var.ab = var.ab)
     od.out <- list(funName = funName, designType = designType,
                    par = par, out = out)
     nlab <- labFun(x = nlab, y = "Level-1 Sample Size: n")
     plab <- labFun(x = plab, y = "Proportion Level-2 Units in Treatment: p")
     varlab <- labFun(x = varlab, y = "Variance")
     vartitle <- labFun(x = vartitle, y = "")

     plot.by <- plotbyFun(x = plot.by, y = list(n = "n", p = "p"))
     nrange <- seq(nlim[1], nlim[2], by = 1)
     prange <- seq(plim[1] + 0.02, plim[2] - 0.02, by = 0.01)
     if (length(plot.by) == 2) figure <- par(mfrow = c (1, 2))
     if (length(plot.by) == 1) figure <- par(mfrow = c (1, 1))
     if (plots) {
       if (!is.null(plot.by$n)) {
         plot.y <- NULL
         for (n in nrange)
           plot.y <- c(plot.y, eval(var.expr))
         graphics::plot(nrange, plot.y,
                        type = "l", lty = 1,
                        xlim = nlim, ylim = varlim,
                        xlab = nlab, ylab = varlab,
                        main = vartitle, col = "black")
         graphics::abline(v = out$n, lty = 2, col = "black")
         n <- out$n
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
         graphics::abline(v = out$p, lty = 2, col = "black")
         p <- out$p
       }
     }
     par(figure)
     od.out <- list(par = par, funName = funName,
                    designType = designType, out = out)
     return(od.out)
    }

  if (test == "joint" | test == "Joint" | test == "JOINT") {
      if (is.null(par$p) & is.null(par$n)) {
        n.of.opt.pars <- 2
        if (verbose) {cat('The ACO algorithm started initilization..', ".\n", sep = "")}
        e.abs <- e # absolute error
        e.rel <- e # relative error
        # initiate parameters
        eval <- 0
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
           # n <- m / ((1 - p) * c1 + p * c1t);
            lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
            lambda.b <- b * sqrt(J * (1 - r2m) /
                                   (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

            (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                    df = J - q.a - 2, lambda.a) +
               pt(qt(sig.level/tside, df = J - q.a - 2),
                  df = J - q.a - 2, lambda.a)) *
              (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                      df = J - q.b - 3, lambda.b) +
                 pt(qt(sig.level/tside, df = J - q.b - 3),
                    df = J - q.b - 3, lambda.b))
          })
        } else {
          pwr <- quote({
        #    n <- m / ((1 - p) * c1 + p * c1t);

            lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
            lambda.b <- b * sqrt(J * (1 - r2m) /
                                   (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

            (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                    df = J - q.a - 2, lambda.a)) *
              (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                      df = J - q.b - 3, lambda.b))
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
            J <- stats::uniroot(function(J) eval(pwr) - power.joint, c(4, 10e4))$root
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
          if (sum((0.001 < o.X[i, 1] & o.X[i, 1] < 0.999),
                  (0.3 < o.X[i, 2]  && o.X[i, 2] < 10000)) == n.of.opt.pars) {

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
          J <- stats::uniroot(function(J) eval(pwr) - power.joint, c(4, 10e4))$root
          m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
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
          # n <- m / ((1 - p) * c1 + p * c1t);

          lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
          lambda.b <- b * sqrt(J * (1 - r2m) /
                                 (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

          (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                  df = J - q.a - 2, lambda.a) +
             pt(qt(sig.level/tside, df = J - q.a - 2),
                df = J - q.a - 2, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                    df = J - q.b - 3, lambda.b) +
               pt(qt(sig.level/tside, df = J - q.b - 3),
                  df = J - q.b - 3, lambda.b))
        })
      } else {
        pwr <- quote({
          #    n <- m / ((1 - p) * c1 + p * c1t);

          lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
          lambda.b <- b * sqrt(J * (1 - r2m) /
                                 (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

          (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                  df = J - q.a - 2, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                    df = J - q.b - 3, lambda.b))
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
          J <- stats::uniroot(function(J) eval(pwr) - power.joint, c(4, 10e4))$root
          m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
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
            J <- stats::uniroot(function(J) eval(pwr) - power.joint, c(4, 10e4))$root
            m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
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
          # n <- m / ((1 - p) * c1 + p * c1t);

          lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
          lambda.b <- b * sqrt(J * (1 - r2m) /
                                 (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

          (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                  df = J - q.a - 2, lambda.a) +
             pt(qt(sig.level/tside, df = J - q.a - 2),
                df = J - q.a - 2, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                    df = J - q.b - 3, lambda.b) +
               pt(qt(sig.level/tside, df = J - q.b - 3),
                  df = J - q.b - 3, lambda.b))
        })
      } else {
        pwr <- quote({
          #    n <- m / ((1 - p) * c1 + p * c1t);

          lambda.a <- a * sqrt((p * (1 - p) *J)/(1 - r2m));
          lambda.b <- b * sqrt(J * (1 - r2m) /
                                 (icc *(1 - r22) + (1 - icc) * (1 - r12) / n));

          (1 - pt(qt(1 - sig.level/tside, df = J- q.a - 2),
                  df = J - q.a - 2, lambda.a)) *
            (1 - pt(qt(1 - sig.level/tside, df = J - q.b - 3),
                    df = J - q.b - 3, lambda.b))
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
          J <- stats::uniroot(function(J) eval(pwr) - power.joint, c(4, 10e4))$root
          m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
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
            J <- stats::uniroot(function(J) eval(pwr) - power.joint, c(5, 10e4))$root
            m <- p * J * (c1t * n + c2t) + (1 - p) * J *(c1*n + c2)
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

