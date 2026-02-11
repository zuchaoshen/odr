#' Optimal sample allocation calculation for two-level MRTs detecting main effects
#'
#' @description The optimal design of two-level
#'     multisite randomized trials (MRTs) detecting main effects is to calculate
#'     the optimal sample allocation that minimize the budget to achieve a
#'     fixed statistical power (e.g., 80%) using the ant colony optimization (ACO)
#'     algorithm. Alternatively, the function can calculate the allocation that
#'     minimizes the variance of a treatment effect under a fixed budget,
#'     which is less precise than the ACO algorithm.
#'     The optimal design parameters include
#'     the level-one sample size per site (\code{n})
#'     and the proportion of level-one unit to be assigned to treatment (\code{p}).
#'     This function solves the optimal \code{n} and/or \code{p}
#'     with and without a constraint.
#'
#' @inheritParams od.4
#' @inheritParams power.2m
#' @inheritParams power.3m
#' @inheritParams od.2.221
#'
#' @param m Total budget, default is the total costs of sampling 60
#'     sites.
#' @param p The proportion of level-1 units within each level 2 unit
#'     to be assigned to treatment.
#' @param c2 The cost of sampling one level-2 unit (site).
#' @param d Standardized effect size, default is 0.1.
#' @param nrange The range of the level-1 sample size per level-2 unit
#'     that used to exclude unreasonable values. Default value is c(2, 10000).
#' @param Jlim The range for solving the root of level-2 sample size
#'     (\code{J}) numerically. Change the default values to a larger range
#'     (e.g., starting with a smaller value) if
#'     f() values at end points are not of opposite sign. For example,
#'     use Jlim = c(1.5, 1e10).
#' @param plot.by Plot the variance by \code{n} and/or \code{p};
#'     default value is plot.by = list(n = "n", p = "p").
#' @param plab The plot label for \code{p},
#'     default value is "Proportion Level-1 Units in Treatment: p".
#' @param verbose Logical; print the values of \code{n}
#'    and \code{p} if TRUE, otherwise not; default value is TRUE.
#' @param q.aco Locality of the ACO search (0,1), default is 0.0001.
#' @param q The number of covariates at level 2. Default is 1.
#' @param d.p The initial sampling domains for p. Default is c(0.1, 0.5).
#' @param d.n The initial sampling domain for n. Default is c(2, 1000).
#' @param aco Logic. If TRUE, the function will use the ant colony
#'     optimization (ACO) algorithm to identify optimal allocations. If FALSE,
#'     the function will use the first-order derivative method to identify
#'     optimal allocations. Default is TRUE.
#' @return
#'     Unconstrained or constrained optimal sample allocation
#'     (\code{n} and \code{p}).
#'     The function also returns the variance of the treatment effect,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.2m
#' @references
#'    Shen, Z., & Kelcey, B. (2022). Optimal sample
#'    allocation in multisite randomized trials.
#'    The Journal of Experimental Education, 90(3), 693-711.
#'    <https://doi.org/10.1080/00220973.2020.1830361>
#'
#' @examples
#' # Unconstrained optimal design #---------
#'   myod1 <- od.2m(icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
#'               c1 = 1, c2 = 10, c1t = 10,
#'               varlim = c(0, 0.005))
#'   myod1$out
#'
#' # Constrained optimal design with p = 0.5 #---------
#'   myod2 <- od.2m(icc = 0.2, omega = 0.02,
#'               r12 = 0.5, r22m = 0.5,
#'               c1 = 1, c2 = 10, c1t = 10,
#'               varlim = c(0, 0.005), p = 0.5)
#'   myod2$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod2)
#'   myre$re # RE = 0.81
#'
#' # Constrained optimal design with n = 5 #---------
#'   myod3 <- od.2m(icc = 0.2, omega = 0.02,
#'               r12 = 0.5, r22m = 0.5, c1 = 1, c2 = 10,
#'               c1t = 10, varlim = c(0, 0.005), n = 5)
#'   myod3$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod3)
#'   myre$re # RE = 0.78
#'
#' # Constrained n and p, no calculation performed #---------
#'   myod4 <- od.2m(icc = 0.2, omega = 0.02, r12 = 0.5, r22m = 0.5,
#'               c1 = 1, c2 = 10, c1t = 10,
#'               varlim = c(0, 0.005), p = 0.5, n = 10)
#'   myod4$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod4)
#'   myre$re # RE = 0.79
#'
od.2m <- function(n = NULL, p = NULL, icc = NULL,
                 r12 = NULL, r22m = NULL,
                 c1 = NULL, c2 = NULL,
                 c1t = NULL, omega = NULL,
                 m = NULL, plots = TRUE, plot.by = NULL,
                 nlim = NULL, plim = NULL, varlim = NULL,
                 nlab = NULL, plab = NULL, varlab = NULL,
                 vartitle = NULL, verbose = TRUE, iter = 100,
                 tol = 1e-10, q = 1, d = 0.1,
                 power = 0.8, aco = TRUE,
                 d.p = c(0.1, 0.5), d.n = c(2, 1000),
                 sig.level = 0.05, two.tailed = TRUE,
                 Jlim = c(2.5, 1e+10),
                 nrange = c(2, 10000),
                 max.value = Inf, max.iter = 300,  e = 1e-10,
                 n.of.ants = 10, n.of.archive = 50, q.aco = 0.0001,
                 xi = 0.5) {
  funName <- "od.2m"
  designType <- "two-level MRTs"
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
  if (!is.null(plot.by) && !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n'))")
  if (!is.numeric(iter) || iter < 2)
    stop("'iter' must be numeric with iter >= 2")
  par <- list(icc = icc,
              r12 = r12, r22m = r22m,
              c1 = c1, c2 = c2, c1t = c1t, omega = omega,
              m = m,
              n = n, p = p, iter = iter,
              power = power, aco = aco, d = d, q = q,
              sig.level = sig.level, two.tailed = two.tailed,
              max.iter = max.iter,
              n.of.ants = n.of.ants, n.of.archive = n.of.archive,
              q.aco = q.aco,
              xi = xi
  )

 if(aco){#ACO = TRUE
    tside <- ifelse(two.tailed == TRUE, 2, 1)
    if (two.tailed) {
      pwr.expr <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J) /
                             (p * (1 - p) * n * omega * (1 - r22m) +
                                (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = J - q - 1),
               df = J - q - 1, lambda) +
          pt(qt(sig.level / tside, df = J - q - 1),
             df = J - q - 1, lambda)
      })
    } else {
      pwr.expr <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J) /
                             (p * (1 - p) * n * omega * (1 - r22m) +
                                (1 - icc) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, J - q - 1),
               df = J - q - 1, lambda)
      })}

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
              eval(pwr.expr) - power, Jlim)$root
            m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
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
                eval(pwr.expr) - power, Jlim)$root
              m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
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
                        archive.design.pars = p.X,
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
            eval(pwr.expr) - power, Jlim)$root
          m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
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
                eval(pwr.expr) - power, Jlim)$root
              m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
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
            eval(pwr.expr) - power, Jlim)$root
          m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
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
                eval(pwr.expr) - power, Jlim)$root
              m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
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
          eval(pwr.expr) - power, Jlim)$root
        m <- J *((1 - p) * c1 * n + p * c1t * n + c2)
        return(list(par = par, funName = funName,
                    designType = designType,
                    out = list(m = m,
                               p = par$p,
                               n = par$n, J = J)))
      }
 } else{ #ACO != TRUE
     if (is.null(par$n)) {
       n.expr <- quote({
         sqrt(((1 - icc) * (1 - r12)*c2) /
                ((p * (1 - p) * omega * (1 - r22m)) *
                ((1 - p) * c1 + p * c1t)))
       })
     } else {
       n.expr <- ({par$n})
     }
     limFun <- function(x, y) {
       if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
     }
     nlim <- limFun(x = nlim, y = c(2, 50))
     plim <- limFun(x = plim, y = c(0, 1))
     varlim <- limFun(x = varlim, y = c(0, 0.05))
     if (is.null(par$p)) {
       p.expr <- quote({
         (omega*(1-r22m)*n*p * (1 - p) + (1 - icc) * (1 - r12))*
           (n*c1t-n*c1)*p*(1-p) -
           (1-2*p)*(1 - icc) * (1 - r12)*((1 - p) * c1 * n +  p * c1t * n+c2)
       })
     }
     if (!is.null(par$n)) {
       if (!is.numeric(n) || n <= 0)
         stop("constrained 'n' must be numeric with n > 0")
     } else {
       n <- sample(2:50, 1)
     }
     if (!is.null(par$p)) {
       if (!is.numeric(par$p) || any(par$p <= 0 | par$p >= 1))
         stop("constrained 'p' must be numeric in (0, 1)")
     } else {
       p <- stats::runif(1, min = 0, max = 1)
     }
     nn <- pp <- NULL
     for (i in 1:iter) {
       if (is.null(par$p)) {
         pp[i] <- stats::uniroot(function(p)
           eval(p.expr), plim)$root
         p <- pp[i]
       } else {
         pp[i] <- par$p
       }
       nn[i] <- eval(n.expr); n <- nn[i]
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
     m <- ifelse(!is.null(m), m, 60 * (p * c1t * n + (1 - p) * c1 * n + c2 ))
     var.expr <- quote({
       (omega * (1 - r22m) * n * p * (1 - p) + (1 - icc) * (1 - r12)) /
         (p * (1 - p) * n * (m / (p * c1t * n + (1 - p) * c1 * n + c2)))
     })
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
  }
