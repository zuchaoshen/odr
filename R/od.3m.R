#' Optimal sample allocation calculation for three-level MRTs detecting main effects
#'
#' @description The optimal design of two-level
#'     multisite randomized trials (MRTs) detecting main effects is to calculate
#'     the optimal sample allocation that minimize the budget to achieve a
#'     fixed statistical power (e.g., 80%) using the ant colony optimization (ACO)
#'     algorithm. Alternatively, the function can calculate the optimal allocation
#'     that minimizes the variance of a treatment effect under a fixed budget,
#'     which is less precise than the ACO algorithm.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n}),
#'     the level-2 sample size per level-3 unit (\code{J}),
#'     and the proportion of level-2 unit to be assigned to treatment (\code{p}).
#'     This function solves the optimal \code{n}, \code{J} and/or \code{p}
#'     with and without constraints.
#'
#' @inheritParams od.2m
#' @inheritParams od.4
#' @inheritParams power.3m
#' @inheritParams power.4m
#' @inheritParams od.2.221
#' @param p The proportion of level-2 units within each level-3 site
#'     to be assigned to treatment.
#' @param m Total budget, default is the total costs of sampling 60
#'     level-3 units.
#' @param c3 The cost of sampling one level-3 unit (site).
#' @param Jrange The range of the level-2 sample size per level-3 unit
#'     that used to exclude unreasonable values. Default value is
#'     c(2, 10000).
#' @param plot.by Plot the variance by \code{n}, \code{J} and/or \code{p};
#'     default value is plot.by = list(n = "n", J = "J", p = "p").
#' @param plab The plot label for \code{p},
#'     default value is "Proportion Level-2 Units in Treatment: p".
#' @param q The number of covariates at level 2. Default is 1.
#' @param d.n The initial sampling domain for n. Default is c(2, 1000).
#' @param d.J The initial sampling domain for J. Default is c(2, 1000).
#' @param ACO Logic. If TRUE, the function will use the ant colony
#'     optimization (ACO) algorithm to identify optimal allocations. If FALSE,
#'     the function will use the first-order derivative method to identify
#'     optimal allocations. Default is TRUE.
#' @param verbose Logical; print the values of \code{n}, \code{J},
#'    and \code{p} if TRUE, otherwise not; default value is TRUE.#'
#' @return
#'     Unconstrained or constrained optimal sample allocation
#'     (\code{n}, \code{J}, and \code{p}).
#'     The function also returns the variance of the treatment effect,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.3m
#'
#' @references
#'   Shen, Z., & Kelcey, B. (2022). Optimal sampling ratios in three-level
#'   multisite experiments. Journal of Research on Educational Effectiveness,
#'   15(1), 130-150.
#'   <https://doi.org/10.1080/19345747.2021.1953200>
#'
#' @examples
#' # Unconstrained optimal design #---------
#'   myod1 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005))
#'   myod1$out # output
#' # Plots by p and J
#'   myod1 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005), plot.by = list(p = 'p', J = 'J'))
#'
#' # Constrained optimal design with p = 0.5 #---------
#'   myod2 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005), p = 0.5)
#'   myod2$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod2)
#'   myre$re # RE = 0.81
#'
#' # Constrained optimal design with n = 5 #---------
#'   myod3 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005), n = 5)
#'   myod3$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod3)
#'   myre$re # RE = 0.89
#'
#' # Constrained n, J and p, no calculation performed #---------
#'   myod4 <- od.3m(icc2 = 0.2, icc3 = 0.1, omega = 0.02,
#'               r12 = 0.5, r22 = 0.5, r32m = 0.5,
#'               c1 = 1, c2 = 5,
#'               c1t = 1, c2t = 200, c3 = 200,
#'               varlim = c(0, 0.005), p = 0.5, n = 15, J = 20)
#'   myod4$out
#' # Relative efficiency (RE)
#'   myre <- re(od = myod1, subod= myod4)
#'   myre$re # RE = 0.75
#'
od.3m <- function(n = NULL, J = NULL, p = NULL, icc2 = NULL, icc3 = NULL,
                 r12 = NULL, r22 = NULL, r32m = NULL,
                 c1 = NULL, c2 = NULL, c3 = NULL,
                 c1t = NULL, c2t = NULL, omega = NULL,
                 m = NULL, plots = TRUE, plot.by = NULL,
                 nlim = NULL, Jlim = NULL, plim = NULL, varlim = NULL,
                 nlab = NULL, Jlab = NULL, plab = NULL, varlab = NULL,
                 Klim = c(6, 1e10), q = 1, d = 0.1,
                 vartitle = NULL,verbose = TRUE, iter = 100, tol = 1e-10,
                 power = 0.8, ACO = TRUE,
                 d.p = c(0.5, 0.9), d.n = c(2, 100), d.J = c(2, 100),
                 sig.level = 0.05, two.tailed = TRUE,
                 nrange = c(2, 1000), Jrange = c(2, 1000),
                 max.value = Inf, max.iter = 300,  e = 1e-10,
                 n.of.ants = 10, n.of.archive = 50, q.aco = 0.0001,
                 xi = 0.5) {
  funName <- "od.3m"
  designType <- "three-level MRTs"
  NumberCheck <- function(x) {!is.null(x) && !is.numeric(x)}
  if (sum(sapply(list(icc2, icc3, r12, r22, r32m,
                      c1, c2, c3, c1t, c2t, omega),
                 function(x) is.null(x))) >= 1)
    stop("All of 'icc2', 'icc3', 'r12', 'r22', 'r32m',
         'c1', 'c2', 'c3', 'c1t', 'c2t', and 'omega' must be specified")
    if (sum(sapply(list(icc2, icc3, r12, r22, r32m, omega), function(x) {
    NumberCheck(x) || any(0 > x | x > 1)
  })) >= 1)
    stop("'icc2', 'icc3', 'r12', 'r22', 'r32m', and 'omega' must be numeric in [0, 1]")
  if (sum(sapply(list(c1, c2, c3, c1t, c2t), function(x) {
    NumberCheck(x) || x < 0})) >= 1)
    stop("'c1', 'c2', 'c3', 'c1t', and 'c2t' must be numeric in [0, inf)")
  if (!is.null(plot.by) && !is.list(plot.by))
    stop("'plot.by' must be in list format (e.g., plot.by = list(n = 'n', J = 'J'))")
  if (!is.numeric(iter) || iter < 2)
    stop("'iter' must be numeric with iter >= 2")
  par <- list(icc2 = icc2, icc3 = icc3, r12 = r12, r22 = r22, r32m = r32m,
              c1 = c1, c2 = c2, c3 = c3,
              c1t =c1t, c2t = c2t,  omega = omega,
              n = n, J = J, p = p, iter = iter,
              tol = tol,
              power = power, ACO = ACO,
              d.p = d.p, d.n = d.n, d.J = d.J,
              sig.level = sig.level, two.tailed = two.tailed,
              nrange = nrange, Jrange = Jrange,
              max.value = max.value, max.iter = max.iter,  e = e,
              n.of.ants = n.of.ants, n.of.archive = n.of.archive, q.aco = q.aco,
              xi = xi)
  if(ACO){ # if ACO == TRUE
    tside <- ifelse(two.tailed == TRUE, 2, 1)
    if (two.tailed == TRUE) {
      pwr.expr <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J * K) /
                             (p * (1 - p) * n * J * omega * (1 - r32m) +
                                n * icc2 * (1 - r22)  +
                                (1 - icc2 - icc3) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = (K - q - 1) ),
               df = (K - q - 1) , lambda) +
          pt(qt(sig.level / tside, df = (K - q - 1)),
             df = (K - q - 1), lambda)
      })
    } else {
      pwr.expr <- quote({
        lambda <- d * sqrt((p * (1 - p) * n * J * K) /
                             (p * (1 - p) * n * J * omega * (1 - r32m) +
                                n * icc2 * (1 - r22)  +
                                (1 - icc2 - icc3) * (1 - r12)));
        1 - pt(qt(1 - sig.level / tside, df = (K - q - 1)),
               df = (K - q - 1), lambda)
      })
    }
      if (is.null(par$p) & is.null(par$n) & is.null(J)){#Unconstrained
        n.of.opt.pars <- 3
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
        n.of.initial <- round(n.of.archive^(1/3), 0)
        n.initial <- seq(from = d.n[1], to = d.n[2], length = n.of.initial)
        p.initial <- seq(from = d.p[1], to = d.p[2], length = n.of.initial)
        J.initial <- seq(from = d.J[1], to = d.J[2], length = n.of.initial)
        n.of.archive <- n.of.initial^3
        nl <- matrix(NA, n.of.archive, n.of.archive-1)
        X <- NULL
        p.X <- NULL
        y <- NULL
        budget <- NULL
        for (n in n.initial){
          for (p in p.initial){
            for (J in J.initial){
              X <- rbind(X, c(p, n, J))
              p.X <- rbind(p.X, c(p, n, J))
              K <- stats::uniroot(function(K)
                eval(pwr.expr) - power, Klim)$root
              m = K*((1 - p) * (c1 * n * J + c2 * J) +
                       p * (c1t * n * J + c2t * J) + c3)
              y <- c(y, 1/m)
              budget <- c(budget, m)
            }
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
          if (sum(apply(dist.mean, 2, stats::sd)) <= e) {
            m = 1/max.y; p = max.X[1]; n = max.X[2]; J = max.X[3];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                        p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "n", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          dist.rank <- pp$gr
          dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
          o.X <- vector()
          o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants,
                                 nl, q.aco, n.of.archive, xi)
          if (length(o.X) == 0) {
            m = 1/max.y; p = max.X[1]; n = max.X[2]; J = max.X[3];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "n", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          X <- NULL
          for (i in 1:n.of.ants){ # exclude unreasonable values
            if (sum((0.001 < o.X[i, 1] & o.X[i, 1] < 0.999),
                    (nrange[1] < o.X[i, 2]  && o.X[i, 2] < nrange[2]),
                    (Jrange[1] < o.X[i, 3]  && o.X[i, 3] < Jrange[2])) == n.of.opt.pars) {
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
              J <- X[j, 3]
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
              K <- stats::uniroot(function(K)
                eval(pwr.expr) - power, Klim)$root
              m = K*((1 - p) * (c1 * n * J + c2 * J) +
                       p * (c1t * n * J + c2t * J) + c3)
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
            m = 1/max.y; p = max.X[1]; n = max.X[2]; J = max.X[3];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "n", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }

          # check if the maximum allowed number of objective function
          # evaluations has not been exceeded
          if (n.iter >= max.iter) {
            m = 1/max.y; p = max.X[1]; n = max.X[2]; J = max.X[3];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "n", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
        }
      } else if (!is.null(par$p) & is.null(par$n) & is.null(J)){#Constrained p
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
        n.of.initial <- round(n.of.archive^(1/2), 0)
        n.initial <- seq(from = d.n[1], to = d.n[2], length = n.of.initial)
        J.initial <- seq(from = d.J[1], to = d.J[2], length = n.of.initial)
        n.of.archive <- n.of.initial^2
        nl <- matrix(NA, n.of.archive, n.of.archive-1)
        X <- NULL
        p.X <- NULL
        y <- NULL
        budget <- NULL
        for (n in n.initial){
          for (J in J.initial){
            X <- rbind(X, c(n, J))
            p.X <- rbind(p.X, c(n, J))
            K <- stats::uniroot(function(K)
              eval(pwr.expr) - power, Klim)$root
            m = K*((1 - p) * (c1 * n * J + c2 * J) +
                     p * (c1t * n * J + c2t * J) + c3)
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
          if (sum(apply(dist.mean, 2, stats::sd)) <= e) {
            m = 1/max.y; p = par$p; n = max.X[1]; J = max.X[2];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("n", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          dist.rank <- pp$gr
          dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
          o.X <- vector()
          o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants,
                                 nl, q.aco, n.of.archive, xi)
          if (length(o.X) == 0) {
            m = 1/max.y; p = par$p; n = max.X[1]; J = max.X[2];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("n", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          X <- NULL
          for (i in 1:n.of.ants){ # exclude unreasonable values
            if (sum((nrange[1] < o.X[i, 1]  && o.X[i, 1] < nrange[2]),
                    (Jrange[1] < o.X[i, 2]  && o.X[i, 2] < Jrange[2])) == n.of.opt.pars) {
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
              J <- X[j, 2]
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
              K <- stats::uniroot(function(K)
                eval(pwr.expr) - power, Klim)$root
              m = K*((1 - p) * (c1 * n * J + c2 * J) +
                       p * (c1t * n * J + c2t * J) + c3)
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
            m = 1/max.y; p = par$p; n = max.X[1]; J = max.X[2];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("n", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }

          # check if the maximum allowed number of objective function
          # evaluations has not been exceeded
          if (n.iter >= max.iter) {
            m = 1/max.y; p = par$p; n = max.X[1]; J = max.X[2];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("n", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
        }
      } else if (is.null(par$p) & !is.null(par$n) & is.null(J)){#Constrained n

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
        n.of.initial <- round(n.of.archive^(1/2), 0)
        p.initial <- seq(from = d.p[1], to = d.p[2], length = n.of.initial)
        J.initial <- seq(from = d.J[1], to = d.J[2], length = n.of.initial)
        n.of.archive <- n.of.initial^2
        nl <- matrix(NA, n.of.archive, n.of.archive-1)
        X <- NULL
        p.X <- NULL
        y <- NULL
        budget <- NULL
        for (p in p.initial){
          for (J in J.initial){
            X <- rbind(X, c(p, J))
            p.X <- rbind(p.X, c(p, J))
            K <- stats::uniroot(function(K)
              eval(pwr.expr) - power, Klim)$root
            m = K*((1 - p) * (c1 * n * J + c2 * J) +
                     p * (c1t * n * J + c2t * J) + c3)
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
          if (sum(apply(dist.mean, 2, stats::sd)) <= e) {
            m = 1/max.y; p = max.X[1]; n = par$n; J = max.X[2];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          dist.rank <- pp$gr
          dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
          o.X <- vector()
          o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants,
                                 nl, q.aco, n.of.archive, xi)
          if (length(o.X) == 0) {
            m = 1/max.y; p = max.X[1]; n = par$n; J = max.X[2];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          X <- NULL
          for (i in 1:n.of.ants){ # exclude unreasonable values
            if (sum((0.001 < o.X[i, 1] & o.X[i, 1] < 0.999),
                    (Jrange[1] < o.X[i, 2]  && o.X[i, 2] < Jrange[2])) == n.of.opt.pars) {
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
              J <- X[j, 2]
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
              K <- stats::uniroot(function(K)
                eval(pwr.expr) - power, Klim)$root
              m = K*((1 - p) * (c1 * n * J + c2 * J) +
                       p * (c1t * n * J + c2t * J) + c3)
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
            m = 1/max.y; p = max.X[1]; n = par$n; J = max.X[2];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }

          # check if the maximum allowed number of objective function
          # evaluations has not been exceeded
          if (n.iter >= max.iter) {
            m = 1/max.y; p = max.X[1]; n = par$n; J = max.X[2];
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "J");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
        }
      } else if (is.null(par$p) & is.null(par$n) & !is.null(J)){
        # Constrained J
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
        n.of.initial <- round(n.of.archive^(1/2), 0)
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
            K <- stats::uniroot(function(K)
              eval(pwr.expr) - power, Klim)$root
            m = K*((1 - p) * (c1 * n * J + c2 * J) +
                     p * (c1t * n * J + c2t * J) + c3)
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
          if (sum(apply(dist.mean, 2, stats::sd)) <= e) {
            m = 1/max.y; p = max.X[1]; n = max.X[2]; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "n");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          dist.rank <- pp$gr
          dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
          o.X <- vector()
          o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants,
                                 nl, q.aco, n.of.archive, xi)
          if (length(o.X) == 0) {
            m = 1/max.y; p = max.X[1]; n = max.X[2]; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "n");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
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
              K <- stats::uniroot(function(K)
                eval(pwr.expr) - power, Klim)$root
              m = K*((1 - p) * (c1 * n * J + c2 * J) +
                       p * (c1t * n * J + c2t * J) + c3)
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
            m = 1/max.y; p = max.X[1]; n = max.X[2]; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "n");
            return(list(archive = pp, archive.design.pars = p.X,
                        archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }

          # check if the maximum allowed number of objective function
          # evaluations has not been exceeded
          if (n.iter >= max.iter) {
            m = 1/max.y; p = max.X[1]; n = max.X[2]; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- c("p", "n");
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
        }
      } else if (!is.null(par$p) & !is.null(par$n) & is.null(J)){
        # Contrained p and n

        n.of.opt.pars <- 1
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
        n.of.initial <- round(n.of.archive, 0)
        J.initial <- seq(from = d.J[1], to = d.J[2], length = n.of.initial)
        n.of.archive <- n.of.initial
        nl <- matrix(NA, n.of.archive, n.of.archive-1)
        X <- NULL
        p.X <- NULL
        y <- NULL
        budget <- NULL
        for (J in J.initial){
          X <- rbind(X, J)
          p.X <- rbind(p.X, J)
          K <- stats::uniroot(function(K)
            eval(pwr.expr) - power, Klim)$root
          m = K*((1 - p) * (c1 * n * J + c2 * J) +
                   p * (c1t * n * J + c2t * J) + c3)
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
          if (sum(apply(dist.mean, 2, stats::sd)) <= e) {
            m = 1/max.y; p = par$p; n = par$n; J = max.X;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "J";
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          dist.rank <- pp$gr
          dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
          o.X <- vector()
          o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants,
                                 nl, q.aco, n.of.archive, xi)
          if (length(o.X) == 0) {
            m = 1/max.y; p = par$p; n = par$n; J = max.X;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "J";
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          X <- NULL
          for (i in 1:n.of.ants){ # exclude unreasonable values
            if (sum(Jrange[1] < o.X[i]  && o.X[i] < Jrange[2]) == n.of.opt.pars) {
              X <- rbind(X, o.X[i])
            }
          }
          if(length(X)>0) {
            p.X <- rbind(p.X, X)
            dim(X) <- c(length(X)/n.of.opt.pars, n.of.opt.pars)

            for (j in 1:dim(X)[1]) {
              # redo power analysis with n.of.ants times for those reasonable
              n.iter <- n.iter + 1
              J <- X[j]
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
              K <- stats::uniroot(function(K)
                eval(pwr.expr) - power, Klim)$root
              m = K*((1 - p) * (c1 * n * J + c2 * J) +
                       p * (c1t * n * J + c2t * J) + c3)
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
            m = 1/max.y; p = par$p; n = par$n; J = max.X;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "J";
            return(list(archive = pp, archive.design.pars = p.X,
                        archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }

          # check if the maximum allowed number of objective function
          # evaluations has not been exceeded
          if (n.iter >= max.iter) {
            m = 1/max.y; p = par$p; n = par$n; J = max.X;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "J";
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
        }
      } else if (!is.null(par$p) & is.null(par$n) & !is.null(J)){
        # Contrained p and J

        n.of.opt.pars <- 1
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
          K <- stats::uniroot(function(K)
            eval(pwr.expr) - power, Klim)$root
          m = K*((1 - p) * (c1 * n * J + c2 * J) +
                   p * (c1t * n * J + c2t * J) + c3)
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
          if (sum(apply(dist.mean, 2, stats::sd)) <= e) {
            m = 1/max.y; p = par$p; n = max.X; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "n";
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          dist.rank <- pp$gr
          dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
          o.X <- vector()
          o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants,
                                 nl, q.aco, n.of.archive, xi)
          if (length(o.X) == 0) {
            m = 1/max.y; p = par$p; n = max.X; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "n";
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          X <- NULL
          for (i in 1:n.of.ants){ # exclude unreasonable values
            if (sum(nrange[1] < o.X[i]  && o.X[i] < nrange[2]) == n.of.opt.pars) {
              X <- rbind(X, o.X[i])
            }
          }
          if(length(X)>0) {
            p.X <- rbind(p.X, X)
            dim(X) <- c(length(X)/n.of.opt.pars, n.of.opt.pars)

            for (j in 1:dim(X)[1]) {
              # redo power analysis with n.of.ants times for those reasonable
              n.iter <- n.iter + 1
              n <- X[j]
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
              K <- stats::uniroot(function(K)
                eval(pwr.expr) - power, Klim)$root
              m = K*((1 - p) * (c1 * n * J + c2 * J) +
                       p * (c1t * n * J + c2t * J) + c3)
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
            m = 1/max.y; p = par$p; n = max.X; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "n";
            return(list(archive = pp, archive.design.pars = p.X,
                        archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }

          # check if the maximum allowed number of objective function
          # evaluations has not been exceeded
          if (n.iter >= max.iter) {
            m = 1/max.y; p = par$p; n = max.X; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "n";
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
        }
      } else if (is.null(par$p) & !is.null(par$n) & !is.null(J)){
        # Contrained n and J

        n.of.opt.pars <- 1
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
          K <- stats::uniroot(function(K)
            eval(pwr.expr) - power, Klim)$root
          m = K*((1 - p) * (c1 * n * J + c2 * J) +
                   p * (c1t * n * J + c2t * J) + c3)
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
          if (sum(apply(dist.mean, 2, stats::sd)) <= e) {
            m = 1/max.y; p = max.X; n = par$n; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "p";
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          dist.rank <- pp$gr
          dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
          o.X <- vector()
          o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants,
                                 nl, q.aco, n.of.archive, xi)
          if (length(o.X) == 0) {
            m = 1/max.y; p = max.X; n = par$n; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "p";
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
          X <- NULL
          for (i in 1:n.of.ants){ # exclude unreasonable values
            if (sum(0.001 < o.X[i]  && o.X[i] < 0.999) == n.of.opt.pars) {
              X <- rbind(X, o.X[i])
            }
          }
          if(length(X)>0) {
            p.X <- rbind(p.X, X)
            dim(X) <- c(length(X)/n.of.opt.pars, n.of.opt.pars)

            for (j in 1:dim(X)[1]) {
              # redo power analysis with n.of.ants times for those reasonable
              n.iter <- n.iter + 1
              p <- X[j]
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
              K <- stats::uniroot(function(K)
                eval(pwr.expr) - power, Klim)$root
              m = K*((1 - p) * (c1 * n * J + c2 * J) +
                       p * (c1t * n * J + c2t * J) + c3)
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
            m = 1/max.y; p = max.X; n = par$n; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "p";
            return(list(archive = pp, archive.design.pars = p.X,
                        archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }

          # check if the maximum allowed number of objective function
          # evaluations has not been exceeded
          if (n.iter >= max.iter) {
            m = 1/max.y; p = max.X; n = par$n; J = par$J;
            K = m /((1 - p) * (c1 * n * J + c2 * J) +
                      p * (c1t * n * J + c2t * J) + c3);
            colnames(p.X) <- "p";
            return(list(archive = pp, archive.design.pars = p.X,
                        n.iter = n.iter, par = par, funName = funName,
                        designType = designType,
                        out = list(m = m, p = p, n = n,
                                   J = J, K = K)))
          }
        }
      } else if (!is.null(par$p) & !is.null(par$n) & !is.null(J)){
        cat("===============================\n",
            "There is no optimization performed
        because all p, n, and J are contrained",
            ".\n===============================\n", sep = "")
        p = par$p; n = par$n; J = par$J;
        K <- stats::uniroot(function(K) eval(pwr.expr) - power, Klim)$root
        m = K* ((1 - p) * (c1 * n * J + c2 * J) +
                  p * (c1t * n * J + c2t * J) + c3);
        return(list(archive = NA, archive.design.pars = NA,
                    n.iter = NA, par = par, funName = funName,
                    designType = designType,
                    out = list(m = m, p = p, n = n,
                               J = J, K = K)))
      }
  } else {
  if (is.null(n)) {
    n.expr <- quote({
      sqrt(((1 - icc2 - icc3) * (1 - r12)) /
        (p * (1 - p) *J * omega * (1 - r32m) + icc2 * (1 - r22)) *
        ((1 - p) * J * c2 + p * J * c2t + c3) /
       ((1 - p) * c1 * J + p * c1t * J))
    })
  } else {
    n.expr <- ({n})
  }
  if (is.null(J)) {
    J.expr <- quote({
      sqrt((n * icc2 * (1 - r22) + (1 - icc2 - icc3) * (1 - r12)) /
             (p * (1 - p) *n * omega * (1 - r32m)) *
             c3 / ((1 - p) * (c1 * n + c2) + p * (c1t * n + c2t)))
    })
  } else {
    J.expr <- ({J})
  }

  limFun <- function(x, y) {
    if (!is.null(x) && length(x) == 2 && is.numeric(x)) {x} else {y}
  }
  nlim <- limFun(x = nlim, y = c(2, 50))
  Jlim <- limFun(x = Jlim, y = c(2, 50))
  plim <- limFun(x = plim, y = c(0, 1))
  varlim <- limFun(x = varlim, y = c(0, 0.05))
  if (is.null(p)) {
    p.expr <- quote({
     (n * J * omega * (1 - r32m) * p * (1 - p) + n * icc2 * (1 - r22) +
        (1 - icc2 - icc3) * (1 - r12)) * (J * (c1t * n + c2t) -
        J * (c1 * n + c2)) * p * (1 - p) -
        (1 - 2 * p) * ((1 - p) * J * (c1 * n + c2) + p * J * (c1t * n + c2t) + c3) *
        (n * icc2 * (1 - r22) + (1 - icc2 - icc3) *(1 - r12))
    })
  }
  if (!is.null(n)) {
    if (!is.numeric(n) || n <= 0)
      stop("constrained 'n' must be numeric with n > 0")
  } else {
    n <- sample(2:50, 1)
  }
  if (!is.null(J)) {
    if (!is.numeric(J) || J <= 0)
      stop("constrained 'J' must be nu meric with J > 0")
  } else {
    J <- sample(2:50, 1)
  }
  if (!is.null(p)) {
    if (!is.numeric(p) || any(p <=0 | p >= 1))
      stop("constrained 'p' must be numeric in (0, 1)")
    p.constr <- p
  } else {
    p.constr <- NULL
    p <- stats::runif(1, min = 0, max = 1)
  }
  nn <- JJ <- pp <- NULL
  for (i in 1:iter) {
    if (is.null(p.constr)) {
      pp[i] <- stats::uniroot(function(p)
        eval(p.expr), plim)$root
      p <- pp[i]
    } else {
      pp[i] <- p
    }
    n <- eval(n.expr); nn[i] <- n
    J <- eval(J.expr); JJ[i] <- J
  }
  if (!is.null(par$n) && !is.null(par$J) && !is.null(par$p)) {
    cat("===============================\n",
        "All of n, J and p are constrained, there is no calculation from other parameters",
        ".\n===============================\n", sep = "")
  }
  if (verbose) {
    if (!is.null(par$n)) {
      cat("The constrained level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    } else {
      cat("The optimal level-1 sample size per level-2 unit (n) is ", n, ".\n", sep = "")
    }
    if (!is.null(par$J)) {
      cat("The constrained level-2 sample size per level-3 unit (J) is ", J, ".\n", sep = "")
    } else {
      cat("The optimal level-2 sample size per level-3 unit (J) is ", J, ".\n", sep = "")
    }
    if (!is.null(par$p)) {
      cat("The constrained proportion of level-2 units in treatment (p) is ", p, ".\n", "\n", sep = "")
    } else {
      cat("The optimal proportion of level-2 units in treatment (p) is ", p, ".\n", "\n" ,sep = "")
    }
  }
  if (nn[iter] - nn[iter-1] <= tol && JJ[iter] - JJ[iter-1] <= tol &&
      pp[iter] - pp[iter-1] <= tol) {
    p <- pp[iter]
    nn <- JJ <- pp <- NULL
  } else {
    cat("===============================\n",
        "The solutions are not converged to specified tolerance,
        please specify a large numer of 'iter' to replace the default value of 100",
        ".\n===============================\n", sep = "")
  }
  m <- ifelse(!is.null(m), m, 60 * (p * (c1t * n * J + c2t * J) +
                                      (1 - p) * (c1 * n * J + c2 * J) + c3))
  var.expr <- quote({
    K <- m / ((1 - p) * (c1 * n * J + c2 * J ) +
                    p * (c1t * n * J + c2t * J) + c3)
    (omega * (1 - r32m) * n * J * p * (1 - p) + icc2 * (1 - r22) * n +
        (1 - icc2 - icc3) * (1 - r12)) / (p * (1 - p) * n * J * K)
  })
  Var <- eval(var.expr)
  par <- c(par, list(m = m))
  out <- list(n = n, J = J, p = p, var = Var)
  od.out <- list(funName = funName, designType = designType,
                 par = par, out = out)
  labFun <- function(x, y) {
    if (!is.null(x) && length(x) == 1 && is.character(x)) {x} else {y}
  }
  nlab <- labFun(x = nlab, y = "Level-1 Sample Size: n")
  Jlab <- labFun(x = Jlab, y = "Level-2 Sample Size: J")
  plab <- labFun(x = plab, y = "Proportion Level-2 Units in Treatment: p")
  varlab <- labFun(x = varlab, y = "Variance")
  vartitle <- labFun(x = vartitle, y = "")
  plotbyFun <- function(x, y) {
    if (!is.null(x) && is.list(x)) {x} else {y}
  }
  plot.by <- plotbyFun(x = plot.by, y = list(n = "n", J = "J", p = "p"))
  nrange <- seq(nlim[1], nlim[2], by = 1)
  Jrange <- seq(Jlim[1], Jlim[2], by = 1)
  prange <- seq(plim[1] + 0.05, plim[2] - 0.05, by = 0.01)
  if (length(plot.by) == 3) figure <- par(mfrow = c (1, 3))
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
    if (!is.null(plot.by$J)) {
      plot.y <- NULL
      for (J in Jrange)
        plot.y <- c(plot.y, eval(var.expr))
      graphics::plot(Jrange, plot.y,
           type = "l", lty = 1,
           xlim = Jlim, ylim = varlim,
           xlab = Jlab, ylab = varlab,
           main = vartitle, col = "black")
      J <- out$J
      graphics::abline(v = J, lty = 2, col = "Blue")
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


