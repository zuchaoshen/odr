#' Optimal sample allocation calculation for single-level randomized controlled
#'     trials (RCTs) investigating mediation effects (1-1-1)
#'
#' @description The optimal design of single-level RCTs
#'     probing mediation effects is to identify the optimal sample
#'     allocation that use the minimum budget to achieve a fixed level of
#'     statistical power. The optimal design parameter is the proportion of
#'     individuals/units to be assigned to the experimental condition.
#'     This function identifies the optimal \code{p}.
#'
#' @inheritParams power.1.111
#' @inheritParams od.1
#' @inheritParams od.2.221
#' @param a The treatment effect on the mediator.
#' @param b The within-treatment correlation between the outcome and
#'     the mediator.
#' @param c The cost of sampling an individual in the control group.
#' @param ct The cost of sampling an individual in the treated group.
#' @param n Total number of individuals in the experimental study, the default
#'     value is NULL.
#' @param nlim The interval/range used to numerically solve for n,
#'     the default values are c(6, 1e7).
#' @param two.tailed Two tailed test, the default value is TRUE.
#' @param sig.level Significance level or type I error rate, default value is 0.05.
#' @param q.a The number of covariates at the mediator model
#'     (except the treatment indicator), the default value is zero.
#' @param q.b The number of covariates in the outcome model (except the treatment
#'     indicator and the mediator), the default value is zero.
#' @param test The type of test will be used to detect mediation effects.
#'     The default is the joint significance test (i.e., test = "joint",
#'     "Joint","JOINT"). Another choice is the Sobel test by
#'     specifying the argument as test = "sobel", "Sobel", or "SOBEL".
#' @param power Statistical power specified, default is .80.
#' @param r.mw The within-treatment correlation between the mediator and the
#'     covariate(s) in the mediator model.
#' @param r.mx The within-treatment correlation between the mediator and the
#'     covariate(s) in the outcome model.
#' @param r.yx The within-treatment correlation between the outcome and the
#'     covariate(s) in the outcome model.
#' @param e Maximum error value used when solution quality used as
#'     the stopping criterion, default is 1e-10.
#' @param max.value Maximal value of optimization when used as
#'     the stopping criterion. Default is infinite.
#' @param d.p The initial sampling domains for p. Default is c(0.10, 0.50).
#' @param max.iter Maximal number of function evaluations when used as
#'     the stopping criterion. Default is 200.
#' @param n.of.archive  Size of the solution archive, default is 20.
#' @param q Locality of the search (0,1), default is 0.0001.
#' @param xi  Convergence pressure (0, Inf), suggested: (0, 1), default is 0.5.
#' @param verbose Print out evaluation process if TRUE, default is TRUE.
#' @param n.of.ants Number of ants used in each iteration after
#'    the initialization stage, the default value is 10.
#' @param tol convergence tolerance.
#'
#' @return
#'     Unconstrained or constrained optimal sample allocation \code{p}).
#'     The function also returns statistical power,
#'     function name, design type,
#'     and parameters used in the calculation.
#'
#' @export od.1.111
#' @examples
#' myod <- od.1.111(a = .3, b = .5, c = 10, ct = 100)
#' myod

od.1.111 <- function(a = NULL, b = NULL,
                     c = NULL, ct = NULL, m = NULL,
                     r.yx = 0, r.mx = 0, r.mw = 0,
                     q.a = 0, q.b = 0,
                     test = "joint",
                     p = NULL, n = NULL,
                     tol = 1e-11, power = 0.80,
                     d.p = c(0.1, 0.5),
                     sig.level = 0.05, two.tailed = TRUE,
                     plim = c(.01, .99),
                     varlim = c(0, 0.001),
                     plab = NULL, varlab = NULL,
                     vartitle = NULL,
                     nlim = c(6, 1e6), verbose = TRUE,
                     max.value = Inf, max.iter = 300,  e = 1e-10,
                     n.of.ants = 10, n.of.archive = 20, q = 0.0001,
                     xi = 0.5
) {
  funName <- "od.1.111"
  designType <- "1-1-1 mediation in single-level RCTs"
  par <- list(a = a, b = b,
              r.yx = r.yx, r.mx = r.mx, r.mw = r.mw,
              c = c, ct =ct,
              n = n, p = p, m = m,
              q.a = q.a, q.b = q.b,
              sig.level = sig.level, two.tailed = two.tailed,
              test = test,
              max.iter = max.iter,
              n.of.ants = n.of.ants, n.of.archive = n.of.archive,
              q = q,
              xi = xi
  )

  if (sum(sapply(list(r.yx, r.mx, r.mw, c, ct),
                 function(x) is.null(x))) >= 1)
    stop("All of 'r.yx', 'r.mx', 'r.mw', 'c', and 'ct'
         must be specified")
  NumberCheck <- function(x) {!is.null(x) & !is.numeric(x)}
  if (sum(sapply(list(r.yx, r.mx, r.mw), function(x) {
    NumberCheck(x) | any(0 > x | x > 1)
  })) >= 1)
    stop("'r.yx', 'r.mx', 'r.mw' must be numeric in [0, 1]")
  if (sum(sapply(list(c, ct), function(x) {
    NumberCheck(x) | x < 0})) >= 1)
    stop("'c', 'ct' must be numeric")
  if (c == 0 & ct == 0 & is.null(par$p))
    stop("when c and ct are both zero, p must be constrained,
         please specify a value  for p")

  labFun <- function(x, y) {
    if (!is.null(x) & length(x) == 1 & is.character(x)) {x} else {y}
  }
  tside <- ifelse(two.tailed == TRUE, 2, 1)

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

      #Power
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
        n <- stats::uniroot(function(n) eval(pwr) - power, nlim)$root
        m <- p*n*ct + (1-p)*n*c
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
      {cat('The ACO algorithm finished initilization of ',
           n.of.archive, ' analyses',".\n", sep = "")}

      while (TRUE) { # the algorithm will stop if one of the criteria is met
        dist.mean <- p.X
        if (sum(apply(dist.mean, 2, stats::sd)) == 0) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
                      designType = designType,
                      out = list(m = 1/max.y, p = max.X, n = par$n)))
        }
        dist.rank <- pp$gr
        dim(dist.mean) <- c(length(pp$v), n.of.opt.pars)
        o.X <- vector()
        o.X <- gen.design.pars(dist.mean, dist.rank, n.of.ants, nl, q, n.of.archive, xi)
        # the algorithm will stop if it converges
        if (length(o.X) == 0) {
          return(list(archive = pp, archive.design.pars = p.X,
                      n.iter = n.iter, par = par, funName = funName,
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

          for (j in 1:dim(X)[1]) {
            # redo power analysis with n.of.ants times for those reasonable

            n.iter <- n.iter + 1
            p <- X[j, 1]
            if (verbose) {cat('Number of tried evaluations is ',
                              n.iter, ".\n", sep = "")}
            n <- stats::uniroot(function(n) eval(pwr) - power, nlim)$root
            m <- p*n*ct + (1-p)*n*c
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
        for (i in 1:n.of.archive) {nl[i,] <- (1:n.of.archive)[1:n.of.archive!=i]}

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
    } else if (!is.null(par$p)) {
      cat("===============================\n",
          "There is no calculation performed
        because p is contrained",
          ".\n===============================\n", sep = "")
      return(list(par = par, funName = funName,
                  designType = designType, test = test,
                  out = c(p = par$p, n = par$n)))
    }
}

