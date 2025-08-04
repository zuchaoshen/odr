#' Using the first-order derivative method to identify the optimal sample
#' allocations for moderation effects in two-level multisite randomized
#'     trials (MRTs)
#' @inheritParams od.2m
#' @inheritParams od.2m.mod
#' @inheritParams plot.power
#' @description The optimal design of two-level
#'     MRTs probing moderation effects
#'     identify the optimal sample allocations.
#'     The optimal design parameters include
#'     the level-1 sample size per level-2 unit (\code{n})
#'     and the proportion of level-1 individuals/units assigned to
#'     the experimental group (\code{p}).
#'     This function solves the optimal \code{n} and/or \code{p}
#'     with and without a constraint using the first-order derivative
#'     method to minimize the variance of the moderation effect
#'     estimator. It includes binary or continuous moderators at
#'     level 2 or level 1.
#' @param m The total cost to plot the variance curve. The default value is
#'     the total cost of sampling 60 sites at the optimal allocation.
#' @export od.2m.only.mod
#' @examples
#' myod <- od.2m.only.mod(icc = .2, r12 = .5, r22m = .5,
#'                        c1 = 10, c1t = 100, c2 = 50, omega = .01)
#' myod$out
od.2m.only.mod <- function (icc = NULL, r12 = NULL, r22m = NULL,
                  c1 = NULL, c1t = NULL, c2 = NULL, omega = .01,
                  Q = .50, n = NULL, p = NULL,
                  m = NULL, iter = 300,
                  binary = TRUE,
                  mod.level = 2,
                  nlim = c(2, 300), plim = c(0.01, 0.99),
                  varlim = c(0, 0.005),
                  by = c("n", "p"),
                  varlab = "Variance",
                  nlab = "Level-One Sample Size (n)",
                  plab = "Proportion (p)",
                  vartitle = ""
                  ){
  funName <- "od.2m.only.mod"
  designType <- "two-level MRTs with moderators"
  par <- list(Q = Q, icc = icc, omega = omega,
              mod.level = mod.level, binary = binary,
              r12 = r12, r22m = r22m, c1 = c1, c2 = c2,
              c1t = c1t, n = n, p = p,
              m = m)
  if(binary){var.mod <- Q*(1-Q)} else {var.mod <- 1}
  if(mod.level==1){
    # variance of a moderation effect estimator
    var.expr <- quote({
      J <- m / ((1 - p) * c1 * n +  p * c1t * n + c2);
      ((1 - icc) * (1 - r12))/(p * (1 - p) * n * J*var.mod) })
    # n expression
    n.expr <- n;
    if(is.null(n)){stop("'n' must be constrained/specified for moderators at
                        level 1")}
    # p expression
    if(is.null(par$p)){
      p.expr <- quote({n*c1t*p^2-(1-p)^2*n*c1+2*p*c2-c2})
    }else{
      p.expr <- p
    }
  }else if(mod.level==2){
    # variance of a moderation effect estimator
    var.expr <- quote({
      J <- m / ((1 - p) * c1 * n +  p * c1t * n + c2);
      (omega*(1-r22m)*p * (1 - p) * n *var.mod +(1 - icc) * (1 - r12))/
        (p * (1 - p) * n * J*var.mod) })
    # n expression
    if(is.null(par$n)){
      n.expr <- quote({
        (1 - icc) * (1 - r12)*c2/(omega*(1-r22m)*p * (1 - p) *var.mod*
                                    ((1 - p) * c1 * n +  p * c1t * n))
      })
    }else{n.expr <- n}
    # p expression
    if(is.null(par$p)){
      p.expr <- quote({
        (omega*(1-r22m)*n*p * (1 - p) *var.mod + (1 - icc) * (1 - r12))*
          (n*c1t-n*c1)*p*(1-p) -
        (1-2*p)*(1 - icc) * (1 - r12)*((1 - p) * c1 * n +  p * c1t * n+c2)
      })
    }else{
      p.expr <- p
    }
  }

  if (is.null(par$n)) {
    n <- sample(2:50, 1)
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
  if(is.null(m)){m <- 60*((1 - p) * c1 * n +  p * c1t * n + c2)}
  Var <- eval(var.expr)
  par <- c(par, list(m = m))
  out <- list(n = n, p = p, var = Var)
  od.out <- list(funName = funName, designType = designType,
                 par = par, out = out)
  nrange <- seq(nlim[1], nlim[2], by = 1)
  prange <- seq(plim[1], plim[2], by = 0.01)
  if (length(by) == 2) figure <- par(mfrow = c (1, 2))
  if (length(by) == 1) figure <- par(mfrow = c (1, 1))
    if ("n" %in% by) {
      plot.y <- NULL
      for (n in nrange)
        plot.y <- c(plot.y, eval(var.expr))
      graphics::plot(nrange, plot.y,
                     type = "l", lty = 1,
                     xlim = nlim, ylim = varlim,
                     xlab = nlab, ylab = varlab,
                     main = vartitle, col = "black")
      n <- out$n
      graphics::abline(v = n, lty = 2, col = "Black")
    }
    if ("p" %in% by) {
      plot.y <- NULL
      for (p in prange)
        plot.y <- c(plot.y, eval(var.expr))
      graphics::plot(prange, plot.y,
                     type = "l", lty = 1,
                     xlim = plim, ylim = varlim,
                     xlab = plab, ylab = varlab,
                     main = vartitle, col = "black")
      p <- out$p
      graphics::abline(v = p, lty = 2, col = "Black")
    }
  par(figure)
  return(od.out)
}

