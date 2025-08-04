#' Plot statistical power curves under a fixed budget across optimal design
#' parameters
#'
#' @description This function plots statistical power curves (for main,
#'     moderation, and/or mediation effects) under a fixed budget
#'     across optimal design parameters.
#' @inheritParams od.2m
#' @param expr Returned objects from an od function (e.g., od.2m, od.2m.mod).
#' @param by Dimensions to plot power curves by the optimal design parameters.
#'     The default value is by all optimal design parameters for a type of design.
#'     For example, default values are by = "p" for single-level designs,
#'     by = c("n", "p") for two-level designs,
#'     and by = c("n", "p", "J") for three-level designs.
#' @param plab Label for the x-axis when the plot is by the optimal design
#'     parameter "p".
#' @param nlab Label for the x-axis when the plot is by the optimal design
#'     parameter "n".
#' @param Jlab Label for the x-axis when the plot is by the optimal design
#'     parameter "J".
#' @param powerlim The limit for plotting power curves.
#' @param powerlab The label for the statistical power.
#' @param legend Logical; present plot legend if TRUE. The default is TRUE.
#' @param plot.title The title of the plot (e.g., plot.title = "Power Curves").
#'      The default is NULL.
#' @export plot.power


plot.power <- function(expr, nlim = c(2, 300), plim = c(0.01, 0.99),
                       powerlim = c(0, 1), plot.title = NULL,
                       m = NULL,
                       by = c("n", "p"), legend = TRUE,
                       nlab = "Level-One Sample Size (n)",
                       plab = "Proportion (p)",
                       Jlab = "Level-Two Sample Size (J)",
                       powerlab = "Statistical Power"){
  if(expr$funName == "od.2m.mod"){
    n = expr$out$n
    p = expr$out$p
    icc <- expr$par$icc
    r12 <- expr$par$r12
    r22m <- expr$par$r22m
    c1 <- expr$par$c1
    c2 <- expr$par$c2
    c1t <- expr$par$c1t
    omega <- expr$par$omega
    gamma <- expr$par$gamma
    q <- expr$par$q
    q.mod <- expr$par$q.mod
    Q <- expr$par$Q
    d <- expr$par$d
    binary <- expr$par$binary
    mod.level <- expr$par$mod.level
    if(is.null(m)){m <- expr$out$m}

  nrange <- seq(nlim[1], nlim[2], by = 1)
  prange <- seq(plim[1], plim[2], by = 0.01)
  if (length(by) == 2) figure <- graphics::par(mfrow = c (1, 2))
  if (length(by) == 1) figure <- graphics::par(mfrow = c (1, 1))
    if ("n" %in% by) {
      power.main <- power.mod <- NULL
      for (n in nrange){
        power.main <- c(power.main, power.2m(m = m, q = q,
                                             p = expr$out$p, icc = icc,
                                             r12 = r12, r22m = r22m,
                                             d = d, c1 = c1, c2 = c2,
                                             omega = omega,
                                             c1t = c1t,
                                             n = n)$out$power);
        power.mod <- c(power.mod, power.2m.mod(expr = expr,
                                               m = m,
                                               gamma = gamma,
                                               n = n)$out$power.mod)
      }
      graphics::plot(nrange, power.mod,
                     type = "l", lty = 1,
                     xlim = nlim, ylim = powerlim,
                     xlab = nlab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$n, lty = 2, col = "black")
      graphics::lines(nrange, power.main, lty = 1, col = "gray")
      if(legend){legend("topright", legend = c("Moderation", "Main"),
             col = c("black", "gray"), lty = 1)}
    }

      if ("p" %in% by) {
        power.main <- power.mod <- NULL
        for (p in prange){
          power.main <- c(power.main, power.2m(m = m, q = q,
                                               n = expr$out$n, icc = icc,
                                               r12 = r12, r22m = r22m,
                                               d = d, c1 = c1, c2 = c2,
                                               c1t = c1t,
                                               omega = omega,
                                               p = p)$out$power)
          power.mod <- c(power.mod, power.2m.mod(expr = expr,
                                                 m = m,
                                                 gamma = gamma,
                                                 p = p)$out$power.mod)
        }
        graphics::plot(prange, power.mod,
                       type = "l", lty = 1,
                       xlim = plim, ylim = powerlim,
                       xlab = plab, ylab = powerlab,
                       main = plot.title, col = "black")
        graphics::abline(v = expr$out$p, lty = 2, col = "black")
        graphics::lines(prange, power.main, lty = 1, col = "gray")
        if(legend){legend("topright", legend = c("Moderation", "Main"),
               col = c("black", "gray"), lty = 1)}
      }

  } else if (expr$funName == "od.2m"){
    n = expr$out$n
    p = expr$out$p
    icc <- expr$par$icc
    r12 <- expr$par$r12
    r22m <- expr$par$r22m
    c1 <- expr$par$c1
    c2 <- expr$par$c2
    c1t <- expr$par$c1t
    omega <- expr$par$omega
    q <- expr$par$q
    d <- expr$par$d
    if(is.null(m)){m <- expr$out$m}

    nrange <- seq(nlim[1], nlim[2], by = 1)
    prange <- seq(plim[1], plim[2], by = 0.01)
    if (length(by) == 2) figure <- graphics::par(mfrow = c (1, 2))
    if (length(by) == 1) figure <- graphics::par(mfrow = c (1, 1))
    if ("n" %in% by) {
      power.main <- NULL
      for (n in nrange){
        power.main <- c(power.main, power.2m(m = m, q = q,
                                             p = expr$out$p, icc = icc,
                                             r12 = r12, r22m = r22m,
                                             d = d, c1 = c1, c2 = c2,
                                             omega = omega,
                                             c1t = c1t,
                                             n = n)$out$power);
      }
      graphics::plot(nrange, power.main,
                     type = "l", lty = 1,
                     xlim = nlim, ylim = powerlim,
                     xlab = nlab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$n, lty = 2, col = "black")
    }

    if ("p" %in% by) {
      power.main <- NULL
      for (p in prange){
        power.main <- c(power.main, power.2m(m = m, q = q,
                                             n = expr$out$n, icc = icc,
                                             r12 = r12, r22m = r22m,
                                             d = d, c1 = c1, c2 = c2,
                                             c1t = c1t,
                                             omega = omega,
                                             p = p)$out$power)
      }
      graphics::plot(prange, power.main,
                     type = "l", lty = 1,
                     xlim = plim, ylim = powerlim,
                     xlab = plab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$p, lty = 2, col = "black")
    }
}

return(figure)
}

