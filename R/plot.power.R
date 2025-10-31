#' Plot statistical power curves under a fixed budget across optimal design
#' parameters
#'
#' @description This function plots statistical power curves (for main,
#'     moderation, and/or mediation effects) under a fixed budget
#'     across optimal design parameters.
#' @inheritParams od.2m
#' @inheritParams od.2
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


plot.power <- function(expr = NULL, nlim = c(2, 300), plim = c(0.01, 0.99),
                       Jlim = c(3, 300),
                       powerlim = c(0, 1), plot.title = NULL,
                       m = NULL, d = NULL,
                       power = .80, q = NULL,
                       by = c("n", "p", "J"), legend = TRUE,
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
  if (length(by) >= 2) figure <- graphics::par(mfrow = c (1, 2))
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
    if (length(by) >= 2) figure <- graphics::par(mfrow = c (1, 2))
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
  } else if (expr$funName == "od.2") {
    n = expr$out$n
    p = expr$out$p
    icc <- expr$par$icc
    r12 <- expr$par$r12
    r22 <- expr$par$r22
    c1 <- expr$par$c1
    c2 <- expr$par$c2
    c1t <- expr$par$c1t
    c2t <- expr$par$c2t
    if(is.null(m)){m <- power.2(expr = expr, d = d, q = q, power = power)$out$m}

    nrange <- seq(nlim[1], nlim[2], by = 1)
    prange <- seq(plim[1], plim[2], by = 0.01)
    if (length(by) >= 2) figure <- graphics::par(mfrow = c (1, 2))
    if (length(by) == 1) figure <- graphics::par(mfrow = c (1, 1))
    if ("n" %in% by) {
      power.main <- NULL
      for (n in nrange){
        power.main <- c(power.main, power.2(m = m, q = q,
                                            p = expr$out$p, icc = icc,
                                            r12 = r12, r22 = r22,
                                            d = d, c1 = c1, c2 = c2,
                                            c1t = c1t, c2t = c2t,
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
        power.main <- c(power.main, power.2(m = m, q = q,
                                             n = expr$out$n, icc = icc,
                                             r12 = r12, r22 = r22,
                                             d = d, c1 = c1, c2 = c2,
                                             c1t = c1t, c2t = c2t,
                                             p = p)$out$power)
      }
      graphics::plot(prange, power.main,
                     type = "l", lty = 1,
                     xlim = plim, ylim = powerlim,
                     xlab = plab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$p, lty = 2, col = "black")
    }

  } else if (expr$funName == "od.3") {
    n = expr$out$n
    p = expr$out$p
    J = expr$out$J
    icc2 <- expr$par$icc2
    icc3 <- expr$par$icc3
    r12 <- expr$par$r12
    r22 <- expr$par$r22
    r32 <- expr$par$r32
    c1 <- expr$par$c1
    c2 <- expr$par$c2
    c3 <- expr$par$c3
    c1t <- expr$par$c1t
    c2t <- expr$par$c2t
    c3t <- expr$par$c3t
    if(is.null(m)){m <- power.3(expr = expr, d = d, q = q, power = power)$out$m}

    nrange <- seq(nlim[1], nlim[2], by = 1)
    Jrange <- seq(Jlim[1], Jlim[2], by = 1)
    prange <- seq(plim[1], plim[2], by = 0.01)
    if (length(by) == 3) figure <- graphics::par(mfrow = c (1, 3))
    if (length(by) == 2) figure <- graphics::par(mfrow = c (1, 2))
    if (length(by) == 1) figure <- graphics::par(mfrow = c (1, 1))
    if ("n" %in% by) {
      power.main <- NULL
      for (n in nrange){
        power.main <- c(power.main, power.3(m = m, q = q, d = d,
                                            icc2 = icc2, icc3 = icc3,
                                            r12 = r12, r22 = r22, r32 = r32,
                                            c1 = c1, c2 = c2, c3 = c3,
                                            c1t = c1t, c2t = c2t, c3t = c3t,
                                            n = n, p = expr$out$p,
                                            J = expr$out$J)$out$power)
      }
      graphics::plot(nrange, power.main,
                     type = "l", lty = 1,
                     xlim = nlim, ylim = powerlim,
                     xlab = nlab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$n, lty = 2, col = "black")
    }

    if ("J" %in% by) {
      power.main <- NULL
      for (J in Jrange){
        power.main <- c(power.main, power.3(m = m, q = q, d = d,
                                            icc2 = icc2, icc3 = icc3,
                                            r12 = r12, r22 = r22, r32 = r32,
                                            c1 = c1, c2 = c2, c3 = c3,
                                            c1t = c1t, c2t = c2t, c3t = c3t,
                                            n = expr$out$n, p = expr$out$p,
                                            J = J)$out$power)
      }
      graphics::plot(Jrange, power.main,
                     type = "l", lty = 1,
                     xlim = Jlim, ylim = powerlim,
                     xlab = Jlab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$J, lty = 2, col = "black")
    }


    if ("p" %in% by) {
      power.main <- NULL
      for (p in prange){
        power.main <- c(power.main, power.3(m = m, q = q, d = d,
                                            icc2 = icc2, icc3 = icc3,
                                            r12 = r12, r22 = r22, r32 = r32,
                                            c1 = c1, c2 = c2, c3 = c3,
                                            c1t = c1t, c2t = c2t, c3t = c3t,
                                            n = expr$out$n, p = p,
                                            J = expr$out$J)$out$power)
      }
      graphics::plot(prange, power.main,
                     type = "l", lty = 1,
                     xlim = plim, ylim = powerlim,
                     xlab = plab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$p, lty = 2, col = "black")
    }



  } else if (expr$funName == "od.3m") {
    n = expr$out$n
    p = expr$out$p
    J = expr$out$J
    icc2 <- expr$par$icc2
    icc3 <- expr$par$icc3
    r12 <- expr$par$r12
    r22 <- expr$par$r22
    r32m <- expr$par$r32m
    c1 <- expr$par$c1
    c2 <- expr$par$c2
    c3 <- expr$par$c3
    c1t <- expr$par$c1t
    c2t <- expr$par$c2t
    omega <- expr$par$omega
    if(is.null(m)){m <- power.3m(expr = expr, d = d, q = q, power = power)$out$m}

    nrange <- seq(nlim[1], nlim[2], by = 1)
    Jrange <- seq(Jlim[1], Jlim[2], by = 1)
    prange <- seq(plim[1], plim[2], by = 0.01)
    if (length(by) == 3) figure <- graphics::par(mfrow = c (1, 3))
    if (length(by) == 2) figure <- graphics::par(mfrow = c (1, 2))
    if (length(by) == 1) figure <- graphics::par(mfrow = c (1, 1))
    if ("n" %in% by) {
      power.main <- NULL
      for (n in nrange){
        power.main <- c(power.main, power.3m(m = m, q = q, d = d,
                                            icc2 = icc2, icc3 = icc3,
                                            r12 = r12, r22 = r22, r32m = r32m,
                                            c1 = c1, c2 = c2, c3 = c3,
                                            c1t = c1t, c2t = c2t, omega = omega,
                                            n = n, p = expr$out$p,
                                            J = expr$out$J)$out$power)
      }
      graphics::plot(nrange, power.main,
                     type = "l", lty = 1,
                     xlim = nlim, ylim = powerlim,
                     xlab = nlab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$n, lty = 2, col = "black")
    }

    if ("J" %in% by) {
      power.main <- NULL
      for (J in Jrange){
        power.main <- c(power.main, power.3m(m = m, q = q, d = d,
                                            icc2 = icc2, icc3 = icc3,
                                            r12 = r12, r22 = r22, r32m = r32m,
                                            c1 = c1, c2 = c2, c3 = c3,
                                            c1t = c1t, c2t = c2t, omega = omega,
                                            n = expr$out$n, p = expr$out$p,
                                            J = J)$out$power)
      }
      graphics::plot(Jrange, power.main,
                     type = "l", lty = 1,
                     xlim = Jlim, ylim = powerlim,
                     xlab = Jlab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$J, lty = 2, col = "black")
    }


    if ("p" %in% by) {
      power.main <- NULL
      for (p in prange){
        power.main <- c(power.main, power.3m(m = m, q = q, d = d,
                                            icc2 = icc2, icc3 = icc3,
                                            r12 = r12, r22 = r22, r32m = r32m,
                                            c1 = c1, c2 = c2, c3 = c3,
                                            c1t = c1t, c2t = c2t, omega = omega,
                                            n = expr$out$n, p = p,
                                            J = expr$out$J)$out$power)
      }
      graphics::plot(prange, power.main,
                     type = "l", lty = 1,
                     xlim = plim, ylim = powerlim,
                     xlab = plab, ylab = powerlab,
                     main = plot.title, col = "black")
      graphics::abline(v = expr$out$p, lty = 2, col = "black")
    }

  } else if (expr$funName == "od.1") {
    p <- expr$out$p
    r12 <- expr$par$r12
    c1 <- expr$par$c1
    c1t <- expr$par$c1t
    if(is.null(m)){m <- power.1(expr = expr, d = d, q = q, power = power)$out$m}
    prange <- seq(plim[1], plim[2], by = 0.01)
    if (length(by) >= 1) figure <- graphics::par(mfrow = c(1, 1))
    if ("p" %in% by) {
      power.main <- NULL
      for (p in prange){
        power.main <- c(power.main, power.1(m = m, q = q,
                                            r12 = r12,
                                            d = d, c1 = c1,
                                            c1t = c1t,
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

