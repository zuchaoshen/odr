#' Optimal Design and Statistical Power for Experimental Studies
#'     Investigating Main, Mediation, and Moderation Effects
#'
#' This package is to help researchers design cost-efficient and well-powered
#'     experimental studies investigating main, mediation, and moderation
#'     effects (and some combinations). Specifically, this package can
#'     (a) identify optimal sample allocations,
#'     (b) compare design efficiency between different sample allocations,
#'     and (c) perform power analyses with and without
#'     accommodating costs and budget.
#'
#' The package covers seven types of experiments aiming to detect
#'     main, moderation, and mediation effects. These experiments are
#'     individual randomized controlled trials (RCTs), two-,
#'     three-, and four-level cluster-randomized trials (CRTs),
#'     two-, three-, and four-level multisite randomized trials (MRTs).
#'     The two categorical functions are
#'     'od' and 'power'. The 'od' function can
#'     calculate the optimal sample allocation with and without constraint for
#'     each type of experiments.
#'     The 'power' function by default can calculate required budget
#'     (and required sample size) for desired
#'     power, minimum detectable effect size (MDES) under a fixed budget,
#'     statistical power under a fixed budget.
#'     The 'power' function can perform conventional power analyses
#'     (e.g., required sample size, power, MDES calculation).
#'     The uniform function 're' (or 'rpe') is to compare
#'     the relative (precision and) efficiency between two designs
#'     with different sample allocations (limited to main effects).
#'
#' @author Zuchao Shen, Benjamin Kelcey
#'
#' Maintainer: Zuchao Shen \href{mailto:zuchao.shen@gmail.com}{zuchao.shen@gmail.com}
#'     (University of Georgia)
#'
#' @docType package
"_PACKAGE"
