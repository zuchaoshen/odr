#' Optimal Design and Statistical Power of Multilevel Randomized Trials
#'
#' This package is to help researchers design cost-efficient experimental studies
#'     assessing main treatment effects
#'     with adequate statistical precision by (a) solving optimal sample allocations,
#'     (b) comparing design precision and efficiency between different sample allocations,
#'     and (c) explicitly accommodating costs and budget in power analyses.
#'
#' The package covers seven types of experiments aiming to detect
#'     main and mediation effects. These experiments are
#'     individual randomized controlled trials (RCTs), two-,
#'     three-, and four-level cluster-randomized trials (CRTs),
#'     two-, three-, and four-level multisite randomized trials (MRTs), and
#'     two-level CRTs investigating mediation effects with group-level mediators.
#'     There are two categorical functions for each type of
#'     experiments and a uniform function for all types of experiments.
#'     The two categorical functions are
#'     'od' and 'power'. The 'od' function can
#'     calculate the optimal sample allocation with and without constraint for
#'     each type of experiments.
#'     The 'power' function by default can calculate required budget
#'     (and required sample size) for desired
#'     power, minimum detectable effect size (MDES) under a fixed budget,
#'     statistical power under a fixed budget.
#'     The 'power' function also can perform conventional power analyses
#'     (e.g., required sample size, power, MDES calculation).
#'     The uniform function 're' (or 'rpe') is to compare
#'     the relative (precision and) efficiency between two designs
#'     with different sample allocations.
#'
#' @author Zuchao Shen, Benjamin Kelcey
#'
#' Maintainer: Zuchao Shen \href{mailto:zuchao.shen@gmail.com}{zuchao.shen@gmail.com}
#'     (Arizona State University)
#'
#' @docType package
"_PACKAGE"
