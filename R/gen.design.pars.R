#' Generate optimal design parameters using ant colony optimization
#'
#' @description This function can generate a set of optimal design parameters
#'     based on given distributions of the rank of optimization target
#'     (or budget).
#'
#' @param dist.mean List of means - coordinates
#' @param dist.rank Rank of the archived values of the objective function(s)
#' @param nl Neighborhood of the search area
#' @param q The locality of the search (0, 1)
#' @param n.of.ants The number of artificial ants in the search.
#' @param xi The convergence pressure (0, Inf)
#' @param n.of.archive The number of the solution archive.
#' @return
#'     Generated optimal design parameter value(s) (i.e., a matrix with n.of.ants
#'     rows and n.of.design.pars columns)
#'
#' @export gen.design.pars
#
#' @references
#'
#'   Socha, K., & Dorigo, M. (2008). Ant colony optimization for
#'   continuous domains. European Journal of Operational Research,
#'   185(3), 1155-1173.
#'
#'   We thank Dr. Krzysztof Socha for providing us the
#'   original code (https://iridia.ulb.ac.be/supp/IridiaSupp2008-001/)
#'   for this function.
gen.design.pars <- function(dist.mean, dist.rank, n.of.ants, nl,
                          q = 0.0001, n.of.archive = 100, xi = 0.50) {
  euc.dist <- function(d) { # Euclidean distance
    return(sqrt(sum(d^2)))
  }
  X <- array(dim = c(n.of.ants, dim(dist.mean)[2]))
  idx <- sample(dim(dist.mean)[1], size = n.of.ants,
    replace = TRUE, prob = stats::dnorm(dist.rank, 1, q*n.of.archive))

  # iterate through the chosen distributions
  for (l in 1:length(idx)) {
    j <- idx[l]
    # rotate the coordinate system
    o.dist.mean <- t(t(dist.mean) - dist.mean[j,])  # translation of origin
    r.dist.mean <- o.dist.mean
    set <- nl[j,]  # set of available neighbors
    vec <- vector()
    for (m in 1:(dim(dist.mean)[2]-1)) {
      dis <- apply(matrix(r.dist.mean[set,m:dim(r.dist.mean)[2]],
        length(set),length(m:dim(r.dist.mean)[2])),1, euc.dist)
      if (sum(dis)==0.0)  return(NULL) # if the distribution have converged
      if (length(set)>1){
        choice <- sample(set,size=1,prob=dis^4)
      } else {
        choice <- set
      }
      vec <- cbind(vec,o.dist.mean[choice,])
      R <- qr.Q(qr(vec), complete=TRUE) # rot. matrix after orthogonalization
      if (det(R)<0) {
        R[,1] <- -R[,1]
      }
      r.dist.mean <-  o.dist.mean %*% R # rotated coordinates
      set <- set[set!=choice]
    }

    dist.sd <- vector()
    for (i in 1:dim(dist.mean)[2]) {
      dist.sd <- c(dist.sd,sum(abs(r.dist.mean[nl[j,],i]-
        r.dist.mean[j,i]))/(n.of.archive-1))
    }
    n.x <- stats::rnorm(dim(dist.mean)[2],r.dist.mean[j,],dist.sd*xi)
    n.x <- R %*% n.x
    n.x <- t(n.x + dist.mean[j,])
    X[l,] <- n.x
  }
  return(X)
}

