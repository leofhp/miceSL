#' Generate random draws from the Dirichlet distribution.
#'
#' @param n Number of random draws to generate.
#' @param alpha Vector of shape parameters.
#'
#' @return A numeric matrix with `n` rows and `length(alpha)` columns, where each row is a Dirichlet random draw.
#'
#' @examples
#' set.seed(42)
#' n <- 10
#' alpha <- c(0.1, 0.4, 0.05, 0.3, 0.15)
#' rdirichlet(n, alpha)
#'
#' @importFrom stats rgamma
#'
#' @export

rdirichlet <- function(n, alpha) {
  if (any(alpha < 0)) stop("'alpha' must not contain negative values")
  k <- length(alpha)
  x <- matrix(rgamma(n * k, shape = alpha), nrow = n, byrow = TRUE)
  return(x / rowSums(x))
}