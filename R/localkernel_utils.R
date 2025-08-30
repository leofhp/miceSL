#' Compute unnormalized kernel values.
#'
#' Helper function used internally by \code{mice.impute.sl}, adapted
#' from \code{superMICE::gaussianKernel}, \code{superMICE::uniformKernel}, and
#' \code{superMICE::triangularKernel} by Laqueur et al. 2022, see
#' <https://github.com/abshev/superMICE>.
#'
#' @param x Numeric vector.
#' @param center Numeric value to center the kernel.
#' @param bw Bandwidth of the kernel.
#' @param kernel One of \code{c("gaussian", "uniform", "triangular")}.
#'
#' @return Kernel values for \code{x} centered at \code{center}.
#'
#' @references
#' Laqueur HS, Shev AB, Kagawa RMC. SuperMICE: An Ensemble Machine Learning Approach to Multiple Imputation by Chained Equations. American Journal of Epidemiology. 2022;191(3):516-525.
#'
#' @importFrom stats dnorm
#'
#' @keywords internal
get_kernel_vals <- function(x, center, bw, kernel) {
  if (kernel == "gaussian") {
    lambda <- bw
  } else if (kernel == "uniform") {
    lambda <- bw / sqrt(1 / 3)
  } else {
    lambda <- bw / sqrt(1 / 6)
  }
  z <- (x - center) / lambda
  if (kernel == "gaussian") {
    out <- dnorm(z) / lambda
  } else if (kernel == "uniform") {
    out <- ((abs(z) <= 1) / lambda) * (1 / 2)
  } else {
    out <- pmax(1 - abs(z), 0) / lambda
  }
  return(out)
}

#' Impute values based on a local normal distribution with a
#' kernel-weighted variance estimate.
#'
#' Helper function used internally by \code{mice.impute.sl}, adapted
#' from \code{superMICE::localImputation} by Laqueur et al. 2022, see
#' <https://github.com/abshev/superMICE>.
#'
#' @param i Index of the missing value to be imputed.
#' @param preds Predictions for missing values.
#' @param y Variable to be imputed.
#' @param delta Missingness indicator of \code{y}.
#' @param bw Numeric vector of kernel bandwidths.
#' @param kernel One of \code{c("gaussian", "uniform", "triangular")}.
#'
#' @return Numeric vector of imputed values.
#'
#' @references
#' Laqueur HS, Shev AB, Kagawa RMC. SuperMICE: An Ensemble Machine Learning Approach to Multiple Imputation by Chained Equations. American Journal of Epidemiology. 2022;191(3):516-525.
#'
#' @importFrom stats rnorm
#'
#' @keywords internal
local_imputation <- function(i, preds, y, delta, bw, kernel) {
  kernel_vals <- get_kernel_vals(preds, preds[delta == 0][i], bw[[i]], kernel)
  weights <- kernel_vals / sum(kernel_vals)
  pi_hat <- sum(kernel_vals * delta) / sum(kernel_vals)
  mu_hat <- sum(weights * delta * y / pi_hat)
  var_hat <- sum(weights * delta * y^2 / pi_hat) - mu_hat^2
  return(rnorm(1, preds[delta == 0][i], sqrt(var_hat)))
}

#' Estimate the jackknife variance from a matrix of kernel values.
#'
#' Helper function used internally by \code{mice.impute.sl}, adapted
#' from \code{superMICE::jackknifeVariance} by Laqueur et al. 2022, see
#' <https://github.com/abshev/superMICE>.
#'
#' @param j Index of the deleted observation in the jackknife procedure.
#' @param kernel_mat Matrix of kernel values centered at \code{j} of
#' shape \code{(n-1)} by \code{k}, where \code{n} is the total number of
#' observations and \code{k} is the number of candidate bandwidths.
#' @param delta Missingness indicator of \code{y}.
#' @param y Variable to be imputed.
#'
#' @return Numeric variance estimate.
#'
#' @references
#' Laqueur HS, Shev AB, Kagawa RMC. SuperMICE: An Ensemble Machine Learning Approach to Multiple Imputation by Chained Equations. American Journal of Epidemiology. 2022;191(3):516-525.
#'
#' @keywords internal
jackknife_var <- function(j, kernel_mat, delta, y) {
  n <- nrow(kernel_mat)
  k <- ncol(kernel_mat)
  weights <- kernel_mat[-j, ] / matrix(colSums(kernel_mat[-j, ]),
                                       nrow = n - 1, ncol = k, byrow = TRUE)
  pi_hat <- colSums(kernel_mat[-j, ] * delta[-j]) / colSums(kernel_mat[-j, ])
  mu_hat <- colSums(weights * delta[-j] * y[-j]) / pi_hat
  mu2_hat <- colSums(weights * delta[-j] * y[-j]^2) / pi_hat
  return(mu2_hat - mu_hat^2)
}

#' Select an optimal bandwidth among a set of candidates using a
#' jackknife procedure.
#'
#' Helper function used internally by \code{mice.impute.sl}, adapted
#' from \code{superMICE::jackknifeBandwidthSelection} by Laqueur et al. 2022,
#' see <https://github.com/abshev/superMICE>.
#'
#' @param i Index of the missing value to be imputed.
#' @param bw_grid Candidate bandwidth values.
#' @param preds Predictions for missing values.
#' @param y Variable to be imputed.
#' @param delta Missingness indicator of \code{y}.
#' @param kernel  One of \code{c("gaussian", "uniform", "triangular")}.
#'
#' @return Optimal bandwidth value.
#'
#' @references
#' Laqueur HS, Shev AB, Kagawa RMC. SuperMICE: An Ensemble Machine Learning Approach to Multiple Imputation by Chained Equations. American Journal of Epidemiology. 2022;191(3):516-525.
#'
#' @keywords internal
select_bandwidth <- function(i, bw_grid, preds, y, delta, kernel) {
  kernel_mat <- do.call(cbind, lapply(
    bw_grid, get_kernel_vals, x = preds, center = preds[i], kernel = kernel
  ))
  n <- nrow(kernel_mat)
  k <- ncol(kernel_mat)
  weights <- kernel_mat / matrix(colSums(kernel_mat),
                                 nrow = n, ncol = k, byrow = TRUE)
  pi_hat <- colSums(kernel_mat * delta) / colSums(kernel_mat)
  mu_hat <- colSums(weights * delta * y) / pi_hat
  mu2_hat <- colSums(weights * delta * y^2) / pi_hat
  var_hat <- mu2_hat - mu_hat^2
  var_hat_jk <- do.call(
    rbind, lapply((1:n)[delta == 1], jackknife_var,
                  kernel_mat = kernel_mat, delta = delta, y = y)
  )
  bias2 <- (var_hat - colMeans(var_hat_jk))^2
  s2 <- colSums((var_hat_jk - matrix(
    colMeans(var_hat_jk), nrow = nrow(var_hat_jk), ncol = k, byrow = TRUE
  ))^2) / (n * (n - 1))
  mse <- bias2 + s2
  return(bw_grid[which.min(mse)])
}