#' Super learning methods for \code{mice}.
#'
#' Unified super learning method for multiple imputation by chained equations (MICE), which includes the local imputation approach of Laqueur et al., 2022 (superMICE): \code{strategy = "localkernel"}, \code{bootstrap = FALSE} and the matching approach of Carpenito & Manjourides, 2022 (MISL): \code{strategy = "matching"}, \code{bootstrap = TRUE} as special cases.
#'
#' @param y Vector to be imputed.
#' @param ry Logical vector of length \code{length(y)} indicating the subset \code{y[ry]} of elements in \code{y} to which the imputation model is fitted. The \code{ry} generally distinguishes the observed (\code{TRUE}) and missing values (\code{FALSE}) in \code{y}.
#' @param x Numeric design matrix with \code{length(y)} rows with predictors for \code{y}. Matrix \code{x} may have no missing values.
#' @param type Integer vector of length \code{length(y)}, as expected by `mice`.
#' @param outcome_type Variable scale of \code{y}. One of \code{c("continuous", "binary", "categorical")} or NULL. In the latter case, \code{y} is assumed to be binary if it is a factor with 2 levels, categorical if is a factor with more than 2 levels, and continuous otherwise.
#' @param continuous_learners Learning algorithms for continuous variables. Any learning algorithms available in the \code{SuperLearner} or \code{sl3} packages may be specified, depending on \code{sl_engine}. Defaults to mean and glm.
#' @param binary_learners Learning algorithms for binary variables. Any learning algorithms available in the \code{SuperLearner} or \code{sl3} packages may be specified, depending on \code{sl_engine}. Defaults to mean and glm.
#' @param categorical_learners Learning algorithms for categorical variables. Any learning algorithms available in the \code{sl3} packages may be specified. Defaults to \code{c("Lrnr_mean", "Lrnr_glmnet")}.
#' @param strategy Imputation strategy for continuous variables. One of \code{c("localkernel", "matching", "dirichlet")}. Defaults to "localkernel", which corresponds to a local imputation strategy with kernel-based variance estimation, introduced as superMICE by Laqueur et al. 2022. 
#' The "matching" strategy applies predictive mean matching to the super learner predictions, introduced as MISL by Carpenito and Manjourides 2022. The "dirichlet" strategy re-weights the super learner ensemble by sampling from a Dirichlet distribution in each iteration, which we term MIDE (Multiple Imputation by Dirichlet-weighted super learning Ensembles).
#' @param bootstrap Logical value indicating whether to train the super learner on a bootstrap sample or the full data in each iteration. Bootstrapping approximates a proper imputation procedure. Defaults to FALSE.
#' @param sl_engine Super learner backend for model training. One of \code{c("SuperLearner", "sl3")}. Must be \code{"sl3"} for categorical \code{outcome_type}. Defaults to \code{"SuperLearner"}.
#' @param cv_folds Number of cross-validation folds to use in super learner training. Defaults to 10.
#' @param kernel If \code{strategy = "localkernel"}, a string specifying the kernel to use in local variance estimation. One of \code{c("gaussian", "uniform", "triangular")}. Defaults to "gaussian".
#' @param bandwidth If \code{strategy = "localkernel"}, a positive numeric value or vector specifying the kernel bandwidth. If a vector of candidate bandwidths is supplied, the jackknife selection procedure of Laqueur et al., 2022 is called. Defaults to 1.
#' @param alpha If \code{strategy = "dirichlet"}, a positive numeric value for the scaling parameter used in sampling weights from the Dirichlet distribution. Defaults to 1.
#' @param ... Additional arguments passed to \code{SuperLearner}, such as learner-specific tuning parameters. Ignored if \code{sl_engine = "sl3"}.
#'
#' @return An imputed vector of `y` values (as expected by `mice`).
#'
#' @references
#' Laqueur HS, Shev AB, Kagawa RMC. SuperMICE: An Ensemble Machine Learning Approach to Multiple Imputation by Chained Equations. American Journal of Epidemiology. 2022;191(3):516-525.
#' @references
#' Carpenito T, Manjourides J. MISL: Multiple imputation by super learning. Statistical Methods in Medical Research. 2022;31(10):1904-1915.
#'
#' @examples
#' # Bivariate data (x, y) with missing values in x based on an MAR mechanism.
#' set.seed(42)
#' N <- 50
#' p <- 0.2
#' x <- rnorm(n = N, mean = 3, sd = 1)
#' y <- 10 + 2 * sin(x^2) - 0.05 * exp(x) - 0.5 * x +
#'   rnorm(n = N, mean = 0, sd = 0.5)
#' x[(runif(N) < (2 * p * pnorm(y, mean = mean(y), sd = sd(y))))] <- NA
#' df <- data.frame(x, y)
#'
#' imps <- mice::mice(
#'   df, m = 5, maxit = 1, method = "sl", strategy = "localkernel",
#'   bootstrap = FALSE, bandwidth = 0.5, cv_folds = 5,
#'   continuous_learners = c("SL.mean", "SL.glm")
#' )
#'
#' @importFrom stats rbinom rnorm runif
#' @importFrom mice mice
#' @import SuperLearner
#'
#' @export

mice.impute.sl <- function(
  y, ry, x, type, outcome_type = NULL,
  continuous_learners = NULL,
  binary_learners = NULL,
  categorical_learners = NULL,
  strategy = c("localkernel", "matching", "dirichlet"),
  bootstrap = FALSE,
  sl_engine = c("SuperLearner", "sl3"),
  cv_folds = 10,
  kernel = c("gaussian", "uniform", "triangular"),
  bandwidth = 1,
  alpha = 1,
  ...
) {
  sl_engine <- match.arg(sl_engine)
  strategy <- match.arg(strategy)
  if (strategy == "localkernel") {
    kernel <- match.arg(kernel)
    if (any(bandwidth <= 0)) {
      stop("'bandwidth' may only contain positive values.")
    }
  }
  if (!is.numeric(cv_folds) || cv_folds < 2) {
    stop("The number of cross-validation folds ('cv_folds') must be at least 2.")
  }

  dots <- list(...)
  if ("wy" %in% names(dots)) dots$wy <- NULL  # remove wy from '...' arguments

  if (is.null(continuous_learners)) {
    continuous_learners <- if (sl_engine == "SuperLearner") c("SL.mean", "SL.glm") else c("Lrnr_mean", "Lrnr_glm")
  } else if (!is.character(continuous_learners) || length(continuous_learners) < 2) {
      stop("'continuous_learners' must be a character vector with at least two learners.")
  }
  if (is.null(binary_learners)) {
    binary_learners <- if (sl_engine == "SuperLearner") c("SL.mean", "SL.glm") else c("Lrnr_mean", "Lrnr_glm")
  } else if (!is.character(binary_learners) || length(binary_learners) < 2) {
    stop("'binary_learners' must be a character vector with at least two learners.")
  }
  if (is.null(categorical_learners)) {
    categorical_learners <- c("Lrnr_mean", "Lrnr_glmnet")
  } else if (!is.character(categorical_learners) || length(categorical_learners) < 2) {
    stop("'categorical_learners' must be a character vector with at least two learners.")
  }

  if (sl_engine == "SuperLearner" && !requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("The 'SuperLearner' package is required. Please install it.")
  }
  if (sl_engine == "sl3" && !requireNamespace("sl3", quietly = TRUE)) {
    stop("The 'sl3' package is required. Please install it.")
  }

  if (is.null(outcome_type)) {
    if (is.factor(y) && nlevels(y) > 2L) {
      outcome_type <- "categorical"
    } else if ((is.factor(y) && nlevels(y) == 2L) ||
                 length(unique(y[!is.na(y)])) == 2) {
      outcome_type <- "binary"
    } else {
      outcome_type <- "continuous"
    }
  }

  if (outcome_type == "binary" && sl_engine == "SuperLearner") {
    stop_message <- "Binary variables must be coded as numeric 0/1 when
      using `sl_engine = 'SuperLearner'`."
    if (is.factor(y)) {
      stop(stop_message)
    } else {
      uniq <- sort(unique(y[!is.na(y)]))
      if (!all(uniq %in% c(0, 1))) {
        stop(stop_message)
      }
    }
  }

  if (bootstrap) {
    boot_idx <- sample(which(ry), replace = TRUE)
    x_boot <- x[boot_idx, , drop = FALSE]
    y_boot <- y[boot_idx]
  } else {
    x_boot <- x[ry, , drop = FALSE]
    y_boot <- y[ry]
  }

  n_mis <- sum(!ry)

  if (sl_engine == "sl3") {
    # sl3 path
    x_df <- as.data.frame(x_boot)
    x_df[[".y"]] <- y_boot

    task <- sl3::sl3_Task$new(data = x_df,
                              covariates = colnames(x),
                              outcome = ".y",
                              outcome_type = outcome_type,
                              folds = cv_folds)

    learners <- switch(outcome_type,
                       "continuous" = continuous_learners,
                       "binary" = binary_learners,
                       "categorical" = categorical_learners)

    learner_list <- lapply(learners, function(lrn) {
      eval(parse(text = paste0("sl3::", lrn, "$new()")))
    })
    stack <- sl3::make_learner(sl3::Stack, learner_list)
    sl <- sl3::Lrnr_sl$new(learners = stack)

    sl_fit <- sl$train(task)

    # Predictions for all rows
    x_all_df <- as.data.frame(x)
    x_all_df[[".y"]] <- y
    task_all <- sl3::sl3_Task$new(data = x_all_df,
                                  covariates = colnames(x),
                                  outcome = ".y",
                                  outcome_type = outcome_type)
    pred <- sl_fit$predict(task_all)

    if (outcome_type == "continuous") {
      if (strategy == "localkernel") {
        delta <- as.numeric(ry)
        missing_indices <- which(!ry)

        if (length(bandwidth) == 1) {
          bw_list <- as.list(rep(bandwidth, times = sum(!ry)))
        } else {
          bw_vec <- sapply(missing_indices, select_bandwidth,
                           bw_grid = bandwidth, preds = pred,
                           y = y, delta = delta,
                           kernel = kernel)
          bw_list <- as.list(bw_vec)
        }

        imp <- sapply(seq_along(missing_indices), local_imputation,
                      preds = pred, y = y, delta = delta,
                      bw = bw_list, kernel = kernel)

      } else if (strategy == "matching") {
        if (bootstrap) {
          # Retrain on all current rows for matching
          x_obs_df <- as.data.frame(x[ry, , drop = FALSE])
          x_obs_df[[".y"]] <- y[ry]

          task_obs <- sl3::sl3_Task$new(data = x_obs_df,
                                        covariates = colnames(x),
                                        outcome = ".y",
                                        outcome_type = outcome_type,
                                        folds = cv_folds)
          sl_fit_obs <- sl$train(task_obs)
          pred_obs <- sl_fit_obs$predict(task_all)[ry]
        } else {
          pred_obs <- pred[ry]
        }

        pred_mis <- pred[!ry]

        imp <- numeric(n_mis)
        for (i in seq_along(pred_mis)) {
          distances <- abs(pred_mis[i] - pred_obs)
          closest <- order(distances)[1:5]
          imp[i] <- sample(y[ry][closest], 1)
        }
      } else if (strategy == "dirichlet") {
        base_learners <- sl_fit$learner_fits
        base_preds <- lapply(base_learners, function(lrnr) {
          lrnr$predict(task_all)[!ry]
        })
        base_preds_mat <- do.call(cbind, base_preds)

        lambda <- sl_fit$metalearner_fit()$coefficients
        lambda[is.na(lambda) | lambda < 0] <- 0
        weights <- rdirichlet(n = n_mis, alpha = alpha * lambda)
        imp <- rowSums(base_preds_mat * weights, na.rm = TRUE)
      }

    } else if (outcome_type == "binary") {
      imp <- rbinom(n_mis, 1, pred[!ry])

    } else if (outcome_type == "categorical") {
      probs <- sl3::unpack_predictions(pred[!ry])
      u <- matrix(runif(nrow(probs) * (ncol(probs) - 1)),
                  nrow = nrow(probs), ncol = ncol(probs) - 1)
      draws <- u > t(apply(probs[, -ncol(probs), drop = FALSE], 1, cumsum))
      idx <- 1 + rowSums(draws)
      levels_y <- levels(y[ry])
      imp <- factor(levels_y[idx], levels = levels_y)
    }

  } else {
    # SuperLearner path
    if (outcome_type == "categorical") {
      stop("Outcome type 'categorical' is currently only supported with sl_engine = 'sl3'.")
    }

    family <- if (outcome_type == "continuous") stats::gaussian() else stats::binomial()

    learners <- switch(outcome_type,
                       "continuous" = continuous_learners,
                       "binary" = binary_learners)

    base_args <- list(
      Y = y_boot, X = as.data.frame(x_boot), family = family,
      SL.library = learners, cvControl = list(V = cv_folds)
    )

    sl_fit <- do.call(SuperLearner::SuperLearner, c(base_args, dots))

    # Predictions for all rows
    pred_obj <- SuperLearner::predict.SuperLearner(
      sl_fit, newdata = as.data.frame(x), onlySL = FALSE
    )
    pred <- as.vector(pred_obj$pred)

    if (outcome_type == "continuous") {
      if (strategy == "localkernel") {
        delta <- as.numeric(ry)
        missing_indices <- which(!ry)

        if (length(bandwidth) == 1) {
          bw_list <- as.list(rep(bandwidth, times = sum(!ry)))
        } else {
          bw_vec <- sapply(missing_indices, select_bandwidth,
                           bw_grid = bandwidth, preds = pred,
                           y = y, delta = delta,
                           kernel = kernel)
          bw_list <- as.list(bw_vec)
        }

        imp <- sapply(seq_along(missing_indices), local_imputation,
                      preds = pred, y = y, delta = delta,
                      bw = bw_list, kernel = kernel)

      } else if (strategy == "matching") {
        if (bootstrap) {
          # Retrain on all current rows for matching
          base_args_obs <- list(
            Y = y[ry], X = as.data.frame(x[ry, , drop = FALSE]),
            family = stats::gaussian(), SL.library = continuous_learners,
            cvControl = list(V = cv_folds)
          )
          sl_fit_obs <- do.call(SuperLearner::SuperLearner, c(base_args_obs, dots))

          pred_obs <- as.vector(SuperLearner::predict.SuperLearner(
            sl_fit_obs, newdata = as.data.frame(x[ry, , drop = FALSE]), onlySL = TRUE
          )$pred)
        } else {
          pred_obs <- pred[ry]
        }

        pred_mis <- pred[!ry]

        imp <- numeric(n_mis)
        for (i in seq_along(pred_mis)) {
          distances <- abs(pred_mis[i] - pred_obs)
          closest <- order(distances)[1:5]
          imp[i] <- sample(y[ry][closest], 1)
        }

      } else if (strategy == "dirichlet") {
        base_preds_mat <- pred_obj$library.predict[!ry, , drop = FALSE]
        lambda <- sl_fit$coef[colnames(base_preds_mat)]
        lambda[is.na(lambda) | lambda < 0] <- 0
        weights <- rdirichlet(n = n_mis, alpha = alpha * lambda)
        imp <- rowSums(base_preds_mat * weights, na.rm = TRUE)
      }

    } else if (outcome_type == "binary") {
      imp <- rbinom(n_mis, 1, pred[!ry])
    }
  }

  return(imp)
}