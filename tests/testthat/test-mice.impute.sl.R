skip_if_not_installed("mice")

#------------------------------- Data Generation -------------------------------
# Generate a small dataframe with MAR missing values in several
# variables on continuous and binary scales.
set.seed(42)
N <- 100
p <- 0.2
x1 <- rnorm(n = N, mean = 3, sd = 1)
x2 <- rnorm(n = N, mean = -1, sd = 0.5)
x3 <- rbinom(n = N, size = 1, prob = 0.3)
y <- 10 + 2 * sin(x1^2) - 0.05 * exp(x2) - 0.5 * x1 +
  rnorm(n = N, mean = 0, sd = 0.5)
x3[(runif(N) < (2 * p * pnorm(x2, mean = mean(x2), sd = sd(x2))))] <- NA
x1[(runif(N) < (2 * p * pnorm(y, mean = mean(y), sd = sd(y))))] <- NA
x2[(runif(N) > (2 * (1 - p) * pnorm(y, mean = mean(y), sd = sd(y))))] <- NA
df <- data.frame(x1, x2, x3, y)

#---------------------------------- Run MICE -----------------------------------
test_that("MICE runs with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = FALSE, bandwidth = 1, cv_folds = 3
    )
  )
  expect_true(mice::is.mids(imps))
})
test_that("MICE runs with sl3", {
  skip_if_not_installed("sl3")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = FALSE, bandwidth = 1, cv_folds = 3, sl_engine = "sl3"
    )
  )
  expect_true(mice::is.mids(imps))
})

#-------------------------------- Binary Factor --------------------------------
# Encode the binary variable as a factor:
x3 <- factor(x3)
df <- data.frame(x1, x2, x3, y)

test_that("Binary factor errors with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  expect_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = FALSE, bandwidth = 1, cv_folds = 3
    )
  )
})
test_that("Binary factor runs with sl3", {
  skip_if_not_installed("sl3")
  expect_no_error(
    suppressWarnings(
      imps <- mice::mice(
        df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
        bootstrap = FALSE, bandwidth = 1, cv_folds = 3, sl_engine = "sl3"
      )
    )
  )
  expect_true(mice::is.mids(imps))
})

#------------------------------ Categorical Factor -----------------------------
# Add a categorical variable with missing data
x4 <- factor(
  sample(c(0, 1, 2), size = N, replace = TRUE, prob = c(0.30, 0.35, 0.35))
)
x4[x2 > -0.5] <- NA
df <- data.frame(x1, x2, x4, y)

test_that("Categorical factor errors with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  expect_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = FALSE, bandwidth = 1, cv_folds = 3
    )
  )
})
test_that("Categorical factor runs with sl3", {
  skip_if_not_installed("sl3")
  expect_no_error(
    suppressWarnings(
      imps <- mice::mice(
        df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
        bootstrap = FALSE, bandwidth = 1, cv_folds = 3, sl_engine = "sl3"
      )
    )
  )
})

#---------------------- Learner set with only one learner ----------------------
test_that("Learner set with only one learner errors", {
  expect_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      continuous_learners = c("SL.glm"), bootstrap = FALSE, bandwidth = 1
    )
  )
})

#------------------------------- Only one cv fold ------------------------------
test_that("Setting only one cv fold errors", {
  expect_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = FALSE, bandwidth = 1, cv_folds = 1
    )
  )
})

#----------------------------- Different strategies ----------------------------
df <- data.frame(x1, x2, y)

# Localkernel with bootstrap
test_that("Localkernel + bootstrap runs with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = TRUE, bandwidth = 1, cv_folds = 3
    )
  )
  expect_true(mice::is.mids(imps))
})
test_that("Localkernel + bootstrap runs with sl3", {
  skip_if_not_installed("sl3")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = TRUE, bandwidth = 1, cv_folds = 3, sl_engine = "sl3"
    )
  )
  expect_true(mice::is.mids(imps))
})

# Matching without bootstrap
test_that("Matching w/o bootstrap runs with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "matching",
      bootstrap = FALSE, cv_folds = 3
    )
  )
  expect_true(mice::is.mids(imps))
})
test_that("Matching w/o bootstrap runs with sl3", {
  skip_if_not_installed("sl3")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "matching",
      bootstrap = FALSE, cv_folds = 3, sl_engine = "sl3"
    )
  )
  expect_true(mice::is.mids(imps))
})


# Matching with bootstrap
test_that("Matching + bootstrap runs with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "matching",
      bootstrap = TRUE, cv_folds = 3
    )
  )
  expect_true(mice::is.mids(imps))
})
test_that("Matching + bootstrap runs with sl3", {
  skip_if_not_installed("sl3")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "matching",
      bootstrap = TRUE, cv_folds = 3, sl_engine = "sl3"
    )
  )
  expect_true(mice::is.mids(imps))
})

# Dirichlet without bootstrap
test_that("Dirichlet w/o bootstrap runs with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "dirichlet",
      bootstrap = FALSE, cv_folds = 3, alpha = 10
    )
  )
  expect_true(mice::is.mids(imps))
})
test_that("Dirichlet w/o bootstrap runs with sl3", {
  skip_if_not_installed("sl3")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "dirichlet",
      bootstrap = FALSE, cv_folds = 3, alpha = 10, sl_engine = "sl3"
    )
  )
  expect_true(mice::is.mids(imps))
})

# Dirichlet with bootstrap
test_that("Dirichlet + bootstrap runs with SuperLearner", {
  skip_if_not_installed("SuperLearner")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "dirichlet",
      bootstrap = TRUE, cv_folds = 3, alpha = 10
    )
  )
  expect_true(mice::is.mids(imps))
})
test_that("Dirichlet + bootstrap runs with sl3", {
  skip_if_not_installed("sl3")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "dirichlet",
      bootstrap = TRUE, cv_folds = 3, alpha = 10, sl_engine = "sl3"
    )
  )
  expect_true(mice::is.mids(imps))
})

#---------------------------- "localkernel" details ----------------------------
test_that("Non-positive bandwidth throws error", {
  expect_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = FALSE, bandwidth = 0, cv_folds = 3
    )
  )
})

test_that("Vector of candidate bandwidths calls selection procedure", {
  skip_if_not_installed("SuperLearner")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = FALSE, bandwidth = c(1, 2), cv_folds = 3
    )
  )
  expect_true(mice::is.mids(imps))
})
test_that("Vector of candidate bandwidths calls selection procedure", {
  skip_if_not_installed("sl3")
  expect_no_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "localkernel",
      bootstrap = FALSE, bandwidth = c(1, 2), cv_folds = 3, sl_engine = "sl3"
    )
  )
  expect_true(mice::is.mids(imps))
})

#----------------------------- "dirichlet" details -----------------------------
test_that("Negative alpha throws error", {
  expect_error(
    imps <- mice::mice(
      df, m = 2, maxit = 2, method = "sl", strategy = "dirichlet",
      bootstrap = FALSE, alpha = -1, cv_folds = 3
    )
  )
})