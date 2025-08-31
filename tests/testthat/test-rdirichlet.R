test_that("Returned matrix by rdirichlet has correct dimensions", {
  set.seed(42)
  n <- 5
  alpha <- c(1, 2, 3)
  out <- rdirichlet(n, alpha)
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(n, length(alpha)))
})

test_that("Rows of rdirichlet output sum to 1", {
  set.seed(42)
  out <- rdirichlet(10, c(1, 1, 1))
  row_sums <- rowSums(out)
  expect_true(all(abs(row_sums - 1) < 1e-8))
})

test_that("Negative alpha throws an error", {
  expect_error(rdirichlet(5, c(-1, 1, 2)), "alpha")
})