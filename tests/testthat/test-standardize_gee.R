test_that("check estimates and standard errors against older package (gee)", {
  require(drgee)

  set.seed(4)
  n <- 300
  ni <- 2
  id <- rep(1:n, each = ni)
  ai <- rep(rnorm(n), each = ni)
  Z <- rnorm(n * ni)
  X <- rnorm(n * ni, mean = ai + Z)
  Y <- rnorm(n * ni, mean = ai + X + Z + 0.1 * X^2)
  dd <- data.frame(id, Z, X, Y)
  fit.std <- standardize_gee(
    formula = Y ~ X + Z + I(X^2),
    link = "identity",
    data = dd,
    values = list(X = seq(-3, 3, 0.5)),
    clusterid = "id"
  )
  expect_equal(unname(fit.std$res_contrast[[1]]$est_table[, 2]), c(
    -2.31169220309552, -2.02480660497755, -1.69617517072297, -1.32579790033181,
    -0.913674793804042, -0.459805851139684, 0.0358089276612683, 0.573169542598816,
    1.15227599367296, 1.7731282808837, 2.43572640423103, 3.14007036371496,
    3.88616015933548
  ), tolerance = 1e-5)
  expect_equal(unname(fit.std$res_contrast[[1]]$est_table[, 3]), c(
    0.189499305124431, 0.153122875180863, 0.124592006350899, 0.10447283361008,
    0.0929753745724416, 0.0893586277257601, 0.0919626188863363, 0.0991487714228512,
    0.110105752426333, 0.12490967174844, 0.144164836355114, 0.168614161078697,
    0.198889498011441
  ), tolerance = 1e-5)
})
