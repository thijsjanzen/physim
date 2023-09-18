context("ddd")

test_that("usage", {
  set.seed(42)

  lambda <- 0.5
  mu <- 0.0
  K <- 100

  found_sizes <- c()
  for (r in 1:100) {
    focal_tree <- DDD::dd_sim(pars = c(lambda, mu, K),
                              age = 7)
    found_sizes[r] <- length(focal_tree$tes$tip.label)
  }

  found_sizes2 <- c()
  for (r in 1:100) {
    focal_tree <- physim::sim_ddd(lambda = lambda,
                                  mu = mu,
                                  K = K,
                                  max_t = 7)
    found_sizes2[r] <- length(focal_tree$tip.label)
  }

  ax <- t.test(found_sizes, found_sizes2)
  testthat::expect_gt(ax$p.value, 0.05)

  found_sizes <- c()
  for (r in 1:100) {
    focal_tree <- DDD::dd_sim(pars = c(lambda * 2, mu, K),
                              age = 7)
    found_sizes[r] <- length(focal_tree$tes$tip.label)
  }

  found_sizes2 <- c()
  for (r in 1:100) {
    focal_tree <- physim::sim_ddd(lambda = lambda * 2,
                                  mu = mu,
                                  K = K,
                                  max_t = 7)
    found_sizes2[r] <- length(focal_tree$tip.label)
  }
  ax <- t.test(found_sizes, found_sizes2)
  testthat::expect_gt(ax$p.value, 0.05)

  # hardcore testing:
  testthat::skip_on_ci()
  testthat::skip_on_cran()
  testthat::skip_on_covr()

  lambda <- 0.5
  mu <- 0.1
  K <- 300 * 1.2
  max_t <- 10

  found <- c()
  num_repl <- 100
  pb <- txtProgressBar(max = num_repl, style = 3)
  for (r in 1:num_repl) {

    tree1 <- physim::sim_ddd(lambda, mu, K, max_t)
    tree2 <- DDD::dd_sim(pars = c(lambda, mu, K), age = max_t)$tes

    s1 <- treestats::calc_all_stats(tree1, normalize = TRUE)
    s2 <- treestats::calc_all_stats(tree2, normalize = TRUE)

    to_add1 <- c(1, as.numeric(unlist(s1)))
    to_add2 <- c(0, as.numeric(unlist(s2)))
    found <- rbind(found, to_add1, to_add2)


    setTxtProgressBar(pb, r)
  }

  colnames(found) <- c("model", names(s1))
  library(tidyverse)
  found <- tibble::as_tibble(found)
  found$model[found$model == 1] <- "physim"
  found$model[found$model == 0] <- "ddd"

  df <- found[, -1]

  vx <- prcomp(df, scale = TRUE)
  for (i in 1:10) {
    a1 <- vx$x[found$model == "ddd", i]
    a2 <- vx$x[found$model == "physim", i]
    ax <- t.test(a1, a2)
    testthat::expect_gt(ax$p.value, 0.05)
  }

})
