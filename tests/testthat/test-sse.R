context("bisse")

test_that("usage", {
  set.seed(42)

  lambda0 <- 0.5
  lambda1 <- 0.7
  mu0 <- 0.0
  mu1 <- 0.0
  q01 <- 0.1
  q10 <- 0.1

  max_t <- 10

  pars <- c(lambda0, lambda1, mu0, mu1, q01, q10)

  found_sizes <- c()
  freq_trait <- c()
  for (r in 1:200) {

    focal_tree <- diversitree::tree.bisse(pars = pars,
                                          max.t = max_t,
                                          x0 = 0,
                                          include.extinct = FALSE)
    found_sizes[r] <- length(focal_tree$tip.label)
    freq_trait[r] <- sum(focal_tree$tip.state) / found_sizes[r]
  }

  found_sizes2 <- c()
  freq_trait2 <- c()
  for (r in 1:1000) {
    focal_tree <- physim::sim_bisse(pars = pars,
                                    crown_age = max_t,
                                    num_species = 100000,
                                    init_state = 0,
                                    cond_num_species = FALSE)
    found_sizes2[r] <- length(focal_tree$phy$tip.label)
    freq_trait2[r] <- sum(focal_tree$traits - 1) / found_sizes2[r]
  }

  ax <- t.test(found_sizes, found_sizes2)
  testthat::expect_gt(ax$p.value, 0.05)

  ax <- t.test(freq_trait, freq_trait2)
  testthat::expect_gt(ax$p.value, 0.05)

  to_plot <- rbind(cbind(found_sizes, "diversitree"),
                   cbind(found_sizes2, "physim"))
  colnames(to_plot) <- c("size", "model")
  to_plot <- as_tibble(to_plot)
  to_plot$size <- as.numeric(to_plot$size)
  ggplot(to_plot, aes(x = size, fill = model)) +
      geom_density(alpha = 0.5)
})
