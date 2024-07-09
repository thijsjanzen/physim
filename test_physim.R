lambda0 <- 0.8
compl_rate <- 1e10
lambda1 <- 0.0
mu0 = 0.0
mu1 = 0.0
crown_age <- 5

local_l <- physim::sim_pbd_cpp(lambda0, mu0, lambda1, mu1,
                               compl_rate, crown_age, 1000, 100)



found <- c()

for (r in 1:100) {
  focal_tree <- physim::pbd_sim_rcpp(pars = c(lambda0,
                                              compl_rate,
                                              lambda1,
                                              mu0,
                                              mu1),
                                     age = crown_age,
                                     upper_species_limit = 1000,
                                     num_tries = 1000)
  found[r] <- treestats::number_of_lineages(focal_tree)
}

mean(found)
found2 <- c()
for (r in 1:100) {
  focal_tree <- ape::rbdtree(birth = lambda0, death = mu0, Tmax = crown_age)
  found2[r] <- treestats::number_of_lineages(focal_tree)
}

found3 <- c()
for (r in 1:100) {
  focal_tree <- PBD::pbd_sim(pars = c(lambda0,
                                      compl_rate,
                                      lambda1,
                                      mu0,
                                      mu1),
                             age = 5)
  found3[r] <- treestats::number_of_lineages(focal_tree$stree_random)
}

require(tidyverse)
to_plot <- rbind(cbind("rcpp", found),
                 cbind("ape", found2),
                 cbind("PBD", found3))
colnames(to_plot) <- c("model", "size")
to_plot <- as_tibble(to_plot)
to_plot$size <- as.numeric(to_plot$size)
ggplot(to_plot, aes(x = size, fill = model)) +
  geom_density(alpha = 0.5)
