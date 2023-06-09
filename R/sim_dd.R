#' function to simulate the linear version of the Diversity-Dependent-Diversification
#' model
#' @param lambda speciation rate
#' @param mu extinction rate
#' @param K level of diversity at which the speciation rate equals the extinction
#' rate
#' @param max_t crown age of the tree
#' @param num_species number of species of the tree
#' @param return_fossils should extinct species be returned in the tree?
#' @param seed seed of the random number generator
#' @description
#' Simulations are performed either conditional on crown age, on the number of
#' extant lineages, or on both, depending on which of max_t and num_species has
#' a value larger than zero.
#' @return phylo object
#' @export
#' @rawNamespace useDynLib(physim, .registration = TRUE)
#' @rawNamespace import(Rcpp)
sim_ddd <- function(lambda,
                    mu,
                    K,
                    max_t = -1,
                    num_species = -1,
                    return_fossils = FALSE,
                    seed = -1) {

  if (is.infinite(max_t) &&
      is.infinite(num_species)) {
    stop("simulations have to be conditioned on either time or number
         of lineages")
  }

  res <- sim_ddd_cpp(lambda, mu, K, max_t, num_species, seed)
  ltable <- res$ltable
  crown_age <- res$crown_age
  if (max_t != -1) crown_age <- max_t
  ltable[, 1] <- crown_age - ltable[, 1]
  notmin1 <- which(ltable[, 4] != -1)
  ltable[notmin1, 4] <- crown_age - c(ltable[notmin1, 4])
  ltable[which(ltable[, 4] == crown_age + 1), 4] <- -1

  phy <- treestats::l_to_phylo(ltable, !return_fossils)
  return(phy)
}
