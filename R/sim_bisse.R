#' function so simulate bisse model
#' @param pars parameters as vector {l0, l1, m0, m1, q01, q10}
#' @param crown_age crown age
#' @param num_species maximum number of species allowed
#' @param init_state initial state, if set to -1, 0 or 1 is randomly drawn
#' @param verbose generate verbose output
#' @param max_tries maximum number of tries to generate tree
#' @param cond_num_species should simulations be conditions on the number of
#' species?
#' @param drop_extinct should extinct species be dropped from the tree?
#' @export
#' @return list
sim_bisse <- function(pars,
                      crown_age = 10,
                      num_species = 300,
                      init_state = -1,
                      verbose = FALSE,
                      max_tries = 1e6,
                      cond_num_species = TRUE,
                      drop_extinct = TRUE) {
  stop("This function is not correct\n")
  res <- bisse_sim_cpp(pars,
                       crown_age,
                       num_species,
                       init_state,
                       verbose,
                       max_tries,
                       cond_num_species)

  if (length(res$traits) < 4) {
    warning("tree went extinct")
    return(list(phy = "ds",
                traits = 0))
  }

  Ltable        <- res$ltable

  initialState  <- res$initial_state
  Ltable[, 1]   <- crown_age - Ltable[, 1] # simulation starts at 0,
  # not at crown age
  notmin1 = which(Ltable[, 4] != -1)
  Ltable[notmin1, 4] = crown_age - c(Ltable[notmin1, 4])
  Ltable[which(Ltable[, 4] == crown_age + 1), 4] = -1

  indices       <- seq(1, length(res$traits), by = 2)
  speciesTraits <- 1 + res$traits[indices]

  phy <- treestats::l_to_phylo(Ltable, drop_extinct = drop_extinct)

  return(list(phy = phy,
              traits = speciesTraits,
              initialState = initialState))
}
