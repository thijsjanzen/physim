#' Function to simulate the protracted speciation process
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b_1,
#' the speciation-initiation rate of good species \cr \code{pars[2]}
#' corresponds to la_1, the speciation-completion rate \cr \code{pars[3]}
#' corresponds to b_2, the speciation-initiation rate of incipient species \cr
#' \code{pars[4]} corresponds to mu_1, the extinction rate of good species \cr
#' \code{pars[5]} corresponds to mu_2, the extinction rate of incipient species
#' \cr
#' @param age Sets the crown age for the simulation
#' @param upper_species_limit stops simulation if the number of good + the
#' number of incipient species exceeds this number
#' @param num_tries number of times to re-start the simulation in case of
#' extinction (all species dead) or overflow (if upper_species_limit is
#' exceeded)
#' @return phy object in case of success, or a string ("extinction" or
#' "overflow",
#' if no tree could be generated within num_tries tries).
#' @author Thijs Janzen
#' @keywords models
#' @export
pbd_sim_rcpp = function(pars,
                        age,
                        upper_species_limit = 100000,
                        num_tries = 100) {
  b1 = pars[1]
  c1 = pars[2]
  b2 = pars[3]
  mu1 = pars[4]
  mu2 = pars[5]

  res <- sim_pbd_cpp(b1,
                     mu1,
                     b2,
                     mu2,
                     c1,
                     age,
                     upper_species_limit,
                     num_tries)

  if (res$status == "success") {
    phy <- try(treestats::l_to_phylo(res$ltable))
    if (class(phy) == "try-error") {
      for (i in 1:nrow(res$ltable)) {
        cat(res$ltable[i, ], "\n")
      }
    }
  } else {
    phy <- res$status
  }

  return(phy)
}



#' Function to simulate the protracted speciation process
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b_1,
#' the speciation-initiation rate of good species \cr \code{pars[2]}
#' corresponds to la_1, the speciation-completion rate \cr \code{pars[3]}
#' corresponds to b_2, the speciation-initiation rate of incipient species \cr
#' \code{pars[4]} corresponds to mu_1, the extinction rate of good species \cr
#' \code{pars[5]} corresponds to mu_2, the extinction rate of incipient species
#' \cr
#' @param age Sets the crown age for the simulation
#' @param min_species minimum number of species in tree
#' @param max_species maximum number of species in tree
#' @param num_tries number of times to re-start the simulation in case of
#' extinction (all species dead) or overflow (if upper_species_limit is
#' exceeded)
#' @return phy object in case of success, or a string ("extinction" or
#' "overflow",
#' if no tree could be generated within num_tries tries).
#' @author Thijs Janzen
#' @keywords models
#' @export
pbd_sim_conditional_rcpp = function(pars,
                        age,
                        min_species = 300,
                        max_species = 300,
                        num_tries = 100) {
  b1 = pars[1]
  c1 = pars[2]
  b2 = pars[3]
  mu1 = pars[4]
  mu2 = pars[5]

  res <- sim_pbd_conditional_cpp(b1,
                                 mu1,
                                 b2,
                                 mu2,
                                 c1,
                                 age,
                                 min_species,
                                 max_species,
                                 num_tries)

  if (res$status == "success") {
    phy <- try(treestats::l_to_phylo(res$ltable))
    if (class(phy) == "try-error") {
      for (i in 1:nrow(res$ltable)) {
        cat(res$ltable[i, ], "\n")
      }
    }
  } else {
    phy <- res$status
  }

  return(phy)
}

