#' Function to simulate the protracted speciation process
#'
#' Simulating the protracted speciation process using the Doob-Gillespie
#' algorithm. This function differs from pbd_sim_cpp that 1) it does not
#' require that the speciation-initiation rate is the same for good and
#' incipient species, and 2) that it simulates the exact protracted speciation
#' process, and not the approximation made by the coalescent point process.
#' This function provides also the conversion to the approximation as output.
#'
#'
#' @param pars Vector of parameters: \cr \cr \code{pars[1]} corresponds to b_1,
#' the speciation-initiation rate of good species \cr \code{pars[2]}
#' corresponds to la_1, the speciation-completion rate \cr \code{pars[3]}
#' corresponds to b_2, the speciation-initiation rate of incipient species \cr
#' \code{pars[4]} corresponds to mu_1, the extinction rate of good species \cr
#' \code{pars[5]} corresponds to mu_2, the extinction rate of incipient species
#' \cr
#' @param age Sets the age for the simulation
#' @return \item{out}{ A list with the following elements: \cr \cr \code{tree}
#' is the tree of extant species in phylo format \cr \code{stree_random} is a
#' tree with one random sample per species in phylo format \cr
#' \code{stree_oldest} is a tree with the oldest sample per species in phylo
#' format \cr \code{stree_youngest} is a tree with the youngest sample per
#' species in phylo format \cr \code{L} is a matrix of all events in the
#' simulation where \cr - the first column is the incipient-level label of a
#' species \cr - the second column is the incipient-level label of the parent
#' of the species \cr - the third column is the time at which a species is born
#' as incipient species\cr - the fourth column is the time of
#' speciation-completion of the species \cr If the fourth element equals -1,
#' then the species is still incipient.  - the fifth column is the time of
#' extinction of the species \cr If the fifth element equals -1, then the
#' species is still extant.  - The sixth column is the species-level label of
#' the species \cr \code{sL_random} is a matrix like L but for
#' \code{stree_random} \cr \code{sL_oldest} is a matrix like L but for
#' \code{stree_oldest} \cr \code{sL_youngest} is a matrix like L but for
#' \code{stree_youngest} \cr \code{igtree.extinct} is the tree in simmap format
#' with incipient and good flags and including extinct species \cr
#' \code{igtree.extant} is the tree in simmap format with incipient and good
#' flags without extinct species \cr \code{recontree} is the reconstructed tree
#' in phylo format, reconstructed using the approximation in Lambert et al.
#' 2014 \cr \code{reconL} is the matrix corresponding to \code{recontree} \cr
#' \code{L0} is a matrix where the crown age is at 0; for internal use only \cr
#' }
#' @author Rampal S. Etienne
#' @keywords models
#' @export
pbd_sim_optim = function(pars,
                         age,
                         max_n = 1e6) {

  local_sample <- function(v) {
    if (length(v) == 1) {
      return(v[1])
    } else {
      return(sample(x = v, size = 1))
    }
  }


  la1 = pars[1]
  la2 = pars[2]
  la3 = pars[3]
  mu1 = pars[4]
  mu2 = pars[5]

  i <- 1
  L1 <- c()
  L2 <- c()
  while (i <= 2) {
    t <- 0
    L <- list()
    if (i == 1) {
      id1 <- 0
      id <- id1 + 1
      Sid1 <- 0
      Sid <- 1
      sg <- id
      si <- NULL
      L[[1]] <- t(as.matrix(c(id, 0, -1e-10, t, -1, 1)))
    }
    if (i == 2) {
      id = id1 + 1
      Sid = Sid1
      sg = NULL
      si = -id
      L[[1]] <- t(as.matrix(c(id, 1, t, -1, -1, 1)))
    }

    Ng = length(sg)
    Ni = length(si)
    probs = c(la1 * Ng, mu1 * Ng, la2 * Ni, la3 * Ni, mu2 * Ni)
    denom = sum(probs)
    probs = probs/denom
    t = t - log(stats::runif(1)) / denom

    while (t <= age) {
      event = sample(1:5,size = 1,prob = probs)
      if (event == 1) {
        parent = local_sample(sg)
        id = id + 1
        new_index <- length(L) + 1
        L[[new_index]] <- c(id,parent,t,-1,-1, L[[abs(parent) - id1]][6])
        si = c(si, -id)
        Ni = Ni + 1
      } else if (event == 2) {
        todie_index <- local_sample(1:length(sg))
        todie <- sg[todie_index]

        L[[todie - id1]][5] <- t
        sg = sg[-todie_index]
        Ng = Ng - 1
      } else if (event == 3) {

        tocomplete_index <- local_sample(1:length(si))
        tocomplete <- abs(si[tocomplete_index])
        L[[tocomplete - id1]][4] <- t
        Sid = Sid + 1
        L[[tocomplete - id1]][6] <- Sid
        sg = c(sg, tocomplete)
        si = si[-tocomplete_index]
        Ng = Ng + 1
        Ni = Ni - 1
      } else if (event == 4) {
        parent = local_sample(si)
        id = id + 1
        new_index <- length(L) + 1
        L[[new_index]] <- c(id,parent,t,-1,-1,L[[abs(parent) - id1]][6])
        si = c(si, -id)
        Ni = Ni + 1
      } else if (event == 5) {
        todie_index <- local_sample(1:length(si))
        todie <- abs(si[todie_index])

        L[[todie - id1]][5] <- t
        si = si[-todie_index]
        Ni = Ni - 1
      }

      probs = c(la1*Ng,mu1*Ng,la2*Ni,la3*Ni,mu2*Ni)
      denom = sum(probs)
      probs = probs/denom
      t = t - log(stats::runif(1)) / denom
    }

    L <- do.call(rbind, L)

    if (i == 1) {
      if ((Ng + Ni) > 0) {
        i = i + 1
        L1 = L
        id1 = id
        Sid1 = Sid
      }
    } else {
      if (i == 2) {
        if (checkgood(L, si, sg) == 1) {
          i = i + 1
          L2 = L
        }
      }
    }
  }

  L = rbind(L1, L2)

  absL = L
  absL[, 2] = abs(L[, 2])

  # Random sampling
  sL_random = sampletree(absL,
                         age,
                         samplemethod = "random")
  stree_random = ape::as.phylo(ape::read.tree(text =
                                                detphy(sL_random, age)))

  return(stree_random)
}
