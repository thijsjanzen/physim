#' analagous function to the DDD function, with intermediate size checking
#' @param pars pars
#' @param age age
#' @param num_lin num_lin
#' @return list
#' @export
dd_sim_R <- function(pars, age, num_lin = -1) {

  L <- c()
  la = pars[1]
  mu = pars[2]
  K = pars[3]
  while (TRUE) {
    done = 0
    while (done == 0) {
      t = rep(0, 1)
      L = matrix(0, 2, 4)
      i = 1
      t[1] = 0
      N = 2
      L[1, 1:4] = c(0, 0, -1, -1)
      L[2, 1:4] = c(0, -1, 2, -1)
      linlist = c(-1, 2)
      newL = 2
      laN = max(0,la - (la - mu) * N[i] / K)
      muN = mu
      denom = (laN + muN) * N[i]
      t[i + 1] = t[i] + stats::rexp(1, denom)
      while (t[i + 1] <= age) {
        i = i + 1
        ranL = DDD::sample2(linlist, 1)
        if ((laN * N[i - 1]/denom) >= stats::runif(1)) {
          N[i] = N[i - 1] + 1
          newL = newL + 1
          L = rbind(L, c(t[i], ranL, sign(ranL) * newL,
                         -1))
          linlist = c(linlist, sign(ranL) * newL)
        }
        else {
          N[i] = N[i - 1] - 1
          L[abs(ranL), 4] = t[i]
          w = which(linlist == ranL)
          linlist = linlist[-w]
          linlist = sort(linlist)
        }
        if (sum(linlist < 0) == 0 | sum(linlist > 0) ==
            0) {
          t[i + 1] = Inf
        }
        else {
          laN = max(0,la - (la - mu) * N[i] / K)
          muN = mu
          denom = (laN + muN) * N[i]
          if (denom == 0) {
            t[i + 1] = Inf
          }
          else {
            t[i + 1] = t[i] + stats::rexp(1, rate = denom)
          }
        }
      }
      if (sum(linlist < 0) == 0 | sum(linlist > 0) == 0) {
        done = 0
      }
      else {
        done = 1
      }
    }
    L[, 1] = age - c(L[, 1])
    notmin1 = which(L[, 4] != -1)
    L[notmin1, 4] = age - c(L[notmin1, 4])
    L[which(L[, 4] == age + 1), 4] = -1

    num_extant_lin <- sum(L[, 4] == -1)
    if (num_extant_lin == num_lin) {
      break
    } else {
  #    cat(num_extant_lin, "\n")
    }
    if (num_lin < 0) break
  }
  tes = treestats::l_to_phylo(L, drop_extinct = TRUE)
  tas = treestats::l_to_phylo(L, drop_extinct = FALSE)
  brts = DDD::L2brts(L, dropextinct = T)
  out = list(tes = tes, tas = tas, L = L, brts = brts)
  return(out)

}
