checkgood = function(L,si,sg,id1)
{
  j = 1
  found = 0
  if(length(sg) > 0)
  {
    found = 1
  } else {
    if(length(si) == 0) { found = 0 } else {
      while(found == 0 & j <= length(si))
      {
        rowinc = which(L[,1] == abs(si[j]))
        parent = L[rowinc,2]
        birth = L[rowinc,3]
        #parent = L[si[j] - id1,2]
        #birth = L[si[j] - id1,3]
        while(found == 0 & parent > 1)
        {
          rowpar = which(L[,1] == parent)
          if(L[rowpar,4] > -1 & L[rowpar,4] < birth)
            #if(L[parent - id1,4] > -1 & L[parent - id1,4] < birth)
          {
            found = 1
          } else
          {
            parent = L[rowpar,2]
            birth = L[rowpar,3]
            #parent = L[parent - id1,2]
            #birth = L[parent - id1,3]
          }
        }
        j = j + 1
      }}}
  invisible(found)
}

detphy = function(L, age, ig = F, dropextinct = T)
{
  #print(L)
  dimL = dim(L)
  if ((dimL[1] == 1))
  {
    linlist = paste("(S-1-1-1:",age,");",sep = "")
  } else {
    L = L[order(L[,1]),1:6]
    if (dropextinct == TRUE)
    {
      sall = which(L[,5] == -1)
      tend = age
    } else {
      sall = which(L[,5] >= -1)
      tend = (L[,5] == -1) * age + (L[,5] > -1) * L[,5]
    }

    linlist = matrix(0,nrow = 1,ncol = 8)
    if (length(sall) == 1)
    {
      linlist[1,] = c(L[sall,],paste("S",paste(L[sall,6],L[sall,6],L[sall,1],sep = "-"),sep = ""),tend)
    } else {
      linlist = cbind(L[sall,],paste("S",paste(L[sall,6],L[sall,6],L[sall,1],sep = "-"),sep = ""),tend)
    }
    done = 0
    while(done == 0)
    {
      j = which.max(linlist[,3])
      daughter = as.numeric(linlist[j,1])
      parent = as.numeric(linlist[j,2])
      parentj = which(linlist[,1] == parent)
      parentinlist = length(parentj)

      if (parentinlist == 1) {
        startedge = as.numeric(linlist[j,3])
        comptime = as.numeric(linlist[j,4])
        endedge = as.numeric(linlist[j,8])
        comptimeparent = as.numeric(linlist[parentj,4])
        endedgeparent = as.numeric(linlist[parentj,8])
        if (ig == FALSE) {
          spec1 = paste(linlist[parentj,7],":",endedgeparent - startedge,sep = "")
          spec2 = paste(linlist[j,7],":",endedge - startedge,sep = "")
        } else {
          if (comptimeparent == -1 | comptimeparent > endedgeparent) {
            comptimeparent = endedgeparent
          }
          itimeparent = max(0,comptimeparent - startedge)
          gtimeparent = endedgeparent - max(comptimeparent,startedge)
          if (itimeparent == 0) {
            spec1 = paste(linlist[parentj,7],":{g,",gtimeparent,":i,",itimeparent,":g,0}",sep = "")
          } else {
            spec1 = paste(linlist[parentj,7],":{g,",gtimeparent,":i,",itimeparent,"}",sep = "")
          }
          if (comptime == -1 | comptime > endedge) {
            comptime = endedge
          }
          itime = comptime - startedge
          gtime = endedge - comptime
          if(itimeparent == 0) {
            spec2 = paste(linlist[j,7],":{g,",gtime,":i,",itime,":g,0}",sep = "")
          } else {
            spec2 = paste(linlist[j,7],":{g,",gtime,":i,",itime,"}",sep = "")
          }
        }
        linlist[parentj,7] = paste("(",spec1,",",spec2,")",sep = "")
        linlist[parentj,8] = linlist[j,3]
        linlist = linlist[-j,]
      } else {
        if (as.numeric(parent) != 0) {
          parentj2 = which(L[,1] == as.numeric(parent))
          comptimeparent2 = L[parentj2,4]

          if(comptimeparent2 > -1 & (comptimeparent2 < as.numeric(linlist[j,3]) | parentj2 <= -1)) {
            linlist[j,4] = L[parentj2,4]
          }
          linlist[j,c(1:3,5)] = L[parentj2,c(1:3,5)]
        }
      }
      if (is.null(nrow(linlist)))    {
        done = 1
        if(ig == FALSE) {
          linlist[7] = paste(linlist[7],":",abs(as.numeric(linlist[3])),";",sep = "")
        } else {
          linlist[7] = paste(linlist[7],";",sep = "")
        }
      } else {
        if(nrow(linlist) == 1)
        {
          done = 1
          if(ig == FALSE)
          {
            linlist[7] = paste("(",linlist[7],":",age,");",sep = "")
          } else {
            linlist[7] = paste(linlist[7],";",sep = "")
          }
        }
      }
    }
  }
  return(linlist[7])
}

#' Internal function to sample a tree from an L table
#' @param L an L-table at the subspecies level,
#'   with number of subspecies of rows
#'   and 6 columns
#' @param age crown age
#' @param samplemethod the sampling method. Can be 'random',
#'   'youngest', 'oldest', 'shortest' or 'longest'
#' @return an L table at the species level
#' @noRd
sampletree = function(L, age, samplemethod = "random")
{
  lenL <- length(L[,1])
  if(samplemethod == "random")
  {
    neworder = DDD::sample2(1:lenL, replace = F)
  }
  if(samplemethod == "youngest")
  {
    neworder = order(L[,3],decreasing = T)
  }
  if (samplemethod == "oldest")
  {
    neworder = order(L[,3],decreasing = F)
  }
  if (samplemethod == "shortest")
  {
    M <- matrix(
      c(rep(-1e10, lenL)),
      nrow = lenL,
      ncol = 1
    )
    # Add column for 'importance', will be sorted on later
    L <- cbind(L, M)
    # Skip 1, because parent 0 is not in L table
    for (i in 2:lenL) {
      # [i,2]: parent
      # [i,6]: species index
      # If parent has different species index than daughter
      if (L[L[i, 2], 6] != L[i, 6]) {
        # If daughter speciation initiation time is greater than parent 'importance'
        if (L[L[i, 2], 7] < L[i, 3]) {
          # Set parent 'importance' to daughter speciation initiation time
          L[L[i, 2], 7] <- L[i, 3]
        }
        # If daughter speciation initiation time is greater than daughter 'importance'
        if (L[i, 7] < L[i, 3]) {
          # Set daughter 'importance' to daughter speciation initiation time
          L[i, 7] <- L[i, 3]
        }
      }
    }
    # Go backwards because we set everyone's parents which goes youngest to oldest, skip 1, because parent 0 is not in L table
    for (i in lenL:2) {
      # If daughter 'importance' is set
      if (L[i, 7] != -1e10) {
        # If parent 'importance' is less than daughter speciation initiation time
        if (L[L[i, 2], 7] < L[i, 3]) {
          # Set parent 'importance' to daughter speciation initiation time
          L[L[i, 2], 7] <- L[i, 3]
        }
      }
    }
    # Skip 1, because parent 0 is not in L table
    for (i in 2:lenL) {
      # If daughter 'importance' is not set
      if (L[i, 7] == -1e10) {
        # If parent 'importance' is set
        if (L[L[i, 2], 7] != -1e10) {
          # If parent 'importance' is greater/equal than daughter speciation initiation time
          if (L[L[i, 2], 7] >= L[i, 3]) {
            # Set daughter 'importance' to daughter speciation initiation time
            L[i, 7] <- L[i, 3]
            # If parent 'importance' is less than daughter speciation initiation time
          } else {
            # Set daughter 'importance' to parent 'importance'
            L[i, 7] <- L[L[i, 2], 7]
          }
        }
      }
    }
    # Order on 'importance'
    neworder <- order(L[, 7], decreasing = TRUE)
    # Remove new column
    L <- L[, -7]
  }
  if (samplemethod == "longest")
  {
    M <- matrix(
      c(rep(-1e10, lenL)),
      nrow = lenL,
      ncol = 1
    )
    # Add column for 'importance', will be sorted on later
    L <- cbind(L, M)
    # Skip 1, because parent 0 is not in L table
    for (i in 2:lenL) {
      # [i,2]: parent
      # [i,6]: species index
      # If parent has different species index than daughter
      if (L[L[i, 2], 6] != L[i, 6]) {
        # If daughter speciation initiation time is greater than parent 'importance'
        if (L[L[i, 2], 7] < L[i, 3]) {
          # Set parent 'importance' to daughter speciation initiation time
          L[L[i, 2], 7] <- L[i, 3]
        }
        # If daughter speciation initiation time is greater than daughter 'importance'
        if (L[i, 7] < L[i, 3]) {
          # Set daughter 'importance' to daughter speciation initiation time
          L[i, 7] <- L[i, 3]
        }
      }
    }
    # Go backwards because we set everyone's parents which goes youngest to oldest, skip 1, because parent 0 is not in L table
    for (i in lenL:2) {
      # If daughter 'importance' is set
      if (L[i, 7] != -1e10) {
        # If parent 'importance' is less than daughter speciation initiation time
        if (L[L[i, 2], 7] < L[i, 3]) {
          # Set parent 'importance' to daughter speciation initiation time
          L[L[i, 2], 7] <- L[i, 3]
        }
      }
    }
    # Skip 1, because parent 0 is not in L table
    for (i in 2:lenL) {
      # If daughter 'importance' is not set
      if (L[i, 7] == -1e10) {
        # If parent 'importance' is set
        if (L[L[i, 2], 7] != -1e10) {
          # If parent 'importance' is greater/equal than daughter speciation initiation time
          if (L[L[i, 2], 7] >= L[i, 3]) {
            # Set daughter 'importance' to daughter speciation initiation time
            L[i, 7] <- L[i, 3]
            # If parent 'importance' is less than daughter speciation initiation time
          } else {
            # Set daughter 'importance' to parent 'importance'
            L[i, 7] <- L[L[i, 2], 7]
          }
        }
      }
    }
    # Order on 'importance'
    neworder <- order(L[, 7], decreasing = FALSE)
    # Remove new column
    L <- L[, -7]
  }
  L2 = L[neworder,]
  ss = NULL
  for(i in 1:lenL)
  {
    if(L2[i,5] == -1)
    {
      if(is.element(L2[i,6],ss) == FALSE)
      {
        ss = c(ss,L2[i,6])
      } else {
        L2[i,5] = age # pseudo extinction just before the present
      }
    }
  }
  L2 = L2[rev(order(L2[,3])),]
  return(L2)
}

pbd_reconstruct = function(L)
{
  L2 = L[order(L[,3]),]
  L3 = L2
  for(i in 1:length(L2[,2]))
  {
    pai = which(abs(L2[,2]) == L2[i,1])
    L3[pai,2] = sign(L2[pai,2]) * i
  }
  orglabs = L3[,1]
  numincspec = length(L3[,2])
  id = 1:numincspec
  L3[,1] = id
  L = L3
  L[1,3] = -1E-10
  if(L[2,3] == 0)
  {
    L[2,3] = 1E-10
  }

  pa = L[,2] # vector of parent species
  ti = L[,3] # vector of speciation-initiation times
  tc = L[,4] # vector of speciation-completion times
  te = L[,5] # vector of extinction times
  sl = L[,6] # vector of species labels
  id2 = id
  tr = NULL ###
  ### print(cbind(id,pa,ti,tc,te,sl))

  # find the branch that went extinct last
  idx1 = which(te == max(te) & te > 0)
  while(length(idx1) != 0)
  {
    # does this extinct branch have offspring?
    # find the offspring who have the extinct branch as parent
    idx2 = rev(which(abs(pa) == idx1))[1]
    if(is.na(idx2))
    {
      # extinct branch does not have offspring
      # extinct branch can be neglected
      ti[idx1] = 0
      tc[idx1] = 0
      te[idx1] = 0
      pa[idx1] = 0
      sl[idx1] = 0
    } else {
      # extinct branch has offspring
      # find the offspring of the offspring of the extinct branch
      idx3 = which(abs(pa) == idx2)
      # was extinct branch good/incipient
      # at time of initiation of offspring?
      if(pa[idx2] > 0)
      {
        # extinct branch was good
        pa[idx3] = idx1
      } else {
        # extinct branch was incipient
        pa[idx3] = sign(pa[idx3]) * idx1
        tc[idx1] = tc[idx2]
      }
      te[idx1] = te[idx2]
      sl[idx1] = sl[idx2] ###
      tr = rbind(tr,c(id2[idx1],id2[idx2])) ##
      id2[idx1] = id2[idx2] ###
      ti[idx2] = 0
      tc[idx2] = 0
      te[idx2] = 0
      pa[idx2] = 0
      sl[idx2] = 0
    }
    #idx1 = rev(which(te > 0))[1]
    idx1 = which(te == max(te) & te > 0)
  }
  ### print(cbind(id,pa,ti,tc,te,sl))
  # eliminate zero rows
  idxs = which(ti != 0)
  diff = (idxs != (1:length(idxs)))
  while(sum(diff) != 0)
  {
    idx1 = (which(diff == 1))[1]
    idx2 = idxs[idx1]
    ti[idx1] = ti[idx2]
    tc[idx1] = tc[idx2]
    te[idx1] = te[idx2]
    pa[idx1] = pa[idx2]
    sl[idx1] = sl[idx2]
    id[idx1] = id[idx2]
    id2[idx1] = id2[idx2]
    ti[idx2] = 0
    tc[idx2] = 0
    te[idx2] = 0
    pa[idx2] = 0
    ### pa[abs(pa) == idx2] = sign(pa[abs(pa) == idx2]) * idx1 ###
    sl[idx2] = 0
    id[idx2] = 0
    id2[idx2] = 0
    idxs = which(ti != 0)
    diff = (idxs != (1:length(idxs)))
  }
  ig = rep(0,length(ti)) # good/incipient flags
  ig[te == -1 & tc != -1] = 1
  ig[te == -1 & tc == -1] = -1
  if(te[1] == -1)
  {
    ig[1] = 1
  }
  zeros = c(which(sl == 0))
  if(length(zeros) > 0)
  {
    id = id[-zeros]
    id2 = id2[-zeros]
    pa = pa[-zeros]
    ti = ti[-zeros]
    te = te[-zeros]
    tc = tc[-zeros]
    sl = sl[-zeros]
    ig = ig[-zeros]
  }
  ### print(cbind(id,pa,ti,tc,te,sl,ig)) ###

  igg = ig # copy of table of good/incipient flags
  ppa = pa # copy of table of parent indices
  its = 0 # index that will run through table
  tt = NULL # table of splitting times
  pp = NULL # table of parent indices
  dd = NULL # table of daughter indices
  sls = NULL # table of species labels
  idxs = which(igg != 0)
  while(idxs[length(idxs)] > 1)
  {
    idx = which.max(ti[idxs])
    di = idxs[idx] # daughter index
    parenti = ppa[di] # parent index (can be negative!)
    pai = which(id == abs(parenti))
    if(igg[pai] == 1 & parenti > 0 & igg[di] == -1)
    {
      #print('1. parent alive, good at event, good at present, daughter inc at present')
      igg[di] = 0
      #igg[pai] = 1 This was already the case
    } else {
      if(igg[pai] == 1 & parenti > 0 & igg[di] == 1)
      {
        #print('2. parent alive, good at event, good at present, daughter good at present')
        igg[di] = 0
        #igg[pai] = 1 This was already the case
        its = its + 1
        tt[its] = ti[di]
        pp[its] = abs(parenti)
        dd[its] = id[di]
        tr = rbind(tr,c(id[di],id2[di])) ##
        sls[its] = sl[di]
      } else {
        if(igg[pai] == 1 & parenti < 0 & igg[di] == -1)
        {
          #print('3. parent alive, inc at event, good at present, daughter inc at present')
          igg[di] = 0
          #igg[pai] = 1 This was already the case
        } else {
          if(igg[pai] == 1 & parenti < 0 & igg[di] == 1)
          {
            #print('4. parent alive, inc at event, good at present, daughter good at present')
            igg[di] = 0
            #igg[pai] = 1 This was already the case
            its = its + 1
            tt[its] = ti[di]
            pp[its] = abs(parenti)
            dd[its] = id[di]
            #dd[its] = id2[di] ##
            tr = rbind(tr,c(id[di],id2[di])) ##
            sls[its] = sl[di]
          } else {
            if(igg[pai] == -1 & parenti < 0 & igg[di] == -1)
            {
              #print('5. parent alive, inc at event, inc at present, daughter inc at present')
              igg[di] = 0
              #igg[pai] = -1 This was already the case
            } else {
              if(igg[pai] == -1 & parenti < 0 & igg[di] == 1)
              {
                #print('6. parent alive, inc at event, inc at present, daughter good at present')
                igg[di] = 0
                igg[pai] = 1
                pp[which(pp == id[di])] = abs(parenti)
                tr = rbind(tr,c(id2[pai],id2[di])) ###
                #id2[pai] = id2[di] ##
                sl[pai] = sl[di] #
              } else {
                if(igg[pai] == 0 & parenti > 0 & igg[di] == -1)
                {
                  #print('7. parent dead, good at event, daughter inc at present')
                  igg[di] = 0
                  igg[pai] = 1
                  tr = rbind(tr,c(id2[pai],id2[di])) ###
                } else {
                  if(igg[pai] == 0 & parenti > 0 & igg[di] == 1)
                  {
                    #print('8. parent dead, good at event, daughter good at present')
                    igg[di] = 0
                    igg[pai] = 1
                    pp[which(pp == id[di])] = abs(parenti)
                    sl[pai] = sl[di]
                    tr = rbind(tr,c(id2[pai],id2[di])) ###
                    ## The daughter keeps her own species label
                  } else {
                    if(igg[pai] == 0 & parenti < 0 & igg[di] == -1)
                    {
                      #print('9. parent dead, inc at event, daughter inc at present')
                      igg[di] = 0
                      igg[pai] = -1
                      tr = rbind(tr,c(id2[pai],id2[di])) ###
                    } else {
                      if(igg[pai] == 0 & parenti < 0 & igg[di] == 1)
                      {
                        #print('10. parent dead, inc at event, daughter good at present')
                        igg[di] = 0
                        igg[pai] = 1
                        pp[which(pp == id[di])] = abs(parenti)
                        tr = rbind(tr,c(id2[pai],id2[di])) ###
                        sl[pai] = sl[di]
                      }
                    }}}}}}}}}
    idxs = which(igg != 0)
  }
  dd = c(1,dd) ##
  pp = c(0,pp) ##
  tt = c(-1e-10,tt) ##
  tt2 = c(0,tt) ##
  te = c(-1,te) ##
  sls = c(sl[1],sls) ##
  dd2 = dd ##
  pp2 = pp ##
  for(i in length(tr[,1]):1)
  {
    dd2[which(dd2 == tr[i,1])] = tr[i,2]
    pp2[which(pp2 == tr[i,1])] = tr[i,2]
  }
  dd = dd2 ##
  pp = pp2 ##
  reconL = cbind(dd,pp,tt,tt,rep(-1,length(dd)),sls,deparse.level = 0)
  ## reconL = rbind(c(1,0,-1e-10,0,-1,1),cbind(dd,pp,tt,tt,rep(-1,length(dd)),sls,deparse.level = 0))
  ### print(reconL) ###
  L = reconL
  L[,1] = orglabs[reconL[,1]]
  L[,2] = c(0,orglabs[reconL[,2]])
  reconL = L; # Leave this semicolon in for teaching
  ### print(reconL) ###
  return(reconL)
}

L2phylo2 = function(L,dropextinct = T)
  # makes a phylogeny out of a matrix with branching times, parent and daughter species, and extinction times
{
  L = L[order(abs(L[,3])),]
  age = L[1,1]
  L[,1] = age - L[,1]
  L[1,1] = -1
  notmin1 = which(L[,4] != -1)
  L[notmin1,4] = age - L[notmin1,4]
  if(dropextinct == T)
  {
    sall = which(L[,4] == -1)
    tend = age
  } else {
    sall = which(L[,4] >= -1)
    tend = (L[,4] == -1) * age + (L[,4] > -1) * L[,4]
  }
  specid = L[,6]
  L = L[,-(4:6)]
  linlist = cbind(L[sall,],paste("S",specid,"-",specid,"-",abs(L[sall,3]),sep = ""),tend)
  done = 0
  while(done == 0)
  {
    #print(linlist)
    j = which.max(linlist[,1])
    daughter = linlist[j,3]
    parent = linlist[j,2]
    parentj = which(parent == linlist[,3])
    parentinlist = length(parentj)
    if(parentinlist == 1)
    {
      spec1 = paste(linlist[parentj,4],":",as.numeric(linlist[parentj,5]) - as.numeric(linlist[j,1]),sep = "")
      spec2 = paste(linlist[j,4],":",as.numeric(linlist[j,5]) - as.numeric(linlist[j,1]),sep = "")
      linlist[parentj,4] = paste("(",spec1,",",spec2,")",sep = "")
      linlist[parentj,5] = linlist[j,1]
      linlist = linlist[-j,]
    } else {
      #linlist[j,1:3] = L[abs(as.numeric(parent)),1:3]
      linlist[j,1:3] = L[which(L[,3] == as.numeric(parent)),1:3]
    }
    if(is.null(nrow(linlist))) { done = 1 }
  }
  linlist[4] = paste(linlist[4],":",linlist[5],";",sep = "")
  phy = ape::read.tree(text = linlist[4])
  tree = ape::as.phylo(phy)
  return(tree)
}


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
#' @param soc Sets whether this age is the stem (1) or crown (2) age
#' @param plotit Sets whether the various trees produced by the function should
#' be plotted or not
#' @param limitsize Sets a maximum to the number of incipient + good species
#' that are created during the simulation; if exceeded, the simulation is
#' aborted and removed.
#' @param add_shortest_and_longest Gives the output of the new samplemethods
#' 'shortest' and 'longest'.
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

pbd_sim_orig = function(pars,age,soc = 2,plotit = FALSE, limitsize = Inf, add_shortest_and_longest = FALSE)
{
  la1 = pars[1]
  la2 = pars[2]
  la3 = pars[3]
  mu1 = pars[4]
  mu2 = pars[5]

  i = 1
  while(i <= soc)
  {
    t = 0
    if(i == 1)
    {
      id1 = 0
      id = id1 + 1
      Sid1 = 0
      Sid = 1
      sg = id
      si = NULL
      L = t(as.matrix(c(id,0,-1e-10,t,-1,1)))
    }
    if(i == 2)
    {
      id = id1 + 1
      Sid = Sid1
      sg = NULL
      si = -id
      L = t(as.matrix(c(id,1,t,-1,-1,1)))
    }

    Ng = length(sg)
    Ni = length(si)
    probs = c(la1*Ng,mu1*Ng,la2*Ni,la3*Ni,mu2*Ni)
    denom = sum(probs)
    probs = probs/denom
    t = t - log(stats::runif(1))/denom

    while(t <= age)
    {
      event = DDD::sample2(1:5,size = 1,prob = probs)
      if(event == 1) # speciation of good species
      {
        parent = as.numeric(DDD::sample2(sg,1))
        id = id + 1
        L = rbind(L,c(id,parent,t,-1,-1,L[abs(parent) - id1,6]))
        si = c(si,-id)
        Ni = Ni + 1
      }
      if(event == 2) # extinction good species
      {
        todie = as.numeric(DDD::sample2(sg,1))
        L[todie - id1,5] = t
        sg = sg[-which(sg == todie)]
        Ng = Ng - 1
      }
      if(event == 3) # completion incipient species
      {
        tocomplete = abs(as.numeric(DDD::sample2(si,1)))
        L[tocomplete - id1,4] = t
        Sid = Sid + 1
        L[tocomplete - id1,6] = Sid
        sg = c(sg,tocomplete)
        si = si[-which(abs(si) == tocomplete)]
        Ng = Ng + 1
        Ni = Ni - 1
      }
      if(event == 4) # speciation incipient species
      {
        parent = as.numeric(DDD::sample2(si,1))
        id = id + 1
        L = rbind(L,c(id,parent,t,-1,-1,L[abs(parent) - id1,6]))
        si = c(si,-id)
        Ni = Ni + 1
      }
      if(event == 5) # extinction incipient species
      {
        todie = abs(as.numeric(DDD::sample2(si,1)))
        L[todie - id1,5] = t
        si = si[-which(abs(si) == todie)]
        Ni = Ni - 1
      }
      if(Ng + Ni > limitsize)
      {
        Ni = 0
        Ng = 0
      }
      probs = c(la1*Ng,mu1*Ng,la2*Ni,la3*Ni,mu2*Ni)
      denom = sum(probs)
      probs = probs/denom
      t = t - log(stats::runif(1))/denom
    }
    if(i == 1)
    {
      if((Ng + Ni) > 0)
      {
        i = i + 1
        L1 = L
        id1 = id
        Sid1 = Sid
        si1 = si
        sg1 = sg
      }
    } else {
      if(i == 2)
      {
        if(checkgood(L,si,sg) == 1)
        {
          i = i + 1
          L2 = L
          si2 = si
          sg2 = sg
        }
      }}
  }
  L = L1
  if(soc == 2)
  {
    L = rbind(L1,L2)
  }
  L0 = L
  absL = L
  absL[,2] = abs(L[,2])
  trans = NULL
  igtree.extinct = phytools::read.simmap(text = detphy(absL,age,ig = T,dropextinct = F))
  igtree.extant = phytools::read.simmap(text = detphy(absL,age,ig = T,dropextinct = T))
  tree = ape::as.phylo(ape::read.tree(text = detphy(absL,age)))

  # Random sampling
  sL_random = sampletree(absL, age, samplemethod = "random")
  stree_random = ape::as.phylo(ape::read.tree(text = detphy(sL_random, age)))
  # Sampling the oldest
  sL_oldest = sampletree(absL, age, samplemethod = "oldest")
  stree_oldest = ape::as.phylo(ape::read.tree(text = detphy(sL_oldest, age)))
  # Sampling the youngest
  sL_youngest = sampletree(absL, age, samplemethod = "youngest")
  stree_youngest = ape::as.phylo(ape::read.tree(text = detphy(sL_youngest, age)))
  # Sampling the shortest
  sL_shortest <- sampletree(absL, age, samplemethod = "shortest")
  stree_shortest <- ape::as.phylo(ape::read.tree(text = detphy(sL_shortest, age)))
  # Sampling the longest
  sL_longest <- sampletree(absL, age, samplemethod = "longest")
  stree_longest <- ape::as.phylo(ape::read.tree(text = detphy(sL_longest, age)))

  sL_random[, 3:5][which(sL_random[, 3:5] == -1)] = age + 1
  sL_random[, 3:5] = age - sL_random[, 3:5]
  sL_random = sL_random[order(sL_random[, 1]), ]
  sL_oldest[, 3:5][which(sL_oldest[, 3:5] == -1)] = age + 1
  sL_oldest[, 3:5] = age - sL_oldest[, 3:5]
  sL_oldest = sL_oldest[order(sL_oldest[, 1]), ]
  sL_youngest[, 3:5][which(sL_youngest[, 3:5] == -1)] = age + 1
  sL_youngest[, 3:5] = age - sL_youngest[, 3:5]
  sL_youngest = sL_youngest[order(sL_youngest[, 1]), ]
  sL_shortest[, 3:5][which(sL_shortest[, 3:5] == -1)] = age + 1
  sL_shortest[, 3:5] = age - sL_shortest[, 3:5]
  sL_shortest = sL_shortest[order(sL_shortest[, 1]), ]
  sL_longest[, 3:5][which(sL_longest[, 3:5] == -1)] = age + 1
  sL_longest[, 3:5] = age - sL_longest[, 3:5]
  sL_longest = sL_longest[order(sL_longest[, 1]), ]

  #sL = sampletree(absL,age)
  #stree = as.phylo(read.tree(text = detphy(sL,age)))
  #sL[,3:5][which(sL[,3:5] == -1)] = age + 1
  #sL[,3:5] = age - sL[,3:5]
  #sL = sL[order(sL[,1]),]

  #rbrts = sort(age - pbd_sim_step2a(L)[[1]])
  reconL = pbd_reconstruct(L0)
  recontree = ape::as.phylo(ape::read.tree(text = detphy(reconL,age)))
  L[,3:5][which(L[,3:5] == -1)] = age + 1
  L[,3:5] = age - L[,3:5]
  L = L[order(L[,1]),]
  #reducL = cbind(L[,3],abs(L[,2:1]),L[,5],L[,4],L[,6])
  #fulltree = L2phylo2(reducL,dropextinct = F)
  reconL[,3:5][which(reconL[,3:5] == -1)] = age + 1
  reconL[,3:5] = age - reconL[,3:5]
  reconL = reconL[order(reconL[,1]),]

  if(plotit == TRUE)
  {
    graphics::par(mfrow = c(3,3))
    cols = stats::setNames(c("gray","black"),c("i","g"))
    phytools::plotSimmap(igtree.extinct,colors = cols)
    phytools::plotSimmap(igtree.extant,colors = cols)
    graphics::plot(tree, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
    graphics::plot(stree_random, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
    graphics::plot(stree_oldest, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
    graphics::plot(stree_youngest, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
    graphics::plot(recontree, edge.width = 2, font = 1, label.offset = 0.1, cex = 1.0)
    graphics::par(mfrow = c(1,1))
  }

  Ltreeslist = list(tree = tree,stree_random = stree_random,stree_oldest = stree_oldest,stree_youngest = stree_youngest,L = L,sL_random = sL_random,sL_oldest = sL_oldest,sL_youngest = sL_youngest,igtree.extinct = igtree.extinct,igtree.extant = igtree.extant,recontree = recontree,reconL = reconL,L0 = L0)
  if (add_shortest_and_longest == TRUE) {
    Ltreeslist <- list(tree = tree,stree_random = stree_random,stree_oldest = stree_oldest,stree_youngest = stree_youngest,L = L,sL_random = sL_random,sL_oldest = sL_oldest,sL_youngest = sL_youngest,igtree.extinct = igtree.extinct,igtree.extant = igtree.extant,recontree = recontree,reconL = reconL,L0 = L0, stree_shortest = stree_shortest, stree_longest = stree_longest, sL_shortest = sL_shortest, sL_longest = sL_longest)
  }
  return(Ltreeslist)
}
