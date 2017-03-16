#---------------------------------------------------------------
# FD.generate
#---------------------------------------------------------------

rFDir <- function(n, alpha, p, tau) {
  D <- length(alpha)
  multin <- matrix(rmultinom(n, 1, p), ncol = D, byrow = TRUE)
  x <- matrix(NA, nrow = n, ncol = D, byrow = TRUE)
  for (i in 1:n) {
    for (j in 1:D) {
      if (multin[i, j] == 1) 
        x[i, j] <- rgamma(1, alpha[j] + tau) else x[i, j] <- rgamma(1, alpha[j])
    }
  }
  somma <- apply(x, 1, sum)
  return(x / as.vector(somma))
}

#---------------------------------------------------------------
# FLEXIBLE DIRICHLET
#---------------------------------------------------------------
# ESTIMATION algorithm: step A=starting values step B=SEM cycle step C=final
# estimation function
#---------------------------------------------------------------

#---------------------------------------------------------------
# step A starting values
#---------------------------------------------------------------

#---------------------------------------------------------------
# A.1 cluster
#---------------------------------------------------------------

#---------------------------------------------------------------
# cluster based on barycenter (median), restricted on certain colums

bar.median.xy.start <- function(data) {
  D <- dim(data)[2]
  # power set of the columns
  lista_vet_col <- vector(mode = "list", length = 2 ^ D)
  for (i in 1:D) for (j in 1:(2^(i - 1))) lista_vet_col[[(2^(i - 1)) + j]] <- c(lista_vet_col[[j]], 
                                                                                i)
  labels <- numeric()
  for (vet.col in lista_vet_col) {
    # only groups consisting of at least two columns are considered
    if (length(vet.col) > 1) {
      bar.m1 <- apply(data, 2, median)
      mat <- t(t(data)/bar.m1)
      # choice of labels (only the subset of columns specified by vet.col is
      # considered)
      etic.bar <- apply(mat[, vet.col], 1, which.max)
      # relabelling according to the positions of the corresponding columns
      etic.bar <- vet.col[etic.bar]
      labels <- cbind(labels, as.vector(etic.bar))
    }
  }
  labels
}

#---------------------------------------------------------------
# cluster based on k-means (k=2...D) with ternary coordinates

ternary.D.start <- function(data, k) {
  D <- dim(data)[2]
  # n <- dim(data)[1] #proviamo a togliere initialitazion of the matrix for the
  # ternary coordinates transformation
  transform.lin <- matrix(0, ncol = D - 1, nrow = D)
  for (i in 1:(D - 1)) transform.lin[i + 1, i] <- sqrt((1 + i)/(2 * i))
  for (i in 1:(D - 2)) for (j in (i + 2):D) transform.lin[j, i] <- 1/sqrt(2 * 
                                                                            i * (1 + i))
  # trasnformation applied to original data
  data.new <- data %*% transform.lin
  # kmeans such that the number of clusters k is random from 2 to D
  labels <- numeric()
  for (k in 2:D) {
    k.agg <- kmeans(data.new, centers = k, iter.max = 100, algorithm = "Hartigan-Wong")
    labels <- cbind(labels, k.agg$cluster)
  }
  labels
}

#-----------------------------------------------------------------
# cluster based on barycenter (median). Only for D=2

bar.median.start <- function(data) {
  bar.m1 <- apply(data, 2, median)
  mat <- t(t(data)/bar.m1)
  etic.bar <- apply(mat, 1, which.max)
  etic.bar
}

#-----------------------------------------------------------------
# cluster based on barycenter (mean). Only for D=2

bar.mean.start <- function(data) {
  bar.m1 <- apply(data, 2, mean)
  mat <- t(t(data)/bar.m1)  #divide each row by barycenter coordinates
  etic.bar <- apply(mat, 1, which.max)
  etic.bar
}

#---------------------------------------------------------------
# cluster based on k-means (k=2) with original coordinates. Only for D=2

orig.start <- function(data) {
  D <- dim(data)[2]
  k.agg <- kmeans(data[, 1:D - 1], centers = 2, iter.max = 100, algorithm = "Hartigan-Wong")
  labels <- k.agg$cluster
  labels
}
#----------------------------------------------------------------


#---------------------------------------------------------------------------
# functions for simulation/application with different initial cluster methods


# FINAL SHORT ternary and median based clustering methods for D >= 3 mean,
# median and k-means for D=2
cluster.f.short <- function(camp, verbose) {
  nsize <- dim(camp)[1]
  Dsize <- dim(camp)[2]
  if (Dsize == 2) {
    labs <- matrix(NA, nrow = nsize, ncol = 3)
    labs[, 1] <- bar.mean.start(camp)
    labs[, 2] <- bar.median.start(camp)
    labs[, 3] <- orig.start(camp)
    if (verbose) 
      cat(c("3/3\n"))
  } else {
    labs <- matrix(NA, nrow = nsize, ncol = (2^Dsize - 2))
    if (verbose) 
      cat(c(1, "/", 2^Dsize - 2, "\n"))
    labs[, 1:(2^Dsize - Dsize - 1)] <- bar.median.xy.start(camp)  #2^D-D-1 labels
    if (verbose) 
      cat(c(2^Dsize - Dsize - 1, "/", 2^Dsize - 2, "\n"))
    labs[, (2^Dsize - Dsize):(2^Dsize - 2)] <- ternary.D.start(camp)  #D-1 labels
    if (verbose) 
      cat(c(2^Dsize - 2, "/", 2^Dsize - 2, "\n"))
  }
  labs
}

#---------------------------------------------------------------------------

#---------------------------------------------------------------
# A.2: relabeling on the basis of mean positions
#---------------------------------------------------------------

#---------------------------------------------------------------------------------
# Funtion which finds the label permutations compatible with the postions of
# maxima. Two main cycles that use w. The first creates vectors given the
# constrains of pos.max. The second fills the gaps (the zeros) with all the
# possible permutations of remaining numbers.

perm.compat <- function(pos.max) {
  D <- length(pos.max)
  p.comp <- matrix(0, nrow = 1, ncol = D)
  assenti <- numeric()
  for (w in 1:D) {
    pos <- which(pos.max == w)
    npos <- length(pos)
    if (npos == 0) 
      assenti <- c(assenti, w) else {
        old.ncomp <- dim(p.comp)[1]
        new.ncomp <- old.ncomp * npos
        p.comp <- matrix(t(p.comp), ncol = D, nrow = new.ncomp, byrow = T)
        for (i in 1:npos) for (j in 1:old.ncomp) p.comp[j + (i - 1) * old.ncomp, 
                                                        pos[i]] <- w
      }
  }
  nassenti <- length(assenti)
  if (nassenti != 0) {
    for (w in 1:nassenti) {
      # print(p.comp)
      nzeri <- nassenti + 1 - w
      old.ncomp <- dim(p.comp)[1]
      pos.zeri <- matrix(0, nrow = old.ncomp, ncol = nzeri)
      for (k in 1:old.ncomp) pos.zeri[k, ] <- which(p.comp[k, ] == 0)
      new.ncomp <- old.ncomp * nzeri
      p.comp <- matrix(t(p.comp), ncol = D, nrow = new.ncomp, byrow = T)
      for (i in 1:nzeri) for (j in 1:old.ncomp) p.comp[j + (i - 1) * old.ncomp, 
                                                       pos.zeri[j, i]] <- assenti[w]
    }
  }
  p.comp
}
#---------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------
# Computation of group means and variances

# ATTENTION: IT CHECKS WHETHER ANY GROUP SIZE IS 0 (Matrices with number of
# rows k<D) OR 1 (NA values in var matrix) RETURNS: group_mean DxD matrix with
# possibly zero rows g_tot and group_var kxD where k=number of not null and
# not unitary size groups group_num Dx1 vector with possibly zero values g_num
# kx1 vector with group sizes != to 0 or 1

group_mean_var <- function(data, etic) {
  D <- dim(data)[2]
  group_num <- numeric(D)
  for (i in 1:D) {
    group_num[i] <- length(which(etic == i))
  }
  # vector of length D possibly containing zeros k<-sum(which(group_num > 1))
  g_num <- group_num[which(group_num > 1)]
  # k <- length(g_num)
  g_mean <- as.matrix(aggregate(data, by = list(etic), FUN = mean)[, -1])  #nrows= number of ni>0
  g_var <- as.matrix(aggregate(data, by = list(etic), FUN = var)[, -1])  #nrows= number of ni>0, possibly NA values if ni=1
  if (any(group_num == 0)) {
    g_new <- matrix(nrow = D, ncol = D)  #DxD (for means) with rows of zeroes when ni=0
    v_new <- matrix(nrow = D, ncol = D)  #DxD (for variances) with rows of zeroes when ni=0
    lab0 <- which(group_num == 0)
    lab1 <- which(group_num != 0)
    g_new[lab0, ] <- rep(0, D)
    g_new[lab1, ] <- g_mean
    v_new[lab0, ] <- rep(0, D)
    v_new[lab1, ] <- g_var
    group_mean <- g_new  #DxD with possible zeroes
    group_var <- v_new  #DxD with possible zeroes
  } else {
    group_mean <- g_mean
    group_var <- g_var
  }
  if (any(group_num == 1)) {
    group_tot <- group_mean[which(group_num > 1), ]
    group_var_tot <- group_var[which(group_num > 1), ]
  }
  if (any(group_num == 1)) 
    return(list(group_mean = group_mean, g_tot = group_tot, group_var = group_var_tot, 
                group_num = group_num, g_num01 = g_num)) else return(list(group_mean = group_mean, g_tot = g_mean, group_var = g_var, 
                                                                          group_num = group_num, g_num01 = g_num))
}
#--------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------
# function for relabeling in the general case (possible ties)

relabel.general <- function(data, etic) {
  ndata <- length(etic)  #n
  centroidi <- group_mean_var(data, etic)$group_mean
  posiz.max <- apply(centroidi, 2, which.max)
  fr <- table(posiz.max)
  D <- length(posiz.max)
  nfr <- length(fr)  #numb. of assigned labels
  if (nfr == D) {
    new.pos <- numeric(ndata)
    for (i in 1:D) new.pos[which(etic == i)] <- which(posiz.max == i)
  } else {
    p.comp <- perm.compat(posiz.max)
    n_p_comp <- dim(p.comp)[1]  #numb. of compatible permutations
    sup.lik <- numeric(n_p_comp)
    for (s in 1:n_p_comp) {
      new.pos <- numeric(ndata)
      for (i in 1:D) new.pos[which(etic == i)] <- which(p.comp[s, ] == i)
      gmv <- group_mean_var(data, new.pos)
      iniz <- iniz_alpha_tau(gmv$group_mean, gmv$g_tot, gmv$group_var, gmv$group_num, 
                             gmv$g_num01)
      camp <- cbind(data, new.pos)
      smv.fdir <- optim(iniz, menologl, data = camp, method = "L-BFGS-B", 
                        hessian = T, lower = rep(0.01, D + 1), upper = rep(Inf, D + 1))
      sup.lik[s] <- -smv.fdir$value
    }
    sup <- which.max(sup.lik)  #position of the permutation which, among the comaptible ones maximizes the lik
    p_comp_chosen <- p.comp[sup, ]
    new.pos <- numeric(ndata)
    for (i in 1:D) new.pos[which(etic == i)] <- which(p_comp_chosen == i)
  }
  new.pos
}
#---------------------------------------------------------------------------------


#-----------------------------------------------------------
# FINAL: relabeling of a group of cluster initializations

labs.new <- function(camp, labs.old, verbose) {
  nsize <- dim(labs.old)[1]
  n.methods <- dim(labs.old)[2]
  if (verbose) {
    tempo <- as.numeric(Sys.time())
    cat(c(1, "/", n.methods, "\n"))
  }
  l.new <- matrix(NA, nrow = nsize, ncol = n.methods)
  for (i in 1:n.methods) {
    if (verbose) 
      if (tempo + 15 < as.numeric(Sys.time())) {
        tempo <- as.numeric(Sys.time())
        cat(c(i, "/", n.methods, "\n"))
      }
    l.new[, i] <- relabel.general(camp, labs.old[, i])
  }
  if (verbose) 
    cat(c(n.methods, "/", n.methods, "\n"))
  l.new
}
#-----------------------------------------------------------


#---------------------------------------------------------------
# A.2b: given a partition it maximizes the loglik wrt alpha and tau
#---------------------------------------------------------------

#-----------------------------------------------------------
# Functions to compute the likelihood of the dirichlet model

dir_moment_match <- function(dati) {
  # D <- ncol(dati) n <- nrow(dati)
  alpha <- apply(dati, 2, mean)
  m2 <- apply(dati * dati, 2, mean)
  s <- (alpha - m2)/(m2 - alpha^2)
  s <- median(s)
  alpha <- alpha * s
  alpha
}

# -logL for Dir
menologv_dir <- function(par, sample) {
  asum <- sum(par)
  n <- dim(sample)[1]
  somme <- apply(log(sample), 2, sum)
  -n * lgamma(asum) + n * sum(lgamma(par)) - sum((par - 1) * somme)
}

logv_dir <- function(par, sample) -menologv_dir(par, sample)

# optim for Dir model
dir.lik <- function(data) {
  D <- dim(data)[2]
  n <- dim(data)[1]
  stime.m <- dir_moment_match(data)
  smv.dir <- optim(stime.m, menologv_dir, sample = data, method = "L-BFGS-B", 
                   hessian = T, lower = rep(0.01, D), upper = rep(Inf, D))
  return(list(smv.a = smv.dir$par, logL.dir = -smv.dir$value, SE.a = sqrt(diag(solve(smv.dir$hessian))), 
              AIC.dir = 2 * D - 2 * -smv.dir$value, BIC.dir = -2 * -smv.dir$value + 
                D * log(n)))
}

#------------------------------------------------------------------------------------------------------------
# alpha and tau initialization given a partition

# ATTENTION: 5 ARGUMENTS (the output of group_mean_var). Initialization of
# relative alphas and tau uses the DxD matrix group_mean initialization of
# alphaplus + tau uses kxD matrices g_tot and group_var THE FUNCTION TO
# ESTIMATE alphaplus + tau USES MEAN OF GROUP ESTIMATES TO GAIN STABILITY IN
# THE CASE OF SMALL VARIANCES

iniz_alpha_tau <- function(m.medie, m.tot, m.var, n.gruppi, n.gruppi01) {
  D <- length(n.gruppi)
  k <- length(n.gruppi01)
  pi <- n.gruppi/sum(n.gruppi)
  di <- n.gruppi01/sum(n.gruppi01)
  m.var <- m.var * (n.gruppi01 - 1)/n.gruppi01
  m_medie_nodiag <- m.medie - diag(diag(m.medie))  #matrice delle medie con zeri sulla diagonale
  if (all(pi != 1)) {
    somme.pesi <- 1 - pi  #D values: for each value the i-th addendum is removed
    stime.alpha <- apply(m_medie_nodiag * pi, 2, sum)/somme.pesi  #weighted means of the elements of each column, apart from the the element on the diagonal
    all <- diag(m.medie) - stime.alpha
    all[all < 0] <- 0
    stima.tau <- mean(all[which(pi > 0)])  #unweighted means of the D estimated of tau (when Pi>0)
    if (k > 1) 
    {
      stima_alpha_tau <- sum((1 - apply(m.tot^2, 1, sum)) * di/apply(m.var, 
                                                                     1, sum)) - 1
    }  # k  <=  D rows of means and var and k weights
    else {
      stima_alpha_tau <- ((1 - sum(m.tot^2))/sum(m.var)) - 1
    }
    # when there is only 1 row of means m.tot and variances m.var
    stima.fin <- c(stime.alpha, stima.tau) * stima_alpha_tau
  } else {
    stima.fin <- c(rep(0, D), which(pi != 0))
  }
  # symbolically alpha = 0 and last element indicates the only label remained
  stima.fin
}
#---------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------
# density function of Dir(alpha_i=alpha + tau e_i) evaluated in x_j=data
fdir_i <- function(par, data, i) {
  D <- length(data)
  alpha <- par[1:D]
  tau <- par[D + 1]
  aplus <- sum(alpha)
  out1 <- lgamma(aplus + tau) + lgamma(alpha[i]) - lgamma(alpha[i] + tau) + 
    tau * log(data[i]) + sum((alpha - 1) * log(data)) - sum(lgamma(alpha))
  exp(out1)
}
#------------------------------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------
# -loglik conditionally on a given partition given by the labels in the last
# (i.e. D + 1) column of data

# IT INCLUDES THE CASES ni=0 AND ni=1

menologl <- function(par, data) {
  D <- dim(data)[2] - 1
  ng <- numeric(D)
  for (i in 1:D) ng[i] <- length(which(data[, D + 1] == i))
  non_null_lab <- which(ng != 0)  #labels with positive frequence
  k <- length(non_null_lab)  #k  <=  D
  out <- numeric()
  for (i in 1:k) {
    X.i <- data[which(data[, D + 1] == non_null_lab[i]), 1:D]
    ni <- ng[non_null_lab[i]]
    fout.i <- numeric()
    if (ni > 1) {
      for (j in 1:ni) fout.i[j] <- fdir_i(par, X.i[j, ], non_null_lab[i])
    } else {
      fout.i <- fdir_i(par, X.i, non_null_lab[i])
    }
    out[i] <- sum(log(unlist(fout.i)))
  }
  -sum(out)
}
#-------------------------------------------------------------------------------------------


#---------------------------------------------------------------
# B: SEM (Given an initial partition, they give a final one according to
# maximization of classified lik)
#---------------------------------------------------------------

# MODIFIED (for large alpha does not manage xij^alphaj)
#---------------------------------------------------
densit_fgn1 <- function(dati, alpha, tau) {
  sa <- sum(alpha)
  D <- length(alpha)
  n <- nrow(dati)
  ci <- lgamma(alpha) - lgamma(alpha + tau)
  c <- lgamma(sa + tau) - sum(lgamma(alpha))
  lcci <- ci + c  #cci <- exp(ci + c)
  x <- matrix(NA, nrow = n, ncol = D)
  for (j in 1:D) {
    for (i in 1:n) {
      x[i, j] <- exp(lcci[j] + sum((alpha - 1) * (log(dati[i, ]))) + tau * 
                       log(dati[i, j]))  #cci[j] * prod(dati[i, ]^(alpha - 1)) * dati[i, j]^tau
    }
  }
  x
}
#---------------------------------------------------

#---------------------------------------------------
log_verosim_fgn1 <- function(dati, alpha, tau, p) {
  x <- densit_fgn1(dati, alpha, tau)
  for (i in 1:ncol(dati)) {
    x[, i] <- x[, i] * p[i]
  }
  x <- apply(x, 1, sum)
  x <- log(x)
  x <- sum(x)
  x
}
#---------------------------------------------------


#--------------------------------------------------------------------------
# nxD matrix containing f_D(x_j;alpha + tau e_i) in jth row and ith column
# (useful for stochastic step3)

fd <- function(data, alpha, tau) {
  sa <- sum(alpha)
  D <- length(alpha)
  n <- nrow(data)
  ci <- lgamma(alpha) - lgamma(alpha + tau)
  c <- lgamma(sa + tau) - sum(lgamma(alpha))
  lcci <- ci + c
  x <- matrix(NA, nrow = n, ncol = D)
  for (j in 1:D) {
    for (i in 1:n) {
      x[i, j] <- exp(lcci[j] + sum((alpha - 1) * (log(data[i, ]))) + tau * 
                       log(data[i, j]))
    }
  }
  x
}
#----------------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------
# Step2: Parameter estimates and value of sup loglik IT INCLUDES THE DIRICHLET
# CASE

step2 <- function(data, etic) {
  n <- dim(data)[1]
  D <- dim(data)[2]
  a0 <- rep(0, D)
  gmv <- group_mean_var(data, etic)
  iniz <- iniz_alpha_tau(gmv$group_mean, gmv$g_tot, gmv$group_var, gmv$group_num, 
                         gmv$g_num01)
  if (iniz[1] != 0) {
    # not dirichlet case
    camp <- cbind(data, etic)
    smv.fdir <- optim(iniz, menologl, data = camp, method = "L-BFGS-B", hessian = T, 
                      lower = rep(0.01, D + 1), upper = rep(Inf, D + 1))
    return(list(sup_loglik = -smv.fdir$value, alpha.est = smv.fdir$par[1:D], 
                p.est = gmv$group_num/n, tau.est = smv.fdir$par[D + 1]))
  } else {
    b0 <- a0
    b0[iniz[D + 1]] <- 1
    return(list(sup_loglik = 0, alpha.est = a0, p.est = b0, tau.est = 0))
  }
}
#-------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------
# Stochastic step 3: given the parameter estimates it creates a partition
# according to SEM

step3.stoch <- function(data, alpha, p, tau) {
  n <- dim(data)[1]
  # D <- dim(data)[2]
  prob.post <- fd(data, alpha, tau)
  prob.post <- t(t(prob.post) * p)  #each i-th column is multiplied by pi
  tot.riga <- apply(prob.post, 1, sum)
  prob.post <- prob.post/tot.riga  #normalization (every element is divedid for the row's total)
  etic.stoch <- numeric(n)
  for (j in 1:n) {
    draw <- rmultinom(1, 1, prob.post[j, ])  #sampled from a Multinom(1, row of prob.post)
    etic.stoch[j] <- which(draw > 0)  #label generation
  }
  etic.stoch
}
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------
# cycle23.stoch for all initializations (identical to cycle23.stoch, but no
# matrix in the output) cycle combining steps 2 (ml) and 3 until convergence
# (SEM)


cycle23.stoch.all <- function(data, etic, iter.max) {
  enne <- dim(data)[1]
  D <- dim(data)[2]
  saved <- numeric()  #it may not reach iter.max (i.e. in Dir case)
  for (s in 1:iter.max) {
    # print(s) TO BE INSERTED AGAIN IN APPLICATIVE CONTEXTS print(s)
    st2 <- step2(data, etic)
    if (st2$sup_loglik != 0) {
      st3 <- step3.stoch(data, st2$alpha.est, st2$p.est, st2$tau.est)
      lik.true <- log_verosim_fgn1(data, st2$alpha.est, st2$tau.est, st2$p.est)
      saved <- rbind(saved, c(st2$alpha.est, st2$p.est, st2$tau.est, lik.true))  #not the sup of the classified lik!!
      etic <- st3
    } else {
      etic <- rep(which(st2$p.est != 0), enne)
      out.dir <- dir.lik(data)
      saved <- rbind(saved, c(out.dir$smv.a, st2$p.est, 0, out.dir$logL.dir))
      break
    }
  }
  n1 <- D + 1
  n2 <- 2 * D
  n3 <- 2 * D + 1
  n4 <- 2 * D + 2
  num.row <- which.max(saved[, n4])
  a.fd <- saved[num.row, 1:D]
  p.fd <- saved[num.row, n1:n2]
  t.fd <- saved[num.row, n3]
  sup.loglik <- saved[num.row, n4]
  return(list(pos.max = num.row, a.fdir = a.fd, p.fdir = p.fd, t.fdir = t.fd, 
              sup.llik = sup.loglik))
}
#----------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------------------
# SEM initialization for a group of clustering methods (ALSO FOR THE CASE D=2)
# Given cluster and relabel, initialize with each method, save initial
# parameter values and suploglik

init.SEM <- function(camp, labels.new, n.iter, verbose) {
  # nsize <- dim(camp)[1]
  D <- dim(camp)[2]
  n.methods <- dim(labels.new)[2]
  nr <- D * 2 + 3  # 1 pos.max,  D alpha,  D pi,  1 tau,  1 suploglik
  init.val <- matrix(NA, nrow = nr, ncol = n.methods)
  for (i in 1:n.methods) {
    if (verbose) 
      cat(c(i, "/", n.methods, "\n"))
    # print(i)
    qq <- cycle23.stoch.all(camp, labels.new[, i], n.iter)
    a.iniz <- qq$a.fdir
    pi.iniz <- qq$p.fdir
    t.iniz <- qq$t.fdir
    p.max <- qq$pos.max
    loglik.iniz <- qq$sup.llik
    # already present in cycle23.stoch.all
    if (t.iniz != 0) {
      # not dirichlet case
      loglik.iniz <- log_verosim_fgn1(camp, a.iniz, t.iniz, pi.iniz)
    } else {
      out.dir <- dir.lik(camp)
      loglik.iniz <- out.dir$logL.dir
    }
    init.val[, i] <- round(c(p.max, a.iniz, pi.iniz, t.iniz, loglik.iniz), 
                           3)
  }
  if (verbose) 
    cat(c(n.methods, "/", n.methods, "\n"))
  init.val
}
#----------------------------------------------------------------------------------------------

#--------------------------------------------------------------------------
# FINAL: best SEM among the ten selected. It prints all positions of best SEMs
# but gives par and loglik of the first only. Remark: With 'which.max' it
# takes the first maximum whereas with 'which(vec == max(vec)) it takes all
# maxima

best.SEM <- function(best.init) {
  D <- (dim(best.init)[1] - 3)/2
  # win <- which.max(best.init[D * 2 + 3, ])
  pos <- which.max(best.init[D * 2 + 3, ])  #position (column n.) of first maximum
  sup.lik <- best.init[D * 2 + 3, pos]
  win <- which(best.init[D * 2 + 3, ] == sup.lik)
  # print(win)
  out <- best.init[, win[1]]
  # print(out)
  out
}
#--------------------------------------------------------------------------

#---------------------------------------------------------------
# STEP C: final estimation function
#---------------------------------------------------------------


#---------------------------------------------------
t_fgn1 <- function(dati, alpha, tau, p) {
  D <- length(alpha)
  n <- nrow(dati)
  x <- densit_fgn1(dati, alpha, tau)
  t <- matrix(NA, nrow = n, ncol = D)
  for (i in 1:D) {
    t[, i] <- x[, i] * p[i]
  }
  somma <- apply(t, 1, sum)
  t <- t/as.vector(somma)
  t
}
#---------------------------------------------------

#---------------------------------------------------
stime_p <- function(dati, alpha, tau, p) {
  t <- t_fgn1(dati, alpha, tau, p)
  stime_p <- apply(t, 2, mean)
  stime_p <- stime_p/sum(stime_p)
  stime_p
}
#---------------------------------------------------

#---------------------------------------------------
menologv_fd <- function(par, data, stime.t) {
  n <- dim(data)[1]
  D <- dim(data)[2]
  alpha <- par[1:D]
  tau <- par[D + 1]
  sa <- sum(alpha) + tau
  t <- stime.t
  somme <- matrix(NA, nrow = n, ncol = 1)
  for (j in 1:n) {
    somme[j, ] <- lgamma(sa) - sum(lgamma(alpha)) + sum(t[j, ] * lgamma(alpha)) - 
      sum(t[j, ] * lgamma(alpha + tau)) + sum(tau * t[j, ] * log(data[j, 
                                                                      ])) + sum((alpha - 1) * log(data[j, ]))
  }
  -sum(somme[, 1])
}
#---------------------------------------------------

#---------------------------------------------------
fit_fd_optim <- function(dati, alpha, tau, p, l, iter) {
  # n <- nrow(dati)
  D <- ncol(dati)
  parametri <- matrix(NA, ncol = (2 * D + 2), nrow = iter)
  niter <- 1
  for (h in 1:iter) {
    old_l <- l
    old_alpha <- alpha
    old_p <- p
    # old_tau <- tau
    stime_dati_t <- t_fgn1(dati, alpha, tau, p)
    p <- stime_p(dati, alpha, tau, p)
    stime.m <- c(alpha, tau)
    stime_hat <- optim(stime.m, menologv_fd, data = dati, stime.t = stime_dati_t, 
                       control = list(maxit = 1), method = "L-BFGS-B", hessian = T, lower = rep(0.01, 
                                                                                                (D + 1)), upper = rep(Inf, (D + 1)))$par
    
    alpha <- stime_hat[-(length(stime_hat))]
    tau <- stime_hat[(length(stime_hat))]
    l <- log_verosim_fgn1(dati, alpha, tau, p)
    parametri[h, ] <- c(l, p, alpha, tau)
    niter <- niter + 1
    if (abs(l - old_l) < 0.001 & max(abs(alpha - old_alpha)) < 0.01 & max(abs(p - 
                                                                              old_p)) < 0.01) {
      # if (abs(l-old_l)<1e-3 & max(abs(alpha-old_alpha))<1e-3 &
      # max(abs(p-old_p))<1e-3) { cat('convergenza')
      break
    }
  }
  # write.table(parametri, 'parametri.csv', sep=';')
  parametri[(niter - 1), ]
}
#---------------------------------------------------


#----------------------------------------------------------------
# FINAL: 'SHORT' FD ESTIMATION (EM starting from the best SEM)

fd.estimation.short <- function(sample, SEM.iter, FIN.iter, verbose) {
  if (verbose) 
    cat("Clustering\n")
  c.f <- cluster.f.short(sample, verbose)
  if (verbose) 
    cat("Labelling\n")
  l.new <- labs.new(sample, labs.old = c.f, verbose)
  if (verbose) 
    cat("Initial SEM\n")
  i.SEM <- init.SEM(sample, l.new, n.iter = SEM.iter, verbose)
  b.SEM <- best.SEM(i.SEM)
  D <- (length(b.SEM) - 3)/2
  if (verbose) 
    cat("Final E-M\n")
  estimates <- fit_fd_optim(sample, b.SEM[2:(D + 1)], b.SEM[(2 * D + 2)], b.SEM[(D + 
                                                                                   2):(2 * D + 1)], b.SEM[(2 * D + 3)], iter = FIN.iter)
  if (verbose) 
    cat("Finished\n")
  estimates
  
}
#------------------------------------------------------------------

#---------------------------------------------------------------
# Density Functions
#---------------------------------------------------------------

dFDir <- function(x, a, p, t) {
  aplus <- sum(a)
  logpart <- lgamma(aplus + t) - sum(lgamma(a)) + sum((a - 1) * log(x)) + log(p) + 
    lgamma(a) - lgamma(a + t) + t * log(x)
  sum(exp(logpart))
}

#--------------------------------------------------------

marginal_univar_dFDir <- function(x, var, a, p, t) {
  aplus <- sum(a)
  p[var] * dbeta(x, a[var] + t, aplus - a[var]) + (1 - p[var]) * dbeta(x, a[var], 
                                                                       aplus - a[var] + t)
}

#------------------------------------------------------------------

#---------------------------------------------------------------
# Graphical representations
#---------------------------------------------------------------


ternary_plot_FD <- function(model, showgrid = T, showdata = T, nlevels = 10) {
  
  data <- model$data
  alpha <- model$alpha
  p <- model$p
  tau <- model$tau
  labs <- colnames(data)
  # D <- dim(data)[2]
  n <- dim(data)[1]
  
  # modified code from compositions package
  usr <- par(list("pty"))
  on.exit(par(usr), add = TRUE)
  par(pty = "s")
  
  s60 <- sin(pi/3)
  c60 <- cos(pi/3)
  s30 <- sin(pi/6)
  c30 <- cos(pi/6)
  att <- seq(0.2, 0.8, by = 0.2)
  label.att <- format(att)
  len.tck <- 0.025
  dist.lab <- 0.03
  dist.axis <- 0.03
  # triangle
  plot(x = c(0, c60, 1, 0), y = c(0, s60, 0, 0), xlim = c(-0.25, 1.25), ylim = c(-0.25, 
                                                                                 s60 + 0.25), main = "FD", type = "n", xlab = "", ylab = "", axes = F)
  lines(x = c(0, c60, 1, 0), y = c(0, s60, 0, 0))
  if (showgrid) {
    # grid
    segments(att, 0, att, -len.tck)
    segments(1 + att * (c60 - 1), att * (s60), 1 + att * (c60 - 1) + len.tck * 
               c30, att * s60 + len.tck * s30)
    segments(c60 + att * -c60, s60 + att * -s60 + len.tck * s30, c60 + att * 
               -c60 + len.tck * -c30, s60 + att * -s60 + len.tck * s30)
    # label
    text(att, -dist.axis, as.graphicsAnnot(label.att), adj = c(0.5, 1))
    text((1 + att * (c60 - 1)) + dist.axis * c30, att * s60 + dist.axis * 
           s30, as.graphicsAnnot(label.att), adj = c(0, 0))
    text(c60 + att * -c60 + len.tck * -c30, s60 + att * -s60 + len.tck * s30, 
         as.graphicsAnnot(label.att), adj = c(1, 0))
    # variables names
    text(-c30 * dist.lab, -s30 * dist.lab, labs[1], adj = c(1, 1))
    text(1 + c30 * dist.lab, -s30 * dist.lab, labs[2], adj = c(0, 1))
    text(c60, s60 + dist.lab, labs[3], adj = c(0.5, 0))
  }
  # end modified code from compositions package
  
  if (showdata) {
    data.ternary <- matrix(NA, ncol = 2, nrow = n)
    for (j in 1:n) data.ternary[j, ] <- c(data[j, 2] + data[j, 3]/2, data[j, 
                                                                          3] * sqrt(3)/2)
    points(data.ternary, pch = 20, cex = 0.25)
  }
  
  aux <- seq(0, 1, by = 0.01)
  points.grid <- expand.grid(x = aux, y = aux)
  npoints <- dim(points.grid)[1]
  points.composit <- matrix(NA, nrow = npoints, ncol = 3)
  points.composit[, 3] <- as.matrix(points.grid[2])/s60
  points.composit[, 2] <- as.matrix(points.grid[1]) - as.matrix(points.grid[2]) * 
    c60/s60
  points.composit[, 1] <- as.matrix(1 - points.composit[, 2] - points.composit[, 
                                                                               3])
  myvalues <- numeric(npoints)
  for (i in 1:npoints) if (any(points.composit[i, ] < 0)) 
    myvalues[i] <- NA else myvalues[i] <- dFDir(points.composit[i, ], alpha, p, tau)
  dim(myvalues) <- c(101, 101)
  
  contour(aux, aux, myvalues, asp = 1, add = T, nlevels = nlevels)
}

#-----------------------------------

zoom_ternary_plot_FD <- function(model, showgrid = T, showdata = T, nlevels = 10) {
  
  data <- model$data
  a <- model$alpha
  p <- model$p
  t <- model$tau
  # labs <- colnames(data)
  D <- dim(data)[2]
  n <- dim(data)[1]
  data.ternary <- matrix(NA, ncol = D - 1, nrow = n)
  for (j in 1:n) data.ternary[j, ] <- c(data[j, 2] + data[j, 3]/2, data[j, 3] * 
                                          sqrt(3)/2)
  
  
  usr <- par(list("pty"))
  on.exit(par(usr), add = TRUE)
  par(pty = "s")
  
  min_X <- min(data.ternary[, 1])
  min_Y <- min(data.ternary[, 2])
  max_X <- max(data.ternary[, 1])
  max_Y <- max(data.ternary[, 2])
  range_X <- max_X - min_X
  range_Y <- max_Y - min_Y
  
  if (showdata) {
    plot(data.ternary, pch = 20, cex = 0.6, axes = showgrid, main = "FD", 
         xlab = "", ylab = "", xlim = c(min_X - range_X * 0.1, max_X + range_X * 
                                          0.1), ylim = c(min_Y - range_Y * 0.1, max_Y + range_Y * 0.1))
  } else {
    plot(x = c(min_X, max_X), y = c(min_Y, max_Y), type = "n", axes = showgrid, 
         main = "FD", xlab = "", ylab = "", xlim = c(min_X - range_X * 0.1, 
                                                     max_X + range_X * 0.1), ylim = c(min_Y - range_Y * 0.1, max_Y + 
                                                                                        range_Y * 0.1))
  }
  
  mm <- mu.FD(a, t)  #group barycenters
  mu.ternary <- matrix(NA, ncol = D - 1, nrow = D)
  for (i in 1:D) {
    mu.ternary[i, ] <- c(mm[i, 2] + mm[i, 3]/2, mm[i, 3] * sqrt(3)/2)
    points(mu.ternary[i, 1], mu.ternary[i, 2], pch = 17, col = 2, cex = 1.2)
  }
  
  # minimap
  lines(x = c(max_X, max_X + range_X * 0.1, max_X + range_X * 0.1/2, max_X), 
        y = c(max_Y, max_Y, max_Y + range_Y * 0.1 * sqrt(3)/2, max_Y), col = 2)
  points(max_X + mu.ternary[, 1] * range_X * 0.1, max_Y + mu.ternary[, 2] * 
           range_Y * 0.1, col = 2, cex = 0.5)
  
  aux_X <- seq(min_X - range_X * 0.1, max_X + range_X * 0.1, by = range_X/85)
  aux_Y <- seq(min_Y - range_Y * 0.1, max_Y + range_Y * 0.1, by = range_Y/85)
  points.grid <- expand.grid(x = aux_X, y = aux_Y)
  c60 <- cos(pi/3)
  s60 <- sin(pi/3)
  npoints <- dim(points.grid)[1]
  points.composit <- matrix(NA, nrow = npoints, ncol = 3)
  points.composit[, 3] <- as.matrix(points.grid[2])/s60
  points.composit[, 2] <- as.matrix(points.grid[1]) - as.matrix(points.grid[2]) * 
    c60/s60
  points.composit[, 1] <- as.matrix(1 - points.composit[, 2] - points.composit[, 
                                                                               3])
  myvalues <- numeric(npoints)
  for (i in 1:npoints) if (any(points.composit[i, ] < 0)) 
    myvalues[i] <- NA else myvalues[i] <- dFDir(points.composit[i, ], a, p, t)
  dim(myvalues) <- c(length(aux_X), length(aux_Y))
  
  contour(aux_X, aux_Y, myvalues, asp = 1, add = T, nlevels = nlevels)
}

#---------------------------------------------

right_triangle_plot_FD <- function(model, var = c(1, 2), showgrid = T, showdata = T, 
                                   nlevels = 10) {
  
  data <- model$data
  alpha <- model$alpha
  p <- model$p
  tau <- model$tau
  labs <- colnames(data)
  var_X <- var[1]
  var_Y <- var[2]
  
  # modified code from compositions package
  usr <- par(list("pty"))
  on.exit(par(usr), add = TRUE)
  par(pty = "s")
  
  att <- seq(0, 1, by = 0.2)
  label.att <- format(att)
  len.tck <- 0.025
  dist.axis <- 0.03
  # triangle
  plot(x = c(0, 1, 0, 0), y = c(0, 0, 1, 0), xlim = c(-0.1, 1.1), ylim = c(-0.1, 
                                                                           1.1), main = "FD", type = "n", xlab = labs[var_X], ylab = labs[var_Y], 
       axes = F)
  lines(x = c(0, 1, 0, 0), y = c(0, 0, 1, 0))
  if (showgrid) {
    # grid
    segments(att, 0, att, -len.tck)
    segments(0, att, -len.tck, att)
    # label
    text(att, -dist.axis, as.graphicsAnnot(label.att), adj = c(0.5, 1))
    text(-dist.axis, att, as.graphicsAnnot(label.att), adj = c(1, 0.5))
  }
  # end modified code from compositions package
  
  if (showdata) {
    points(data[, var], pch = 20, cex = 0.25)
  }
  
  aux <- seq(0, 1, by = 0.01)
  points.grid <- expand.grid(x = aux, y = aux)
  npoints <- dim(points.grid)[1]
  points.composit <- matrix(NA, nrow = npoints, ncol = 3)
  points.composit[, var_X] <- as.matrix(points.grid[, 1])
  points.composit[, var_Y] <- as.matrix(points.grid[, 2])
  points.composit[, -var] <- 1 - points.composit[, var_X] - points.composit[, 
                                                                            var_Y]
  myvalues <- numeric(npoints)
  for (i in 1:npoints) if (any(points.composit[i, ] < 0)) 
    myvalues[i] <- NA else myvalues[i] <- dFDir(points.composit[i, ], alpha, p, tau)
  dim(myvalues) <- c(101, 101)
  
  contour(aux, aux, myvalues, asp = 1, add = T, nlevels = nlevels)
}

#-----------------------------------

zoom_right_triangle_plot_FD <- function(model, var = c(1, 2), showgrid = T, showdata = T, 
                                        nlevels = 10) {
  
  data <- model$data
  a <- model$alpha
  p <- model$p
  t <- model$tau
  labs <- colnames(data)
  # D <- dim(data)[2] n <- dim(data)[1]
  var_X <- var[1]
  var_Y <- var[2]
  usr <- par(list("pty"))
  on.exit(par(usr), add = TRUE)
  par(pty = "s")
  
  min_X <- min(data[, var_X])
  min_Y <- min(data[, var_Y])
  max_X <- max(data[, var_X])
  max_Y <- max(data[, var_Y])
  range_X <- max_X - min_X
  range_Y <- max_Y - min_Y
  
  if (showdata) {
    plot(data[, var], pch = 20, cex = 0.6, axes = showgrid, main = "FD", xlab = labs[var_X], 
         ylab = labs[var_Y], xlim = c(min_X - range_X * 0.1, max_X + range_X * 
                                        0.1), ylim = c(min_Y - range_Y * 0.1, max_Y + range_Y * 0.1))
  } else {
    plot(x = c(min_X, max_X), y = c(min_Y, max_Y), type = "n", axes = showgrid, 
         main = "FD", xlab = labs[var_X], ylab = labs[var_Y], xlim = c(min_X - 
                                                                         range_X * 0.1, max_X + range_X * 0.1), ylim = c(min_Y - range_Y * 
                                                                                                                           0.1, max_Y + range_Y * 0.1))
  }
  
  mm <- mu.FD(a, t)  #group barycenters
  points(mm[, var], pch = 17, col = 2, cex = 1.2)
  
  # minimap
  lines(x = c(max_X, max_X + range_X * 0.1, max_X, max_X), y = c(max_Y, max_Y, 
                                                                 max_Y + range_Y * 0.1, max_Y), col = 2)
  points(max_X + mm[, var_X] * range_X * 0.1, max_Y + mm[, var_Y] * range_Y * 
           0.1, col = 2, cex = 0.5)
  
  aux_X <- seq(min_X - range_X * 0.1, max_X + range_X * 0.1, by = range_X/85)
  aux_Y <- seq(min_Y - range_Y * 0.1, max_Y + range_Y * 0.1, by = range_Y/85)
  points.grid <- expand.grid(x = aux_X, y = aux_Y)
  npoints <- dim(points.grid)[1]
  points.composit <- matrix(NA, nrow = npoints, ncol = 3)
  points.composit[, var_X] <- as.matrix(points.grid[, 1])
  points.composit[, var_Y] <- as.matrix(points.grid[, 2])
  points.composit[, -var] <- 1 - points.composit[, var_X]
  -points.composit[, var_Y]
  myvalues <- numeric(npoints)
  for (i in 1:npoints) if (any(points.composit[i, ] < 0)) 
    myvalues[i] <- NA else myvalues[i] <- dFDir(points.composit[i, ], a, p, t)
  dim(myvalues) <- c(length(aux_X), length(aux_Y))
  
  contour(aux_X, aux_Y, myvalues, asp = 1, add = T, nlevels = nlevels)
}


#------------------------------------------------

marginal_plot_FD <- function(model, var, zoomed = T, showdata = T, showgrid = T) {
  x <- model
  data <- x$data
  a <- x$alpha
  p <- x$p
  t <- x$tau
  labs <- colnames(data)
  # D <- dim(data)[2] n <- dim(data)[1]
  min_X <- min(data[, var])
  max_X <- max(data[, var])
  range_X <- max_X - min_X
  max_Y <- max(c(hist(data[, var], plot = F)$density, marginal_univar_dFDir(seq(0, 
                                                                                1, 0.001), var, a, p, t)))
  
  if (zoomed) 
    curve(marginal_univar_dFDir(x, var = var, a = a, p = p, t = t), min_X - 
            range_X/10, max_X + range_X/10, ylim = c(0, max_Y * 1.1), xlab = labs[var], 
          ylab = "density", lty = 1, axes = showgrid)  #main='FD',
  else curve(marginal_univar_dFDir(x, var = var, a = a, p = p, t = t), 0, 1, 
             n = 1001, ylim = c(0, max_Y * 1.1), xlab = labs[var], ylab = "density", 
             lty = 1, lwd = 2, axes = showgrid)  #main='FD',
  if (showdata) 
    hist(data[, var], freq = F, add = TRUE)
}

#---------------------------------------------------------------
# Barycenters of clusters
#---------------------------------------------------------------

# FD.barycenters

mu.FD <- function(a, t) {
  D <- length(a)
  a.rel <- a/sum(a)
  t.rel <- t/sum(a)
  id <- diag(D)
  mu <- matrix(NA, nrow = D, ncol = D)
  for (i in 1:D) {
    mu[i, ] <- a.rel/(1 + t.rel) + t.rel * id[i, ]/(1 + t.rel)
  }
  mu
}

#---------------------------------------------------------------
# Clusters distances
#---------------------------------------------------------------

# FD.clusterdistances

Kappa <- function(tau, alpha) {
  tau * (digamma(alpha + tau) - digamma(alpha))
}

#---------------------------------------------------------------
# Marginals from FDir model
#---------------------------------------------------------------

### E(X_i)
E.Fdir <- function(alpha, p, tau) {
  (alpha + tau * p)/(sum(alpha) + tau)
}

### V(X_i)
V.Fdir <- function(alpha, p, tau) {
  aplus <- sum(alpha)
  medie <- E.Fdir(alpha, p, tau)
  medie * (1 - medie)/(aplus + tau + 1) + (tau^2 * p * (1 - p))/((aplus + tau) * 
                                                                   (aplus + tau + 1))
}

### Cov(X_i, X_r)
Cov.Fdir <- function(alpha, p, tau) {
  D <- length(alpha)
  medie <- E.Fdir(alpha, p, tau)
  aplus <- sum(alpha)
  covarianze <- matrix(NA, ncol = D, nrow = D)
  varianze <- V.Fdir(alpha, p, tau)
  for (i in 1:D) {
    for (j in 1:D) {
      covarianze[i, j] <- ifelse(i != j, -(medie[i] * medie[j])/(aplus + 
                                                                   tau + 1) - (tau^2 * p[i] * p[j])/((aplus + tau) * (aplus + tau + 
                                                                                                                        1)), varianze[i])
    }
  }
  covarianze
}

### Cor(X_i, X_r)
Cor.Fdir <- function(alpha, p, tau) {
  D <- length(alpha)
  varianze <- V.Fdir(alpha, p, tau)
  covarianze <- Cov.Fdir(alpha, p, tau)
  correlazioni <- matrix(NA, ncol = D, nrow = D)
  for (i in 1:D) {
    for (j in 1:D) {
      correlazioni[i, j] <- covarianze[i, j]/sqrt(varianze[i] * varianze[j])
    }
  }
  correlazioni
}

#---------------------------------------------------------------
# AIC & BIC
#---------------------------------------------------------------

AIC.Fdir <- function(model) {
  data <- model$data
  logl <- model$logL
  D <- dim(data)[2]
  n <- dim(data)[1]
  npar <- 2 * D
  AIC.Fdir <- 2 * npar - 2 * logl
  BIC.Fdir <- -2 * logl + npar * log(n)
  return(list(AIC = AIC.Fdir, BIC = BIC.Fdir))
}

#---------------------------------------------------------------
# CONDITIONAL BOOTSTRAP
#---------------------------------------------------------------

#---------------------------------------------------
prob.cond <- function(alpha, p, tau, data) {
  D <- dim(data)[2]
  n <- dim(data)[1]
  lcost <- lgamma(alpha) - lgamma(alpha + tau)
  numeratore <- matrix(NA, ncol = D, nrow = n)
  for (i in 1:n) numeratore[i, ] <- exp(log(p) + lcost + tau * log(data[i, ]))
  denominatore <- apply(numeratore, 1, sum)
  output <- numeratore/denominatore
  output
}
#---------------------------------------------------

#---------------------------------------------------
Z.multi <- function(alpha, p, tau, data) {
  D <- dim(data)[2]
  n <- dim(data)[1]
  mat.p <- prob.cond(alpha, p, tau, data)
  zeta.star <- matrix(NA, ncol = D, nrow = n)
  for (j in 1:n) zeta.star[j, ] <- rmultinom(1, 1, mat.p[j, ])
  zeta_star_sum <- apply(zeta.star, 2, sum)
  return(list(Z.mat = zeta.star, Z.sum = zeta_star_sum))
}
#---------------------------------------------------

#---------------------------------------------------
# MMatrix of cji
CC <- function(alpha, p, tau, data) {
  D <- dim(data)[2]
  n <- dim(data)[1]
  lcost <- lgamma(alpha) - lgamma(alpha + tau)
  ris <- matrix(NA, ncol = D, nrow = n)
  for (j in 1:n) ris[j, ] <- exp(lcost + tau * log(data[j, ]))
  ris
}
#---------------------------------------------------

#---------------------------------------------------
# matrix of cji divided for cj.=sum_i (pi * cji) (weighted sum, by row)
CC.rel <- function(mat, p) {
  D <- dim(mat)[2]
  n <- dim(mat)[1]
  mat.pond <- matrix(NA, ncol = D, nrow = n)
  for (j in 1:n) mat.pond[j, ] <- mat[j, ] * p
  cj. <- apply(mat.pond, 1, sum)
  risultato <- mat/cj.
  risultato
}
#---------------------------------------------------

#---------------------------------------------------
# Matrix assigning at position if the total of the column minus the element ai
# position ij.
CC.meno <- function(mat) {
  D <- dim(mat)[2]
  n <- dim(mat)[1]
  ris <- matrix(NA, ncol = D, nrow = n)
  c.i <- apply(mat, 2, sum)
  for (j in 1:n) {
    for (i in 1:D) {
      ris[j, i] <- c.i[i] - mat[j, i]
    }
  }
  ris
}
#---------------------------------------------------

#---------------------------------------------------
# Matrix DxD with double sums on J and on j'\neq j of cji rel mat1 is CC.rel
# and mat2 is CC.meno #check
CC.sum <- function(mat1, mat2) {
  D <- dim(mat1)[2]
  ris <- matrix(NA, ncol = D, nrow = D)
  for (i in 1:D) {
    for (h in 1:D) {
      ris[i, h] <- mat1[, i] %*% mat2[, h]
    }
  }
  ris
}
#---------------------------------------------------

#---------------------------------------------------
# matrix A (inf osservata di dim (D-1)x(D-1) relativa ai pi)
expectedA <- function(alpha, p, tau, data) {
  CC.A <- CC(alpha, p, tau, data)
  CC_rel_A <- CC.rel(CC.A, p)
  CC_meno_A <- CC.meno(CC_rel_A)
  CC_sum_A <- CC.sum(CC_rel_A, CC_meno_A)
  D <- dim(CC_sum_A)[1]
  ris <- matrix(NA, ncol = (D - 1), nrow = (D - 1))
  for (i in 1:(D - 1)) {
    for (h in 1:(D - 1)) {
      ris[i, h] <- CC_sum_A[i, D] + CC_sum_A[h, D]
      -CC_sum_A[i, h] - CC_sum_A[D, D]
    }
  }
  ris
}
#---------------------------------------------------

#---------------------------------------------------
# matrix B (inf osservata di dim (D-1)x(D) relativa ai pi e agli ai)
expectedB <- function(mya, myp, myt, data) {
  D <- length(mya)
  myn <- dim(data)[1]
  # sum in j e j' /neq j of the products of the relative Cij
  CC.B <- CC(mya, myp, myt, data)
  CC_rel_B <- CC.rel(CC.B, myp)
  CC_meno_B <- CC.meno(CC_rel_B)
  CC_sum_B <- CC.sum(CC_rel_B, CC_meno_B)
  # first vector of values for the computation of B varying i
  B.V1 <- colSums(CC_rel_B)
  # second vector of values for the computation of B varying h
  mya.tot <- sum(mya)
  B.V2 <- myn * digamma(mya) - myn * digamma(mya.tot + myt) - colSums(log(data))
  # third vector of values for the computation of B varying h
  B.V3 <- digamma(mya + myt) - digamma(mya)
  ris <- matrix(NA, ncol = (D), nrow = (D - 1))
  B.M1 <- matrix(NA, ncol = (D), nrow = (D - 1))
  for (i in 1:(D - 1)) {
    for (h in 1:(D)) {
      B.M1[i, h] <- myp[h] * (CC_sum_B[i, h] - CC_sum_B[D, h])
      # special cases when h=i and h=D
      if (h == i) 
        B.M1[i, h] <- B.M1[i, h] + B.V1[h]
      if (h == D) 
        B.M1[i, h] <- B.M1[i, h] - B.V1[h]
      # final results
      ris[i, h] <- (B.V1[i] - B.V1[D]) * (B.V2[h]) + B.M1[i, h] * B.V3[h]
    }
  }
  ris
}
#---------------------------------------------------

#---------------------------------------------------
# matrix D (inf osservata di dim (D-1)x(1) relativa ai pi al tau)
expectedD <- function(mya, myp, myt, data) {
  D <- length(mya)
  myn <- dim(data)[1]
  # sum in j e j' /neq j of the products of the relative Cij
  CC.D <- CC(mya, myp, myt, data)
  CC_rel_D <- CC.rel(CC.D, myp)
  CC_meno_D <- CC.meno(CC_rel_D)
  CC_sum_D <- CC.sum(CC_rel_D, CC_meno_D)
  # first vector of values for the computation of D varying i
  mya.tot <- sum(mya)
  D.V1 <- rowSums(t(CC_rel_D) * (-myn * digamma(mya.tot + myt) + digamma(mya + 
                                                                           myt) - t(log(data))))
  # second vector of values for the computation of D varying h
  D.V2 <- myp * digamma(mya + myt)
  D.M1 <- t(t(CC_sum_D[1:(D - 1), , drop = F]) - CC_sum_D[D, ])
  D.M2 <- CC_meno_D[, 1:(D - 1), drop = F] - CC_meno_D[, D]
  D.M3 <- t(t(-log(data)) * myp * t(CC_rel_D))
  # final vector
  ris <- D.V1[1:(D - 1)] - D.V1[D] + colSums(D.V2 * t(D.M1))
  +colSums(D.M2 * rowSums(D.M3))
  as.matrix(ris, ncol = 1)
}
#---------------------------------------------------

#---------------------------------------------------
matC.1 <- function(data, alpha, tau, zeta.s) {
  n <- dim(data)[1]
  D <- length(alpha)
  a.sum <- sum(alpha)
  diagonale <- zeta.s * (trigamma(alpha + tau) - trigamma(alpha)) + n * trigamma(alpha)
  matrix(-n * trigamma(a.sum + tau), ncol = D, nrow = D) + diag(diagonale, D)
}
#---------------------------------------------------

#---------------------------------------------------
matE.1 <- function(data, alpha, tau, zeta.s) {
  n <- dim(data)[1]
  # D <- length(alpha)#proviamo
  a.sum <- sum(alpha)
  as.matrix(-n * trigamma(a.sum + tau) + zeta.s * (trigamma(alpha + tau)))
}
#---------------------------------------------------

#---------------------------------------------------
matF.1 <- function(data, alpha, tau, zeta.s) {
  n <- dim(data)[1]
  # D <- length(alpha)#proviamo
  a.sum <- sum(alpha)
  as.matrix(-n * trigamma(a.sum + tau) + sum(zeta.s * (trigamma(alpha + tau))))
}
#---------------------------------------------------

#---------------------------------------------------
matC.2 <- function(data, alpha, tau, zeta.s) {
  n <- dim(data)[1]
  D <- length(alpha)
  a.sum <- sum(alpha)
  data <- log(data)
  x.tilde <- apply(data, 2, sum)
  matC.2 <- matrix(NA, ncol = D, nrow = D)
  for (i in 1:D) {
    for (h in 1:D) {
      matC.2[i, h] <- (n * digamma(a.sum + tau) - n * digamma(alpha[i]) + 
                         x.tilde[i] - zeta.s[i] * (digamma(alpha[i] + tau) - digamma(alpha[i]))) * 
        (n * digamma(a.sum + tau) - n * digamma(alpha[h]) + x.tilde[h] - 
           zeta.s[h] * (digamma(alpha[h] + tau) - digamma(alpha[h])))
    }
  }
  matC.2
}
#---------------------------------------------------

#---------------------------------------------------
matC <- function(data, alpha, tau, zeta.s) {
  matC.1(data, alpha, tau, zeta.s) - matC.2(data, alpha, tau, zeta.s)
}
#---------------------------------------------------

#---------------------------------------------------
matE.2 <- function(data, prob, alpha, tau, zeta.mat) {
  n <- dim(data)[1]
  # D <- length(alpha) #proviamo
  a.sum <- sum(alpha)
  zeta.s <- apply(zeta.mat, 2, sum)
  zeta.m <- zeta.mat
  data <- log(data)
  x.tilde <- apply(data, 2, sum)
  as.matrix((n * digamma(a.sum + tau) - n * digamma(alpha) + x.tilde - zeta.s * 
               (digamma(alpha + tau) - digamma(alpha))) * (n * digamma(a.sum + tau) - 
                                                             sum(zeta.s * digamma(alpha + tau)) + sum(zeta.m * data)))
}
#---------------------------------------------------

#---------------------------------------------------
matE <- function(data, prob, alpha, tau, zeta.s, zeta.mat) {
  matE.1(data, alpha, tau, zeta.s) - matE.2(data, prob, alpha, tau, zeta.mat)
}
#---------------------------------------------------

#---------------------------------------------------
matF.2 <- function(data, alpha, tau, zeta.s, zeta.mat) {
  n <- dim(data)[1]
  data <- log(data)
  a.sum <- sum(alpha)
  (n * digamma(a.sum + tau) - sum(zeta.s * digamma(alpha + tau)) + sum(zeta.mat * 
                                                                         data))^2
}
#---------------------------------------------------

#---------------------------------------------------
matF <- function(data, alpha, tau, zeta.s, zeta.mat) {
  matF.1(data, alpha, tau, zeta.s) - matF.2(data, alpha, tau, zeta.s, zeta.mat)
}
#---------------------------------------------------

#---------------------------------------------------------------------------------------
### SPEEDED BOOTSTRAP
matI.zeri <- function(data, prob, alpha, tau, zeta.s, zeta.mat) {
  D <- dim(data)[2]
  mat_I <- matrix(NA, ncol = 2 * D, nrow = 2 * D)
  mat_I[1:(D - 1), 1:(D - 1)] <- 0
  mat_I[1:(D - 1), D:(2 * D - 1)] <- 0
  mat_I[1:(D - 1), (2 * D)] <- 0
  mat_I[D:(2 * D - 1), 1:(D - 1)] <- 0
  mat_I[D:(2 * D - 1), D:(2 * D - 1)] <- matC(data, alpha, tau, zeta.s)
  mat_I[D:(2 * D - 1), (2 * D)] <- matE(data, prob, alpha, tau, zeta.s, zeta.mat)
  mat_I[(2 * D), 1:(D - 1)] <- 0
  mat_I[(2 * D), D:(2 * D - 1)] <- t(matE(data, prob, alpha, tau, zeta.s, zeta.mat))
  mat_I[(2 * D), (2 * D)] <- matF(data, alpha, tau, zeta.s, zeta.mat)
  mat_I
}

boot.FD <- function(alpha, p, tau, data, B = 500) {
  matrix_I <- 0
  for (b in 1:B) {
    output <- Z.multi(alpha, p, tau, data)
    zeta.sum <- output$Z.sum
    zeta.matrix <- output$Z.mat
    matrice_I <- matI.zeri(data, prob = p, alpha = alpha, tau = tau, zeta.s = zeta.sum, 
                           zeta.mat = zeta.matrix)
    matrix_I <- matrix_I + matrice_I
  }
  matrix_I <- matrix_I/B
  D <- dim(data)[2]
  matrix_I[1:(D - 1), 1:(D - 1)] <- expectedA(alpha, p, tau, data)
  matrix_I[1:(D - 1), D:(2 * D - 1)] <- expectedB(alpha, p, tau, data)
  matrix_I[1:(D - 1), (2 * D)] <- expectedD(alpha, p, tau, data)
  matrix_I[D:(2 * D - 1), 1:(D - 1)] <- t(expectedB(alpha, p, tau, data))
  matrix_I[(2 * D), 1:(D - 1)] <- t(expectedD(alpha, p, tau, data))
  matrix_I.inv <- solve(matrix_I)
  S_E_smv <- sqrt(abs(diag(matrix_I.inv)))
  p_D <- sqrt(abs(sum(matrix_I.inv[1:(D - 1), 1:(D - 1)])))
  output <- t(c(S_E_smv[D:(2 * D - 1)], S_E_smv[1:(D - 1)], p_D, S_E_smv[(2 * 
                                                                            D)]))
  colnames(output) <- c(rep("a", D), rep("p", D), "t")
  output
}