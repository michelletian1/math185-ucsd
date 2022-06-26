#' ## Part A
#' 
load("./smokers.rda")
value = smokers$heart_rate
group = smokers$smoking_group
sum.ks.test <- function(x, g, B=1e4){
  lev <- unique(g)
  lst <- lapply( seq_along(lev), function(i) x[group == lev[i]] )
  names(lst)<-lev
  f <- function(i, j){
    stat_obs = ks.test(i, j)$statistic
    z = c(i, j)
    m = length(i)
    stat_test = numeric(B)
    for (b in 1:B) {
      z_perm = sample(z) # permute the combined sample
      i_perm = z_perm[1:m]
      j_perm = z_perm[-(1:m)]
      stat_test[b] = ks.test(i_perm, j_perm)$statistic
    }
    pval = (sum(abs(stat_test) >= abs(stat_obs))+1)/(B+1)
    return(pval)
  }
  res <- lapply(lst, function(x) lapply(lst, function(y) f(x, y)))
  res<-unlist(res)
  res <- matrix(res, nrow = length(lst), ncol = length(lst), byrow = T)
  row.names(res) <- colnames(res) <- names(lst)
  return(res)
}
pairwise_pval_sum = sum.ks.test(value, group)
#' 
#' ## Part B
#'
load("./smokers.rda")
value = smokers$heart_rate
group = smokers$smoking_group
max.ks.test <- function(x, g, B=1e4){
  lev <- unique(g)
  lst <- lapply( seq_along(lev), function(i) x[group == lev[i]] )
  names(lst)<-lev
  f <- function(i, j){
    stat_obs = ks.test(i, j, alternative='greater')$statistic
    z = c(i, j)
    m = length(i)
    stat_test = numeric(B)
    for (b in 1:B) {
      z_perm = sample(z) # permute the combined sample
      i_perm = z_perm[1:m]
      j_perm = z_perm[-(1:m)]
      stat_test[b] = ks.test(i_perm, j_perm, alternative='greater')$statistic
    }
    pval = (sum(abs(stat_test) >= abs(stat_obs))+1)/(B+1)
    return(pval)
  }
  res <- lapply(lst, function(x) lapply(lst, function(y) f(x, y)))
  res<-unlist(res)
  res <- matrix(res, nrow = length(lst), ncol = length(lst), byrow = T)
  row.names(res) <- colnames(res) <- names(lst)
  return(res)
}
pairwise_pval_max = max.ks.test(value, group)
#'
#' ## Part C
#' compare results :
pairwise_pval_max
pairwise_pval_sum
#' The non and the light have the same p-value and we failed to reject the null hypothesis. The moderate also failed to reject null hypothesis and for the heavy we successfully reject the null hypothesis.
#' For sum moderate we can see it display more higher outcomes than the max same goes for the heavy.
#' <BR><BR>
#'
#' # Problem 2
#' 
I = c(5, 10, 20, 50)
J = c(2, 5, 10)
M = 1e4
k = 1
pval_matrix = matrix(0, length(I) * length(J), M)
for (i in I){
  for (j in J){
    dat = matrix(rnorm(i*j, mean=0, sd=1), nrow=i, ncol=j)
    for (m in 1:M){
      dat_permute = t(apply(dat, MARGIN=1, FUN=sample))
      pval_matrix[k, m] = pval_matrix[k, m] + friedman.test(dat_permute)$p.value
    }
    k = k+1
  }
}
for (k in 1:12){
  d <- density(pval_matrix[k, ]) # returns the density data
  plot(d) # plots the results 
}
#' If we increase the j, the plot will become smoother or will be more stable. If we increase the I the distribution will become more smoother or more uniform. If we decrease the value of I and J the plot will become much shorter and will not become uniform. 
#' 
#' <BR><BR>
#' 
#' Problem 3
#' 
gene.pval <- function(y) {
  row_length = nrow(y)
  col_length = ncol(y)
  k = row_length / 2
  p_values = rep(0, col_length)
  for (j in 1:col_length){
    p_values[j] = t.test(y[1:k, j], y[-1:-k, j])$p.value
  }
  return(p_values)
}
m = 1e3
k = 10
for (m0 in c(m-0, m-10, m-100)){
  mean_g1 = 0
  mean_g2_t = mean_g1 + 0
  deltas = seq(1, 2.8, 0.2)
  j = 1
  
  FWEPs.bon = numeric(length(deltas))
  FDPs.bon = numeric(length(deltas))
  FNPs.bon = numeric(length(deltas))
  
  FWEPs.holm = numeric(length(deltas))
  FDPs.holm = numeric(length(deltas))
  FNPs.holm = numeric(length(deltas))
  
  FWEPs.hoch = numeric(length(deltas))
  FDPs.hoch = numeric(length(deltas))
  FNPs.hoch = numeric(length(deltas))
  
  FWEPs.bh = numeric(length(deltas))
  FDPs.bh = numeric(length(deltas))
  FNPs.bh = numeric(length(deltas))
  
  FWEPs.by = numeric(length(deltas))
  FDPs.by = numeric(length(deltas))
  FNPs.by = numeric(length(deltas))
  for (delta in deltas){
    C = 100
    FWEP.bon = numeric(C)
    FDP.bon = numeric(C)
    FNP.bon = numeric(C)
    
    FWEP.holm = numeric(C)
    FDP.holm = numeric(C)
    FNP.holm = numeric(C)
    
    FWEP.hoch = numeric(C)
    FDP.hoch = numeric(C)
    FNP.hoch = numeric(C)
    
    FWEP.bh = numeric(C)
    FDP.bh = numeric(C)
    FNP.bh = numeric(C)
    
    FWEP.by = numeric(C)
    FDP.by = numeric(C)
    FNP.by = numeric(C)
    for (i in 1:C){
      mean_g2_f = mean_g1 + delta
      g1 = matrix(rnorm(k*m, mean=mean_g1, sd=1), k, m)
      g2_t = matrix(rnorm(k*m0, mean=mean_g2_t, sd=sqrt(2)), k, m0)
      g2_f = matrix(rnorm(k*(m-m0), mean=mean_g2_f, sd=sqrt(2)), k, (m-m0))
      g2 = cbind(g2_t, g2_f)
      dat = rbind(g1, g2)
      pval = gene.pval(dat)
      
      pval.bon = p.adjust(pval, "bon")
      reject.bon = (pval.bon <= 0.10)
      R.bon = sum(reject.bon)
      T1.bon = sum((reject.bon[1:900] == rep(TRUE, 900)))
      T2.bon = sum(tail(reject.bon, 100) == rep(FALSE, 100))
      FWEP.bon[i] = T1.bon > 1
      FDP.bon[i] = T1.bon / max(R.bon, 1)
      FNP.bon[i] = T2.bon / max(m-m0, 1)
      
      pval.holm = p.adjust(pval, "holm")
      reject.holm = (pval.holm <= 0.10)
      R.holm = sum(reject.holm)
      T1.holm = sum((reject.holm[1:900] == rep(TRUE, 900)))
      T2.holm = sum(tail(reject.holm, 100) == rep(FALSE, 100))
      FWEP.holm[i] = T1.holm > 1
      FDP.holm[i] = T1.holm / max(R.holm, 1)
      FNP.holm[i] = T2.holm / max(m-m0, 1)
      
      pval.hoch = p.adjust(pval, "hoch") 
      reject.hoch = (pval.hoch <= 0.10)
      R.hoch = sum(reject.hoch) 
      T1.hoch = sum((reject.hoch[1:900] == rep(TRUE, 900)))
      T2.hoch = sum(tail(reject.hoch, 100) == rep(FALSE, 100))
      FWEP.hoch[i] = T1.hoch > 1
      FDP.hoch[i] = T1.hoch / max(R.hoch, 1)
      FNP.hoch[i] = T2.hoch / max(m-m0, 1)
      
      pval.bh = p.adjust(pval, "BH")
      reject.bh = (pval.bh <= 0.10)
      R.bh = sum(reject.bh) 
      T1.bh = sum((reject.bh[1:900] == rep(TRUE, 900)))
      T2.bh = sum(tail(reject.bh, 100) == rep(FALSE, 100))
      FWEP.bh[i] = T1.bh > 1
      FDP.bh[i] = T1.bh / max(R.bh, 1)
      FNP.bh[i] = T2.bh / max(m-m0, 1)
      
      pval.by = p.adjust(pval, "BY")
      reject.by = (pval.by <= 0.10)
      R.by = sum(reject.by)
      T1.by = sum((reject.by[1:900] == rep(TRUE, 900)))
      T2.by = sum(tail(reject.by, 100) == rep(FALSE, 100))
      FWEP.by[i] = T1.by > 1
      FDP.by[i] = T1.by / max(R.by, 1)
      FNP.by[i] = T2.by / max(m-m0, 1)
    }
    FWEPs.bon[j] = mean(FWEP.bon)
    FDPs.bon[j] = mean(FDP.bon)
    FNPs.bon[j] = mean(FNP.bon)
    
    FWEPs.holm[j] = mean(FWEP.holm)
    FDPs.holm[j] = mean(FDP.holm)
    FNPs.holm[j] = mean(FNP.holm)
    
    FWEPs.hoch[j] = mean(FWEP.hoch)
    FDPs.hoch[j] = mean(FDP.hoch)
    FNPs.hoch[j] = mean(FNP.hoch)
    
    FWEPs.bh[j] = mean(FWEP.bh)
    FDPs.bh[j] = mean(FDP.bh)
    FNPs.bh[j] = mean(FNP.bh)
    
    FWEPs.by[j] = mean(FWEP.by)
    FDPs.by[j] = mean(FDP.by)
    FNPs.by[j] = mean(FNP.by)
    
    j = j+1
  }
  par(mar = c(5, 5, 1, 10))
  matplot(deltas, cbind(FWEPs.bon, FWEPs.holm, FWEPs.hoch, FWEPs.bh, FWEPs.by), type = "l", lwd=3, ylab="FWEP fraction")
  legend('right', xpd = TRUE, inset = c(-0.4, 0), c('bonferroni', 'holm', 'hochberg', 'benjamini-hochberg', 'benjamini-yekutieli'), lty=1:6, col=1:6, lwd=3, bg='white')
  
  matplot(deltas, cbind(FDPs.bon, FDPs.holm, FDPs.hoch, FDPs.bh, FDPs.by), type = "l", lwd=3, ylab="FDP fraction")
  legend('right', xpd = TRUE, inset = c(-0.4, 0), c('bonferroni', 'holm', 'hochberg', 'benjamini-hochberg', 'benjamini-yekutieli'), lty=1:6, col=1:6, lwd=3, bg='white')
  
  matplot(deltas, cbind(FNPs.bon, FNPs.holm, FNPs.hoch, FNPs.bh, FNPs.by), type = "l", lwd=3, ylab="FNP fraction")
  legend('right', xpd = TRUE, inset = c(-0.4, 0), c('bonferroni', 'holm', 'hochberg', 'benjamini-hochberg', 'benjamini-yekutieli'), lty=1:6, col=1:6, lwd=3, bg='white')
}
#' 
#' <BR><BR>
#' lower<x<upper
