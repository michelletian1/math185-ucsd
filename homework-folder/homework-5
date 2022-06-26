#' # Problem 1
#' 
dat = read.table("./father_son_Pearson.txt", header = TRUE)
str(dat)
attach(dat)

ks.indep.stat <- function(x, y){
  cdf_x = sort(x)
  cdf_y = sort(y)
  
  i = 1
  j = 1
  d = 0.0
  fnx = 0.0
  fny = 0.0
  lx = length(cdf_x)
  ly = length(cdf_y)
  while (i < lx && j < ly){
    dx = cdf_x[i]
    dy = cdf_y[j]
    if (dx <= dy){
      i = i + 1
      fnx = i / lx
    }
    if (dy <= dx){
      j = j + 1
      fny = j / ly
    }
    dist = abs(fny - fnx)
    if (dist > d){
      d = dist
    }
  }
  return(dist)
}
ks.indep.test <- function(x, y, B=1e4){
  stat_obs = ks.indep.stat(x, y)
  z = c(x, y)
  m = length(x)
  stat_test = numeric(B)
  for (b in 1:B) {
    z_perm = sample(z) # permute the combined sample
    x_perm = z_perm[1:m]
    y_perm = z_perm[-(1:m)]
    stat_test[b] = ks.indep.stat(x_perm, y_perm)
  }
  pval = (sum(abs(stat_test) >= abs(stat_obs))+1)/(B+1)
  return(pval)
}
ks.indep.test(dat$Father, dat$Son)
#' <BR><BR>
#' 
#' # Problem 2
#' 
fpath_c = paste0(getwd(), "/curry.txt")
curry_load = read.table(fpath_c, header = FALSE, sep = "", fill = TRUE)
fpath_d = paste0(getwd(), "/durant.txt")
durant_load = read.table(fpath_d, header = FALSE, sep = "", fill = TRUE)

durant_mat = durant_load[,c(1:3, 10)]
curry_mat = curry_load[,c(1:3, 10)]

both_mat = rbind(curry_mat, durant_mat)
common_mat = matrix(0, nrow = (79+62), ncol = 2)
common_mat[1:79, 1] = curry_mat[, 4]
common_mat[(79+1):(79+62), 2] = durant_mat[, 4]

remove_ind = numeric(0) # the row indices to remove due to same game
for (j in 1:79){ # curry has 79 games 
  for (i in (79+1):(79+62)){ # durant has 62 games 
    if ( all( both_mat[j, 1:3] == both_mat[i, 1:3]) ){
      remove_ind = c(remove_ind, i)
      common_mat[j, 2] = common_mat[i, 2]
    }
  }
}
common_mat = common_mat[-remove_ind, ]

colnames(common_mat) = c('curry', 'durant')
dat = data.frame( common_mat ) 
attach(dat)

#' All:
boxplot(dat, main="All Games of Curry and Durant 2016-2017", xlab = "player", ylab = "score")

plot(curry, durant, pch=16, main="All Scores During 2016-2017", xlab="curry", ylab="durant")

cor.test(curry, durant, method='spearman')
ks.indep.test(curry, durant, B=1e4)

#' From the boxplot we can see that durant have more points rather than steph curry, while the median of curry and durant is not that far away.
#' 
#' Without 0 in either:
curry_0 = curry[(curry!=0 & durant!=0)]
durant_0 = durant[(curry!=0 & durant!=0)]

boxplot(data.frame(curry_0, durant_0), main="All Games of Curry and Durant 2016-2017", xlab = "player", ylab = "score")

plot(curry_0, durant_0, pch=16, main="All Scores During 2016-2017", xlab="curry", ylab="durant")

cor.test(curry_0, durant_0, method='spearman')
ks.indep.test(curry_0, durant_0, B=1e4)

#'The zero effects the results of the graph. If we don't include the zero we can see that the gap of points between Kevin Durant and Steph Curry is smaller. 
#' 
#' <BR><BR>
