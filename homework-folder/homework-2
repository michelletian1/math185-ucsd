#' # Problem 1
#' Part A
permChisqTest <- function(tab, B=1e4){
  tab_margins <- addmargins(tab)
  n_rows <- nrow(tab_margins)
  n_cols <- ncol(tab_margins)
  totals <- tab_margins['Sum', 'Sum']
  probabilities <- rep(1/length(tab), length(tab))
  expected <- totals * probabilities
  simulations <- matrix(sample.int((n_rows-1) * (n_cols-1), B * totals, replace = TRUE), nrow = totals)
  chi <- apply(simulations, 2L, function(x, E, k) {
    sum((table(factor(x, levels = 1L:k)) - E)^2 / E)}, E = array(expected), k = (n_rows-1) * (n_cols-1))
  statistic <- sum((tab - expected)^2/expected)
  print(chi)
  p_value = (1 + sum(chi >= statistic)) / (B + 1)
  p_value
}

#' Part B
school_admission <- read.table(url("https://media.pearsoncmg.com/aw/aw_deveaux_stats_2/activstats/text/Ch03_Magnet_schools.txt"),
                               sep='\t',
                               header=TRUE)
ethnicity_admission <- table(school_admission$Ethnicity, school_admission$Admission.Decision)
ethnicity_admission
begin_time <- proc.time()
permChisqTest(ethnicity_admission, B=1e6)
proc.time() - begin_time
begin_time <- proc.time()
chisq.test(ethnicity_admission, simulate.p.value=TRUE, B=1e6)
proc.time() - begin_time

#' <BR><BR>
#'
#' # Problem 2
#' We create the normal distribution function and histogram
problem_2_norm <- function(n){
  set.seed(2)
  x <- rnorm(n, mean=5, sd=1)
  median_x <- median(x)
  medians <- replicate(n = 1000, expr = {
    x_i <- sample(x, length(x), replace = T)
    median(x_i)
  })
  hist(x, prob=TRUE)
  lines(density(x))
  statistic = sum((medians - median_x)^2/median_x)
  statistic
}
problem_2_norm(10)
problem_2_norm(100)
problem_2_norm(10000)

problem_2_unif <- function(n){
  set.seed(2)
  x <- runif(n, min=0, max=10)
  median_x <- median(x)
  medians <- replicate(n = 1000, expr = {
    x_i <- sample(x, length(x), replace = T)
    median(x_i)
  })
  hist(x, prob=TRUE)
  lines(density(x))
  statistic = sum((medians - median_x)^2/median_x)
  statistic
}
problem_2_unif(10)
problem_2_unif(100)
problem_2_unif(10000)

#' <BR><BR>
