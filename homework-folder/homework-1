#' # Problem 1
#'
#' ## Part A
#' From the response we learned that the sum is 101 not 100 this happened because a lot of times
#' when we are dealing with percentages on table, we round up the number. The result have margin error of 3, which indicates that each of the answer 
#' 
#' ## Part B
plot1 <- c(0.47, 0.52, 0.02)
barplot(plot1, xlab = 'Responses', 
        ylab = 'Percentages',
        main = 'Bar Chart for responses',
        names.arg = c('Federal Goverment', 'Each State', 'Unsure/No Answer'),
        col = c('Blue', 'Red', 'Black')
        )

pie(plot1,
    label = plot1,
    main = 'Pie chart for responses',
    col = c('Blue', 'Green', 'Yellow'))

# adding legend for the pie chart
legend('topright',
       c('Fedeal Goverment', 'Each State', 'Unsure/ No Answer'),
       fill = c('Blue', 'Green', 'Yellow')
       )

#' ## Part C
#' Null hypothesis : The response has no prefrence between Federal Goverment and Each State
#' Alternative Hypothesis : The response favor Each State better than Federal Goverment.
observed_test_statistic <- abs(0.5 - 0.52)
simulated_proportion <- rmultinom(n = 10000, size = 1644, prob = c(.5, .5))
simulated_test_statistic <- abs(simulated_proportion[2,] / 1644 - 0.5)
hist(simulated_test_statistic)
abline(v=observed_test_statistic, col='red')

#' now compute p-value
p_val <- mean(simulated_test_statistic >= observed_test_statistic)
print(p_val)
#' Using confidence level of 0.05 we didn’t have enough evidence to reject the null hypothesis
#' 
#' <BR><BR>
#'
#' # Problem 2
#' test statistic is asymptotically X^2 distributed. If there's an unusual regularity of asssumptions it is said tne d(y) is to asymptotically X^2 distributed with S - 1 as their degrees of freedom.
#' Hence, it is accrute when the H0 holds true. 
#' 
#' Drawing histogram
a <- 1
n <- 10
b <- 5
x <- 100
range <- seq(a, b, 0.01)
y <- dunif(range, a,b)
samples <- ceiling(runif(n=10, min=1, max=b+1)) - 1
setdiffn <- ceiling(runif(n = 100, min = 1, max = 100 * 5)) - 1
par(mfrow=c(3,1))
hist(range, y , ylim = c(5, 10, 20, 50))
hist(range, y, ylim = c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000))
hist(range)
#'
#'For the first historgram it didn't have any bin as it show the samples of raw data
#'
#'For the second histogram it shows the range as the bin of 10
#'
#'The third shows range as the bin is set of 100
#'
#' <BR><BR>
#'
#' # Problem 3
#'
#' ## Part A
school_admission <- read.table(url("https://media.pearsoncmg.com/aw/aw_deveaux_stats_2/activstats/text/Ch03_Magnet_schools.txt"),
                               sep='\t',
                               header=TRUE)
ethnicity_admission <- table(school_admission$Ethnicity, school_admission$Admission.Decision)
ethnicity_admission

#' now produce the barplot for evey ethnicity
barplot(t(ethnicity_admission), col=c('red', 'green', 'blue'), main='Admission by ethnicity')
legend('topleft', colnames(ethnicity_admission), fill=c('red', 'green', 'blue'))

#' ## Part B
#' Null hypothesis : There is no correlation between ethnicity and admission as the both are independent from each other
#' 
#' Alternative hypothesis : There is a correlation between ethnicity and admission as they both depend on each other to produce an outcome
#' 
#' Making the observed table :
ethnicity_admission
#' Making the expected table :
chisq.test(ethnicity_admission)$expected
chisq.test(ethnicity_admission)

#' From the test we learned that we failed to reject the null hypothesis as the p value is smaller than 0.05. Hence, there is no correlation between ethnicity and admission. 

#' <BR><BR>
