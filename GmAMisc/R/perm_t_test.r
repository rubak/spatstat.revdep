#' R function for permutation-based t-test
#'
#' The function allows to perform a permutation-based t-test to compare two independent groups. The
#' test's results are graphically displayed within the returned chart.
#'
#' A permutation t-test proves useful when the assumption of 'regular' t-test are not met. In
#' particular, when the two groups being compared show a very skewed distribution, and when the
#' sample sizes are very unbalanced.\cr
#'
#' "The permutation test is useful even if we plan to use the two-sample t test. Rather than relying
#' on Normal quantile plots of the two samples and the central limit theorem, we can directly check
#' the Normality of the sampling distribution by looking at the permutation distribution.
#' Permutation tests provide a “gold standard” for assessing two-sample t tests. If the two P-values
#' differ considerably, it usually indicates that the conditions for the two-sample t don’t hold for
#' these data. Because permutation tests give accurate P-values even when the sampling distribution
#' is skewed, they are often used when accuracy is very important." (Moore, McCabe, Craig,
#' "Introduction to the Practice of Statistics", New York: W. H. Freeman and Company, 2009).\cr
#'
#' @param data Dataframe containing the data.
#' @param format It takes "long" if the data are arranged in two columns, with the left-hand one
#'   containing the values, and the righ-hand one containing a grouping variable; it takes "short"
#'   if the values of the two groups being compared are stored in two different adjacent columns.
#' @param sample1.lab Label for the first sample being tested (default: smpl 1).
#' @param sample2.lab Label for the first sample being tested (default: smpl 2).
#' @param B Desired number of permutations (set at 999 by default).
#'
#' @return The frequency histogram returned by the function displays the distribution of the
#' permuted mean difference between the two samples; a solid dot indicates the observed mean
#' difference, while an hollow dot represents the mean of the permuted differences.
#' Two dashed blue lines indicates the 0.025 and 0.975 percentile of the permuted
#' distribution. A rug plot at the bottom histgram indicates the individual permuted mean
#' differences. At the bottom of the chart, some information are displayed.
#' In particular, the observed mean difference and the permuted p-values are reported.
#' In the last row, the result of the regular (parametric) t-test (both assuming and
#' not assuming equal variances) is reported to allow users to compare the outcome of
#' these different versions of the test.
#'
#' @keywords perm.t.test
#'
#' @export
#'
#' @importFrom utils unstack
#' @importFrom stats quantile
#' @importFrom graphics hist polygon
#'
#' @examples
#' #load the 'resample' package which stores a toy dataset
#' library(resample)
#'
#' #load the 'Verizon' dataset
#' data("Verizon")
#'
#' #performs the permutation-based t-test using 199 permutations
#' perm.t.test(Verizon, format="long", B=199)
#'
perm.t.test <- function (data,format,sample1.lab=NULL,sample2.lab=NULL,B=999){
  options(scipen=999)

  if (format=="long") {
    unstacked.data <- utils::unstack(data)
    sample1 <- unstacked.data[[1]]
    sample2 <- unstacked.data[[2]]
  } else {
    sample1 <- data[,1]
    sample2 <- data[,2]
  }

  #if the parameters for samples' label are NULL, assign default labels, otherwise do nothing
  if (is.null(sample1.lab)==TRUE) {
    sample1.lab <- "smpl 1"
    sample2.lab <- "smpl 2"
  } else {}

  #get some statistics for the two samples
  n1 <- length(sample1)
  n2 <- length(sample2)
  mean1 <- round(mean(sample1), 2)
  mean2 <- round(mean(sample2),2)
  error1 <- qnorm(0.975)*sd(sample1)/sqrt(n1)
  error2 <- qnorm(0.975)*sd(sample2)/sqrt(n2)
  sample1_lci <- round(mean1 - error1,2)
  sample1_uci <- round(mean1 + error1,2)
  sample2_lci <- round(mean2 - error2,2)
  sample2_uci <- round(mean2 + error2,2)

  #get regular t-test results (equal variance)
  p.equal.var <- round(t.test(sample1, sample2, var.equal=TRUE)$p.value, 4)

  #get regular t-test results (unequal variance)
  p.unequal.var <- round(t.test(sample1, sample2, var.equal=FALSE)$p.value, 4)

  #start permutation procedures
  pooledData <- c(sample1, sample2)
  size.sample1 <- length(sample1)
  size.sample2 <- length(sample2)
  size.pooled <- size.sample1 + size.sample2
  nIter <- B
  meanDiff <- numeric(nIter + 1)
  meanDiff[1] <- round(mean1 - mean2, digits=2)

  #set the progress bar to be used inside the loop
  pb <- txtProgressBar(min = 0, max = B, style = 3)

  for(i in 2:B){
    index <- sample(1:size.pooled, size=size.sample1, replace=F)
    sample1.perm <- pooledData[index]
    sample2.perm <- pooledData[-index]
    meanDiff[i] <- mean(sample1.perm) - mean(sample2.perm)
    setTxtProgressBar(pb, i)
  }

  p.lowertail <- (1 + sum (meanDiff[-1] < meanDiff[1])) / (1 + B)
  p.uppertail <- (1 + sum (meanDiff[-1]  > meanDiff[1])) / (1 + B)
  two.sided.p <- 2 * min(p.lowertail, p.uppertail)


  graphics::hist(meanDiff, main=paste0("Distribution of permuted mean differences\n(number of permutations: ", B, ")"),
                 xlab="",
                 sub=paste0(sample1.lab, " (n: ", n1,") (95% CI lower bound., mean, 95% CI upper bound.): ", sample1_lci, ", ", mean1, ", ", sample1_uci, "\n", sample2.lab, " (n: ", n2,") (95% CI lower bound., mean, 95% CI upper bound.): ", sample2_lci, ", ", mean2, ", ", sample2_uci,"\nobs. mean diff. (solid dot): ", meanDiff[1],"; perm. p.value mean ", sample1.lab, " < ", sample2.lab, ": ", round(p.lowertail, 4), "; perm. p.value mean ", sample1.lab, " > ", sample2.lab, ": ", round(p.uppertail,4), "; perm. p.value (2-sided): ", round(two.sided.p,4),"\nregular t-test p-values (2-sided): ", round(p.equal.var,4)," (equal variance); ",round(p.unequal.var,4), " (unequal variance)"),
                 cex.main=0.85,
                 cex.sub=0.70)

  rug(meanDiff, col="#0000FF")

  abline(v=stats::quantile(meanDiff, 0.025), lty=2, col="blue")
  abline(v=stats::quantile(meanDiff, 0.975), lty=2, col="blue")

  points(x=meanDiff[1], y=0, pch=20, col = "black")
  points(x=mean(meanDiff[-1]), y=0, pch=1, col="black")
}
