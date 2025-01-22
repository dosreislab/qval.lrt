#' Calculate q-values from LRT statistics
#'
#' @param lrt numeric, vector of LRT values
#' @param df numeric, degrees of freedom
#' @param one.sided logical, whether the test is one-sided
#'
#' @details
#' P-values are calculated from the \code{lrt} values using a chi-squared
#' distribution with \code{df} degrees of freedom. The q-values, which control the
#' false-discovery rate, are then calculated using the \code{qvalue} function from
#' the \code{qvalue} package.
#'
#' If \code{one.sided = FALSE}, q-values are calculated using Benjamini and
#' Hochber (1995, BH) and Storey and Tibshirani (2003, ST) methods in the usual
#' way.
#'
#' If \code{one.sided = TRUE}, p-values are first divided by two (Self and
#' Liang, 1987), and the BH and ST q-values are calculated. In the case of ST,
#' the proportion of true null tests is estimated as \eqn{\pi_0 = \min(1,
#' 2n_0/N)}, where \eqn{n_0} is the number of tests with LRT equal to zero
#' (p-values equal to one) and \eqn{N} is the total number of tests, and then
#' used to obtain the q-values. An alternative robust form of the ST method
#' (STR) is also calculated by first removing all LRT values equal to zero, and
#' then calculating the q-values on the remaining unhalved p-values.
#'
#' When an LRT test is one-sided, half of true null tests are expected to
#' generate LRT values equal to zero (Self and Liang, 1987), and Storey and
#' Tibhirani's smoother estimate of \eqn{\pi_0} on the histogram of p-values at
#' 1 cannot be used, and a corrected estimate is needed. If the one-sided test
#' is misspecified (see examples), so that the alternative model contains true
#' parameter values in the disallowed parameter area of the test, then the
#' proportion of LRT values equal to zero can be much larger than 50\%. In this
#' case, the STR q-values are robust, while ST q-values are more conservative.
#'
#' @references
#' Benjamini, Y. and Hochberg, Y. 1995, Controlling the false discovery rate: a
#' practical and powerful approach to multiple testing. Journal of the Royal
#' Statistical Society: Series B (Methodological) 57:289–300.
#'
#' Self, S. G. and Liang, K.-Y. 1987, Asymptotic properties of maximum
#' likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. Journal of the American Statistical Association 82:605–610.
#'
#' Storey, J. D. and Tibshirani, R. 2003, Statistical significance for
#' genomewide studies. Proceedings of the National Academy of Sciences
#' 100:9440–9445.
#'
#' @return
#' A list with elements \code{qvals.bh}, \code{qvals.st}, and \code{pvals}, with
#' the q-values estimated under the BH and ST methods (as output by
#' \code{qvalue::qvalue}) respectively, and the p-values. If \code{one.sided =
#' TRUE}, an additional element, \code{qvals.stR}, with a robust version of ST
#' is also returned. If \code{one.sided = TRUE}, the returned p-values are
#' halved.
#'
#' @examples
#' par(mfcol=c(2,4));
#' par(mai=c(.3,.1,0,0), omi=c(.1,.2,.1,.1)) # b, l, t, r
#' pch=19; cex=.5; ylim=c(0, 5); div <- 50
#'
#' # H0: mu = 0; H1: mu != 0
#' # *** H0 is always true: ***
#' N <- 1e4L; p0 <- 1; n <- 1e2L
#' x0 <- matrix(rnorm(round(n * N), 0, 1), ncol=n)
#' lrt <- apply(x0, 1, lrt.norm)
#' qvals <- qval.lrt(lrt, df=1, one.sided=FALSE)
#' pvals <- qvals$pvals
#' i <- rank(pvals, ties="max")
#' # plots
#' hist(pvals, n=div, ylim=ylim, main="", prob=TRUE); abline(1, 0, lty=2)
#' plot(i/N, pvals, pch=pch, cex=cex); abline(0, 1)
#' all.equal(p.adjust(qvals$pvals, "BH"), qvals$qvals.bh$qvalues)
#'
#' # H0: mu = 0; H1: mu != 0
#' # *** H0 is not always true: ***
#' p0 <- .9; p1 <- 1 - p0; mu <- .5
#' n0 <- round(p0 * N); n1 <- round(p1 * N)
#' # simulate a mixture of positive and negative tests:
#' x0 <- matrix(rnorm(round(n * n0), 0, 1), ncol=n)
#' x1 <- matrix(rnorm(round(n * n1), mu, 1), ncol=n)
#' lrt <- apply(rbind(x1, x0), 1, lrt.norm)
#' qvals <- qval.lrt(lrt, df=1, one.sided=FALSE)
#' pvals <- qvals$pvals
#' i <- rank(pvals, ties="max")
#' hist(pvals, n=div, ylim=ylim, main="", prob=TRUE)
#' abline(1, 0, lty=2); abline(qvals$qvals.st$pi0, 0, lty=2, col="red")
#' plot(i/N, pvals, pch=pch, cex=cex); abline(0, 1)
#'
#' # Observed FDR:
#' 1 - sum(qvals$qvals.bh$qvalues[1:n1] < .05) / sum(qvals$qvals.bh$qvalues < .05)
#' 1 - sum(qvals$qvals.st$qvalues[1:n1] < .05) / sum(qvals$qvals.st$qvalues < .05)
#'
#' # H0: mu = 0; H1: mu > 0
#' # *** H0 is not always true: (test is one-sided)***
#' lrt <- apply(rbind(x1, x0), 1, lrt.norm.onesided)
#' qvals <- qval.lrt(lrt, df=1, one.sided=TRUE)
#' pvals <- qvals$pvals
#' i <- rank(pvals, ties="max")
#' hist(pvals, n=div, ylim=ylim, main="", prob=TRUE)
#' abline(1, 0, lty=2); abline(qvals$qvals.st$pi0, 0, lty=2, col="red")
#' plot(i/N, pvals, pch=pch, cex=cex); abline(0, 1)
#'
#' # Observed FDR:
#' 1 - sum(qvals$qvals.bh$qvalues[1:n1] < .05) / sum(qvals$qvals.bh$qvalues < .05)
#' 1 - sum(qvals$qvals.st$qvalues[1:n1] < .05) / sum(qvals$qvals.st$qvalues < .05)
#' 1 - sum(qvals$qvals.stR$qvalues[1:n1] < .05) / sum(qvals$qvals.stR$qvalues < .05)
#'
#' # H0: mu = 0; H1: mu > 0
#' # *** H0 is not always true: (test is one-sided) ***
#' # *** H0 is misspecified in some tests ***
#' p0 <- p0. <- .45; p1 <- 1 - p0 - p0.; mu. <- -.1
#' n0 <- round(p0 * N)
#' n0. <- round(p0. * N)
#' n1 <- round(p1 * N)
#' x0 <- matrix(rnorm(n0 * n, 0, 1), ncol=n)
#' x0. <- matrix(rnorm(n0. * n, mu., 1), ncol=n)
#' x1 <- matrix(rnorm(n1 * n, mu, 1), ncol=n)
#' lrt <- apply(rbind(x1, x0, x0.), 1, lrt.norm.onesided)
#' qvals <- qval.lrt(lrt, df=1, one.sided=TRUE)
#' pvals <- qvals$pvals
#' i <- rank(pvals, ties="max")
#' hist(pvals, n=div, ylim=ylim, main="", prob=TRUE)
#' plot(i/N, pvals, pch=pch, cex=cex); abline(0, 1)
#'
#' # Observed FDR:
#' 1 - sum(qvals$qvals.bh$qvalues[1:n1] < .05) / sum(qvals$qvals.bh$qvalues < .05)
#' 1 - sum(qvals$qvals.st$qvalues[1:n1] < .05) / sum(qvals$qvals.st$qvalues < .05)
#' 1 - sum(qvals$qvals.stR$qvalues[1:n1] < .05) / sum(qvals$qvals.stR$qvalues < .05)
#'
#' @export
qval.lrt <- function(lrt, df, one.sided) {
    pvals <- pchisq(q=lrt, df=df, lower.tail = FALSE)
    if (!one.sided) {
        qvals.bh <- qvalue::qvalue(pvals, pi0 = 1)
        qvals.st <- qvalue::qvalue(pvals)
        obj <- list(pvals=pvals, qvals.bh=qvals.bh, qvals.st=qvals.st)
    }
    if (one.sided) {
        p1 <- pvals == 1
        n0 <- sum(p1)
        N <- length(pvals)
        pi0 = min(1, 2 * n0 / N)
        pvals2 <- pvals; pvals2[!p1] <- pvals2[!p1] / 2
        qvals.bh <- qvalue::qvalue(pvals2, pi0 = 1)
        qvals.st <- qvalue::qvalue(pvals2, pi0 = pi0)
        qvals.stR <- qvalue::qvalue(pvals[!p1])
        obj <- list(pvals=pvals2, qvals.bh=qvals.bh, qvals.st=qvals.st,
                    qvals.stR=qvals.stR)
    }
    return (obj)
}
