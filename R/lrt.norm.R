# likelihood ratio test of the normal distribution with known variance
# returns the p-value of the test
# H0: mu = a, H1: mu != a
#' @export
lrt.norm <- function(x, mu=0, sd=1) {
    m. <- mean(x)
    dn0 <- sum( dnorm(x, mean=mu, sd=sd, log=TRUE) )
    dn1 <- sum( dnorm(x, mean=m., sd=sd, log=TRUE) )
    return( 2 * (dn1 - dn0) )
}

# One-sided version of the test
# returns the (uncorrected) p-values
# H0: mu = a, H1: mu > a
#' @export
lrt.norm.onesided <- function(x, mu=0, sd=1) {
    m. <- mean(x)
    dn0 <- sum( dnorm(x, mean=mu, sd=sd, log=TRUE) )
    if (m. > mu) {
        dn1 <- sum( dnorm(x, mean=m., sd=sd, log=TRUE) )
    } else {
        dn1 <- dn0
    }
    return( 2 * (dn1 - dn0) )
}
