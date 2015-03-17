#' The multivariate gamma function.
#' 
#' Compute background parameters (overall and group means, within- and between group variances and covariances matrices) for the two level likelihood ratio.
#' @param a a non-negative numeric.
#' @param p a non-negative numeric.
#' @param log a logical value. If \code{TRUE}, multivariate gamma function is given as \code{log}. Default is \code{log=TRUE}.
#' @return A numeric value.
#' @keywords internal
#' @author Alexandre Thiery 
MultivariateGamma <- function(a, p, log = TRUE) 
{
    res = sum(lgamma( a + (-p-1:p)/2 ) ) + log(pi)*(p*(p-1)/4) 
    if(log)
        return(res)
    return(exp(res))
}
