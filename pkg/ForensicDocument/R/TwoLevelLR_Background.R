#' Compute background parameters for the two level likelihood ratio.
#' 
#' Compute background parameters (overall and group means, within- and between group variances and covariances matrices) for the two level likelihood ratio.
#' @param data a \eqn{n x p} numeric matrix, with \eqn{p \geq 2}. The background data containing \eqn{n} measurements on \eqn{p} variables.
#' @param fac a factor of length \eqn{p}, indicating the 'population' of each measurement.
#' @return A list containing the overall and group means (\code{all.means} and \code{group.means}), within- and between group variances and covariances matrices (\code{W} and \code{B}).
#' @seealso \code{\link{TwoLevelLR}}
#' @references Bozza S, Taroni F, Marquis R and Schmittbuhl M (2008). "Probabilistic evaluation of handwriting evidence: likelihood ratio for authorship." \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{57} (3), pp. 329-341.
#' @author Silvia Bozza \cr Alexandre Thiery 
#' @export
TwoLevelLR_Background = function(data, fac)
{
    if(missing(data)) stop("data is missing.")
    if(missing(fac)) stop("fac is missing.")
    
    # data
    if(!is.matrix(data)) stop("data must be a matrix.")
    if(nrow(data) == 0) stop("no rows (i.e. measurements) in data")
    if(ncol(data) < 2) stop("number of columns in data (i.e. variables) is < 2.")
    if(sum(is.na(data)) > 0) stop("data contains 'NA' value.")
    
    # fac
    if(!is.factor(fac)) stop("fac must be a factor")
    if(length(fac) != nrow(data)) stop("data and fac dimensions don't match.")
    ## remove factor (or population) that have only one or zero measurement
    fac_count = sapply(levels(fac), function(x) sum(fac == x))
    fac_subset = (fac %in% levels(fac)[fac_count < 2])
    if(sum(fac_subset) > 0)
    {
        data = subset(data, !fac_subset)
        fac = subset(fac, !fac_subset)
        warning(sprintf("%d population(s) have been removed: the number of measurements in those are < 2:\n%s.",
                        sum(fac_subset), paste(levels(fac)[fac_count < 2], collapse = ", ")))
    }
    fac = droplevels(fac)
    if(nrow(data) == 0) stop("no rows (i.e. measurements) in data")
    
    lev = levels(fac)
    nlev = nlevels(fac)
    p = ncol(data)
   
    
    # overall population mean
    all.means = colMeans(data)    
    # which rows refer to repeated measurements on the same item
    group.index = lapply(lev, function(x) which(fac == x))
    # assign values to the group means for the population
    group.means = t(sapply(group.index, function(x) colMeans(data[x,])))    
    
    S = sapply(1:nlev, function(x) {
        as.numeric((group.means[x,] - all.means) %*% t(group.means[x,] - all.means)* length(group.index[[x]]))
    })
    Sw = sapply(1:nlev, function(x) {
        rowSums(apply(data[group.index[[x]],], 1, function(y) {
            as.numeric((y-group.means[x,])%*%t(y-group.means[x,]))
        }))
    })
    S = matrix(rowSums(S), ncol = p)
    Sw = matrix(rowSums(Sw), ncol = p)
    
    # convert matrix Sw to matrix U
    W = Sw / (nrow(data) - nlev)    
    # convert matrix S to matrix C
    B = S/(nlev-1)
    
    result = list(group.means,all.means,W,B)
    names(result) = c('group.means','all.means','W','B')
    return(result)
}
