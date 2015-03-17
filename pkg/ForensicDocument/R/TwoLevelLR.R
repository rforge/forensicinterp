#' Two Level Multivariate Likelihood ratio
#' 
#' Computes a two level multivariate likelihood ratio.
#' 
#' @details This methodology has been developed by \emph{Bozza et al. (2008)} for the assesment of evidence through the derivation of a likelihood ratio for multivariate data. It allows to take into account the correlation between variables, and the non-constant variability within sources.
#' 
#' In the context of handwritten expertise suppose that: (i) an anonymous letter (\emph{i.e.} the \emph{questioned document}) is available for comparative analysis, and (ii) written material from a suspect is selected for comparative purposes (\emph{i.e.} the \emph{reference document}. For the compuation of the likelihood ratio, we consider the following propositions of interest:
#' \itemize{
#'  \item{\eqn{H_p}}{the suspect is the author of the questioned document;}
#'  \item{\eqn{H_d}}{the suspect is not the author of the questioned document - a random person wrote the document.}
#' }
#' 
#' @param data1 measurements from the 'reference' material.
#' @param data2 measurements from the 'questioned' material.
#' @param background background parameters for the overall and group means, within- and between group variances and covariances matrice.
#' @param n.iter number of MCMC iterations. Default is \code{10000}.
#' @param n.burnin number of burn-in iterations. Default is \code{1000}.
#' @param nw degrees of freedom for the inverse Wishart distribution. Considering p variables in the data, nw must be > 2*p+4.
#' @return The value of the log likelihood ratio (\eqn{log(LR)}).
#' @seealso \code{\link{TwoLevelLR_Background}}
#' @references Bozza S, Taroni F, Marquis R and Schmittbuhl M (2008). "Probabilistic evaluation of handwriting evidence: likelihood ratio for authorship." \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{57} (3), pp. 329-341.
#' @author Silvia Bozza \cr Alexandre Thiery 
#' @examples 
#' \dontrun{
#' ## Example where Hp is true. 
#' ## That is, when the suspect is the author of the questioned document
#' 
#' ## DATASET: We use the `characterO` dataset. It contains the extracted Fourier (`n.fourier = 4`) parameters from 554 handwritten character loops, written by 11 writers. For more information, see `?characterO`.
#' data(characterO)
#' 
#' ## PARAMETERS:
#' # data1: (reference document), the first 23 characters of writer 1
#' # data1: (questioned document), the last 23 characters of writer 1
#' # background: data from the remaining 10 writers
#' # n.iter: 11000
#' # n.burnin: 1000
#' # nw: 50
#' 
#' 
#' # data1, data2
#' data_reference = subset(characterO$measurements[,-1], subset = (characterO$info$writer == 1))[1:23,]
#' data_questioned = subset(characterO$measurements[,-1], subset = (characterO$info$writer == 1))[-(1:23),]
#' 
#' # background
#' subset = characterO$info$writer != 1
#' data_back = subset(characterO$measurements[,-1], subset = subset)
#' background = TwoLevelLR_Background(data_back, fac = as.factor(characterO$info$writer[subset]))
#' 
#' # others
#' n.iter = 11000
#' n.burnin = 1000
#' nw  = 50
#' 
#' # compute LLR
#' LLR = TwoLevelLR(data1 = data_reference,  
#'                  data2 = data_questioned,
#'                  background = background, 
#'                  n.iter = n.iter, n.burnin = n.burnin,
#'                  nw = nw)
#' LLR
#' }
#' @export
TwoLevelLR = function(data1, data2, background,
                      n.iter = 1100, n.burnin = 100, nw)
{
    if(missing(data1)) stop("data1 is missing.")
    if(missing(data2)) stop("data2 is missing.")
    if(missing(background)) stop("background is missing.")
    if(missing(nw)) stop("nw (degrees of freedom) is missing.")
    
    
    # data1
    if(!is.matrix(data1)) stop("data1 must be a matrix.")
    if(nrow(data1) == 0) stop("no rows in data1.")
    if(sum(is.na(data1)) > 0) stop("data1 contains 'NA' value.")
    # data2
    if(!is.matrix(data1)) stop("data2 must be a matrix!")
    if(nrow(data2) == 0) stop("no rows in data2")
    if(sum(is.na(data2)) > 0) stop("data2 contains 'NA' value.")
    # data1 & data2
    if(ncol(data1) != ncol(data2)) stop("data1 and data2 columns don't match.")
    
    # n.iter, ...
    if(!is.numeric(n.iter)) stop("n.iter is not numeric.")
    if(!is.numeric(n.burnin)) stop("n.burnin is not numeric.")
    
    n.iter <- round(n.iter,0)
    n.burnin <- round(n.burnin)
    
    if(n.iter < 0) stop("n.iter is negative.")
    if(n.iter == 0) stop("n.iter is set to 0.")
    if(n.burnin < 0) stop("n.burnin is negative.")
    
    if(n.iter - n.burnin <= 0) stop("n.burnin > n.iter.")
    
    # p (number of variables)
    p = ncol(data1)
    
    # nw
    if(!is.numeric(nw)) stop("nw is not numeric.")
    nw <- round(nw,0)
    if(nw <= 2*p+4) stop("nw must be > 2*p+4.")
    
    
    QR = rbind(data1, data2) 
    WB = background
    
    
    ## get parameters
    # within and between covariance
    W = WB$W
    B = WB$B
    # groups mean
    group.means = WB$group.means
    # overall mean
    mu = t(WB$all.means)
    #
    invB=solve(B)
    invBmu=invB%*%t(mu)
    logdetB=log(det(B))
    # theta
    theta.1=matrix(mvtnorm::rmvnorm(1,mu,B),nrow=1)
    theta.2=theta.1
    
    # U 
    U=W*(nw-2*p-2)
    logdetU=log(det(U))
    # W
    W.1= MCMCpack::riwish(nw,U)
    IW.1=solve(W.1)
    IW.2=IW.1
    
    ## compute
    L_hp = TwoLevelLR_numerator( QR, p, nw, n.iter, n.burnin, mu, B, logdetB, invB, invBmu, U, logdetU, IW.1)
    L_hd = TwoLevelLR_denominator(data1, data2, p, nw, n.iter, n.burnin, mu, B, logdetB, invB, invBmu, U, logdetU, IW.1, IW.2)
    
    LR = L_hp - L_hd
    LR
}