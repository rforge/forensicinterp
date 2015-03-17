#' @rdname TwoLevelLR_func
#' @name TwoLevelLR_func
#' @aliases TwoLevelLR_func
#' @aliases TwoLevelLR_denominator
#' @title Functions for the Two Level Multivariate Likelihood ratio numerator and denominator.
#' @seealso \code{\link{TwoLevelLR}}
#' @references Bozza S, Taroni F, Marquis R and Schmittbuhl M (2008). "Probabilistic evaluation of handwriting evidence: likelihood ratio for authorship." \emph{Journal of the Royal Statistical Society: Series C (Applied Statistics)}, \strong{57} (3), pp. 329-341.
#' @keywords internal
#' @author Silvia Bozza \cr Alexandre Thiery 
#' @import mvtnorm
#' @import MCMCpack
TwoLevelLR_denominator =function(data.1, data.2, p, nw, n.iter, n.burnin, mu, B, logdetB, invB, invBmu, U, logdetU, IW.1, IW.2) {
    
    ### checks ###
    data.1 = as.matrix(data.1)
    data.2 = as.matrix(data.2)
    
    ### COMPUTATION ###
    # number of samples
    nr.1=dim(data.1)[1]
    nr.2=dim(data.2)[1]
    # nw* (Inverse-Wishart degrees of freedom)
    nwstar.1=nr.1 + nw       
    nwstar.2=nr.2 + nw
    # ӯ1 and ӯ2
    bary.1= colMeans(data.1)
    bary.2= colMeans(data.2)
    # S1 = sum_i ( sum_j (y1_ij - ӯ1)(y1_ij - ӯ1)' )
    # S2 = sum_i ( sum_j (y2_ij - ӯ2)(y2_ij - ӯ2)' )
    S.1 = matrix(rowSums(apply(data.1, 1, function(x) as.numeric((x-bary.1)%*%t(x-bary.1)))), ncol = p) 
    S.2 = matrix(rowSums(apply(data.2, 1, function(x) as.numeric((x-bary.2)%*%t(x-bary.2)))), ncol = p) 
    
    ### SAMPLE n.iter ψ1 = (Ѳ1, W1) and ψ2 = (Ѳ2, W2)
    IW.temp.1 = IW.1
    IW.temp.2 = IW.2
    res1 = lapply(1:n.iter, function(i) { 
        # sample B1* and μ1*
        Bstar.1  = solve(invB+nr.1*IW.temp.1)
        mustar.1 = Bstar.1%*%(invBmu+nr.1*IW.temp.1%*%(bary.1))
        # sample Ѳ1
        theta.1  = matrix(mvtnorm::rmvnorm(1,mustar.1,Bstar.1),nrow=1)   
        # sample B2* and μ2*
        Bstar.2  = solve(invB+nr.2*IW.temp.2)
        mustar.2 = Bstar.2%*%(invBmu+nr.2*IW.temp.2%*%(bary.2))
        # sample Ѳ2
        theta.2  = matrix(mvtnorm::rmvnorm(1,mustar.2,Bstar.2),nrow=1)   
        
        # sample U1*
        Ustar.1 = nr.1*t(theta.1-bary.1)%*%(theta.1-bary.1)+U+S.1        
        # sample W1 and its inverse
        W.1     = MCMCpack::riwish(nwstar.1,Ustar.1) 
        IW.temp.1 <<- solve(W.1)
        # sample U2*
        Ustar.2 = nr.2*t(theta.2-bary.2)%*%(theta.2-bary.2)+U+S.2    
        # sample W2 and its inverse
        W.2     = MCMCpack::riwish(nwstar.2,Ustar.2)     
        IW.temp.2 <<- solve(W.2)
        
        list(theta.1, IW.temp.1, theta.2, IW.temp.2)
    })
    
    
    ### ESTIMATE THE (y1,y2) DENSITY
    ### l.f(y | ψ* , H2) = l.f(y1 | ψ1* , H2) + l.f(y2 | ψ2* , H2)
    logf = sapply(res1, function(x) {
        -(p*nr.1/2)*log(2*pi)+(nr.1/2)*log(det(x[[2]]))- # l.f(y1 | ψ1* , H2)
            sum(apply(data.1, 1, function(y) (y-x[[1]]) %*% x[[2]] %*% t(y-x[[1]]))/2) +
            -(p*nr.2/2)*log(2*pi)+(nr.2/2)*log(det(x[[4]]))-  # l.f(y2 | ψ2* , H2)
            sum(apply(data.2, 1, function(y) (y-x[[3]]) %*% x[[4]] %*% t(y-x[[3]]))/2)
    }) 
    ## index to maximize l.f(y | ψ* , H2)
    log.max = which(logf == max(logf))[1]
    ## Ѳ1*, W1*
    theta.1.star = res1[[log.max]][[1]]
    IW.1.star    = res1[[log.max]][[2]]
    ## Ѳ2*, W2*
    theta.2.star = res1[[log.max]][[3]]
    IW.2.star    = res1[[log.max]][[4]]
    ### l.f(y | ψ* , H2)
    logf.star  = logf[log.max]
    
    
    ### ESTIMATE THE PRIOR ORDINATES ( π(ψ* | y) )
    ### prior independence leads to posterior independence of π(Ѳi* | y) and π(Wi* | y)
    ### l.π(ψ* | y) = l.π(Ѳ1* | y) + l.π(W1* | y) ) + l.π(Ѳ2* | y) + l.π(W2* | y) )
    res2 = sapply((n.burnin):n.iter, function(x) {
        
        ## compute B1*(g) and μ1*(g)
        Bstar.1.g=solve(invB+nr.1*res1[[x]][[2]])
        mustar.1.g=Bstar.1.g%*%(invBmu+nr.1*res1[[x]][[2]]%*%(bary.1))         
        ##  π(Ѳ1* | y1, μ1*(g), B1*(g)) -- in log
        temp1 = mvtnorm::dmvnorm(theta.1.star, mustar.1.g,Bstar.1.g,log=TRUE)
        
        ## compute B2*(g) and μ2*(g)
        Bstar.2.g=solve(invB+nr.2*res1[[x]][[4]])
        mustar.2.g=Bstar.2.g%*%(invBmu+nr.2*res1[[x]][[4]]%*%(bary.2))          
        ##  π(Ѳ1* | y2, μ2*(g), B2*(g)) -- in log
        temp2 = mvtnorm::dmvnorm(theta.2.star, mustar.2.g,Bstar.2.g,log=TRUE)
        
        
        ## compute U1*(g)
        Ustar.1.g=nr.1* t(res1[[x]][[1]]-bary.1) %*% (res1[[x]][[1]]-bary.1)+U+S.1 
        ## π(W1* | y1, nw1*, U1*(g)) -- in log, without constant part
        temp3 = (nwstar.1-p-1)/2*(log(det(Ustar.1.g)))- sum(diag(IW.1.star%*%Ustar.1.g))/2             
        
        ## compute U2*(g)
        Ustar.2.g=nr.2* t(res1[[x]][[3]]-bary.2) %*% (res1[[x]][[3]]-bary.2)+U+S.2   
        ## π(W2* | y2, nw2*, U2*(g)) -- in log, without constant part
        temp4 = (nwstar.2-p-1)/2*(log(det(Ustar.2.g)))- sum(diag(IW.2.star%*%Ustar.2.g))/2       
        
        c(temp1, temp2, temp3, temp4)
    })
    n.iter=n.iter-n.burnin
    
    ## π(Ѳi* | yi) -- in log
    lpihat.theta1.star = log(sum(exp(res2[1,]))/n.iter)
    lpihat.theta2.star = log(sum(exp(res2[2,]))/n.iter)
    
    ## π(Wi* | yi) -- in log
    const1 = log(det(IW.1.star))*(nwstar.1)/2 - log(2)*(nwstar.1-p-1)*(p/2) - MultivariateGamma(nwstar.1/2, p)
    const2 = log(det(IW.2.star))*(nwstar.2)/2 - log(2)*(nwstar.2-p-1)*(p/2) - MultivariateGamma(nwstar.2/2, p)
    lpihat.W1.star = log(sum(exp(res2[3,]+ const1))/n.iter)
    lpihat.W2.star = log(sum(exp(res2[4,]+ const2))/n.iter)
    
    ###  l.π(ψ* | y) = l.π(Ѳ1* | y1)       +  l.π(W1* | y1)  + l.π(Ѳ2* | y2)       +  l.π(W2* | y2) )
    lpihat.psi.star  = lpihat.theta1.star + lpihat.W1.star + lpihat.theta2.star + lpihat.W2.star
    
    
    ### ESTIMATE THE POSTERIOR ORDINATES### 
    ### π(ψ* | y, H1) = π(Ѳ1* | y1, H2) x π(W1* | y1, H2) x π(Ѳ2* | y2, H2) x π(W2* | y2, H2)
    ## π(Ѳi* | yi, H1) -- in log
    lpi.theta1.star = mvtnorm::dmvnorm(x = theta.1.star, mean =  mu, sigma = solve(invB), log = TRUE)
    lpi.theta2.star = mvtnorm::dmvnorm(x = theta.2.star, mean =  mu, sigma = solve(invB), log = TRUE)
    ## (Wi* | yi, H1) -- in log
    const1 = log(det(IW.1.star))*(nw)/2 - log(2)*(nw-p-1)*(p/2) - MultivariateGamma(nw/2, p)
    const2 = log(det(IW.2.star))*(nw)/2 - log(2)*(nw-p-1)*(p/2) - MultivariateGamma(nw/2, p)
    lpi.W1.star = const1 +  (nw-p-1)/2*(logdetU)- sum(diag(IW.1.star%*%U))/2
    lpi.W2.star = const2 +  (nw-p-1)/2*(logdetU)- sum(diag(IW.2.star%*%U))/2
    ## π(ψ* | y, H1) -- in log
    ## l.π(ψ* | H1) = l.π(Ѳ1* | y1, H2) + l.π(W1* | y1, H2) + l.π(Ѳ2* | y2, H2) + l.π(W2* | y2, H2)
    lpi.psi.star    = lpi.theta1.star   + lpi.W1.star       + lpi.theta2.star   + lpi.W2.star
    
    ### COMPUTE marginal density
    ### l.m(y | H1) = l.f(y | ψ* , H1) + π(ψ* | H1)   - π(ψ* | y, H1)
    lmhatHp.den     = logf.star        + lpi.psi.star - lpihat.psi.star
    return(lmhatHp.den)
    
}
