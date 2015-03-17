#' @rdname TwoLevelLR_func
#' @aliases TwoLevelLR_numerator
#' @keywords internal
TwoLevelLR_numerator =function( data, p, nw, n.iter, n.burnin, mu, B, logdetB, invB, invBmu, U, logdetU, IW) {
    
    ### checks ###
    data = as.matrix(data)
    
    ### COMPUTATION ###
    # number of samples
    nr=nrow(data)
    # nw* (Inverse-Wishart degrees of freedom)
    nwstar=nr+nw
    # ӯ
    bary =  colMeans(data)
    # S = sum_i ( sum_j (y_ij - ӯ)(y_ij - ӯ)' )
    S = matrix(rowSums(apply(data, 1, function(x) as.numeric((x-bary)%*%t(x-bary)))),  ncol = p) 
    
    
    ### SAMPLE n.iter ψ = (Ѳ, W)
    IW.temp = IW
    res1 = lapply(1:n.iter, function(i) {
        # sample B* and μ*
        Bstar  = solve(invB+nr*IW.temp)
        mustar = Bstar%*%(invBmu+nr*IW.temp%*%(bary))
        # sample Ѳ
        theta  = matrix(mvtnorm::rmvnorm(1,mustar,Bstar),nrow=1)   
        
        # sample U*
        Ustar = nr * t(theta-bary) %*% (theta-bary) + U + S    
        # sample W and its inverse
        W     = MCMCpack::riwish(nwstar,Ustar)  
        IW.temp <<- solve(W)
        
        list(theta, IW.temp)
    })
    
    
    ### ESTIMATE THE y DENSITY
    ### l.f(y | ψ* , H1)
    logf = sapply(res1, function(x) {
        logf=-(p*nr/2)*log(2*pi)+(nr/2)*log(det(x[[2]]))   
        logf-sum(apply(data, 1, function(y) (y-x[[1]]) %*% x[[2]] %*% t(y-x[[1]]))/2)
    })
    ## index to maximize l.f(y | ψ* , H1)
    log.max = which(logf == max(logf))[1]
    ## Ѳ*, W*
    theta.star = res1[[log.max]][[1]]
    IW.star    = res1[[log.max]][[2]]
    ### l.f(y | ψ* , H1)
    logf.star  = logf[log.max]
    
    
    ### ESTIMATE THE PRIOR ORDINATES ( π(ψ* | y) )
    ### prior independence leads to posterior independence of π(Ѳ* | y) and π(W* | y)
    ### l.π(ψ* | y) = l.π(Ѳ* | y) + l.π(W* | y) )
    res2 = sapply((n.burnin+1):n.iter, function(g) {
        
        ## compute B*(g) and μ*(g)
        Bstar.g = solve(invB+nr*res1[[g]][[2]])
        mustar.g = Bstar.g%*%(invBmu+nr*res1[[g]][[2]]%*%(bary))
        ##  π(Ѳ* | y, μ*(g), B*(g)) -- in log
        temp1 = mvtnorm::dmvnorm( x = theta.star, mean = mustar.g, sigma = Bstar.g, log=TRUE)
        
        ## compute U*(g)
        Ustar.g = nr*t(res1[[g]][[1]]-bary) %*% (res1[[g]][[1]]-bary)+ U + S    
        ## π(W* | y, nw*, U*(g)) -- in log, without constant part
        temp2 = (nwstar-p-1)/2*(log(det(Ustar.g)))- sum(diag(IW.star%*%Ustar.g))/2
        
        c(temp1, temp2)
    })
    n.iter=n.iter-n.burnin
    
    ## π(Ѳ* | y) -- in log
    lpihat.theta.star = log(sum(exp(res2[1,]))/n.iter)
    
    ## π(W* | y) -- in log
    # with log(prod(gamma( (nwstar -p-1:p)/2) )*pi^(p*(p-1)/4))
    #    = sum(log( gamma( (nwstar -p-1:p)/2 ) )) + log( pi^(p*(p-1)/4) )
    const = log(det(IW.star))*(nwstar)/2 - log(2)*(nwstar-p-1)*(p/2) - MultivariateGamma(nwstar/2, p)
    lpihat.W.star = log(sum(exp(res2[2,]+ const))/n.iter)
    
    ###  l.π(ψ* | y) = l.π(Ѳ* | y)       + l.π(W* | y) )
    lpihat.psi.star = lpihat.theta.star + lpihat.W.star
    
    ### ESTIMATE THE POSTERIOR ORDINATES
    ### π(ψ* | H1) = π(Ѳ* | y, H1) x π(W* | y, H1)
    ## π(Ѳ* | H1) -- in log
    lpi.theta.star=dmvnorm(x = theta.star, mu, solve(invB), log = TRUE)
    ## (W* | y, H1) -- in log
    const = log(det(IW.star))*(nw)/2 - log(2)*(nw-p-1)*(p/2) - MultivariateGamma(nw/2, p)
    lpi.W.star = const +  (nw-p-1)/2*(logdetU) - sum(diag(IW.star %*% U))/2
    ## π(ψ* | y, H1) -- in log
    ### l.π(ψ* | H1) = l.π(Ѳ* | y, H1)  x l.π(W* | y, H1)
    lpi.psi.star     = lpi.theta.star   + lpi.W.star
    
    ### COMPUTE marginal density
    ### l.m(y | H1) = l.f(y | ψ* , H1) + π(ψ* | H1)   - π(ψ* | y, H1)
    lmhatHp.num     = logf.star        + lpi.psi.star - lpihat.psi.star
    return(lmhatHp.num)
}
