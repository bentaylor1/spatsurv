##' mcmcPriors function
##'
##' A function to define priors for mcmc
##'
##' @param betaprior prior for beta, the covariate effects
##' @param omegaprior prior for omega, the parameters of the baseline hazard
##' @param etaprior prior for eta, the parameters of the latent field
##' @param call function to evaluate the log-prior e.g. logindepGaussianprior
##' @param derivative function to evaluate the first and second derivatives of the prior 
##' @return an onject of class mcmcPriors
##' @seealso \link{survspat},
##' @export

mcmcPriors <- function(betaprior=NULL,omegaprior=NULL,etaprior=NULL,call=NULL,derivative=NULL){
    retlist <- list()
    retlist$betaprior <- betaprior
    retlist$omegaprior <- omegaprior
    retlist$etaprior <- etaprior
    retlist$call <- call
    retlist$derivative <- derivative
    class(retlist) <- "mcmcPriors"
    return(retlist)
}



##' betapriorGauss function
##'
##' A function to define Gaussian priors for beta.
##'
##' @param mean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param sd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @return an object of class "betapriorGauss"
##' @export

betapriorGauss <- function(mean,sd){
    retlist <- list()
    retlist$mean <- mean
    retlist$sd <- sd
    class(retlist) <- "betapriorGauss"
    return(retlist) 
}



##' omegapriorGauss function
##'
##' A function to define Gaussian priors for omega.
##'
##' @param mean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param sd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @return an object of class "omegapriorGauss"
##' @export

omegapriorGauss <- function(mean,sd){
    retlist <- list()
    retlist$mean <- mean
    retlist$sd <- sd
    class(retlist) <- "omegapriorGauss"
    return(retlist) 
}



##' etapriorGauss function
##'
##' A function to define Gaussian priors for eta.
##'
##' @param mean the prior mean, a vector of length 1 or more. 1 implies a common mean.
##' @param sd the prior standard deviation, a vector of length 1 or more. 1 implies a common standard deviation.
##' @return an object of class "etapriorGauss"
##' @export

etapriorGauss <- function(mean,sd){
    retlist <- list()
    retlist$mean <- mean
    retlist$sd <- sd
    class(retlist) <- "etapriorGauss"
    return(retlist) 
}




##' indepGaussianprior function
##'
##' A function to evaluate the prior for independent normals
##'
##' @param beta parameter beta at which prior is to be evaluated 
##' @param omega parameter omega at which prior is to be evaluated
##' @param eta parameter eta at which prior is to be evaluated
##' @param priors an object of class mcmcPriors, see ?mcmcPriors
##' @return the log prior
##' @export

indepGaussianprior <- function(beta=NULL,omega=NULL,eta=NULL,priors){
    
    lp <- 0
    if(!is.null(priors$betaprior)){
        lp <- lp + sum(dnorm(beta,priors$betaprior$mean,priors$betaprior$sd,log=TRUE))
    }

    if(!is.null(priors$omegaprior)){
        lp <- lp + sum(dnorm(omega,priors$omegaprior$mean,priors$omegaprior$sd,log=TRUE))
    }
    
    if(!is.null(priors$etaprior)){
        lp <- lp + sum(dnorm(eta,priors$etaprior$mean,priors$etaprior$sd,log=TRUE))
    }
    
    return(lp)
}


##' derivindepGaussianprior function
##'
##' A function to compute the first and second derivatives of the log-density assuming independent Gaussian priors for each of the parameters.
##'
##' @param beta a vector, the parameter beta
##' @param omega a vector, the parameter omega
##' @param eta a vector, the parameter eta 
##' @param priors an object of class 'mcmcPrior', see ?mcmcPrior
##' @return ...
##' @export

derivindepGaussianprior <- function(beta=NULL,omega=NULL,eta=NULL,priors){
    deriv1 <- c((-1/priors$betaprior$sd^2)*(beta-priors$betaprior$mean),(-1/priors$omegaprior$sd^2)*(omega-priors$omegaprior$mean),(-1/priors$etaprior$sd^2)*(eta-priors$etaprior$mean))
    sdbeta <- priors$betaprior$sd
    sdomega <- priors$omegaprior$sd
    sdeta <- priors$etaprior$sd
    if (length(priors$betaprior$sd)<length(beta)){
        sdbeta <- rep(priors$betaprior$sd,length(beta))
    }
    if (length(priors$omegaprior$sd)<length(omega)){
        sdomega <- rep(priors$omegaprior$sd,length(omega))
    }
    if (length(priors$etaprior$sd)<length(eta)){
        sdeta <- rep(priors$etaprior$sd,length(eta))
    }
    deriv2 <- c(-1/sdbeta^2,-1/sdomega^2,-1/sdeta^2)
    return(list(deriv1=deriv1,deriv2=deriv2))
}