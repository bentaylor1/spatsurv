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


##' logindepnormalprior function
##'
##' A function to evaluate the log prior for independent normals
##'
##' @param beta parameter beta at which prior is to be evaluated 
##' @param omega parameter omega at which prior is to be evaluated
##' @param betapriormean prior mean for beta 
##' @param betapriorsd prior standard deviation for beta
##' @param omegapriormean prior mean fpr omega 
##' @param omegapriorsd prior standard deviation for omega 
##' @return the log prior
##' @export

logindepnormalprior <- function(beta,omega,betapriormean,betapriorsd,omegapriormean,omegapriorsd){
    return(sum(dnorm(beta,betapriormean,betapriorsd,log=TRUE))+sum(dnorm(omega,omegapriormean,omegapriorsd,log=TRUE)))
}



##' logindepGaussianprior function
##'
##' A function to evaluate the log prior for independent normals
##'
##' @param beta parameter beta at which prior is to be evaluated 
##' @param omega parameter omega at which prior is to be evaluated
##' @param eta parameter eta at which prior is to be evaluated
##' @param priors an object of class mcmcPriors, see ?mcmcPriors
##' @return the log prior
##' @export

logindepGaussianprior <- function(beta=NULL,omega=NULL,eta=NULL,priors){
    
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