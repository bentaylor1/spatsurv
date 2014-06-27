##' distinfo.weibull function
##'
##' A function to return the internal transformation function (and its inverse) for the baseline hazard type. E.g. for an Exponential baseline hazard, we work with the log rate, so log is the transformation function. 
##' 
##' @return the transformation and inverse transformation, jacobian and hessian
##' @export

distinfo.weibull <- function(){
    retlist <- list()
    retlist$npars <- 2
    retlist$parnames <- c("alpha","lambda")
    retlist$trans <- log
    retlist$itrans <- exp
    retlist$jacobian <- exp
    retlist$hessian <- list(exp,exp)
    return(retlist)
}





#################################################################################
# Weibull survival model
#################################################################################


##' basehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

basehazard.weibull <- function(pars){
    fun <- function(t){
        return(pars[2]*pars[1]*t^(pars[1]-1)) # in this case alpha=pars[1], lambda=pars[2] 
    }
    return(fun)  
}



##' gradbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

gradbasehazard.weibull <- function(pars){
    fun <- function(t){
        return(t^(pars[1]-1)*cbind(pars[2]*(1+pars[1]*log(t)),pars[1])) # in this case alpha=pars[1], lambda=pars[2]
    }
    return(fun)
    
}



##' hessbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

hessbasehazard.weibull <- function(pars){
    funfun <- function(t,pars){
        m <- matrix(0,2,2) # note m[2,2]=0 i.e. d2h_0/dlambda^2 = 0
        m[1,2] <- m[2,1] <- t^(pars[1]-1)*(1+pars[1]*log(t))
        m[1,1] <- pars[2]*t^(pars[1]-1)*log(t)*(2+pars[1]*log(t))
        return(as.vector(m)) # in this case alpha=pars[1], lambda=pars[2]
    }
    
    fun <- function(t){
        return(lapply(t,funfun,pars=pars))
    }
    return(fun)
    
}



##' cumbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

cumbasehazard.weibull <- function(pars){
    fun <- function(t){
        return(pars[2]*t^(pars[1])) # in this case alpha=pars[1], lambda=pars[2]
    }
    return(fun) 
}




##' gradcumbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

gradcumbasehazard.weibull <- function(pars){
    fun <- function(t){
        return(t^(pars[1])*cbind(pars[2]*log(t),1)) # in this case alpha=pars[1], lambda=pars[2]
    }
    return(fun)    
}




##' hesscumbasehazard.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export
hesscumbasehazard.weibull <- function(pars){
    funfun <- function(t,pars){
        m <- matrix(0,2,2) # note m[2,2]=0 i.e. d2H_0/dlambda^2 = 0
        other <- log(t)*t^pars[1]
        m[1,2] <- m[2,1] <- other 
        m[1,1] <- pars[2]*other*log(t)
        return(m) # in this case alpha=pars[1], lambda=pars[2]
    }
    
    fun <- function(t){
        return(lapply(t,funfun,pars=pars))
    }
    return(fun)
}



#################################################################################



##' densityquantile.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @param other X
##' @return ...
##' @export

densityquantile.weibull <- function(pars,other){
    fun <- function(probs){
        return((-log(1-probs)/(pars[2]*other$expXbetaplusY))^(1/pars[1]))
    }
    return(fun)    
}