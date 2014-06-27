##' distinfo.exp function
##'
##' A function to return the internal transformation function (and its inverse) for the baseline hazard type. E.g. for an Exponential baseline hazard, we work with the log rate, so log is the transformation function. 
##' 
##' @return the transformation and inverse transformation, jacobian and hessian
##' @export

distinfo.exp <- function(){
    retlist <- list()
    retlist$npars <- 1
    retlist$parnames <- "rate"
    retlist$trans <- log
    retlist$itrans <- exp
    retlist$jacobian <- exp
    retlist$hessian <- list(exp)
    return(retlist)
}




#################################################################################
# exponential survival model
#################################################################################



##' basehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

basehazard.exp <- function(pars){
    fun <- function(t){
        return(rep(pars,length(t))) # in this case pars is a 1-vector, the rate    
    }
    return(fun)
}



##' gradbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

gradbasehazard.exp <- function(pars){
    fun <- function(t){
        return(rep(1,length(t))) # in this case pars is a 1-vector, the rate
    }
    return(fun)
}


##' hessbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

hessbasehazard.exp <- function(pars){
    fun <- function(t){
        return(as.list(rep(0,length(t)))) # in this case pars is a 1-vector, the rate
    }
    return(fun)
}



##' cumbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

cumbasehazard.exp <- function(pars){
    fun <- function(t){
        return(pars*t) # in this case pars is a 1-vector, the rate
    }
    return(fun)  
}



##' gradcumbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

gradcumbasehazard.exp <- function(pars){
    fun <- function(t){
        return(t) # in this case pars is a 1-vector, the rate
    }
    return(fun)    
}



##' hesscumbasehazard.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

hesscumbasehazard.exp <- function(pars){
    fun <- function(t){
        return(as.list(rep(0,length(t)))) # in this case pars is a 1-vector, the rate
    }
    return(fun)
}





#################################################################################



##' densityquantile.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @param other X
##' @return ...
##' @export

densityquantile.exp <- function(pars,other){
    fun <- function(probs){
        return(-log(1-probs)/(pars*other$expXbetaplusY))
    }
    return(fun)    
}