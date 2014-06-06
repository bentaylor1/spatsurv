##' print.mcmcspatsurv function
##'
##' A function to print summary tables from an MCMC run 
##'
##' @method print mcmcspatsurv
##' @param x an object inheriting class mcmcspatsurv 
##' @param probs vector of quantiles to return
##' @param digits see help file for format
##' @param scientific see help file for format
##' @param ... additional arguments 
##' @return summary tables
##' @export

print.mcmcspatsurv <- function(x,probs=c(0.5,0.025,0.975),digits = 3, scientific = -3,...){

    quant <- quantile(x,probs)

    cat("\n")
    cat("Fixed Effects:\n\n")    
    print(quant$betaquant,digits=digits,scientific=scientific)
    cat("\n")
    
    cat("Baseline Hazard Parameters:\n\n")
    print(quant$omegaquant,digits=digits,scientific=scientific)
    cat("\n")

    cat("Spatial Covariance Parameters:\n\n")   
    print(quant$etaquant,digits=digits,scientific=scientific)
    cat("\n")
} 

##' quantile.mcmcspatsurv function
##'
##' A function to extract quantiles of the parameters from an mcmc run
##'
##' @method quantile mcmcspatsurv
##' @param x an object of class mcmc spatsurv 
##' @param probs vector of probabilities
##' @param ... other arguments to be passed to the function
##' @return ...
##' @export

quantile.mcmcspatsurv <- function(x,probs,...){
    m1 <- t(apply(x$betasamp,2,quantile,probs=probs))
    m2 <- t(apply(x$omegasamp,2,quantile,probs=probs))
    m3 <- t(apply(x$etasamp,2,quantile,probs=probs))
    return(list(betaquant=m1,omegaquant=m2,etaquant=m3))
}



##' summary.mcmcspatsurv function
##'
##' A function to print summary tables from an MCMC run 
##'
##' @method summary mcmcspatsurv
##' @param object an object inheriting class mcmcspatsurv 
##' @param probs vector of quantiles to return
##' @param ... additional arguments 
##' @return summary tables
##' @export

summary.mcmcspatsurv <- function(object,probs=c(0.5,0.025,0.975),...){
    quant <- quantile(object,probs)
    return(rbind(quant$betaquant,quant$omegaquant,quant$etaquant))
} 


##' frailtylag1 function
##'
##' A function to produce and return the lag 1 autocorrelation for each of the spatially correlated frailty chains
##'
##' @param object an object of class mcmcspatsurv 
##' @param ... other arguments to be passed to the plot function 
##' @return the lag 1 autocorrelation for each of the spatially correlated frailty chains 
##' @export

frailtylag1 <- function(object,...){
    lag1acf <- apply(object$Ysamp,2,function(x){acf(x,plot=FALSE)$acf[2]})
    plot(lag1acf,xlab="Frailty Index",ylab="Lag 1 Autocorrelation",xlim=c(-1,1),...)
    return(lag1acf)
}




##' spatialpars function
##'
##' A function to return the mcmc chains for the spatial covariance function parameters
##'
##' @param x an object of class mcmcspatsurv
##' @return the mcmc chains
##' @export

spatialpars <- function(x){
    return(x$etasamp)
}



##' hazardpars function
##'
##' A function to return the mcmc chains for the hazard function parameters
##'
##' @param x an object of class mcmcspatsurv
##' @return the mcmc chains
##' @export

hazardpars <- function(x){
    return(x$omegasamp)
}



##' fixedpars function
##'
##' A function to return the mcmc chains for the covariate effects
##'
##' @param x an object of class mcmcspatsurv
##' @return the mcmc chains
##' @export

fixedpars <- function(x){
    return(x$etasamp)
}



##' randompars function
##'
##' A function to return the mcmc chains for the spatially correlated frailties
##'
##' @param x an object of class mcmcspatsurv
##' @return the mcmc chains
##' @export

randompars <- function(x){
    return(x$etasamp)
}



##' baselinehazard function
##'
##' A function to 
##'
##' @param x X 
##' @param t X 
##' @param n X
##' @param probs X
##' @param plot X 
##' @return ...
##' @export

baselinehazard <- function(x,t=NULL,n=100,probs=c(0.025,0.5,0.975),plot=TRUE){
    # extract and transform onto appropriate scale    
    transfun <- get(paste("transformestimates.",x$dist,sep=""))
    if(ncol(x$omegasamp)==1){    
        omegasamp <- matrix(apply(x$omegasamp,1,transfun))
    }
    else{
        omegasamp <- t(apply(x$omegasamp,1,transfun))
    }    
    
    if(is.null(t)){
        t <- seq(0,max(x$mlmod$data$Y[,"time"]),length.out=n)
    }
    
    fun <- function(pars){
        f <- get(paste("basehaz.",x$dist,sep=""))(pars)
        return(f(t))
    } 
    
    samp <- t(apply(omegasamp,1,fun))   

    toreturn <- t(apply(samp,2,quantile,probs=probs))
    
    rownames(toreturn) <- t 
    
    if(plot){
        if(length(probs)==3){
            matplot(toreturn,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="time",ylab="Baseline Hazard")
            legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))
        }
        else{
            matplot(toreturn,type="l",xlab="time",ylab="Baseline Hazard")
        }
    }    
    
    return(toreturn)
}


##' basehaz.exp function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

basehaz.exp <- function(pars){
    f <- function(t){
        return(pars)
    }
    return(f)
}



##' basehaz.weibull function
##'
##' A function to 
##'
##' @param pars X 
##' @return ...
##' @export

basehaz.weibull <- function(pars){
    alpha <- pars[1]
    lambda <- pars[2]
    f <- function(t){
        return(lambda*alpha*t^(alpha-1))
    }
    return(f)
}



##' hazard_exp function
##'
##' A function to compute the hazard function for an individual where the baseline hazard comes from an exponential survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the hazard function for the individual
##' @export

hazard_exp <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    rate <- inputs$omega
    
    f <- function(t){
        return(expXbeta_plus_Y*rate)
    }
    return(f)      
}




##' survival_exp function
##'
##' A function to compute the survival function for an individual where the baseline hazard comes from an exponential survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the survival function for the individual
##' @export
survival_exp <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    rate <- inputs$omega
    
    f <- function(t){
        return(exp(-expXbeta_plus_Y*rate*t))
    }
    return(f)      
}



##' density_exp function
##'
##' A function to compute the density function for an individual where the baseline hazard comes from an exponential survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the density function for the individual
##' @export

density_exp <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    rate <- inputs$omega
    
    f <- function(t){
        return(expXbeta_plus_Y*rate*exp(-expXbeta_plus_Y*rate*t))
    }
    return(f)      
}




##' hazard_weibull function
##'
##' A function to compute the hazard function for an individual where the baseline hazard comes from a Weibull survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the hazard function for the individual
##' @export

hazard_weibull <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    alpha <- inputs$omega[1]
    lambda <- inputs$omega[2] 
    
    f <- function(t){
        return(expXbeta_plus_Y*lambda*alpha*t^(alpha-1))
    }
    return(f)    
}



##' survival_weibull function
##'
##' A function to compute the survival function for an individual where the baseline hazard comes from a Weibull survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the survival function for the individual
##' @export

survival_weibull <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    alpha <- inputs$omega[1]
    lambda <- inputs$omega[2]   
    
    f <- function(t){
        return(exp(-expXbeta_plus_Y*lambda*t^alpha))
    }
    return(f)    
}



##' density_weibull function
##'
##' A function to compute the density function for an individual where the baseline hazard comes from a Weibull survival model
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the density function for the individual
##' @export

density_weibull <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    alpha <- inputs$omega[1]
    lambda <- inputs$omega[2]
    
    f <- function(t){
        return(expXbeta_plus_Y*lambda*alpha*t^(alpha-1)*exp(-expXbeta_plus_Y*lambda*t^alpha))
    }
    return(f)    
}


##' predict.mcmcspatsurv function
##'
##' A function to 
##'
##' @method predict mcmcspatsurv
##' @param object an object of class mcmcspatsurv
##' @param newdata X
##' @param type can be "density", "hazard" or "survival". Default is "density".
##' @param t X
##' @param n X
##' @param probs X
##' @param ... other arguments 
##' @return ...
##' @export

predict.mcmcspatsurv <- function(object,newdata=NULL,type="density",t=NULL,n=110,probs=c(0.025,0.5,0.975),...){

    if(is.null(newdata)){
        newdata <- object$mlmod$data$X
    }
    nobs <- nrow(newdata)
    
    if(is.null(t)){
        t <- seq(0,max(object$mlmod$data$Y[,"time"]),length.out=n)
    }
    
    nits <- nrow(object$Ysamp)

    predictedtime <- c()
    predictmat <- matrix(NA,nobs,length(probs))    
    
    pb <- txtProgressBar(min = 1, max = nobs)    
    
    for(i in 1:nobs){
        dat <- matrix(NA,nits,n)
        inputs <- list()
        inputs$X <- newdata[i,]
        for (j in 1:nits){
            inputs$Y <- object$Ysamp[j,i]
            inputs$beta <- object$betasamp[j,]
            inputs$omega <- object$omegasamp[j,]
            fun <- get(paste(type,"_",object$dist,sep=""))(inputs)
            dat[j,] <- fun(t)
            
            #if(type=="density"){ 
            #    expectedvalue <- function(t){return(t*fun(t))}
            #    predictedtime[j] <- integrate(expectedvalue,lower=0,upper=Inf)$value
            #}       
        }
        
        #if(type=="density"){
        #    predictmat[i,] <- quantile(predictedtime,probs=probs)
        #}

        #toplot <- t(apply(dat,2,quantile,probs=probs))        
        #matplot(t,toplot,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="time",ylab=type,main=paste("Individual",i))
        #legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))        
        
        setTxtProgressBar(pb,i)       
    }
    close(pb)
   
    browser()
    
}