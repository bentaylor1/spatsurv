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
    cat("Fixed Effects:\n")    
    print(quant$betaquant,digits=digits,scientific=scientific)
    cat("\n")
    
    cat("Baseline Hazard Parameters:\n")
    print(quant$omegaquant,digits=digits,scientific=scientific)
    cat("\n")

    cat("Spatial Covariance Parameters:\n")   
    print(quant$etaquant,digits=digits,scientific=scientific)
    cat("\n")
    
    cat("Deviance Information Criterion: ",x$DIC,"\n")
    cat("\n")
    
    cat("MCMC Details:\n")
    m <- matrix(c(x$mcmc.control$nits,x$mcmc.control$burn,x$mcmc.control$thin),3,1)
    rownames(m) <- c("Number of Iterations","Burnin length","Thining")
    colnames(m) <- ""
    print(m,digits=digits,scientific=scientific)
    cat("\n")
    
    cat("Running Time:\n")
    print(x$time.taken)
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

quantile.mcmcspatsurv <- function(x,probs=c(0.025,0.5,0.975),...){
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



##' vcov.mcmcspatsurv function
##'
##' A function to 
##'
##' @method vcov mcmcspatsurv
##' @param object X 
##' @param ... X 
##' @return ...
##' @export

vcov.mcmcspatsurv <- function(object,...){
    return(cov(cbind(object$betasamp,object$omegasamp,object$etasamp)))
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
    plot(lag1acf,xlab="Frailty Index",ylab="Lag 1 Autocorrelation",ylim=c(-1,1),...)
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

    omegasamp <- x$omegasamp   
    
    if(is.null(t)){
        t <- seq(0,max(x$survivaldata[,"time"]),length.out=n)
    }
    
    fun <- function(pars){
        f <- get(paste("basehazard.",x$dist,sep=""))(pars)
        return(f(t))
    } 
    
    samp <- t(apply(omegasamp,1,fun))   

    toreturn <- t(apply(samp,2,quantile,probs=probs))
    
    rownames(toreturn) <- t 
    
    if(plot){
        if(length(probs)==3){
            matplot(t,toreturn,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="time",ylab="Baseline Hazard")
            legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))
        }
        else{
            matplot(t,toreturn,type="l",xlab="time",ylab="Baseline Hazard")
        }
    }    
    
    return(list(t=t,qts=toreturn))
}






##' hazard_PP function
##'
##' A function to compute the hazard function 
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the hazard function for the individual
##' @export

hazard_PP <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    f <- function(t){
        h <- get(paste("basehazard.",inputs$dist,sep=""))(inputs$omega)
        return(expXbeta_plus_Y*h(t))
    }
    return(f)      
}




##' survival_PP function
##'
##' A function to compute the survival function 
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the survival function for the individual
##' @export
survival_PP <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    f <- function(t){
        H <- get(paste("cumbasehazard.",inputs$dist,sep=""))(inputs$omega)
        return(exp(-expXbeta_plus_Y*H(t)))
    }
    return(f)      
}



##' density_PP function
##'
##' A function to compute the density function 
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return the density function for the individual
##' @export

density_PP <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    f <- function(t){
        h <- get(paste("basehazard.",inputs$dist,sep=""))(inputs$omega)    
        H <- get(paste("cumbasehazard.",inputs$dist,sep=""))(inputs$omega)
        return(expXbeta_plus_Y*h(t)*exp(-expXbeta_plus_Y*H(t)))
    }
    return(f)      
}



##' Et_PP function
##'
##' A function to compute the expected survival time 
##'
##' @param inputs X 
##' @return ...
##' @export

Et_PP <- function(inputs){
    f <- density_PP(inputs)
    expect <- function(t){
        return(t*f(t))
    }
    int <- integrate(expect,0,Inf)$value
    return(int) 
}



##' densityquantile_PP function
##'
##' A function to compute quantiles of the density function 
##'
##' @param inputs inputs for the function including the model matrix, frailties, fixed effects and the parameters of the baseline hazard derived from this model
##' @return quantiles of the density function for the individual
##' @export

densityquantile_PP <- function(inputs){
    X <- inputs$X
    Y <- inputs$Y
    beta <- inputs$beta
    
    Xbeta <- sum(X*beta)
    expXbeta <- exp(Xbeta)
    expXbeta_plus_Y <- expXbeta*exp(Y)
    
    dq <- get(paste("densityquantile.",inputs$dist,sep=""))(inputs$omega,other=list(expXbetaplusY=expXbeta_plus_Y)) 
    
    return(dq)    
}





##' predict.mcmcspatsurv function
##'
##' A function to 
##'
##' @method predict mcmcspatsurv
##' @param object an object of class mcmcspatsurv
##' @param type can be "densityquantile","density", "hazard" or "survival". Default is "densityquantile".
##' @param newdata X
##' @param t X
##' @param n X
##' @param indx X
##' @param probs X
##' @param plot X
##' @param pause X
##' @param ... other arguments 
##' @return ...
##' @export

predict.mcmcspatsurv <- function(object,type="densityquantile",newdata=NULL,t=NULL,n=110,indx=NULL,probs=c(0.025,0.5,0.975),plot=TRUE,pause=TRUE,...){

    if(is.null(newdata)){
        newdata <- object$X
    }

    if(is.null(indx)){
        indx <- 1:nrow(newdata)
    }    
    
    nobs <- length(indx) # nrow(newdata)
    
    if(is.null(t)){
        t <- seq(0,max(object$survivaldata[,"time"]),length.out=n)
    }

    predictmat <- NULL    
    if(type=="densityquantile"){
        t <- probs
        n <- length(t)
        pb <- txtProgressBar(min = 0, max = nobs)
        predictmat <- matrix(NA,nobs,length(probs))
    }
    
    if(type=="Et"){
        n <- nobs
        pb <- txtProgressBar(min = 0, max = nobs)
    }
    
    
    
    nits <- nrow(object$Ysamp)
    
    omegasamp <- object$omegasamp    
    dat <- matrix(NA,nits,n)
    for(i in indx){
        if(type!="Et"){
            dat <- matrix(NA,nits,n)
        }
        inputs <- list()
        inputs$X <- newdata[i,]
        inputs$dist <- object$dist
        for (j in 1:nits){
            inputs$Y <- object$Ysamp[j,i]
            inputs$beta <- object$betasamp[j,]
            inputs$omega <- omegasamp[j,]
            fun <- get(paste(type,"_PP",sep=""))(inputs)
            if(type!="Et"){
                dat[j,] <- fun(t)
            }
            else{
                dat[j,i] <- fun
            }      
        }

        if(type!="densityquantile" & type!="Et"){          
            toplot <- t(apply(dat,2,quantile,probs=probs))        
            matplot(t,toplot,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="time",ylab=type,main=paste("Individual",i))
            legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))
            if(pause){
                cat("[press [enter] to continue]")
                ob <- scan(n=1,quiet=TRUE)
            }
        }
        else{        
            setTxtProgressBar(pb,i)
        }
        
        if(!is.null(predictmat)){
            if(type!="Et"){
                predictmat[i,] <- colMeans(dat)
            }
        }     
    }
    
    if(type=="densityquantile" & type=="Et"){
        close(pb)
    }    
    
    if(type=="Et"){
        predictmat <- colMeans(dat)
        attr(predictmat,"empirical") <- t(dat) 
    }
    
    if(length(indx==1)&(type=="hazard"|type=="survival"|type=="density")){
        predictmat <- toplot
    }
       
    
    return(list(t=t,predict=predictmat))
    
}


## plot.mcmcspatsurv function
##
## A function to produce diagnostic plots for objects of class mcmcspatsurv
##
## @method plot mcmcspatsurv
## @param x an object of class mcmcspatsurv
## @param n number of time points to consider
## @param pr optional predictions, if they've already been computed, must be type="densityquantile", see ?predict.mcmcspatsurv. If NULL, the predicted median survival time is used.
## @param alpha significance level, default is 0.05. 
## @param ... other arguments  
## @return produces some diagnostic plots (currently only one diagnostic plot...)
## @export

#plot.mcmcspatsurv <- function(x,n=1000,pr=NULL,alpha=0.05,...){
#
#    warning("This function is currently being de-bugged.",immediate.=TRUE)
#
#    survdat <- getsurvdata(x)
#    times <- survdat[,"time"]
#    cens <- survdat[,"status"]
#    maxt <- max(times)
#    tm <- seq(min(times),maxt,length.out=n)
#    
#    if(is.null(pr)){
#        cat("To save time you can use the pr argument to feed in precomputed predictions using the predict function with type='Et', see ?predict.mcmcspatsurv.\n")
#        pr <- predict(x,type="Et")
#    }
#    
#    p <- rep(0,n)
#    np <- rep(0,n)
#    o <- rep(0,n)
#    no <- rep(0,n)
#    for(i in 1:n){
#        temptimes <- times
#        temptimes[times<=tm[i] & cens==0] <- NA
#        p[i] <- sum(pr>tm[i])
#        np[i] <- length(times)
#        o[i] <- sum(temptimes>tm[i],na.rm=TRUE)
#        no[i] <- sum(!is.na(temptimes))
#    }
#    
#    rr <- (p/np)/(o/no)
#    logrr <- log(rr)
#    
#    selogrr <- sqrt(1/p-1/np+1/o-1/no)
#    z <- qnorm(1-alpha/2)
#    lower <- exp(logrr - z*selogrr)
#    upper <- exp(logrr + z*selogrr)
#    
#    lower[is.na(lower)|is.infinite(lower)] <- NA
#    upper[is.na(upper)|is.infinite(upper)] <- NA
#    
#    tmax <- tm[min(which(is.na(lower))[1],which(is.na(upper))[1])]   
#    
#    plot(NULL,xlim=c(0,tmax),ylim=range(c(lower,upper),na.rm=TRUE),xlab="Time",ylab="propn. predicted to survive / propn. observed to survive")
#    lines(tm,lower,col="red",lty="dashed")
#    lines(tm,upper,col="red",lty="dashed")
#    lines(tm,rr)
#    abline(h=1,col="blue")
#    return(list(pr=pr,tm=tm,rr=rr,selogrr=selogrr,lower=lower,upper=upper))
#}




##' priorposterior function
##'
##' A function to 
##'
##' @param x X 
##' @param breaks X 
##' @param ylab X 
##' @param main X 
##' @param pause X 
##' @param ... X 
##' @return ...
##' @export

priorposterior <- function(x,breaks=30,ylab="Density",main="",pause=TRUE,...){
    nbeta <- ncol(x$betasamp)
    nomega <- ncol(x$omegasamp)
    neta <- ncol(x$etasamp)
    
    for(i in 1:nbeta){
        h <- hist(x$betasamp[,i],ylab=ylab,xlab=colnames(x$betasamp)[i],breaks=breaks,freq=FALSE,main=main,...)
        xrg <- range(h$breaks)
        r <- seq(xrg[1],xrg[2],length.out=1000)
        if(length(x$priors$betaprior$mean)==1){
            lines(r,dnorm(r,mean=x$priors$betaprior$mean,sd=x$priors$betaprior$sd),col="red",lwd=2)
        }
        else{
            lines(r,dnorm(r,mean=x$priors$betaprior$mean[i],sd=x$priors$betaprior$sd[i]),col="red",lwd=2)
        }
        
        if(pause){
            cat("[press [enter] to continue]")
            scan(n=1,quiet=TRUE)
        }
    }
    
      
    samp <- x$omegasamp
    samp <- x$omegatrans(samp)
    for(i in 1:nomega){        
        h <- hist(samp[,i],ylab=ylab,xlab=paste("Transformed",colnames(x$omegasamp)[i]),breaks=breaks,freq=FALSE,main=main,...)
        xrg <- range(h$breaks)
        r <- seq(xrg[1],xrg[2],length.out=1000)
        if(length(x$priors$betaprior$mean)==1){
            lines(r,dnorm(r,mean=x$priors$omegaprior$mean,sd=x$priors$omegaprior$sd),col="red",lwd=2)
        }
        else{
            lines(r,dnorm(r,mean=x$priors$omegaprior$mean[i],sd=x$priors$omegaprior$sd[i]),col="red",lwd=2)
        }
        
        if(pause){
            cat("[press [enter] to continue]")
            scan(n=1,quiet=TRUE)
        }
    }
    
    
    
    for(i in 1:neta){
        samp <- x$cov.model$trans[[i]](x$etasamp[,i])
        h <- hist(samp,ylab=ylab,xlab=paste("Transformed",colnames(x$etasamp)[i]),breaks=breaks,freq=FALSE,main=main,...)
        xrg <- range(h$breaks)
        r <- seq(xrg[1],xrg[2],length.out=1000)
        if(length(x$priors$etaprior$mean)==1){
            lines(r,dnorm(r,mean=x$priors$etaprior$mean,sd=x$priors$etaprior$sd),col="red",lwd=2)
        }
        else{
            lines(r,dnorm(r,mean=x$priors$etaprior$mean[i],sd=x$priors$etaprior$sd[i]),col="red",lwd=2)
        }
        
        if(pause){
            cat("[press [enter] to continue]")
            scan(n=1,quiet=TRUE)
        }
    }
}



##' posteriorcov function
##'
##' A function to produce a plot of the posterior covariance function.
##'
##' @param x an object of class mcmcspatsurv 
##' @param probs vector of probabilities to be fed to quantile function  
##' @param rmax  maximum distance in space to compute this distance up to
##' @param n the number of points at which to evaluate the posterior covariance.
##' @param plot whether to plot the result
##' @param ... other arguments to be passed to matplot function 
##' @return produces a plot of the posterior spatial covariance function.
##' @export

posteriorcov <- function(x,probs=c(0.025,0.5,0.975),rmax=NULL,n=100,plot=TRUE,...){
    nr <- nrow(x$etasamp)
    nc <- ncol(x$etasamp)
    pars <- matrix(NA,nr,nc)
    for(i in 1:nc){
        pars[,i] <- x$cov.model$trans[[i]](x$etasamp[,i]) # transform (e.g. to log-scal)
    }
    rmaxx <- 0.25*sum(apply(bbox(x$data),1,diff))/2 # approx 1/4 of mean length of observation window
    if(!is.null(rmax)){
        rmaxx <- rmax
    }
    
    r <- seq(0,rmaxx,length.out=n)
    covs <- t(apply(pars,1,function(pp){x$cov.model$eval(r,pars=pp)}))
    qts <- t(apply(covs,2,quantile,probs=probs))
    rownames(qts) <- r
    
    if(plot){
        if(length(probs)==3){
            matplot(r,qts,type="l",col=c("purple","black","blue"),lty=c("dashed","solid","dashed"),xlab="Distance",ylab="Covariance")
            legend("topright",lty=c("dashed","solid","dashed"),col=rev(c("purple","black","blue")),legend=rev(probs))
        }
        else{
            matplot(r,qts,type="l",xlab="Distance",ylab="Covariance")
        }
    }
    
    return(list(r=r,qts=qts)) 
}



## spatialpredict function
##
## A function to 
##
## @param object X 
## @param xgrid X 
## @param ygrid X 
## @param cellwidth X 
## @return ...
## @export

#spatialpredict <- function(object,xgrid=NULL,ygrid=NULL,cellwidth=NULL){
#
#    if(is.null(xgrid)&is.null(ygrid)){
#        if(is.null(cellwidth)){
#            stop("You must specify cellwidth or xgrid and ygrid")
#        }
#        grid <- FFTgrid(spatialdata=object$data,cellwidth=cellwidth,ext=1)
#        xgrid <- grid$mcens
#        ygrid <- grid$ncens  
#    }
#
#    gr <- as.matrix(expand.grid(xgrid,ygrid))
#    crds <- rbind(gr,coordinates(object$data))
#    
#    n <- nrow(object$X)
#    N <- nrow(gr)    
#    mx <- nrow(crds)
#    
#    u <- as.vector(as.matrix(dist(crds)))
#    
#    nits <- nrow(object$Ysamp)
#    
#    etasamp <- sapply(1:length(object$cov.model$trans),function(i){object$cov.model$trans[[i]](object$etasamp[,i])})
#    
#    predY <- NULL
#    
#    Y <- colMeans(object$Ysamp)
#    eta <- colMeans(object$etasamp)
#    
#    mu_22 <- -object$cov.model$itrans[[object$control$sigmaidx]](eta[object$control$sigmaidx])^2/2                    
#    
#    etapars <- sapply(1:cov.model$npar,function(i){cov.model$itrans[[i]](eta[i])})
#    sigma <- matrix(EvalCov(cov.model=object$cov.model,u=u,parameters=etapars),mx,mx)
#   
#    sigma_12 <- sigma[1:N,(N+1):mx]
#    sigma_22inv <- solve(sigma[(N+1):mx,(N+1):mx])
# 
#    predY <- matrix(sigma_12%*%sigma_22inv%*%(Y-mu_22),length(xgrid),length(ygrid))
#    
#    return(list(xgrid=xgrid,ygrid=ygrid,predY=predY))
#}


##' makegreycale function
##'
##' A function to 
##'
##' @param v X 
##' @return ...
##' @export

makegreyscale <- function(v){
    v <- v-min(v,na.rm=TRUE) # smallest v now zero
    v <- v / max(v,na.rm=TRUE)
    return(grey(1-v))
}



##' MCE function
##'
##' A function to 
##'
##' @param object X 
##' @param fun X 
##' @return ...
##' @export

MCE <- function(object,fun){
    nits <- nrow(object$betasamp)
    
    result <- lapply(1:nits,function(i){fun(beta=object$betasamp[i,],omega=object$omegasamp[i,],eta=object$etasamp[i,],Y=object$Ysamp[i,])})
    
    return((1/nits)*Reduce('+',result))
}




##' hazardexceedance function
##'
##' A function to 
##'
##' @param threshold X 
##' @param direction X 
##' @return ...
##' @export

hazardexceedance <- function(threshold,direction="upper"){
    fun <- function(beta,omega,eta,Y){
        EY <- exp(Y)
        d <- length(Y)
        len <- length(threshold)
        A <- matrix(NA,len,d)
        
        for(i in 1:len){
            if(direction=="upper"){
                A[i,] <- as.numeric(EY>threshold[i])
            }
            else{
                A[i,] <- as.numeric(EY<threshold[i])
            }
        }
        return(A)        
    }
    attr(fun,"threshold") <- threshold
    attr(fun,"direction") <- direction
    return(fun)
}
