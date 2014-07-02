library(spatsurv)

set.seed(10)

dat <- simsurv(mcmc.control=mcmcpars(nits=100,burn=10,thin=11),dist="exp")

X <- as.data.frame(dat$X) # covariates

survtimes <- dat$survtimes
n <- length(survtimes)
censtimes <- runif(n,min(survtimes),max(survtimes))                                    
survdat <- gencens(survtimes,censtimes)                                    

plot(survfit(survdat~1))
                
dat1 <- simsurv(mcmc.control=mcmcpars(nits=100,burn=10,thin=11),dist="weibull",omega=c(1,0.5))                
