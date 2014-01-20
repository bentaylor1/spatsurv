\name{survspat}
\alias{survspat}
\title{survspat function}
\usage{
survspat(formula, data, dist, covmodel, mcmc.control, priors, nn = NULL)
}
\arguments{
  \item{formula}{see ?flexsurvreg}

  \item{data}{a SpatialPointsDataFrame object}

  \item{dist}{choice of distribution function for baseline
  hazard options are: "exp"}

  \item{covmodel}{an object of class covmodel, see
  ?covmodel}

  \item{mcmc.control}{mcmc control parameters, see
  ?mcmcpars}

  \item{priors}{an object of class Priors, see ?mcmcPriors}

  \item{nn}{If set to null (the default,) computations are
  performed with the full covariance matrix, if nn is set
  to a positive integer, this is the number of nearest
  neighbours to use in forming the sparse approximation to
  the covariance.}
}
\value{
the mcmc output
}
\description{
A function to run a Bayesian analysis on right censored
survial data assuming a proportional hazards model with
baseline hazard derived from the exponential model.
}

