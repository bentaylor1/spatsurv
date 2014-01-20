\name{proposalvariance.weibull}
\alias{proposalvariance.weibull}
\title{proposalvariance.weibull function}
\usage{
proposalvariance.weibull(X, delta, tm, betahat, omegahat, Yhat, priors,
  covmodel, u, sparse, ninfo)
}
\arguments{
  \item{X}{design matrix}

  \item{delta}{censoring indicator}

  \item{tm}{observed times}

  \item{betahat}{estimate of model parameter beta}

  \item{omegahat}{estimate of model parameter omega}

  \item{Yhat}{estimate of the latent field}

  \item{priors}{the priors, an object of class
  'mcmcPriors', see ?mcmcPriors}

  \item{covmodel}{an object of class 'covmodel', see
  ?covmodel}

  \item{u}{vector of distances between points}

  \item{sparse}{logical: whether to use a full, or a sparse
  representation of the covariance structure.}

  \item{ninfo}{neighbourhood information for constructing
  sparse matrix}
}
\value{
estimates of eta, gamma and a proposal variance matrixc for
use in the MALA algorithm
}
\description{
A function to compute an approximate proposal variance for
the MALA algorithm for a parametric proportional hazards
model with baseline hazard derived from the Weibull
distribution.
}

