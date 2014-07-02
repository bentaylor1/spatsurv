\name{indepGaussianprior}
\alias{indepGaussianprior}
\title{indepGaussianprior function}
\usage{
indepGaussianprior(beta = NULL, omega = NULL, eta = NULL, priors)
}
\arguments{
  \item{beta}{parameter beta at which prior is to be
  evaluated}

  \item{omega}{parameter omega at which prior is to be
  evaluated}

  \item{eta}{parameter eta at which prior is to be
  evaluated}

  \item{priors}{an object of class mcmcPriors, see
  ?mcmcPriors}
}
\value{
the log prior
}
\description{
A function to evaluate the prior for independent normals
}

