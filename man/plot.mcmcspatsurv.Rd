\name{plot.mcmcspatsurv}
\alias{plot.mcmcspatsurv}
\title{plot.mcmcspatsurv function}
\usage{
\method{plot}{mcmcspatsurv}(x, n = 1000, pr = NULL, alpha = 0.05, ...)
}
\arguments{
  \item{x}{an object of class mcmcspatsurv}

  \item{n}{number of time points to consider}

  \item{pr}{optional predictions, if they've already been
  computed, must be type="densityquantile", see
  ?predict.mcmcspatsurv. If NULL, the predicted median
  survival time is used.}

  \item{alpha}{significance level, default is 0.05.}

  \item{...}{other arguments}
}
\value{
produces some diagnostic plots (currently only one
diagnostic plot...)
}
\description{
A function to produce diagnostic plots for objects of class
mcmcspatsurv
}

