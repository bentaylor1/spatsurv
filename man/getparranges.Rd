% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatsurvMisc.R
\name{getparranges}
\alias{getparranges}
\title{getparranges function}
\usage{
getparranges(priors, leneta, mult = 1.96)
}
\arguments{
\item{priors}{an object of class mcmcPriors}

\item{leneta}{the length of eta passed to the function}

\item{mult}{defaults to 1.96 so the grid formed will be mean plus/minus 1.96 times the standard deviation}
}
\value{
an appropriate range used to calibrate the MCMC: the mean of the prior for eta plus/minus 1.96 times the standard deviation
}
\description{
A function to extract parameter ranges for creating a grid on which to evaluate the log-posterior, used in calibrating the MCMC. This function
is not intended for general use.
}

