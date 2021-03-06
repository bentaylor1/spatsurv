% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/outputSummary.R
\name{Summarise}
\alias{Summarise}
\title{Summarise function}
\usage{
Summarise(obj, digits = 3, scientific = -3, inclIntercept = FALSE,
  printmode = "LaTeX", displaymode = "console", ...)
}
\arguments{
\item{obj}{an object produced by a call to lgcpPredictSpatialPlusPars, lgcpPredictAggregateSpatialPlusPars, lgcpPredictSpatioTemporalPlusPars orlgcpPredictMultitypeSpatialPlusPars}

\item{digits}{see the option "digits" in ?format}

\item{scientific}{see the option "scientific" in ?format}

\item{inclIntercept}{logical: whether to summarise the intercept term, default is FALSE.}

\item{printmode}{the format of the text to return, can be 'LaTeX' (the default) or 'text' for plain text.}

\item{displaymode}{default is 'console' alternative is 'rstudio'}

\item{...}{other arguments passed to the function "format"}
}
\value{
A text summary, that can be pasted into a LaTeX document and later edited.
}
\description{
A function to completely summarise the output of an object of class mcmcspatsurv.
}

