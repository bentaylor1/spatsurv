% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatgrid.R
\name{getOptCellwidth}
\alias{getOptCellwidth}
\title{getOptCellwidth function}
\usage{
getOptCellwidth(dat, cellwidth, ext = 2, plot = TRUE, boundingbox = NULL)
}
\arguments{
\item{dat}{any spatial data object whose bounding box can be computed using the function bbox.}

\item{cellwidth}{an initial suggested cellwidth}

\item{ext}{the extension parameter for the FFT transform, set to 2 by default}

\item{plot}{whether to plot the grid and data to illustrate the optimal grid}

\item{boundingbox}{optional bounding box over which to construct computational grid, supplied as an object on which the function 'bbox' returns the bounding box}
}
\value{
the optimum cell width
}
\description{
A function to compute an optimal cellwidth close to an initial suggestion. This maximises the efficiency of the
MCMC algorithm when in the control argument of the function survspat, the option gridded is set to TRUE
}

