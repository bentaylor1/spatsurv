\name{gencens}
\alias{gencens}
\title{gencens function}
\usage{
gencens(survtimes, censtimes, type = "right")
}
\arguments{
  \item{survtimes}{a vector of survival times}

  \item{censtimes}{a vector of censoring times for left or
  right censored data, 2-column matrix of censoring times
  for interval censoring (number of rows equal to the
  number of observations).}

  \item{type}{the type of censoring to generate can be
  'right' (default), 'left' or 'interval'}
}
\value{
a named list containing 'obstimes', the observed time of
the event; and 'censored', the censoring indicator which is
equal to 1 if the event is observed and 0 otherwise.
}
\description{
A function to generate observed times given a vector of
true survival times and a vector of censoring times. Used
in the simulation of survival data
}

