\name{artsim}
\alias{artsim}
\title{Simulation of stationary ARTFIMA
}
\description{
Simulation of stationary ARTFIMA or bootstrap a fitted model.
}
\usage{
artsim(n = 100, d = 0, lambda = 0, phi = numeric(0), 
    theta = numeric(0), mean = 0, sigma2 = 1, obj = NULL)
}

\arguments{
   \item{n}{
length of time series
}
  \item{d}{
artfima difference parameter, real value greater than zero. 
If d=0, ARFIMA model is used.
}
  \item{lambda}{
lambda artfima temper decay parameter
}
  \item{phi}{
AR coefficients
}
  \item{theta}{
MA coefficients
}
  \item{mean}{
mean of series
}
  \item{sigma2}{
innovation variance
}
 \item{obj}{
output from artfima(). If obj is not output from artfima() 
then the other arguments are used to determine the time series parameters, 
except for the series length n.
}
}

\value{
vector of length n, the simulated time series
}
\references{
TBA
}
\author{
A. I. McLeod, aimcleod@uwo.ca
}

\examples{
z <- artsim(5000, d=5/6, lambda=0.045)
var(z)
tacvfARTFIMA(d=5/6, lambda=0.045, maxlag=1)[1]
}

\keyword{ ts }
\keyword{ datagen }