\name{sdfartfima}
\alias{sdfartfima}

\title{Computation of theoretical spectral density function (SDF)}

\description{
Computes the theoetical SDF at the Fourier frequencies for a time series of
length n. Used for Whittle MLE. Assumes model parameters are valid
for a stationary process (no checking done).
}
\usage{
sdfartfima(n, d, lambda, phi = numeric(0), theta = numeric(0))
}

\arguments{
  \item{n}{length of time series}
  \item{d}{ARTFIMA difference parameter, any real value. When d=numeric(0), 
  reduces to ARMA and lambda is ignored.}
  \item{lambda}{ARTFIMA tempered decay parameter. When lambda=numeric(0), 
      reduces to ARFIMA} 
  \item{phi}{AR coefficients}
  \item{theta}{MA coefficients, Box-Jenkins definition}
}

\value{vector of length floor(n/2) containing the values of the SDF at the 
  Fourier frequencies, 1/n, ..., floor(n/2)/n.
}

\references{TBA}

\author{A. I. McLeod, aimcleod@uwo.ca}

\seealso{
\code{\link{tacvfARTFIMA}}
}


\examples{
sdfartfima(100, d=1, lambda=0.1, phi=0.5)
}

\keyword{ ts }

