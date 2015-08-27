\name{tacvfARTFIMA}
\alias{tacvfARTFIMA}

\title{Autocovariance function of ARTFIMA}

\description{Theoretical autocovariance function of ARTFIMA model}

\usage{
tacvfARTFIMA(d = numeric(0), lambda = numeric(0), phi = numeric(0), 
      theta = numeric(0), maxlag, sigma2 = 1)
}

\arguments{
  \item{d}{ARTFIMA difference parameter, any real value. When d=0, reduces to 
      ARMA and lambda is ignored.}
  \item{lambda}{ARTFIMA tempered decay parameter. When lambda=0, 
      reduces to ARFIMA} 
  \item{phi}{AR coefficients}
  \item{theta}{MA coefficients, Box-Jenkins definition}
  \item{maxlag}{maxlag+1 lags computed corresponding to 0,1,...,maxlag}
  \item{sigma2}{innovation variance}
}

\details{TBA}

\value{vector of length maxlag+1 of the specified autocovariances}

\references{TBA}

\author{A. I. McLeod, aimcleod@uwo.ca}

\seealso{\code{\link{ARMAacf}}}


\examples{
tacvfARTFIMA(d=1, lambda=0.1, maxlag=25)
}

\keyword{ts}

