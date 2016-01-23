\name{seARTFIMA}
\alias{seARTFIMA}

\title{
Standard errors of ARTFIMA based on Fisher Information
}
\description{
Estimated standard errors are obtained from the large sample expected
Fisher information matrix
}
\usage{
seARTFIMA(d, lambda, n)
}

\arguments{
  \item{d}{
ARTFIMA d parameter
}
  \item{lambda}{
ARTFIMA lambda parameter
}
  \item{n}{
length of time series
}
}
\value{
vector with the estimated standard errors
}
\author{
A. I. McLeod
}


\seealso{
\code{\link{informationMatrixARTFIMA}}
}
\examples{
seARTFIMA(d=5/8, lambda=0.06, n=1000)
}
\keyword{ ts }
