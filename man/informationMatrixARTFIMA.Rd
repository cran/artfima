\name{informationMatrixARTFIMA}
\alias{informationMatrixARTFIMA}
\title{
Information matrix for ARTFIMA
}
\description{
The information matrix for the lambda and d in ARTFIMA model.
}
\usage{
informationMatrixARTFIMA(d, lambda)
}

\arguments{
  \item{d}{
d parameter
}
  \item{lambda}{
lambda parameter
}
}

\value{
The 2-by-2 asymptotic information matrix per observation
}

\author{
A. I. McLeod
}

\seealso{
\code{\link{artfima}}, 
\code{\link{seARTFIMA}}
}

\examples{
informationMatrixARTFIMA(d=5/6, lambda=0.06)
}

\keyword{ ts }
