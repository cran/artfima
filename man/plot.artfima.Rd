\name{plot.artfima}
\alias{plot.artfima}

\title{Plot Method for "arfima" Object}

\description{
Plots the observed periodogram and the fitted spectral density function.
}

\usage{\method{plot}{artfima}(x, type, ...)}

\arguments{
  \item{x}{ object of class "artfima" }
  \item{type}{ "log" or "log-log" plot }
  \item{...}{ optional arguments }
}

\value{None}

\references{TBA}

\author{A. I. McLeod, aimcleod@uwo.ca}

\seealso{\code{\link{artfima}}}

\examples{
z <- artsim(n=500, d=5/6, lambda=0.045)
ans <- artfima(z)
plot(ans)
plot(ans, type="log-log")
}

\keyword{ts}

