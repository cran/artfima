\name{artfima}
\alias{artfima}

\title{MLE for ARTFIMA model}

\description{
Maximum likelihood estimation of the ARTFIMA model as well as the edge cases
ARIMA and ARFIMA. Exact MLE and Whittle approximate MLE are implemented.
}
\usage{
artfima(z, glp = c("ARTFIMA", "ARFIMA", "ARIMA"), armaOrder = c(0, 0, 0), 
  constant = TRUE, seQ = TRUE, likAlg = c("Exact", "Whittle"), blueQ = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{time series data}
  \item{glp}{general linear process type: ARTFIMA, ARFIMA or ARMA.}
  \item{armaOrder}{c(p,D,q), where p is the AR order, D is the regular difference
  parameter and q is the MA order.}
  \item{constant}{When D is zero, this is the mean. Otherwise it is the mean
  of the differenced series.}
  \item{seQ}{TRUE, standard errors of the glp parameters are estimated using
  the observed Hessian}
  \item{likAlg}{either "Exact" or "Whittle"}
  \item{blueQ}{TRUE, exact MLE for the constant parameter is done. Otherwise
  the sample mean is used.}
}

\details{TBA}

\value{A lengthy list is produced. A terse summary is provided by the associated
print method.
}

\references{TBA}

\author{A. I. McLeod, aimcleod@uwo.ca}

\examples{
artfima(Ogden, likAlg = "Whittle")
}

\keyword{ts}