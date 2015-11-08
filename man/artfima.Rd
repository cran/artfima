\name{artfima}
\alias{artfima}

\title{MLE for ARTFIMA model}

\description{
Maximum likelihood estimation of the ARTFIMA model as well as the edge cases
ARIMA and ARFIMA. Exact MLE and Whittle approximate MLE are implemented.
}
\usage{
artfima(z, glp = c("ARTFIMA", "ARFIMA", "ARIMA"), armaOrder = c(0, 0, 0), 
  constant = TRUE, seQ = TRUE, likAlg = c("Exact", "Whittle"), blueQ = FALSE, 
  fixd = NULL)
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
  \item{fixd}{only used with ARTFIMA, default setting fixd=NULL means the MLE
  for the parameter d is obtained other if fixed=d0, where d0 is a numeric
  value in the interval (-2, 2) the d parameter in ARTFIMA is fixed at this 
  value while the remaining parameters are estimated.}
}

\details{TBA}

\value{A lengthy list is produced. A terse summary is provided by the associated
print method.
}

\references{TBA}

\author{A. I. McLeod, aimcleod@uwo.ca}

\examples{
artfima(Ogden, likAlg = "Whittle")
#
set.seed(323277)
d0 <- 5/6
z <- artSim(n=1000, lambda=0.05, d=5/6)
artfima(z, fixd=d0) #constrained, exact
artfima(z, likAlg="Whittle", fixd=d0) #constrained, Whittle
}

\keyword{ts}