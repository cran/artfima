\name{artfima}
\alias{artfima}

\title{MLE for ARTFIMA model}

\description{
Maximum likelihood estimation of the ARTFIMA model as well as the edge cases
ARIMA and ARFIMA. Exact MLE and Whittle approximate MLE are implemented.
}
\usage{
artfima(z, glp = c("ARTFIMA", "ARFIMA", "ARIMA"), arimaOrder = c(0, 0, 0), 
  likAlg = c("exact", "Whittle"),  fixd = NULL)
}

\arguments{
  \item{z}{time series data}
  \item{glp}{general linear process type: ARTFIMA, ARFIMA or ARMA.}
  \item{arimaOrder}{c(p,D,q), where p is the AR order, D is the regular difference
  parameter and q is the MA order.}
  \item{likAlg}{either "exact" or "Whittle"}
  \item{fixd}{only used with ARTFIMA, default setting fixd=NULL means the MLE
  for the parameter d is obtained other if fixed=d0, where d0 is a numeric
  value in the interval (-2, 2) the d parameter in ARTFIMA is fixed at this 
  value while the remaining parameters are estimated.}
}

\details{The ARFIMA and ARIMA are subsets or edge-cases
of the ARTFIMA model. The likelihood and probability density function
for these models is defined by the multivariate normal distribution.
The log-likelihood, AIC and BIC are comparable across models.
When the Whittle MLE algorithm is used, the final log-likelihood
is obtained by plugging this estimates into the exact log-likelihood.
}

\value{A lengthy list is produced. A terse summary is provided by the associated
print method.
}

\references{TBA}

\author{A. I. McLeod, aimcleod@uwo.ca}

\seealso{
\code{\link{bestModels}}
}

\examples{
artfima(Nile)
artfima(Nile, likAlg = "exact")
}

\keyword{ts}