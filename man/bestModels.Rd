\name{bestModels}
\alias{bestModels}

\title{Best BIC Models}

\description{ARIMA(p,0,q), ARFIMA(p,0,q) and ARTFIMA(p,0,q) models are fit for
various p=0,1,..., and q=0,1,... and the best models according to the BIC
criterion are selected.}

\usage{
bestModels(z, p = 2, q = 2, nbest = 5, likAlg = c("exact", "Whittle"), 
    use_parallel = TRUE)
}

\arguments{
  \item{z}{
time series data
}
  \item{p}{
maximum order of AR component
}
  \item{q}{
maximum order of AR component
}

  \item{nbest}{
number of models in selection
}
  \item{likAlg}{
likelihood method to use
}
  \item{use_parallel}{
By default the R parallel package is used to speed up the computations. 
The speedup is modest - about 20 to 30 percent less time - see example below. 
Requires that you have at least
3 cores. Set use_parallel=FALSE if you don't want to use parallel.
}
}
\details{
The models are ranked according to the AIC/BIC criterion and the best ones are 
shown.
}
\value{
The plausiblity score or the actual AIC/BIC values.
}
\author{
A.I. McLeod
}

\seealso{
\code{\link{best_glp_models}}
}

\examples{
\dontrun{#comparing times w/o parallelization
z<-tseg(1000, "BJARMA11")
#comparing times using parallel and not using
startTime <- proc.time()[3]
bestModels(z)
proc.time()[3]-startTime #about 132 sec
startTime <- proc.time()[3]
bestModels(z, use_parallel=FALSE)
proc.time()[3]-startTime #about 167 sec
#Remark: The correct model is ARIMA(1,0,1) but BIC tends to overfit. Usually
#the ARIMA(1,0,1) is in the top 5 and it is usually the most parsimonious 
#model in this group. If we use a stronger penalty such as BICq, q=1/4,
#the ARIMA(1,0,1) would be selected as the top model more often. 
}
}

\keyword{ ts }
