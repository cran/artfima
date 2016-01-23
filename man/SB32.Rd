\name{SB32}
\alias{SB32}
\docType{data}
\title{
Turbulent flow data from Station SB32
}
\description{
Turbulent flow water time series, Lake Huron, druing 2009-2010.
}
\usage{data("SB32")}
\format{
  The format is:
 num [1:5374] -2.2 1.4 -0.6 -0.4 -1.5 -2.6 -0.9 0.5 -0.9 1.5 ...
}
\details{
See paper by Meerschaert, Sabzikar, Phanikumar and Zeleke (2014).
}
\references{
M.M. Meerschaert, Farzad Sabzikar, M.S. Phanikumar, and A. Zeleke, 
Tempered fractional time series model for turbulence in geophysical flows, 
Journal of Statistical Mechanics: Theory and Experiment, 
Vol. 2014 p. P09023 (13 pp.) doi:10.1088/1742-5468/2014/09/P09023.
}

\examples{
%Figure 2 in Meerschaert, Sabzikar, Phanikumar and Zeleke (2014)
data(SB32)
ans <- spec.pgram(SB32, fast=FALSE, detrend=FALSE, plot=FALSE, taper=0)
lsdf <- log((ans$spec)/(2*pi))
lfre <- log(ans$freq) 
plot(lfre, lsdf, xlab="log frequency", ylab="log spectral density", 
     col=gray(0, 0.7))
sdf <- sdfartfima(n=length(SB32), d=5/6, lambda=0.045)
lines(lfre, log(sdf), col="red", lwd=3)
}
\keyword{datasets}
