\name{NileMin}
\alias{NileMin}
\docType{data}
\title{ Nile Annual Minima, 622 AD to 1284 AD}
\description{Annual Minimum flow of Nile River. See below for details.}
\usage{data(NileMin)}
\format{
  The format is:
 Time-Series [1:663] from 622 to 1284: 11.57 10.88 11.69 11.69  9.84 ...
 - attr(*, "title")= "Nile River minima series"
}
\details{
  The minimum annual level of the Nile has been recorded over many
  centuries and was given by Toussoun (1925). The data over the period
  622 AD to 1284 AD is considered more homogenous and reliable than
  the full dataset and has been analyzed by Beran (1994) and Percival and
  Walden (2000).  
  The full dataset is available  StatLib Datasets hipel-mcleod
  archive -- file: Minimum. 
}
\source{
Toussoun, O. (1925).
Memoire sur l'Histoire du Nil.
In Memoires a l'Institut d'Egypte, 18, 366-404.
}

\references{
Beran, J. (1994).
Statistics for Long-Memory Processes.
Chapman and Hall, New York.

Percival, D.B. and Walden, A.T. (2000)
Wavelet Methods for Time Series Analysis.
Cambridge University Press.
}
\examples{
#Example 1
#Compute Hurst's K estimate of H
data(NileMin)
plot(NileMin)


}
\keyword{datasets}
