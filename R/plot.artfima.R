plot.artfima <- function(x, type=c("log", "log-log"), ...) {
  stopifnot("artfima" %in% class(x))
  type <- match.arg(type)
  p <- x$arimaOrder[1]
  q <- x$arimaOrder[3]
  sigmaSq <- x$sigmaSq
  mainTitle <- paste0(x$glp, "(", p, ",", 0, ",", q, ")" )
  if ("artfmap" %in% class(x)) {
    Ip <- x$w
    n <- 2*length(Ip)
  } else {
    w <- x$w
    n <- length(w)
    Ip <- Periodogram(w)
  }
  fr <- (1/n)*(1:length(Ip))
  if (type=="log") {
    plot(fr, log(Ip), xlab="frequency", ylab="log power", 
         main=mainTitle, type="l", col=rgb(0,0,1,0.5),
         sub="logged periodogram and spectral density")
  } else {
    plot(log(fr), log(Ip), xlab="log frequency", ylab="log power", 
         main=mainTitle, type="l", col=rgb(0,0,1,0.5),
         sub="logged periodogram and spectral density")   
  }
  nplot <- max(n, 400) #0.5*nplot = number of ordinates
  y <- sigmaSq*artfimaSDF(n=nplot, obj=x, plot="none")
  fr <- (1/nplot)*(1:length(y))
  if (type=="log") {
    lines(fr, log(y), type="l", lwd=2, col="red")
  } else {
    lines(log(fr), log(y), type="l", lwd=2, col="red")
  }
}
