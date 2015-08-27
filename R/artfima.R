artfima <-
function(z, glp=c("ARTFIMA", "ARFIMA", "ARIMA"), armaOrder=c(0,0,0), 
                   constant=TRUE, seQ=TRUE, likAlg=c("Exact","Whittle"), 
                   blueQ=FALSE) {
  alg <- 1 # default optimization is "L-BFGS-B"
  glp <- match.arg(glp)
  likAlg <- match.arg(likAlg)
  p <- armaOrder[1]
  d0 <- armaOrder[2]
  q <- armaOrder[3]
  glpOrder <-switch(glp, "ARTFIMA"=2, "ARFIMA"=1, "ARIMA"=0)
  stopifnot(all(armaOrder>=0)) 
  is.wholenumber <- function(x) abs(x - round(x)) < .Machine$double.eps^0.5
  stopifnot(is.wholenumber(p), is.wholenumber(d0), is.wholenumber(q))
#
  w <- if(d0>0) diff(z, differences=d0) else z
  if(constant) mnw <- mean(w) else mnw<-0
  varw <- var(w)
  w <- w-mnw
  n <- length(w)
#initialization
  nbeta <- p+q+glpOrder
  binit <- numeric(nbeta)
#Whittle method is fast because this is only done once
  if (likAlg=="Whittle")
    Ip <- (spec.pgram(w, fast=FALSE, detrend=FALSE, plot=FALSE, taper=0)$spec)/(2*pi)  
  nullModelLoglikelihood <- (-n/2)*log(sum(w^2)/n)
  entropyPenalty <-  switch(likAlg,
      Exact=-nullModelLoglikelihood,
      Whittle=sum(w^2)
  )
  entropyPenalty <- entropyPenalty+2*abs(entropyPenalty)
#negative log-likelihood=entropy
#begin Entropy. beta - pacf parameterization
#note - ***parameters passed using lexical scoping***
#also tacvf r is exported to arfima environment
  count <- 0
  Entropy<-function(beta) {
#in the optimization, put lambda, d, phi, theta
    phi<-theta<-lambda <- d <- numeric(0)
    count <<- count+1
    if (glpOrder==2) {
      lambda <- beta[1]
      d <- beta[2]
    } else {
      if (glpOrder==1) d <- beta[1]
        }
    if(p>0) phi <- PacfToAR(beta[(1+glpOrder):(p+glpOrder)])
    if(q>0) theta <- PacfToAR(beta[(p+glpOrder+1):(p+q+glpOrder)]) 
    if (likAlg=="Exact") {
      r <- tacvfARTFIMA(d=d, lambda=lambda, phi = phi, theta = theta, 
                        maxlag = n-1) 
#needed in some cases, eg. NileMin with p=1, q=3, glp="FGN"
      if (any(is.na(r))) {
        negLL <- NA
      } else {
        negLL <- try(-DLLoglikelihood(r, w), silent=TRUE)
      }
      negLL <- ifelse(is.numeric(negLL), negLL, entropyPenalty)
    } else {
      fp <- sdfartfima(n=n, d=d, lambda=lambda, phi=phi, theta=theta)
      negLL <- 2*mean(Ip/fp)
    }
#cat("\n ***********iter = ", iter, fill=TRUE)
#cat("count = ", count, fill=TRUE)
#cat("negLL=", negLL, fill=TRUE)
#cat("beta=",beta,fill=TRUE)
#cat("r[1:10] = ", r[1:10], fill=TRUE)
#cat("w[1:10] = ", w[1:10], fill=TRUE)
  negLL
  }#end Entropy
#
#lower and upper limits with "L-BFGS-B"
  lambdaLo <- 0.001
  lambdaHi <- 3
  dHi <- 2
  dfHi <- 0.49
  if (glp=="ARTFIMA") {
    blo<- c(lambdaLo, -dHi, rep(-0.99,p+q))
    bhi<- c(lambdaHi, dHi,  rep(0.99,p+q))
  } else {
      if (glp=="ARFIMA") {
        blo<- c(-dfHi, rep(-0.99,p+q))
      } else { #ARMA
          blo<- rep(-0.99,p+q) 
        }
      bhi <- -blo
  }
  
#while mean not converged############################################################ 
  w0 <- w
  etol <- maxIter <- 1
  if(blueQ) maxIter <- 5
  iter <- meanMLE <- 0
  while(etol> 1e-06 && iter<maxIter){
    iter<-iter+1
#trace=6 for full output
    ans<-optim(par=binit, fn=Entropy, method="L-BFGS-B",
         lower=blo, upper=bhi, control=list(trace=0), hessian=seQ)
    if(ans$convergence != 0) {#convergence problem. Use Nelder-Mead with penalty function
      alg<-2
      ans<-optim(par=binit, fn=Entropy, method="Nelder-Mead", hessian=seQ)
      if(ans$convergence != 0) {#convergence problem. Use SANN with penalty function
        alg<-3
       ans<-optim(par=binit, fn=Entropy, method="SANN", hessian=seQ)
      }
    }
    negLL <- ans$value
    bHat <- ans$par
    binit <- bHat
    lambdaHat <- dHat <- phiHat <- thetaHat <- numeric(0)
    onBoundary <- FALSE
    if (glpOrder > 0) dHat <- bHat[glpOrder]
    if (glpOrder == 1 && abs(dHat) > dfHi) onBoundary <- TRUE
    if (glpOrder==2) {
      dHat <- bHat[glpOrder]
      lambdaHat <-bHat[1]
      if (lambdaHat>=lambdaHi || lambdaHat <= lambdaLo) onBoundary <- TRUE
      if (abs(dHat) >= dHi) onBoundary <- TRUE
    }
    if (p > 0) phiHat <- PacfToAR(bHat[(1+glpOrder):(p+glpOrder)])
    if (q > 0) thetaHat <- PacfToAR(bHat[(p+1+glpOrder):(p+q+glpOrder)])
    rHat <- tacvfARTFIMA(d=dHat, lambda=lambdaHat, 
                         phi = phiHat, theta = thetaHat, maxlag = n-1)
    if(maxIter > 1) {#if MaxIt==0, sample mean is used
      meanMLEPrev <- meanMLE
      meanMLE <- TrenchMean(rHat, w0)
      w <- w0-meanMLE
      etol <- abs(meanMLE-meanMLEPrev)/(abs(meanMLE)+0.01)
    }
  }
  convergence <- ans$convergence
#end while###########################################################################  
#since Entropy is used, Hessian is pd
  if(seQ&&likAlg=="Exact") {
    Hinv<-try(solve(ans$hessian),silent=TRUE)
    if(!all(is.numeric(Hinv))) Hinv <- NA
    if(!any(is.na(Hinv))){
      sebHat <- suppressWarnings(sqrt(diag(Hinv)))} else {
          sebHat<-rep(NA,length(bHat))}
      } else {#either not requested or Whittle used
          sebHat <- rep(NA,length(bHat))
      }
  rhoHat <- rHat[-1]/rHat[1]
  seMean <- sqrt(var(w)/n*(1+2*sum(1-(1:(n-1))/n*rhoHat^2))/n)
  constantHat <- mnw+meanMLE
  ansExact <- exactLoglikelihood(rHat, w)
  LL <- ansExact$LL
  sigmaSq <- ansExact$sigmaSq
  snr <- (varw-sigmaSq)/sigmaSq
  res <- DLResiduals(rHat, w)
  out<-list(dHat=dHat, lambdaHat=lambdaHat, phiHat=phiHat, thetaHat=thetaHat, 
            constant=constantHat, seMean=seMean, se=sebHat, n=n, 
            sigmaSq=sigmaSq, snr=snr, likAlg=likAlg, convergence=convergence, 
            blueQ=blueQ, LL=LL, algorithm=alg, constantQ=constant, glp=glp, 
            armaOrder=armaOrder, glpOrder=glpOrder, tacvf=rHat, res=res,
            nullModelLogLik=nullModelLoglikelihood, onBoundary=onBoundary)
 class(out) <- "artfima"
 out
}
