#Source: bestmodels.R
#includes functions: best_glp_models, best_models_serial, best_models_parll,
# KEval
best_glp_models <- function(z, glp = c("ARTFIMA", "ARFIMA", "ARIMA"), p=2, q=2, 
                            likAlg = c("exact", "Whittle")) {
  stopifnot(p>0 || q>0)
  stopifnot(is.numeric(z))
  stopifnot(length(z)>=50)
  ps <- 0:p
  qs <- 0:q
  LL_model <- matrix(numeric(0), nrow=p+1, ncol=q+1)
  bicPenalty <- aicPenalty <- LL_model
  rn <- paste0(paste0("AR(", ps),")")
  cn <- paste0(paste0("MA(", qs),")")
  rownames(LL_model) <- rn
  colnames(LL_model) <- cn
  startTime <- proc.time()[3]
  glpOrder <-switch(glp, "ARTFIMA"=2, "ARFIMA"=1, "ARIMA"=0)
  for (i in 1:(p+1)) {
    for (j in 1:(q+1)) {
      LL_model[i,j] <- artfima(z, glp=glp, arimaOrder=c(ps[i], 0, qs[j]),
                               likAlg=likAlg)$LL
      aicPenalty[i,j] <- 2*(ps[i]+qs[j]+glpOrder+1)
      bicPenalty[i,j] <- log(length(z))*(ps[i]+qs[j]+glpOrder+1)
    }
  }
  totTime <- proc.time()[3] - startTime
  #
  aic_model <- (-2*LL_model) + aicPenalty
  indaic <- min(aic_model)==aic_model
  pOpt = ps[as.logical(apply(indaic, 1, sum))]
  qOpt = qs[as.logical(apply(indaic, 2, sum))]
  bestaicModel <- paste0(glp, "(", pOpt, ",", qOpt, ")" )
  bestaic <- min(aic_model)
  #
  bic_model <- (-2*LL_model) + bicPenalty
  indbic <- min(bic_model)==bic_model
  pOpt = ps[as.logical(apply(indbic, 1, sum))]
  qOpt = qs[as.logical(apply(indbic, 2, sum))]
  bestbicModel <- paste0(glp, "(", pOpt, ",", qOpt, ")" )
  bestbic <- min(bic_model)
  ans <- list(LL=LL_model, artfima_time = totTime, 
              aic=list(bestaic=bestaic, bestaicModel=bestaicModel, 
                       aic=aic_model),
              bic=list(bestbic=bestbic, bestbicModel=bestbicModel, 
                       bic=bic_model)
  )
  ans
}


best_models_serial <- function(z, p=2, q=2, nbest=5, 
            likAlg = c("exact", "Whittle"), plausibilityQ=TRUE) {
  outARIMA <- best_glp_models(z, glp="ARIMA", p=p, q=q, likAlg=likAlg)
  outARFIMA <- best_glp_models(z, glp="ARFIMA", p=p, q=q, likAlg=likAlg)
  outARTFIMA <- best_glp_models(z, glp="ARTFIMA", p=p, q=q, likAlg=likAlg)
  #
  aics <- c(as.vector(outARIMA$aic$aic), as.vector(outARFIMA$aic$aic),
            as.vector(outARTFIMA$aic$aic))
  bics <- c(as.vector(outARIMA$bic$bic), as.vector(outARFIMA$bic$bic),
            as.vector(outARTFIMA$bic$bic))
  p1 <- nrow(outARIMA$bic$bic)
  q1 <- ncol(outARIMA$bic$bic)
  armaLab1 <- as.vector(outer(0:(p1-1), 0:(q1-1), FUN=function(x,y) 
    paste0("ARIMA(", x, ",0,", y, ")")))
  p1 <- nrow(outARFIMA$bic$bic)
  q1 <- ncol(outARFIMA$bic$bic)
  armaLab2 <- as.vector(outer(0:(p1-1), 0:(q1-1), FUN=function(x,y) 
    paste0("ARFIMA(", x, ",0,", y, ")")))
  p1 <- nrow(outARTFIMA$bic$bic)
  q1 <- ncol(outARTFIMA$bic$bic)
  armaLab3 <- as.vector(outer(0:(p1-1), 0:(q1-1), FUN=function(x,y) 
    paste0("ARTFIMA(", x, ",0,", y, ")")))
  armaLab <- c(armaLab1, armaLab2, armaLab3)
  ind <- order(aics)
  aics <- aics[ind]
  armaLab <- armaLab[ind]
  names(aics) <- armaLab
  aics <- aics[1:nbest]
  armaLab <- c(armaLab1, armaLab2, armaLab3)
  ind <- order(bics)
  bics <- bics[ind]
  armaLab <- armaLab[ind]
  names(bics) <- armaLab
  bics <- bics[1:nbest]
  ans <- list(aic=aics, bic=bics)
  class(ans) <- "bestmodels"
  ans 
}


best_models_parll <- function(z, p=2, q=2, nbest=5, 
            likAlg = c("exact", "Whittle")) {
  KEval <- function(k, p, q) {#scoping used to pass variables
    ans <- NULL
    if (k==1) {
      ans <- best_glp_models(z, glp="ARIMA", p=p, q=q, likAlg=likAlg)
      ans <- c(ans$aic$aic, ans$bic$bic)
      return(ans)
    } 
    if (k==2) {
      ans <- best_glp_models(z, glp="ARFIMA", p=p, q=q, likAlg=likAlg)
      ans <- c(ans$aic$aic, ans$bic$bic)
      return(ans)
    } 
    if (k==3) {
      ans <- best_glp_models(z, glp="ARTFIMA", p=p, q=q, likAlg=likAlg)
      ans <- c(ans$aic$aic, ans$bic$bic)
      return(ans)
    } 
    ans
  }
  cl <- makeCluster(spec=3, type="PSOCK")
  #Export variables
  clusterExport(cl, list("z", "likAlg"), envir=environment())
  clusterExport(cl, list("best_glp_models", "likAlg"), envir=environment())
  clusterEvalQ(cl, library("artfima"))
  #clusterExport(cl, list("z", "likAlg"))
  #clusterExport(cl, list("best_glp_models", "likAlg"))
  #clusterEvalQ(cl, library("artfima"))
  out <- parSapply(cl, 1:3, KEval, p=p, q=q)
  stopCluster(cl)
  armaLab <- as.vector(outer(0:p, 0:q, FUN=function(x,y) 
    paste0("(", x, ",0,", y, ")")))
  ind <- 1:(nrow(out)/2)
  outaic = out[ind,]
  outbic = out[-ind,]
  #
  dimnames(outaic) <- list(armaLab, c("ARIMA", "ARFIMA", "ARTFIMA"))
  aic <- as.vector(outaic)
  aic_names <- c(paste0(colnames(outaic)[1], rownames(outaic)),
                 paste0(colnames(outaic)[2], rownames(outaic)),
                 paste0(colnames(outaic)[3], rownames(outaic)))
  names(aic) <- aic_names
  bestAic <- sort(aic)[1:nbest]
  #
  dimnames(outbic) <- list(armaLab, c("ARIMA", "ARFIMA", "ARTFIMA"))
  bic <- as.vector(outbic)
  bic_names <- c(paste0(colnames(outbic)[1], rownames(outbic)),
                 paste0(colnames(outbic)[2], rownames(outbic)),
                 paste0(colnames(outbic)[3], rownames(outbic)))
  names(bic) <- bic_names
  bestBic <- sort(bic)[1:nbest]
  #
  ans <- list(aic=bestAic, bic=bestBic)
  class(ans) <- "bestmodels"
  ans
}

bestModels <- function(z, p=2, q=2, nbest=5, likAlg = c("exact", "Whittle"), 
                        use_parallel=TRUE) {
  stopifnot(is.numeric(z))
  stopifnot(length(z)>=50)
  if (use_parallel) {
    ans <- best_models_parll(z, p, q, nbest, likAlg)
  } else {
    ans <- best_models_serial(z, p, q, nbest, likAlg)  
  }
  ans
}
