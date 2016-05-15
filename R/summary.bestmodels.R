#summary.bestmodels
summary.bestmodels <- function(object, ...) {
  cnames <- c("best", "2nd best", "3rd best", "4th best", "5th best")
  aicN <- names(object$aic)
  aicX <- formatC(object$aic, digits=2, format="f")
  aicP <- formatC(exp(0.5*(object$aic[1] - object$aic)), digits=3, format="f")
  bicN <- names(object$bic)
  bicX <- formatC(object$bic, digits=2, format="f")
  bicP <- formatC(exp(0.5*(object$bic[1] - object$bic)), digits=3, format="f")
  m <- matrix(c(aicN,aicX,aicP,bicN,bicX,bicP), byrow=TRUE, nrow=6)
  dimnames(m) <- list(c("AIC models","AIC","p(AIC)","BIC models","BIC","p(BIC)"),
                      cnames)
  as.data.frame.matrix(m)
}


