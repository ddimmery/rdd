#' Summarizing Regression Discontinuity Designs
#' 
#' \code{summary} method for class \code{"RD"}
#' 
#' @method summary RD
#' @param object an object of class \code{"RD"}, usually a result of a call to \code{\link{RDestimate}}
#' @param ... unused
#' @include RDestimate.R
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>

summary.RD<-function(object,...){
  cat("\n")
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  cat("Type:\n")
  cat(object$type,"\n\n")
  
  cat("Bandwidth:\n")
  cat(object$bw,"\n\n")
  
  #If the model wasn't included in the output, we need to get it
  mod<-FALSE
  if("model" %in% names(object$call)) mod<-object$call$model
  if(!mod){
    object$call$model<-TRUE
    object$call$verbose<-FALSE
    object<-eval(object$call)
  }
  
  cat("Observations within bandwidth:\n")
  if(object$type=="sharp") 
    cat(length(residuals(object$model)),"\n\n")
  else 
    cat(length(residuals(object$model$iv)),"\n\n")
  
  cat("Estimate:\n")
  #Need to get this to give at least as much as stata does in fuzzy designs
  p<-if(object$p<.0001) "<0.0001" else sprintf("%.4f",object$p)
  stars<-if (object$p<0.001) "***" else if(object$p<0.01) "**" else if(object$p<0.05) "*" else if(object$p<0.1) "." else ""
  out<-matrix(c(sprintf("%0.3f",object$est),sprintf("%0.3f",object$se),sprintf("%.3f",object$z),p,stars),nrow=1)
  rownames(out)<-"LATE"
  colnames(out)<-c("Estimate","Std. Error","z value","Pr(>|z|)","")
  print(out,digits=4,quote=FALSE)
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1\n\n")
  if(object$type=="sharp") {
    p<-pf(summary(object$model)$fstatistic[[1]],
          summary(object$model)$fstatistic[[2]],
          summary(object$model)$fstatistic[[3]])
    p<-2*min(p,1-p)
    p<-if(object$p<.00001) "<0.00001" else sprintf("%.4f",p)
    cat("F-statistic:",
    summary(object$model)$fstatistic[[1]],
    "on",
    summary(object$model)$fstatistic[[2]],
    "and",
    summary(object$model)$fstatistic[[3]],
    "DoF, p-value:",
    p,"\n\n")

  } else {
    p<-summary(object$model$iv)$waldtest[2]
    p<-if(object$p<.00001) "<0.00001" else sprintf("%.4f",p)
    cat("F-statistic:",
        summary(object$model$iv)$waldtest[1],
        "on",
        summary(object$model$iv)$waldtest[3],
        summary(object$model$iv)$waldtest[4],
        "DoF, p-value:",
        p,"\n\n")
  }
}