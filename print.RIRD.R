print.RIRD<-function(x,digits=max(3, getOption("digits") - 3),...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  xo<-x
  x<-x[-length(x)] # Get rid of the call which we don't care about now

  cat("Output:\n")
  
  ostr<-"diff.mean"
  # put in the string to get the right statistic, ie diff.mean
  o<-NULL

  
  for(out in x) {
    est<-out[[ostr]]
    bw<-out$bw
    ln<-c(est,bw)
    o<-rbind(o,ln)
    rownames(o)[nrow(o)]<-out$name
  }
  colnames(o) <- c("Est.","bw")

  print.default(o,quote=FALSE,digits=digits)
  invisible(xo)
  
}
