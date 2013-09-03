summary.RIRD<-function(x,Y,X,digits=max(3, getOption("digits") - 3),
                       method="binomial",...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  x<-x[-length(x)] # Get rid of the call which we don't care about now

  cat("Output:\n")
  
  ostr<-"diff.mean"
  # put in the string to get the right statistic, ie diff.mean
  o<-NULL
  for(out in x) {
    prob<-out$prob
    est<-out[[ostr]]
    bw<-out$bw
    # ci.null<-out$t.test$conf.int
    BN.p<-out$meansBN.pval
    FM.p<-out$meansFM.pval

    ci.inv<-getRICI(Y=Y,X=X,bw=bw,method=method,...)

    ln<-c(est,bw,prob,BN.p,FM.p,ci.inv$ci,ci.inv$alpha)
    o<-rbind(o,ln)
    rownames(o)[nrow(o)]<-out$name 
  }
  colnames(o)<-c("Est.","bw","Pr(Treatment)","p (BN RI)","p (FM RI)",
                 "CI lower","CI upper",
                 "alpha")

  print.default(o,quote=FALSE,digits=digits)
  invisible(o)
  
}
