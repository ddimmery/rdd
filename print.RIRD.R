print.RIRD<-function(x,Y,X,digits=max(3, getOption("digits") - 3),...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Output:\n")
  
  ostr<-"diff.mean"
  # put in the string to get the right statistic, ie diff.mean
  o<-matrix(NA,nrow=length(x)+1,ncol=10)
  colnames(o)
  for(out in x) {
    prob<-x$prob
    est<-x[[ostr]]
    bw<-x$bw
    ci.null<-x$t.test$conf.int
    BN.p<-x$meansBN.pval
    FM.p<-x$meansFM.pval
    ci.inv<-getCI(Y,X,bw)

    ln<-c(est,bw,prob,BN.p,FM.p,ci.null,ci.inv$ci,ci.inv$alpha)
    format(ln,digits = digits)
    cat("\n")
  }
  invisible(x)
  
}
