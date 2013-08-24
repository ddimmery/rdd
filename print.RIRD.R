print.RIRD<-function(x,digits=max(3, getOption("digits") - 3),...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
      "\n\n", sep = "")
  
  cat("Coefficients:\n")
  print.default(format(x$est,digits = digits),print.gap=2,quote=FALSE)
  cat("\n")
  invisible(x)
  
}