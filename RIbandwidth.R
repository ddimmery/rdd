#' Randomization-Inference Bandwidth Selection
#' 
#' \code{RIbandwidth} calculates the Randomization-Inference bandwidth
#' for Randomization-Inference Regression discontinuity designs.
#' 
#' @param X a numerical vector which is the running variable
#' @param Y a numerical vector which is the outcome variable. This is typically a single covariate one wishes to ensure balance on.
#' @param cutpoint the cutpoint
#' @param verbose logical flag indicating whether to print more information to the terminal. Default is \code{FALSE}.
#' @param alpha numeric, alpha-level of significance to use for rejection region
#' @param max.sim numeric, when exact randomization test is not feasible, how many simulations to perform
#' @param w.test numeric vector, window widths to test. Assumes symmetric window.
#' @param statistic
#' @return The bandwidth
#' @references Titiunik etal
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>

RIbandwidth <- function(Y,X,cutpoint=NULL,verbose=FALSE,method=c("binomial","fixed.margins"),
                        statistic=c("means","ks","wilcox"),w.test=NULL,alpha=.05,max.sim=10000) {
  sub<-complete.cases(X)&complete.cases(Y)
  X <- X[sub]
  Y <- Y[sub]
  Nx<-length(X)
  Ny<-length(Y)
  binary<- sum(Y==1 | Y==0) == length(Y)
  
  if(Nx!=Ny)
    stop("Running and outcome variable must be of equal length")
  if(is.null(cutpoint)) {
    cutpoint<-0
    if(verbose) cat("Using default cutpoint of zero.\n")
  } else {
    if(! (typeof(cutpoint) %in% c("integer","double")))
      stop("Cutpoint must be of a numeric type")
  }

  # Sort to the right
  indxr = order(X[X>=cutpoint], decreasing=FALSE)
  Yr = Y[X>=cutpoint][indxr]
  rr = X[X>=cutpoint][indxr]
  maxr = length(Yr)
  
  # Sort to the left
  indxl = order(X[X<c], decreasing=TRUE)
  Yl = Y[X<cutpoint][indxl]
  rl = X[X<cutpoint][indxl]
  maxl = length(Yl)
  mx = min(maxr, maxl)
  
  #Need progress bar
  function(w) {
    wr <- cutpoint + w
    wl <- cutpoint - w
    ir <- rr <= wr
    il <- rl <= wl
    nr <- sum(ir)
    nl <- sum(il)
    n <- nr + nl
    if( nl == 0 | nr ==0) return(NULL)
    ri<-RIestimate(y0 = Yl[il], Y1 = Yr[ir], M = M, method = method, statistic = statistic)
    
    out<-c(nl,nr,wl,wr)
    
  }
  #Begin RI window selection
  
  return(bw)
}