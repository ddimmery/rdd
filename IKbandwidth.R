#' Imbens-Kalyanaraman Optimal Bandwidth Calculation
#' 
#' \code{IKbandwidth} calculates the Imbens-Kalyanaraman optimal bandwidth
#' for local linear regression in Regression discontinuity designs.
#' 
#' @param X a numerical vector which is the running variable
#' @param Y a numerical vector which is the outcome variable
#' @param cutpoint the cutpoint
#' @param verbose logical flag indicating whether to print more information to the terminal. Default is \code{FALSE}.
#' @param kernel string indicating which kernel to use. Options are \code{"triangular"} (default and recommended), 
#' \code{"rectangular"}, \code{"epanechnikov"}, \code{"quartic"}, \code{"triweight"}, \code{"tricube"}, 
#' \code{"gaussian"}, and \code{"cosine"}.
#' @return The optimal bandwidth
#' @references Imbens, Guido and Karthik Kalyanaraman. (2009) "Optimal Bandwidth Choice for the regression discontinuity estimator," \emph{NBER Working Paper Series}. 14726. \url{http://www.nber.org/papers/w14726}
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>

IKbandwidth <-function (X,Y,cutpoint=NULL,verbose=FALSE,kernel="triangular") {
  #Implementation of Imbens-Kalyanaraman optimal bandwidth
  # for regression discontinuity
  if(length(X)!=length(Y))
    stop("Running and outcome variable must be of equal length")
  Nx<-length(X[complete.cases(X)])
  Ny<-length(Y[complete.cases(Y)])
  if(is.null(cutpoint)) {
    cutpoint<-0
    if(verbose) cat("Using default cutpoint of zero.\n")
  } else {
    if(! (typeof(cutpoint) %in% c("integer","double")))
      stop("Cutpoint must be of a numeric type")
  }
  #Now we should be ready to start
  #Pilot bandwidth
  h1<-1.84*sd(X)*Nx^(-1/5)
  left<-X>=(cutpoint-h1) & X<cutpoint
  right<-X>=cutpoint & X<=(cutpoint+h1)
  Nl<-sum(left)
  Nr<-sum(right)
  fbarx<-(Nl+Nr)/(2*Nx*h1)
  varY<-(sum((Y[left]-mean(Y[left],na.rm=TRUE))^2,na.rm=TRUE)+sum((Y[right]-mean(Y[right],na.rm=TRUE))^2,na.rm=TRUE)) / (sum(c(left,right)))
  left<-X<cutpoint
  right<-X>=cutpoint
  Nl<-sum(left)
  Nr<-sum(right)
  medXl<-median(X[left])
  medXr<-median(X[right])
  cX<-X-cutpoint
  if(!any(c(X[left]>medXl,X[right]<medXr)) )
    stop("Insufficient data in vicinity of the cutpoint to calculate bandwidth.")
  #Model a cubic within the pilot bandwidth
  mod<-lm(Y~I(X>=cutpoint)+poly(cX,3,raw=TRUE),subset=(X>=medXl&X<=medXr))
  m3<-6*coef(mod)[5]
  if(m3^2<.01) warning("Low curvature near the cutpoint. May be an unstable estimate.")
  #New bandwidth estimate
  h2l<-3.56*(Nl^(-1/7))*(varY/(fbarx*m3^2))^(1/7)
  h2r<-3.56*(Nr^(-1/7))*(varY/(fbarx*m3^2))^(1/7)
  left<-(X>=(cutpoint-h2l)) & (X<cutpoint)
  right<-(X>=cutpoint) & (X<= (cutpoint+h2r))
  Nl<-sum(left)
  Nr<-sum(right)
  if(Nl==0 | Nr==0)
    stop("Insufficient data in vicinity of the cutpoint to calculate bandwidth.")
  #Estimate quadratics for curvature estimation
  mod<-lm(Y~poly(cX,2,raw=TRUE),subset=right)
  m2r<-2*coef(mod)[3]
  mod<-lm(Y~poly(cX,2,raw=TRUE),subset=left)
  m2l<-2*coef(mod)[3]
  rl<-720*varY/(Nl*(h2l^4))
  rr<-720*varY/(Nr*(h2r^4))
  #Which kernel are we using?
  # Method for finding these available in I--K p. 6
  if(kernel=="triangular") {
    ck<-3.43754
  } else if (kernel=="rectangular") {
    ck<-2.70192 #5.4?
  } else if(kernel=="epanechnikov") {
    ck<-3.1999
  } else if(kernel=="quartic" | kernel=="biweight") {
    ck<-3.65362
  } else if(kernel=="triweight") {
    ck<-4.06065
  } else if(kernel=="tricube") {
    ck<-3.68765
  } else if(kernel=="gaussian") {
    ck<-1.25864
  } else if(kernel=="cosine") {
    ck<-3.25869
  } else {
    stop("Unrecognized kernel.") 
  }
  #And there's our optimal bandwidth
  optbw<-ck*(2*varY/(fbarx*((m2r-m2l)^2+rr+rl)))^(1/5)*(Nx^(-1/5))
  left<-(X>=(cutpoint-optbw)) & (X<cutpoint)
  right<-(X>=cutpoint) & (X<= (cutpoint+optbw))
  if(all(c(left,right)%in%c(0,NA)))
    stop("Insufficient data in the calculated bandwidth.")
  names(optbw)<-NULL
  if(verbose) cat("Imbens-Kalyanamaran Optimal Bandwidth: ",sprintf("%.3f",optbw),"\n")
  return(optbw)
  }
