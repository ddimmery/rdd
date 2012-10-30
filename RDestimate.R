#' Regression Discontinuity Estimation
#' 
#' \code{RDestimate} supports both sharp and fuzzy RDD
#' utilizing the \pkg{AER} package for 2SLS regression
#' under the fuzzy design. Local linear regressions are performed
#' to either side of the cutpoint using the Imbens-Kalyanaraman
#' optimal bandwidth calculation, \code{\link{IKbandwidth}}.
#' 
#' @param formula the formula of the RDD. This is supplied in the
#' format of \code{y ~ x} for a simple sharp RDD, or \code{y ~ x | c1 + c2}
#' for a sharp RDD with two covariates. Fuzzy RDD may be specified as
#' \code{y ~ x + z} where \code{x} is the running variable, and 
#' \code{z} is the endogenous treatment variable. Covariates are then included in the 
#' same manner as in a sharp RDD.
#' @param data an optional data frame
#' @param subset an optional vector specifying a subset of observations to be used
#' @param cutpoint the cutpoint. If omitted, it is assumed to be 0.
#' @param bw the bandwidth. If omitted, it is calculated using the Imbens-Kalyanaraman method.
#' @param kernel a string specifying the kernel to be used in the local linear fitting. 
#' \code{"triangular"} kernel is the default and is the "correct" theoretical kernel to be used for 
#' edge estimation as in RDD (Lee and Lemieux 2010). Other options are \code{"rectangular"}, 
#' \code{"epanechnikov"}, \code{"quartic"}, 
#' \code{"triweight"}, \code{"tricube"}, \code{"gaussian"} and \code{"cosine"}.
#' @param se.type this specifies the robust SE calculation method to use. Options are,
#' as in \code{\link{vcovHC}}, \code{"HC3"}, \code{"const"}, \code{"HC"}, \code{"HC0"}, 
#' \code{"HC1"}, \code{"HC2"}, \code{"HC4"}, \code{"HC4m"}, \code{"HC5"}. This option 
#' is overriden by \code{cluster}.
#' @param cluster an optional vector specifying clusters within which the errors are assumed
#' to be correlated. This will result in reporting cluster robust SEs. This option overrides
#' anything specified in \code{se.type}. It is suggested that data with a discrete running 
#' variable be clustered by each unique value of the running variable (Lee and Card 2008).
#' @param verbose will provide some additional information printed to the terminal.
#' @param frame logical. If \code{TRUE}, the data frame used in model fitting will be returned.
#' @param model logical. If \code{TRUE}, the model object will be returned.
#' @return \code{RDestimate} returns an object of \link{class} "\code{RD}".
#' The functions \code{summary} and \code{plot} are used to obtain and print a summary and plot of
#' the estimated regression discontinuity. The object of class \code{RD} is a list 
#' containing the following components:
#' \item{type}{a string denoting either \code{"sharp"} or \code{"fuzzy"} RDD.}
#' \item{est}{the estimate of the discontinuity in the outcome under a sharp design, 
#' or the Wald estimator in the fuzzy design}
#' \item{se}{the standard error}
#' \item{z}{the z statistic}
#' \item{p}{the p value}
#' \item{ci}{the vector of the 95% confidence interval, \code{c("CI Lower Bound","CI Upper Bound")}}
#' \item{bw}{the bandwidth}
#' \item{call}{the matched call}
#' \item{na.action}{the observations removed from fitting due to missingness}
#' \item{model}{(if requested) For a sharp design, the \code{lm} object is returned.
#' For a fuzzy design, a list is returned with two elements: \code{firststage}, the first stage \code{lm}
#' object, and \code{iv}, the \code{ivreg} object.}
#' \item{frame}{(if requested) Returns the model frame used in fitting.}
#' @seealso \code{\link{summary.RD}}, \code{\link{plot.RD}}, \code{\link{DCdensity}} \code{\link{IKbandwidth}}, \code{\link{kernelwts}}, \code{\link{vcovHC}}, 
#' \code{\link{ivreg}}, \code{\link{lm}}
#' @references Lee, David and Thomas Lemieux. (2010) "Regression Discontinuity Designs in Economics," \emph{Journal of Economic Literature}. 48(2): 281-355. \url{http://www.aeaweb.org/articles.php?doi=10.1257/jel.48.2.281}
#' @references Imbens, Guido and Thomas Lemieux. (2010) "Regression discontinuity designs: A guide to practice," \emph{Journal of Econometrics}. 142(2): 615-635. \url{http://dx.doi.org/10.1016/j.jeconom.2007.05.001}
#' @references Lee, David and David Card. (2010) "Regression discontinuity inference with specification error," \emph{Journal of Econometrics}. 142(2): 655-674. \url{http://dx.doi.org/10.1016/j.jeconom.2007.05.003}
#' @references Angrist, Joshua and Jorn-Steffen Pischke. (2009) \emph{Mostly Harmless Econometrics}. Princeton: Princeton University Press.
#' @import AER sandwich lmtest
#' @include IKbandwidth.R kernelwts.R
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>
#' @examples
#' x<-runif(1000,-1,1)
#' cov<-rnorm(1000)
#' y<-3+2*x+3*cov+10*(x>=0)+rnorm(1000)
#' RDestimate(y~x)
#' # Efficiency gains can be made by including covariates
#' RDestimate(y~x|cov)

RDestimate<-function(formula, data, subset=NULL, cutpoint=NULL, bw=NULL, kernel="triangular", se.type="HC1", cluster=NULL, verbose=FALSE, model=FALSE, frame=FALSE) {
  call<-match.call()
  if(missing(data)) data<-environment(formula)
  formula<-as.Formula(formula)
  X<-model.frame(formula,rhs=1,lhs=0,data=data,na.action=na.pass)[[1]]
  Y<-model.frame(formula,rhs=0,lhs=NULL,data=data,na.action=na.pass)[[1]]
  if(!is.null(subset)){
    X<-X[subset]
    Y<-Y[subset]
  }
  if (!is.null(cluster)) {
    robust.se <- function(model, cluster){
      M <- length(unique(cluster))
      N <- length(cluster)           
      K <- model$rank
      dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
      uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
      rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
      rcse.se <- coeftest(model, rcse.cov)
      return(list(rcse.cov, rcse.se))
    }
  }
  na.ok<-complete.cases(X)&complete.cases(Y)
  if(length(all.vars(formula(formula,rhs=1,lhs=F)))>1){
   type<-"fuzzy" 
   Z<-model.frame(formula,rhs=1,lhs=0,data=data,na.action=na.pass)[[2]]
   if(!is.null(subset)) Z<-Z[subset]
   na.ok<-na.ok&complete.cases(Z)
   if(length(all.vars(formula(formula,rhs=1,lhs=F)))>2)
     stop("Invalid formula. Read ?RDestimate for proper syntax")
  } else {
   type="sharp" 
  }
  covs<-NULL
  if(length(formula)[2]>1){
    covs<-model.frame(formula,rhs=2,lhs=0,data=data)
    if(!is.null(subset)) covs<-subset(covs,subset)
    na.ok<-na.ok&complete.cases(covs)
    covs<-subset(covs,na.ok)
  }
  X<-X[na.ok]
  Y<-Y[na.ok]
  if(type=="fuzzy") Z<-Z[na.ok]

  if(is.null(cutpoint)) {
    cutpoint<-0
    if(verbose) cat("No cutpoint provided. Using default cutpoint of zero.\n")
  }
  if(frame) {
    if(type=="sharp") {
      if (!is.null(covs))
        dat.out<-data.frame(X,Y,covs)
      else
        dat.out<-data.frame(X,Y)
    } else {
      if (!is.null(covs))
        dat.out<-data.frame(X,Y,Z,covs)
      else
        dat.out<-data.frame(X,Y,Z)
    }
  }
  if(is.null(bw)) bw<-IKbandwidth(X=X,Y=Y,cutpoint=cutpoint,kernel=kernel, verbose=verbose)
  #Subset to within the bandwidth, except for when using gaussian weighting
  sub<-X>=(cutpoint-bw) & X<=(cutpoint+bw)
  
  if(kernel=="gaussian") 
    sub<-TRUE
  Y<-Y[sub]
  X<-X[sub]

  w<-kernelwts(X,cutpoint,bw,kernel=kernel)
               
  Xl<-(X<cutpoint)*X
  Xr<-(X>=cutpoint)*X
  Tr<-as.integer(X>=cutpoint)
  if(type=="sharp"){
    if(verbose) {
      cat("Running Sharp RD\n")
      cat("Running variable:",all.vars(formula(formula,rhs=1,lhs=F))[1],"\n")
      cat("Outcome variable:",all.vars(formula(formula,rhs=F,lhs=1))[1],"\n")
      if(!is.null(covs)) cat("Covariates:",paste(names(covs),collapse=", "),"\n")
    }
    if(!is.null(covs)) {
      covs<-subset(covs,sub)
      data<-data.frame(Y,Tr,Xl,Xr,covs,w)
      form<-as.formula(paste("Y~Tr+Xl+Xr+",paste("Tr*",names(covs),collapse="+",sep=""),sep=""))
    } else {
      data<-data.frame(Y,Tr,Xl,Xr,w)
      form<-as.formula(Y~Tr+Xl+Xr)
    }

    mod<-lm(form,weights=w,data=data)
    est<-coef(mod)["Tr"]
    names(est)<-"LATE Estimate"
    if (is.null(cluster)) {
      se<-coeftest(mod,vcovHC(mod,type=se.type))[2,2]
    } else {
      se<-robust.se(mod,cluster)[[2]][2,2]
    }
    names(se)<-"LATE SE"
    z<-est/se
    names(z)<-"z statistic"
    p<-2*pnorm(abs(z),lower.tail=F)
    names(p)<-"P>|Z|"
    ci<-c(est-qnorm(.975)*se,est+qnorm(.975)*se)
    names(ci)<-c("CI Lower Bound","CI Upper Bound")
    o<-list(type=type,est=est,se=se,z=z,p=p,ci=ci,bw=bw,call=call,na.action=(which(na.ok==FALSE)))
    if(model) o$model=mod
    if(frame) o$frame=dat.out
    class(o)<-"RD"
    return(o)
  } else {
    if(verbose){
      cat("Running Fuzzy RD\n")
      #CLEAN UP
      cat("Running variable:",all.vars(formula(formula,rhs=1,lhs=F))[1],"\n")
      cat("Outcome variable:",all.vars(formula(formula,rhs=F,lhs=1))[1],"\n")
      cat("Treatment variable:",all.vars(formula(formula,rhs=1,lhs=F))[2],"\n")
      if(!is.null(covs)) cat("Covariates:",paste(names(covs),collapse=", "),"\n")
    }
    Z<-Z[sub]
    if(!is.null(covs)) {
      covs<-subset(covs,sub)
      data<-data.frame(Y,Tr,Xl,Xr,Z,covs,w)
      form<-as.Formula(paste(
              "Y~Z+Xl+Xr+",paste(names(covs),collapse="+"),
              "|Tr+Xl+Xr+",paste(names(covs),collapse="+"),sep=""))
      form1<-as.Formula(paste("Z~Tr+Xl+Xr+",paste("Tr*",names(covs),collapse="+",sep="")))
    } else {
      data<-data.frame(Y,Tr,Xl,Xr,Z,w)
      form<-as.Formula(Y~Z+Xl+Xr|Tr+Xl+Xr)
      form1<-as.formula(Z~Tr+Xl+Xr)
    }
    mod1<-lm(form1,weights=w,data=data)
    mod<-ivreg(form,weights=w,data=data)
    if(verbose==TRUE) {
      cat("First stage:\n")
      print(summary(mod1))
      cat("IV-RD\n")
      print(summary(mod))
    }
    est<-coef(mod)["Z"]
    names(est)<-"LATE Estimate"
    if (is.null(cluster)) {
      se<-coeftest(mod,vcovHC(mod,type=se.type))[2,2]
    } else {
      se<-robust.se(mod,cluster)[[2]][2,2]
    }
    names(se)<-"LATE SE"
    z<-est/se
    names(z)<-"z statistic"
    p<-2*pnorm(abs(z),lower.tail=F)
    names(p)<-"P>|Z|"
    ci<-c(est-qnorm(.975)*se,est+qnorm(.975)*se)
    names(ci)<-c("CI Lower Bound","CI Upper Bound")
    o<-list(type=type,est=est,se=se,z=z,p=p,ci=ci,bw=bw,call=call,na.action=(which(na.ok==FALSE)))
    if(model) o$model=list(firststage=mod1,iv=mod)
    if(frame) o$frame=dat.out
    class(o)<-"RD"
    return(o)
  }
}