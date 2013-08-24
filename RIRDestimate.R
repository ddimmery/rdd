RIRDestimate<-function(formula,w.test,cutpoint,bws=NULL,w.q=c(0,.05),method=c("binomial","fixed.margins"),
                       statistic=c("means","ks","wilcox"),...) {
  call<-match.call()
  if(missing(data)) data<-environment(formula)
  formula<-as.Formula(formula)
  X<-model.frame(formula,rhs=1,lhs=0,data=data,na.action=na.pass)[[1]]
  Y<-model.frame(formula,rhs=0,lhs=NULL,data=data,na.action=na.pass)[[1]]
  if(!is.null(subset)){
    X<-X[subset]
    Y<-Y[subset]
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
    covs<-model.frame(formula,rhs=2,lhs=0,data=data,na.action=na.pass)
    if(!is.null(subset)) covs<-subset(covs,subset)
    covs<-subset(covs,na.ok)
  }
  X<-X[na.ok]-cutpoint
  Y<-Y[na.ok]
  if(type=="fuzzy") Z<-as.double(Z[na.ok])
  
  doBN <- "binomial"%in%method
  # Do exact fixed margins inference
  doFM <- "fixed.margins"%in%method
  
  # Do difference in means tests
  doMeans <- "means" %in% statistic
  # Do Kolmogorov-Smirnov tests
  doKS <- "ks" %in% statistic
  # Do Wilcoxon tests
  doWX <- "wilcox" %in% statistic
  
  estit<-function(bw) {
    Y0<-Y[X < 0 & X >= -bw]
    Y1<-Y[X >= 0 & X <= bw]
    RIestimate(Y0,Y1,method=method,statistic=statistic,...)
  }
  
  if(is.null(bws)) {
    wins<-apply(covs,2,RIbandwidth,w.test=w.test,method=method,statistic=statistic,...) # Rewrite the next bit just directly using the doKS etc stuff
    if(missing(cutpoint)) cutpoint=mean(wins[[1]][1:2])
    getwins<-function(x,type){
      wins<-sapply(x,function(x) return(x[[type]][1:2])) #"window", liberal or conservative bands
      wins<-abs(wins[1,]-wins[2,])/2
      return(wins)
    }
    winwin<-sapply(wins,getwins,type="window")
    libwin<-sapply(wins,getwins,type="liberal")
    conwin<-sapply(wins,getwins,type="conservative")
    out<-list(lapply(winwin,estit),liberal=lapply(libwin,estit),conservative=lapply(conwin,estit))
    )
  } else {
    out<-lapply(bws,estit)
  }
  class(out)<-"RIRD"
  return(out)
}