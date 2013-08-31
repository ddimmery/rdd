RIRDestimate<-function(
                       formula,
                       w.test,
                       cutpoint=NULL,
                       bws=NULL,
                       w.q=c(0,.05),
                       method=c("binomial","fixed.margins"),
                       statistic=c("means","ks","wilcox"),
                       data,
                       subset=NULL,
                       ...) {
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
  
  if(is.null(cutpoint)) {
    cutpoint<-0
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
  
  estit<-function(bw,tau=0) {
    Y0<-Y[X < 0 & X >= -bw]
    Y1<-Y[X >= 0 & X <= bw]
    RIestimate(Y0,Y1,method=method,statistic=statistic,tau=tau,...)
  }
  
  if(is.null(bws) & !is.null(covs)) {
    wins<-apply(X=covs,MARGIN=2,FUN=RIbandwidth,w.test=w.test,method=method,statistic=statistic,R=X,...) # Rewrite the next bit just directly using the doKS etc stuff
    return(wins)

    getwins<-function(x,win.type,stat,ri.type){
      wins<-x[[stat]][[ri.type]][[win.type]][1:2] 
      #"window", liberal or conservative bands
      #Currently the bounds are symmetric, but we wont assume that here, and take average
      wins<-{abs(wins[1])+abs(wins[2])}/2
      names(wins)<-NULL
      return(wins)
    }
    
    libwin<-NULL
    conwin<-NULL

    if(doKS) {
      libwinKS<-sapply(wins,getwins,
                       win.type="liberal",stat="ks",ri.type="window")
      conwinKS<-sapply(wins,getwins,
                       win.type="conservative",stat="ks",ri.type="window")
      libwin<-c(libwin,libwinKS)
      conwin<-c(conwin,conwinKS)
      if(doBN) {
        libwinKSBN<-sapply(wins,getwins,
                           win.type="liberal",stat="ks",ri.type="BN")
        conwinKSBN<-sapply(wins,getwins,
                           win.type="conservative",stat="ks",ri.type="BN")
        libwin<-c(libwin,libwinKSBN)
        conwin<-c(conwin,conwinKSBN)
      }
      
      if(doFM) {
        libwinKSFM<-sapply(wins,getwins,
                           win.type="liberal",stat="ks",ri.type="FM")
        conwinKSFM<-sapply(wins,getwins,
                           win.type="conservative",stat="ks",ri.type="FM")
        libwin<-c(libwin,libwinKSFM)
        conwin<-c(conwin,conwinKSFM)
      }
    }
    
    if(doMeans) {
      libwinMn<-sapply(wins,getwins,
                       win.type="liberal",stat="means",ri.type="window")
      conwinMn<-sapply(wins,getwins,
                       win.type="conservative",stat="means",ri.type="window")
      libwin<-c(libwin,libwinMn)
      conwin<-c(conwin,conwinMn)

      if(doBN) {
        libwinMnBN<-sapply(wins,getwins,
                           win.type="liberal",stat="means",ri.type="BN")
        conwinMnBN<-sapply(wins,getwins,
                           win.type="conservative",stat="means",ri.type="BN")
        libwin<-c(libwin,libwinMnBN)
        conwin<-c(conwin,conwinMnBN)
      }

      if(doFM) {
        libwinMnFM<-sapply(wins,getwins,
                           win.type="liberal",stat="means",ri.type="FM")
        conwinMnFM<-sapply(wins,getwins,
                           win.type="conservative",stat="means",ri.type="FM") 
        libwin<-c(libwin,libwinMnFM)
        conwin<-c(conwin,conwinMnFM)
      }
    }
    
    if(doWX) {
      libwinWX<-sapply(wins,getwins,
                       win.type="liberal",stat="wilcox",ri.type="window")
      conwinWX<-sapply(wins,getwins,
                       win.type="conservative",stat="wilcox",ri.type="window")
      libwin<-c(libwin,libwinWX)
      conwin<-c(conwin,conwinWX)


      if(doBN) {
        libwinWXBN<-sapply(wins,getwins,
                           win.type="liberal",stat="wilcox",ri.type="BN")
        conwinWXBN<-sapply(wins,getwins,
                           win.type="conservative",stat="wilcox",ri.type="BN")
      libwin<-c(libwin,libwinWXBN)
      conwin<-c(conwin,conwinWXBN)

      }
      
      if(doBN) {
        libwinWXFM<-sapply(wins,getwins,
                           win.type="liberal",stat="wilcox",ri.type="FM")
        conwinWXFM<-sapply(wins,getwins,
                           win.type="conservative",stat="wilcox",ri.type="FM")
        libwin<-c(libwin,libwinWXFM)
        conwin<-c(conwin,conwinWXFM)
      }

    }
    
    out<-list(liberal=lapply(libwin,estit),conservative=lapply(conwin,estit))
  } else if (!is.null(bws)) {
    out<-lapply(bws,estit)
  } else {
    stop("No bandwidth or means of calculating bandwidth found.")
  }

  class(out)<-"RIRD"
  return(out)
}
