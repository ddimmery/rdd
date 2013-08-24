#' Perform randomization inference on the given data
#' 
#' \code{RIestimate} performs randomization inference on the given data. It can assume
#' either a binomial distribution for treatment or fixed margins.
#' 
#' @param Y0 a numerical vector which contains the control sample
#' @param Y1 a numerical vector which contains the treated sample
#' @param M the maximum number of simulations to perform (will do exact inference if the necessary permutations are less than this value)
#' @param tau the constant null effect to test against. default is the sharp null.
#' @param verbose logical flag indicating whether to print more information to the terminal. Default is \code{FALSE}.
#' @param alpha numeric, alpha-level of significance to use for rejection region
#' @param method randomization method, one or more of binomial of fixed.margins
#' @param statistic statistics to test; can be one or more of "means" (for difference in means), "ks" (for Kolmogorov-Smirnov) or "wilcox" (Wilcoxon rank sum)
#' @return The estimates
#' @references Titiunik etal
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>

RIestimate <- function(Y0,Y1,M,tau=0,verbose=FALSE,alpha=0.05,method=c("binomial","fixed.margins"),statistic=c("means","ks","wilcox")) {
  n1 <- length(Y1)
  n0 <- length(Y0)
  n <- n1 + n0
  pr <- n1/n
  Y1 <- Y1-tau
  Y <- c(Y1, Y0)
  Tr <- c(rep(1,n1), rep(0,n0))
  
  # Do exact binomial inference
  doBN <- "binomial" %in% method
  exactBN <- 2^n<M & doBN
  # Do exact fixed margins inference
  doFM <- "fixed.margins"%in%method
  exactFM <- choose(n,n1)<M & doFM
  # is the outcome dichotomous
  isBinary <- sum(Y==1 | Y==0) == n
  
  # Do difference in means tests
  doMeans <- "means" %in% statistic
  # Do Kolmogorov-Smirnov tests
  doKS <- "ks" %in% statistic
  # Do Wilcoxon tests
  doWX <- "wilcox" %in% statistic
  
  if( !(doMeans|doKS|doWX)) stop("Must select at least one statistic to calculate")
  if( !doMeans & isBinary) stop("Can only calculate difference in means statistics with binary variable.")
  cat("\n")
  if(verbose) {
    cat("Number treated: ",n1,"\n")
    cat("Number untreated: ",n0,"\n")
    cat("Binomial randomization inference: ", ifelse(doBN,"YES","NO"),"\n")
    cat("Exact binomial RI: ", ifelse(exactBN,"YES","NO"),"\n")
    cat("Fixed Margins randomization inference: ", ifelse(doFM,"YES","NO"),"\n")
    cat("Exact FM RI: ", ifelse(exactFM,"YES","NO"),"\n")
    cat("Difference of Means: ", ifelse(doMeans,"YES","NO"),"\n")
    cat("KS Test: ", ifelse(doKS,"YES","NO"),"\n")
    cat("Wilcoxon Test: ", ifelse(doWX,"YES","NO"),"\n")
    cat("Binary variable: ",ifelse(isBinary,"YES","NO"),"\n")
  }
  
  #Initialize output list
  out <- list()
  out$prob <- pr
  
  if(doMeans) {
    mean0<-mean(Y1) - mean(Y0)
    out$diff.mean = mean0
    out$t.test <- try(t.test(Y1, Y0, alternative = "two.sided", mu = 0, paired = FALSE, var.equal = FALSE, conf.level = (1-alpha)))
    if(is(out$t.test, "try-error")) out$t.test = list(estimate=c(mean(Y1), mean(Y0)), p.value=NA)
  }
  
  if(!isBinary) {
    if(doKS) {
      out$ks.test <- ks.test(Y1,Y0,alternative="two.sided")
      ks0 <- out$ks.test$statistic
      out$ks.stat <- ks0
    }
    if(doWX) {
      out$wilcox.test<-wilcox.test(Y1, Y0, alternative = c("two.sided"), mu = 0,paired = FALSE, exact = NULL, correct = TRUE, conf.int = FALSE)
      wx0 <- sum(rank(Y,ties.method="average")[Tr=1])
      out$wilcox.stat<-wx0
    }
  }
  
  if(exactFM) {
    #Perform exact inference when feasible
    TrFM <- combn(1:n,n1)
    ncomb <- ncol(TrFM)
    stats <- apply(TrFM,2,RIstatsFM,Y = Y, doKS = doKS, doWX = doWX, doMeans = doMeans, binary = isBinary)
    
    # Get data
    if(doMeans) {
      mns<-stats[1,]
      out$meansFM.pval<-sum(abs(mns)>=abs(mean0))/ncomb
    }
    if(!isBinary){
      if(doKS) {
        kss<-stats[doMeans+1,]
        out$ksFM.pval<-sum(kss >= ks0)/ncomb
      }
      if(doWX) {
        wxs<-stats[doMeans+doKS+1,]
        out$wxFM.pval<-sum(abs(wxs - mean(wxs)) >= abs(wx0 - mean(wxs)))/ncomb
      }
    }
  } else if (doFM) {
    # Otherwise, perform M simulations
    simFM<-function() {
      TrFM<-sample(Tr,replace=FALSE)
      return(RIstatsFM(TrFM,Y = Y, doKS = doKS, doWX = doWX, doMeans = doMeans, binary = isBinary))
    }
    stats<-replicate(M,simFM())
    
    if(doMeans) {
      mns<-stats[1,]
      out$meansFM.pval<-sum(abs(mns)>=abs(mean0))/M
    }
    if(!isBinary){
      if(doKS) {
        kss<-stats[doMeans+1,]
        out$ksFM.pval<-sum(kss >= ks0)/M
      }
      if(doWX) {
        wxs<-stats[doMeans+doKS+1,]
        out$wxFM.pval<-sum(abs(wxs - mean(wxs)) >= abs(wx0 - mean(wxs)))/M
      }
    }
  }
  
  if(exactBN) {
    TrBN<-matrix(NA,nrow=n,ncol=(2^n)-2)
    indx<-1:n
    p=NULL
    col = 1
    for (q in 1:(n-1)) {
      CC = combn(1:n,q)
      ncomb = ncol(CC)
      for(u in 1:ncomb) {
        TrBN[CC[,u],col] = 1
        TrBN[indx[-CC[,u]],col] = 0
        col = col + 1
      }
      prob = pr^q * ((1-pr)^(n-q))
      p = c(p,rep(prob,ncomb))
      p = p/sum(p) # normalize so they sum to one -- remember I am excluding two possible outcomes (all Tr, all Co)
    }
    stats<- apply(TrBN,2,RIstatsBN, Y=Y, doKS = doKS, doWX = doWX, doMeans = doMeans, binary = isBinary)
    
    if(doMeans) {
      mns<-stats[1,]
      out$meansBN.pval<-sum(p[abs(mns)>=abs(mean0)])
    }
    if(!isBinary){
      if(doKS) {
        kss<-stats[doMeans+1,]
        out$ksBN.pval<-sum(p[kss >= ks0])
      }
      if(doWX) {
        wxs<-stats[doMeans+doKS+1,]
        out$wxBN <- mean(wxs)
        out$wxBN.pval<-sum(p[abs(wxs - out$wxBN) >= abs(wx0 - out$wxBN)])
      }
    }
  } else if(doBN) {
    indxr<-matrix(runif(n*M),nrow=n)
    TrBN<-indxr<=pr
    inullBN<-apply(TrBN,2,function(x) if(sum(x)==n | sum(x)==0) return(1) else return(0))
    out$nullBN<-sum(inullBN)
    if(out$nullBN>0) TrBN<-TrBN[,!inullBN]
    stats<-apply(TrBN,2,RIstatsBN,Y=Y, doKS = doKS, doWX = doWX, doMeans = doMeans, binary = isBinary)
    if(doMeans) {
      mns<-stats[1,]
      out$meansBN.pval<-sum(abs(mns)>=abs(mean0))/M
    }
    if(!isBinary){
      if(doKS) {
        kss<-stats[doMeans+1,]
        out$ksBN.pval<- sum(kss>=ks0)/M
      }
      if(doWX) {
        wxs<-stats[doMeans+doKS+1,]
        out$wxBN <- mean(wxs)
        out$wxBN.pval<-sum(abs(wxs - out$wxBN)>=abs(wx0 - out$wxBN))/M
      }
    }
    
  }
  
  if(doWX) out$wx <- (1/2 * n1 * (n+1))
  if(doFM) out$exactFM <- exactFM
  if(doBN) out$exactBN <- exactBN
  return(out)
}