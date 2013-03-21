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

RIbandwidth <- function(Y, X, w.test, cutpoint=NULL, verbose=FALSE, method=c("binomial","fixed.margins"),
                        statistic=c("means","ks","wilcox"), alpha=.05, max.sim=10000, diag.out=FALSE) {
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
  
  doBN <- "binomial"%in%method
  # Do exact fixed margins inference
  doFM <- "fixed.margins"%in%method

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
  
  #Need progress bar
  test.window<-function(w) {
    wr <- cutpoint + w
    wl <- cutpoint - w
    ir <- rr <= wr
    il <- rl <= wl
    nr <- sum(ir)
    nl <- sum(il)
    n <- nr + nl
    if( nl == 0 | nr ==0) return(NULL)
    ri<-RIestimate(y0 = Yl[il], Y1 = Yr[ir], M = M, method = method, statistic = statistic, alpha = alpha)
    
    diag<-c(Obsl=nl, Obsr=nr, Windowl=wl, Windowr=wr, exactFM=ri$exactFM, exactBN=ri$exactBN)
    diag<-c(diag,t.pval=ri$t.test$p.value,meanl=ri$t.test$estimate[1],meanr=ri$t.test$estimate[2])
    diag<-c(diag,meansBN.pval=ri$meansBN.pval,meansFM.pval=ri$meansFM.pval,nullBN=ri$nullBN)
    
    if(!isBinary) {
      diag<-c(diag,wilcox.pval=ri$wilcox.test$p.value,ks.pval=ri$ks.test$p.value)
      diag<-c(diag,wxBN.pval=ri$wxBN.pval,wxFM.pval=ri$wxFM.pval)
      diag<-c(diag,ksBN.pval=ri$ksBN.pval,ksFM.pval=ri$ksFM.pval)
    }
    return(diag)
  }
  #Begin RI window selection
  
  out<-sapply(w.test, test.window)
  
  ret<-list()
  if(diag.out) ret$diag <- out
  
  if(doMeans) { #Retrieve window for diff of means
    ret$means<-list()
    ret$means$window <- getWindow(pval=out["t.pval",], alpha = alpha, rr = out["Windowr",], 
                                  rl = out["Windowl",], obsr=out["Obsr",] , obsl=out["Obsr",])
    
    if(doBN) ret$means$BN <- getWindow(pval=out["meansBN.pval",], alpha = alpha, rr = out["Windowr",], 
                                       rl = out["Windowl",], obsr=out["Obsr",] , obsl=out["Obsr",])
    
    if(doFM) ret$means$FM <- getWindow(pval=out["meansFM.pval",], alpha = alpha, rr = out["Windowr",], 
                                       rl = out["Windowl",], obsr=out["Obsr",] , obsl=out["Obsr",])
  }
  
  if(doKS) { # For KS
    ret$ks <- list()
    ret$ks$window <- getWindow(pval=out["ks.pval",], alpha = alpha, rr = out["Windowr",], 
                       rl = out["Windowl",], obsr=out["Obsr",] , obsl=out["Obsr",])
    
    if(doBN) ret$ks$BN <- getWindow(pval=out["ksBN.pval",], alpha = alpha, rr = out["Windowr",], 
                                       rl = out["Windowl",], obsr=out["Obsr",] , obsl=out["Obsr",])
    
    if(doFM) ret$ks$FM <- getWindow(pval=out["ksFM.pval",], alpha = alpha, rr = out["Windowr",], 
                                       rl = out["Windowl",], obsr=out["Obsr",] , obsl=out["Obsr",])
  }
  if(doWX) { # For Wilcoxon
    ret$wilcox<-list()
    ret$wilcox$window <- getWindow(pval=out["wilcox.pval",], alpha = alpha, rr = out["Windowr",], 
                               rl = out["Windowl",], obsr=out["Obsr",] , obsl=out["Obsr",])
    
    if(doBN) ret$wilcox$BN <- getWindow(pval=out["wxBN.pval",], alpha = alpha, rr = out["Windowr",], 
                                    rl = out["Windowl",], obsr=out["Obsr",] , obsl=out["Obsr",])
    
    if(doFM) ret$wilcox$FM <- getWindow(pval=out["wxFM.pval",], alpha = alpha, rr = out["Windowr",], 
                                    rl = out["Windowl",], obsr=out["Obsr",] , obsl=out["Obsr",])
  }
  
  return(ret)
}