# getCI will take the output of RIestimate and return a conf interval
#  which assumes constant treatment effects
#  ** This should be folded out into an external function **
#  ** Call it when summary.RIRD is called **
  getCI <- function(Y,X,bw,alpha=.05,method,statistic="means",prec=.0001,...) {
    #So the data resides in x$t.test$p.value
    #Want to make this equal to alpha/2 & 1-alpha/2 by
    #changing tau
    #Start out only doing it for diff in means
    Y0<-Y[X < 0 & X >= -bw]
    Y1<-Y[X >= 0 & X <= bw]
#Take method and make ri.type one of FM or BN
    if(method=="binomial") 
      ri.type<-"BN"
    else if(method=="fixed.margins")
      ri.type<-"FM"
    else
      stop("RI type must be either 'binomial' or 'fixed.margins'.")

    ri.init<-RIestimate(Y0,Y1,
                     method=method,statistic=statistic,tau=0,...)
    linit<-ri.init$t.test$conf.int[1]
    uinit<-ri.init$t.test$conf.int[2]

    esttau<-ri.init$diff.mean

    sddelt <- sum(abs(esttau-c(linit,uinit)))/2

    tau<-c(sign(linit)*sddelt/2+linit,linit)
    print(tau)
    print(sddelt)
    print(linit)
    print(uinit)
    pstr<-paste0(statistic,ri.type,".pval")
    pcalc<-c(NA,NA) 
    pcalc[1]<-RIestimate(Y0,Y1,
                         method=method,statistic=statistic,tau=tau[1],...)[pstr][[1]]
    pcalc[2]<-RIestimate(Y0,Y1,
                         method=method,statistic=statistic,tau=tau[2],...)[pstr][[1]]
    p.target<-alpha/2 #Start looking for lower bound

    while(TRUE) {
      
      pdif<-pcalc[2]-pcalc[1] # How did p change
      tdif<-tau[2]-tau[1]
      slope<-pdif/tdif
      ord<-order(abs(pcalc-p.target)) #smallest first
      tau<-tau[ord]
      pcalc<-pcalc[ord]
      
      print(tau)
      print(pcalc)
      #In case a p value is stuck in no man's land
      if(any(abs(pcalc)<1e-8)) {
        ip <- which(abs(pcalc)<1e-8)[1]
        tau[ip]<-rnorm(1,linit,sddelt/2)
        pcalc[ip] <- RIestimate(Y0,Y1,
                               method=method,
                               statistic=statistic,
                               tau=tau[ip],...)[pstr][[1]]
        next
      }

      if(abs(pcalc[2]-p.target)<prec) break
      if(abs(tdif)<prec/10) break
      if(abs(pdif)<prec/10) {
        tau[2] <- rnorm(1,tau[2],sddelt/2)
        pcalc[2] <- RIestimate(Y0,Y1,
                               method=method,
                               statistic=statistic,
                               tau=tau[2],...)[pstr][[1]]
        next
      }
      tau[2]<-(p.target-pcalc[2]+slope*tau[2])/slope
      pcalc[2]<-RIestimate(Y0,Y1,
                           method=method,
                           statistic=statistic,
                           tau=tau[2],...)[pstr][[1]]
    } 
    ltau <- tau[2]
    lp <- pcalc[2]

    tau<-c(sign(uinit)*sddelt/2+uinit,uinit)

    pcalc[1]<-RIestimate(Y0,Y1,
                         method=method,statistic=statistic,tau=tau[1],...)[pstr][[1]]
    pcalc[2]<-RIestimate(Y0,Y1,
                         method=method,statistic=statistic,tau=tau[2],...)[pstr][[1]]

    p.target<-1-p.target #And now upper bound
    
    
    while(TRUE) {
      
      pdif<-pcalc[2]-pcalc[1] # How did p change
      tdif<-tau[2]-tau[1]
      slope<-pdif/tdif
      ord<-order(abs(pcalc-p.target)) #smallest first
      tau<-tau[ord]
      pcalc<-pcalc[ord]
      
      print(tau)
      print(pcalc)
      #In case a p value is stuck in no man's land
      if(any(abs(pcalc)<1e-8)) {
        ip<-which(abs(pcalc)<1e-8)[1]
        tau[ip]<-rnorm(1,uinit,sddelt/2)
        pcalc[ip] <- RIestimate(Y0,Y1,
                               method=method,
                               statistic=statistic,
                               tau=tau[ip],...)[pstr][[1]]
        next
      }

      if(abs(pcalc[2]-p.target)<prec) break
      if(abs(tdif)<prec/10) break
      if(abs(pdif)<prec/10) {
        tau[2] <- rnorm(1,tau[2],sddelt)
        pcalc[2] <- RIestimate(Y0,Y1,
                               method=method,
                               statistic=statistic,
                               tau=tau[2],...)[pstr][[1]]
        next
      }
      tau[2]<-(p.target-pcalc[2]+slope*tau[2])/slope
      pcalc[2]<-RIestimate(Y0,Y1,
                           method=method,
                           statistic=statistic,
                           tau=tau[2],...)[pstr][[1]]
    } 

    utau <- tau[2]
    up <- pcalc[2]

    alpha <- 1 - up + lp
    ci <- c(lwr=ltau,upr=utau)
    return(list(alpha = alpha, ci = ci))
  }

