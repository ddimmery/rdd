getRIWindow <- function(pval, alpha, rr, rl, obsr, obsl) { #Need to rewrite this to only use a bandwidth, not right and left window points
  if(sum(pval>=alpha)==0) {
    pval.min = 0
    wr.i = length(rr)
    wr = rr[wr.i]
    wl = rl[wr.i]
    wr.con = wr
    wl.con = wl
    Nr = obsr[wr.i]
    Nl = obsl[wr.i]  
    Nr.con = Nr
    Nl.con = Nl  
  } else {
    pval.min = min(pval[pval >= alpha]) # this will be window limit assuming monotonicity, otherwise will be too liberal
    wr.i    = min(which(pval >= alpha)) # first time that p-value is alpha or more
    wr = rr[wr.i]
    wl = rl[wr.i]
    Nr = obsr[wr.i]
    Nl = obsl[wr.i]
    ## conservative window limit: wr such that p-value is alpha or more and never falls below alpha again
    cum = cumsum(pval >= alpha)
    test = cum[2:length(cum)] - cum[1:(length(cum)-1)]
    testi = test == 0 
    if(!any(testi)) 
      wr.i.con <- length(testi)+2
    else
      wr.i.con = max(which(testi))+2
    
    if(wr.i.con > length(rr)) wr.i.con = length(rr)
    wr.con = rr[wr.i.con]
    wl.con = rl[wr.i.con]  
    Nr.con = obsr[wr.i.con]
    Nl.con = obsl[wr.i.con]  
  }
  return(list(liberal=c(wr=wr, wl=wl, Nr=Nr, Nl=Nl),conservative=c(wr=wr.con, wl= wl.con, Nr=Nr.con, Nl=Nl.con)))
}
