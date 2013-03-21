getWindow <- function(pval, alpha, rr, rl, obsr, obsl) {
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
    ## conservative window limit: wr such that p-value is 0.2 or more and never falls below 0.2 again
    cum = cumsum(pval >= alpha)
    test = cum[2:length(cum)] - cum[1:(length(cum)-1)]
    wr.i.con = max(which(test == 0))+2
    if(wr.i.con > length(rr)) wr.i.con = length(rr)
    wr.con = rr[wr.i.con]
    wl.con = rl[wr.i.con]  
    Nr.con = obsr[wr.i.con]
    Nl.con = obsl[wr.i.con]  
  }
  return(c(wr=wr, wl=wl, wr.con=wr.con, wl.con= wl.con, Nr=Nr, Nl=Nl, Nr.con=Nr.con, Nl.con=Nl.con))
}