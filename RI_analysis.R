setwd("~/Dropbox/rdd-dev/RI/")

require("foreign") #For Stata data
require("dummies") #Create dummies

d<-read.dta("PoliticalDynastiesData.dta")

d$housesenate <- ifelse(d$careerpath==3,1,0)

#Set up location and year indicators
reg <- data.frame(dummy("region",data=d))
dec <- data.frame(dummy("decade",data=d))
reg<-reg[,-1]
dec<-dec[,-c(1,10,20,21,22)]
#region 1 and 1880 as basis
nreg<-names(reg)
ndec<-names(dec)
d<-cbind(d,reg)
d<-cbind(d,dec)

stt <- data.frame(dummy("state",data=d))
yr <- data.frame(dummy("year",data=d))
stt<-stt[,-10]
yr<-yr[,-53]
# #PA and 1892 as basis
nstt<-names(stt[,-10])
nyr<-names(yr[,-53])
d<-cbind(d,stt)
d<-cbind(d,yr)
rm(list=c("stt","yr"))
rm(list=c("dec","reg"))

covs<-c(nreg,ndec,"democrat","republican","female","collegeatt",
        "outsider","preanyoffice","age","military","farmer",
        "lawyer","business")
sub<-d$year==d$yearenter & d$realbirthyear<=1910 & d$diedinoffice==0 & d$prerelative==0

clus <- d$state

OUT <- cbind(d$democrat,d$female,d$collegeatt,d$outsider,d$preanyoffice,d$age,d$military,d$farmer,d$lawyer,d$business,d$deathage)
# short names
OUTnm <- cbind("Democrat","Female","College","Outsider","Prev Office","Age","Military","Farmer","Lawyer","Business","Death Age")

alpha = 0.05
M = 10000
step = 0.005
start.length=0.005
verbose=TRUE
wmax = NULL
c = 0 
filepath = "./output/"
filename = paste("01-chooseW-balance-fixedW",sep="")
w.length = seq(0.5, 30, by=0.25)

for(v in 1:length(OUTnm)) {
  
  cat("Starting analysis of covariate", OUTnm[v], "\n")
  t0 = proc.time()
  
  dat <- data.frame(OUT[,v], d$postrelative,d$marginvote)
  dim(d)
  dim(dat)
  nona <- complete.cases(dat)
  binary <- (sum(dat[nona,1]== 0 | dat[nona,1]== 1) == length(dat[nona,1]))
  
  # check if missing values make the number of obs go below 10 and 10 in first window
  cat("Checking number of observations in first window \n")
  Y = dat[nona,1]
  Tr=dat[nona,2]  
  r = dat[nona,3]
  ii = (r>= -0.5 & r<= 0.5)
  Y = Y[ii]
  ntr= length(Y[Tr[ii]==1])
  nco= length(Y[Tr[ii]==0])  
  cat("There are", ntr, "obs in Tr group and", nco, "obs in Co group \n\n")
  
  # do balance tests
  #Y= dat[nona,1]; Tr=dat[nona,2]; r=dat[nona,3];cov.name=OUTnm[v]; cov.binary=binary 
  balout = select.balance.window.fixedW(Y= dat[nona,1], Tr=dat[nona,2], r=dat[nona,3], c=c, verbose=verbose, w.length=w.length,
                                        cov.name=OUTnm[v], cov.binary=binary, alpha=alpha, M = M, 
                                        file.name = filename, file.path = filepath, subsample.name = NULL)
  
  cat("Total time for covariate", OUTnm[v],":", (proc.time()[3]- t0[3])/60, " minutes \n\n")
  cat("---------------------------------------------------------------------------------------\n")
}
