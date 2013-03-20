RIstatsFM<-function(Tr,Y, doKS = TRUE, doWX = TRUE, doMeans = TRUE, binary=NULL) {
  if(is.null(binary)) binary <- sum(Y==1 | Y==0) == length(Y)
  o<-c()
  if(doMeans) o<-c(o,mn=mean(Y[Tr]) - mean(Y[-Tr]))
  if(!binary) {
    if(doKS) o<-c(o,ks=ks.test(Y[Tr],Y[-Tr],alternative=="two.sided")$statistic)
    if(doWX) o<-c(o,wx=sum(rank(Y,ties.method="average")[Tr]))
  }
  return(o)
}