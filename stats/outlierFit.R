library(stats4)
library(fitdistrplus)
library(extraDistr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
countsFile=args[1]
maxBeta=as.integer(args[2])
outname=str_c(args[3],'bisbeeFit.csv')
if (length(args)>3){
  sampleFile=args[4]
}
data=read.csv(countsFile)
iso1_idx=which(endsWith(colnames(data),"iso1"))
iso2_idx=which(endsWith(colnames(data),"iso2"))
sample_names=str_replace(colnames(data)[iso1_idx],"_iso1","")

if (exists("sampleFile")){
  incSamples=read.table(sampleFile)
  sampleIdx=match(incSamples[,1],sample_names)
}else{
  sampleIdx=1:length(sample_names)
}

sampleCount=length(sampleIdx)
eventCount=nrow(data)
print(paste('sampleCount:',sampleCount))
print(paste('eventCount:',eventCount))

iso1=sapply(data[,iso1_idx],as.numeric)
iso2=sapply(data[,iso2_idx],as.numeric)
#print(head(iso1))
info_idx=seq(1,which(colnames(data)=="event_jid"))
print(info_idx)
cnames<-c(colnames(data)[info_idx],"alpha","beta","depth_min","depth_q05","depth_q25","depth_median","depth_q75","depth_q95","depth_max","psi_min","psi_q05","psi_q25","psi_median","psi_q75","psi_q95","psi_max","log-likelihood")
depthQ=t(apply(iso1+iso2,1,quantile,c(0,0.05,0.25,0.5,0.75,0.95,1)))
psi=iso1/(iso1+iso2)
meanPSI=rowMeans(psi,na.rm = TRUE)
sumDepth=rowSums(iso1+iso2)
hIdx=which(meanPSI>0.5)
mCounts=iso1
mCounts[hIdx,]=iso2[hIdx,]
mle_bb_mode0<-function(x){
  nLL<-function(v1,v2) -sum(dbbinom(mCounts[x,],iso1[x,]+iso2[x,],(1-1/maxBeta)/(1+exp(v1))+1/maxBeta,(maxBeta-1)/(1+exp(v2))+1,log=TRUE))
  mleRes<-try(mle(nLL,start=list(v1=0,v2=0),method = "Nelder-Mead"))
  if(class(mleRes)=="try-error"){
    mleRes<-try(mle(nLL,start=list(v1=0,v2=0),method = "BFGS"))
    if(class(mleRes)=="try-error"){
      mleRes<-NULL
    }
  }
  return(mleRes)  
}
a<-rep(1,length(meanPSI))
b<-sapply(sumDepth+1,function(x) min(x,maxBeta))
ll<-rep(NA,length(meanPSI))
fIdx<-which(rowSums(mCounts)>0)
if(length(fIdx)>0){
  param<-sapply(fIdx,mle_bb_mode0)
  idx=which(!sapply(param,is.null))
  if(length(idx)>0){
    a[fIdx[idx]]<-sapply(idx,function(x) (1-1/maxBeta)/(1+exp(param[[x]]@coef['v1']))+1/maxBeta)
    b[fIdx[idx]]<-sapply(idx,function(x) (maxBeta-1)/(1+exp(param[[x]]@coef['v2']))+1)
    ll[fIdx[idx]]<-sapply(idx,function(x) param[[x]]@min)
  }
}
alpha=a
alpha[hIdx]=b[hIdx]
beta=b
beta[hIdx]=a[hIdx]
if(length(fIdx)>0){
  nIdx=which(sapply(param,is.null))
  alpha[fIdx[nIdx]]=NA
  beta[fIdx[nIdx]]=NA 
}

psiQ=t(apply(iso1/(iso1+iso2),1,quantile,c(0,0.05,0.25,0.5,0.75,0.95,1),na.rm=TRUE))
write.table(cbind(data[,info_idx],alpha,beta,depthQ,psiQ,ll),outname,col.names=cnames,row.names=FALSE,append = FALSE,sep=",")

