library(stats4)
library(fitdistrplus)
library(extraDistr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
countsFile=args[1]
sample_file=args[2]
outname=str_c(args[3],".bisbeeDiff.csv")
maxW=as.integer(args[4])
if (length(args)>4){
  groups=c(args[5],args[6])
}

data=read.csv(countsFile)
iso1_idx=which(endsWith(colnames(data),"iso1"))
iso2_idx=which(endsWith(colnames(data),"iso2"))
sample_names=str_replace(colnames(data)[iso1_idx],"_iso1","")

sample_table=read.table(sample_file)
if(exists('groups')){
  sample_table$V2=factor(sample_table$V2,levels=groups)
}else{
  groups=levels(sample_table$V2) 
}

g1idx=sample_names %in% sample_table$V1[unclass(sample_table$V2)==1]
g2idx=sample_names %in% sample_table$V1[unclass(sample_table$V2)==2]
sidx=g1idx | g2idx
print(sample_names[sidx])

iso1=sapply(data[,iso1_idx],as.numeric)
iso2=sapply(data[,iso2_idx],as.numeric)

info_idx=seq(1,iso1_idx[1]-1)
cnames<-c(colnames(data)[info_idx],'ll_ratio',paste('psi',groups[1],sample_names[g1idx],sep='_'),paste('psi',groups[2],sample_names[g2idx],sep='_'),paste('depth',groups[1],sample_names[g1idx],sep='_'),paste('depth',groups[2],sample_names[g2idx],sep='_'),paste('fitPSI',groups[1],sep="_"),paste('fitPSI',groups[2],sep="_"),'fitPSI_all','fitW_2group','fitW_all')
initParam<-function(iso1,iso2){
  psi=iso1/(iso1+iso2)
  idx=which(apply(iso1,1,max)>0 & apply(iso2,1,max)>0 & rowSums(is.finite(psi))>1)
  betaParam<-sapply(idx,function(x) fitdist(psi[x,is.finite(psi[x,])],'beta',method="mme"))
  initPSI=rowMeans(psi,na.rm=T)
  if(length(idx)>0){
    initPSI[idx]<-sapply(1:ncol(betaParam),function(x) betaParam['estimate',x][[1]][['shape1']]/(betaParam['estimate',x][[1]][['shape1']]+betaParam['estimate',x][[1]][['shape2']]))
  }
  initPSI[is.na(initPSI)]=0.5
  initPSI[initPSI<1E-5]=1E-5
  initPSI[initPSI>1-1E-5]=1-1E-5
  initW=rowMeans(log(iso1+iso2+1))
  if(length(idx)>0){
    initW[idx]=sapply(1:ncol(betaParam),function(x) betaParam['estimate',x][[1]][['shape1']]+betaParam['estimate',x][[1]][['shape2']])
  }
  initW[initW<2.5 | is.na(initW)]=2.5
  return(data.frame(initPSI=initPSI,initW=initW))
}
initParam1<-initParam(iso1[,g1idx],iso2[,g1idx])
initParam2<-initParam(iso1[,g2idx],iso2[,g2idx])
initPSI1<-initParam1$initPSI
initPSI2<-initParam2$initPSI
initW<-rowMeans(cbind(initParam1$initW,initParam2$initW),na.rm=T)
initParamA<-initParam(iso1[,sidx],iso2[,sidx])
initPSIa<-initParamA$initPSI
initWa<-initParamA$initW
#maxW<-apply(iso1[,(g1idx | g2idx)]+iso2[,(g1idx | g2idx)],1,max)+2
initW[initW>=maxW-1]=maxW-1
mle_bb_2group<-function(x){
  nLL<-function(v1,v2,v3) -sum(dbbinom(iso1[x,g1idx],iso1[x,g1idx]+iso2[x,g1idx],(maxW/(1+exp(v3))+2)*exp(v1)/(1+exp(v1)),(maxW/(1+exp(v3))+2)*(1-exp(v1)/(1+exp(v1))),log=TRUE))-sum(dbbinom(iso1[x,g2idx],iso1[x,g2idx]+iso2[x,g2idx],(maxW/(1+exp(v3))+2)*exp(v2)/(1+exp(v2)),(maxW/(1+exp(v3))+2)*(1-exp(v2)/(1+exp(v2))),log=TRUE))
  mleRes<-try(mle(nLL,start=list(v1=log(initPSI1[x]/(1-initPSI1[x])),v2=log(initPSI2[x]/(1-initPSI2[x])),v3=log((maxW-initW[x])/initW[x])-2),method = "BFGS"))
  if(class(mleRes)=="try-error"){
    mleRes<-try(mle(nLL,start=list(v1=log(initPSI1[x]/(1-initPSI1[x])),v2=log(initPSI2[x]/(1-initPSI2[x])),v3=log((maxW-initW[x])/initW[x])-2),method = "Nelder-Mead"))
    if(class(mleRes)=="try-error"){
      mleRes<-NULL
    }
  }
  return(mleRes)
}
param<-sapply(1:nrow(iso1),mle_bb_2group)
idx=which(!sapply(param,is.null))
fitPSI1<-rep(NA,length(param))
fitPSI2<-rep(NA,length(param))
fitW<-rep(NA,length(param))
ll_2group<-rep(NA,length(param))
if(length(idx)>0){
  fitPSI1[idx]<-exp(sapply(idx,function(x) param[[x]]@coef['v1']))/(1+exp(sapply(idx,function(x) param[[x]]@coef['v1'])))
  fitPSI2[idx]<-exp(sapply(idx,function(x) param[[x]]@coef['v2']))/(1+exp(sapply(idx,function(x) param[[x]]@coef['v2'])))
  fitW[idx]<-sapply(idx,function(x) maxW/(1+exp(param[[x]]@coef['v3']))+2)
  ll_2group[idx]<-sapply(idx,function(x) param[[x]]@min)
}
initWa[initWa>=maxW-1]=maxW-1
mle_bb_1group<-function(x){
  nLL<-function(v,v3) -sum(dbbinom(iso1[x,sidx],iso1[x,sidx]+iso2[x,sidx],(maxW/(1+exp(v3))+2)*exp(v)/(1+exp(v)),(maxW/(1+exp(v3))+2)*(1-exp(v)/(1+exp(v))),log=TRUE))
  mleRes<-try(mle(nLL,start=list(v=log(initPSIa[x]/(1-initPSIa[x])),v3=log((maxW-initWa[x])/initWa[x])-2),method = "BFGS"))
  if(class(mleRes)=="try-error"){
    mleRes<-try(mle(nLL,start=list(v=log(initPSIa[x]/(1-initPSIa[x])),v3=log((maxW-initWa[x])/initWa[x])-2),method = "Nelder-Mead"))
    if(class(mleRes)=="try-error"){
      mleRes<-NULL
    }
  }
  return(mleRes)  
}
#  mle_bb_1group<-function(x){
#    nLL<-function(psi,lnW) -sum(dbbinom(iso1[x,sidx],iso1[x,sidx]+iso2[x,sidx],exp(lnW)*psi,exp(lnW)*(1-psi),log=TRUE))
#    return(tryCatch(mle(nLL,start=list(psi=initPSIa[x],lnW=depthLNmean[x]),upper=c(1-1E-5,log(depthWmean[x]+1)),lower=c(1E-5,0)),error=function(e) NULL))
#  }
param<-sapply(1:nrow(iso1),mle_bb_1group)
idx=which(!sapply(param,is.null))
fitPSIa<-rep(NA,length(param))
fitWa<-rep(NA,length(param))
ll_1group<-rep(NA,length(param))
if(length(idx)>0){
  fitPSIa[idx]<-exp(sapply(idx,function(x) param[[x]]@coef['v']))/(1+exp(sapply(idx,function(x) param[[x]]@coef['v'])))
  fitWa[idx]<-sapply(idx,function(x) maxW/(1+exp(param[[x]]@coef['v3']))+2)
  ll_1group[idx]<-sapply(idx,function(x) param[[x]]@min)
}
psi=iso1/(iso1+iso2)
res<-cbind(data[,info_idx],ll_1group-ll_2group,psi[,g1idx],psi[,g2idx],iso1[,g1idx]+iso2[,g1idx],iso1[,g2idx]+iso2[,g2idx],fitPSI1,fitPSI2,fitPSIa,fitW,fitWa)
write.table(res,outname,col.names=cnames,row.names=FALSE,append = FALSE,sep=",")


