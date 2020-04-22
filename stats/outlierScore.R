library(extraDistr)
library(stringr)

args = commandArgs(trailingOnly=TRUE)
fitFile=args[1]
countsFile=args[2]
outname=args[3]

psiStats=read.table(fitFile,sep="\t",header = TRUE,stringsAsFactors = FALSE)

data=read.csv(countsFile)
iso1_idx=which(endsWith(colnames(data),"iso1"))
iso2_idx=which(endsWith(colnames(data),"iso2"))
sample_names=str_replace(colnames(data)[iso1_idx],"_iso1","")

sampleCount=length(sample_names)
eventCount=nrow(data)

print(paste('sampleCount:',sampleCount))
print(paste('eventCount:',eventCount))

iso1=sapply(data[,iso1_idx],as.numeric)
iso2=sapply(data[,iso2_idx],as.numeric)
info_idx=seq(1,iso1_idx[1]-1)

idx=match(data$event_jid,psiStats$event_jid)
print(idx[which(!is.na(idx))])
alpha=psiStats[idx,"alpha"]
beta=psiStats[idx,"beta"]
scores=matrix(NA,nrow = length(alpha),ncol = sampleCount)
colnames(scores)<-sample_names
idx1=which(alpha>=beta)
idx2=which(alpha<beta)
if(length(idx1)>0){
  scores[idx1,]=t(sapply(idx1, function(x) log(pmax(pbbinom(iso1[x,],iso1[x,]+iso2[x,],alpha[x],beta[x],log.p=FALSE),.Machine$double.xmin))))
}
if(length(idx2)>0){
  scores[idx2,]=t(sapply(idx2,function(x) -log(pmax(pbbinom(iso2[x,],iso1[x,]+iso2[x,],beta[x],alpha[x],log.p=FALSE),.Machine$double.xmin))))
}

write.table(cbind(data[,info_idx],scores),outname,row.names=FALSE,append = FALSE,sep="\t")
