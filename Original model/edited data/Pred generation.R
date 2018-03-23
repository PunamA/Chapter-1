rm(list=ls(all=TRUE))
set.seed(123)
library('mvtnorm')

#get data
dat=read.csv('original data/unsampled_Mal.csv',as.is=T)

covs=c('Dist_HF','Dist_UC','Elevation','NDVI','Dist_water','Dist_Roads','LST_N','NTL')

hist(dat$NTL)
dat$NTL=log(dat$NTL+1)

#standardize covariates for each survey
nsurv=max(dat$survey)
for (i in 1:nsurv){
  cond=dat$survey==i
  dat1=dat[cond,]

  #standardize village level covariates
  nome=paste('Original model/edited data/stats',i,'.csv',sep='')
  stats=read.csv(nome,as.is=T)
  rownames(stats)=stats$X
  

  for (j in 1:length(covs)){
    tipo=covs[j]
    dat1[,tipo]=(dat1[,tipo]-stats['media1',tipo])/stats['sd1',tipo]
  }
  
  nomes=paste('Original model/edited data/pred',i,'.csv',sep='')
  dat2=dat1[,c('Latitude','Longitude',covs)]
  print(i)
  print(apply(dat2,2,mean))
  print(apply(dat2,2,sd))
  write.csv(dat2,nomes,row.names=F)
}

