rm(list=ls(all=TRUE))
set.seed(123)
library('mvtnorm')

#get data
dat=read.csv('original data/fake data.csv',as.is=T)  #for testing
#dat=read.csv('original data/Mal_dat.csv',as.is=T)

ind=which(colnames(dat)%in%c('ComID',"ChdID",'Mal',"Age",'Latitude','Longitude','survey','village'))
covs=colnames(dat[-ind])

cor1=cor(dat[,covs]) #strong correlation NDVI & LST_D
cond=cor1 > -0.6 & cor1 < 0.6
cor1[cond]=NA
diag(cor1)=NA
cor1

#remove strongly correlated covariates: Slope, Prec, Rain_all, LST_D
covs=c('Dist_HF','Dist_UC','Elevation', 'Aspect',"pop",'NDVI','Dist_water','Dist_Roads','LST_N','NTL')

#tranform NTL to log form
hist(dat$NTL)
dat$NTL=log(dat$NTL+1)

#standardize covariates for each survey
nsurv=max(dat$survey)
for (i in 1:nsurv){
  cond=dat$survey==i
  dat1=dat[cond,]
  indiv=dat1[,c('village','Latitude','Longitude','Age','Mal')]
  
  #standardize village level covariates
  covmat=unique(dat1[,c('village',covs)])
  media1=colMeans(covmat[,covs])
  sd1=apply(covmat[,covs],2,sd)
  for (j in 1:length(covs)){
    tipo=covs[j]
    covmat[,tipo]=(covmat[,tipo]-media1[tipo])/sd1[tipo]
  }
  print(c(i,nrow(covmat)))
  
  fim=merge(indiv,covmat,all=T); 
  print(c(i,nrow(indiv),nrow(fim)))
  
  #standardize age
  age.m=mean(fim$Age)
  age.sd=sd(fim$Age)
  fim$Age=(fim$Age-age.m)/age.sd
  
  media1=c(media1,age.m); names(media1)[length(media1)]='Age'
  sd1=c(sd1,age.sd);      names(sd1)[length(sd1)]='Age'
  
  nomes=paste('Original model/edited data/',c('stats','AP'),i,'.csv',sep='')
  stats=rbind(media1,sd1)
  write.csv(stats,nomes[1])
  write.csv(fim,nomes[2],row.names=F)
}

