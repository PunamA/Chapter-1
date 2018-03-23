rm(list=ls(all=TRUE))
set.seed(18)
library('mvtnorm')

#for loop for each individual survey
for (jjj in 1:6){
nome=paste('Original model/edited data/AP',jjj,'.csv',sep='')
dat=read.csv(nome,as.is=T)

#bring in gibbs functions
source('Original model/Scripts/gibbs functions district.R')

nomes=c('Latitude','Longitude','village','Dist_HF','Dist_UC','Elevation',
        'NDVI','Dist_water','Dist_Roads','LST_N','NTL')
clust=unique(round(dat[,sort(nomes)],4))
clust=clust[order(clust$village),]
nclust=nrow(clust); nclust

#get dmat
dmat=create.dmat(dat)
dtd=t(dmat)%*%dmat

#get wmat and gamma
ind=which(colnames(clust)%in%c('Latitude','Longitude','village'))
wmat=data.matrix(cbind(1,clust[,-ind]))
gamma=rep(0,ncol(wmat))

#get xmat and betas
nomes.x=c('Age')
xmat=data.matrix(dat[,nomes.x])
betas=rep(0,ncol(xmat))
xtx=t(xmat)%*%xmat

#distance matrix
dist1=as.matrix(dist(clust[,c('Latitude','Longitude')]))
nvil=nrow(dist1)

#discretize rho to improve mixing
cor1=seq(from=0.05,to=0.95,by=0.01)
rho.vals=-quantile(dist1,0.1)/log(cor1)
ldet.vals=rep(NA,length(cor1))
Sigma.vals=invSigma.vals=list()
for (i in 1:length(rho.vals)){
  Sigma.vals[[i]]=exp(-dist1/rho.vals[i])
  invSigma.vals[[i]]=solve(Sigma.vals[[i]])
  ldet.vals[i]=determinant(Sigma.vals[[i]],logarithm = T)$modulus[[1]]
}

#get initial values
ind.rho=3
rho=rho.vals[ind.rho]
alpha=rmvnorm(1,mean=rep(0,ncol(dist1)),sigma=Sigma.vals[[ind.rho]])
sig2=2

#get z
z=rep(NA,nrow(dat))
cond=!is.na(dat$Mal) & dat$Mal;    z[cond]=1
cond=!is.na(dat$Mal) & dat$Mal==0; z[cond]=-1
cond=is.na(dat$Mal);                    z[cond]=-0.5;

#useful things for gibbs sampler
ngibbs=10000
vec.alpha=matrix(NA,ngibbs,length(alpha))
vec.outros=matrix(NA,ngibbs,2) #rho+sig2
vec.gamma=matrix(NA,ngibbs,length(gamma))
vec.betas=matrix(NA,ngibbs,length(betas))
param=list(alpha=t(alpha),rho=rho,gamma=gamma,betas=betas,
           invSigma=invSigma.vals[[ind.rho]],z=z,sig2=sig2)

options(warn=2)

#Gibbs Sampler
for (i in 1:ngibbs){
  print(i)
  param$alpha=update.alpha(param)#t(alpha.true)#

  tmp=update.rho(param)
  param$rho=tmp$rho
  param$invSigma=tmp$invSigma
  
  param$sig2=update.sig2(param)
  param$gamma=update.gamma(param)#gamma.true#
  param$betas=update.betas(param)#betas.true#
  
  #sample states
  param$z=update.z(param)#z.true#

  #store results
  vec.alpha[i,]=param$alpha
  vec.outros[i,]=c(param$rho,param$sig2) 
  vec.gamma[i,]=param$gamma
  vec.betas[i,]=param$betas
}
#change colnames to recognise parameters
colnames(vec.outros)=c('rho','sig2')
colnames(vec.gamma)=colnames(wmat)
colnames(vec.betas)=colnames(xmat)

#thinning to half
seq1=floor(seq(from=(ngibbs/2),to=ngibbs,length.out=5000))
nomes=paste('Original model/model results/',c('outros','gamma','betas','alpha'),jjj,'.csv',sep='')
write.csv(vec.outros[seq1,],nomes[1],row.names=F)
write.csv(vec.gamma[seq1,],nomes[2],row.names=F)
write.csv(vec.betas[seq1,],nomes[3],row.names=F)
write.csv(vec.alpha[seq1,],nomes[4],row.names=F)
}
