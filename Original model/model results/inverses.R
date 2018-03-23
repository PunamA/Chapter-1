rm(list=ls(all=TRUE))
set.seed(18)
library('mvtnorm')

for (jjj in 1:6){
  nome=paste('Original model/edited data/AP',jjj,'.csv',sep='')
  dat=read.csv(nome,as.is=T)
  
  nomes=c('Latitude','Longitude','village','Dist_HF','Dist_UC','Elevation',
          'NDVI','Dist_water','Dist_Roads','LST_N','NTL')
  clust=unique(round(dat[,sort(nomes)],4))
  clust=clust[order(clust$village),]
  nclust=nrow(clust); nclust
  
  #distance matrix
  dist1=as.matrix(dist(clust[,c('Latitude','Longitude')]))
  nvil=nrow(dist1)
  
  #discretize rho to improve mixing
  cor1=seq(from=0.05,to=0.95,by=0.01)
  rho.vals=-quantile(dist1,0.1)/log(cor1)
  nrho=length(rho.vals)
  
  #create inverses
  invmat=matrix(NA,nrho,nclust*nclust)
  for (i in 1:nrho){
    invmat[i,]=solve(exp(-dist1/rho.vals[i]))
  }
  
  nome=paste('Original model/model results/inverse',jjj,'.csv',sep='')
  write.csv(invmat,nome,row.names=F)
}
