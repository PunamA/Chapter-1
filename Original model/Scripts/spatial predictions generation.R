rm(list=ls(all=TRUE))
set.seed(1)
nsurvey=6
source('Original model/Scripts/spatial predictions functions.R')

options(warn=2)
for (i in 1:nsurvey){ 
  print(c('AP',i))
  #-------------------------------------------------
  #OLD STUFF
  nome=paste('Original model/edited data/AP',i,'.csv',sep='')
  dat=read.csv(nome,as.is=T)
  
  
  #create old and new datasets fo coordinates
  nomes=c('Latitude','Longitude','village','Dist_HF','Dist_UC','Elevation',
          'NDVI','Dist_water','Dist_Roads','LST_N','NTL')
  
  clust=unique(round(dat[,sort(nomes)],4))
  clust=clust[order(clust$village),]
  
  ind=which(colnames(clust)%in%c('Latitude','Longitude','village'))
  wmat.old=data.matrix(cbind(1,clust[,-ind]))
  coord.old=clust[,c('Latitude','Longitude')]
  
  #distance matrix
  dist1=as.matrix(dist(coord.old))

  #discretize rho to improve mixing
  cor1=seq(from=0.05,to=0.95,by=0.01)
  rho.vals=-quantile(dist1,0.1)/log(cor1)
  
  #-------------------------------------------------
  #NEW STUFF
  nome=paste('Original model/edited data/pred',i,'.csv',sep='')
  dat.new=read.csv(nome,as.is=T)
  coord.new=dat.new[,c('Latitude','Longitude')]
  ind=which(colnames(coord.new)%in%c('Latitude','Longitude'))
  tmp=dat.new[,-ind]
  tmp1=tmp[,sort(colnames(tmp))]
  wmat.new=data.matrix(cbind(1,tmp1))
  
  #-------------------------------------------------
  #GET PARAMETERS
  nomes=paste('Original model/model results/',c('alpha','gamma','inverse','outros'),i,'.csv',sep='')
  
  alpha=data.matrix(read.csv(nomes[1],as.is=T))
  gamma=data.matrix(read.csv(nomes[2],as.is=T))
  inverse=data.matrix(read.csv(nomes[3],as.is=T))
  outros=read.csv(nomes[4],as.is=T)

  # colnames(wmat.new); colnames(wmat.old); colnames(gamma) #CHECK THIS
  #-------------------------------------------------
  #MAKE SPATIAL PREDICTIONS
  fim=make.spat.pred(coord.old=coord.old,
                     coord.new=coord.new,
                     outros=outros,
                     invSigma=inverse,
                     gamma=gamma,
                     wmat.old=wmat.old,
                     wmat.new=wmat.new,
                     alpha=alpha,
                     rho.vals=round(rho.vals,5))      
  
  nome=paste('Original model/model predictions/pred',i,'.csv',sep='')
  write.csv(fim,nome,row.names=F)
}
