make.spat.pred=function(coord.old,coord.new,outros,invSigma,gamma,wmat.old,wmat.new,alpha,rho.vals){
  nsamp=nrow(coord.new)
  nsim=nrow(gamma)
  nvil=nrow(coord.old)
  
  coord.new$prev=coord.new$lo=coord.new$hi=coord.new$missing=NA
  for(k in 1:nsamp){
    print(k)
    dist=sqrt(((coord.new$Longitude[k]-coord.old$Longitude)^2)+((coord.new$Latitude[k]-coord.old$Latitude)^2))
    media.nova=wmat.new[k,]%*%t(gamma)
      
    alpha.new=rep(NA,nsim)
    for (j in 1:nsim){
      rho=round(outros$rho[j],5)
      ind=which(rho==rho.vals)
      inv.Sigma.oo1=matrix(invSigma[ind,],nvil,nvil)
      
      gamma1=t(t(gamma[j,]))
      err=alpha[j,]-(wmat.old%*%gamma1)

      Sigmaon=exp(-dist/rho)
      mu=media.nova[j] + Sigmaon%*%inv.Sigma.oo1%*%err
      Sigma=outros[j,'sig2']*(1-Sigmaon%*%inv.Sigma.oo1%*%t(t(Sigmaon)))
      
      if (Sigma>0) alpha.new[j]=rnorm(1,mean=mu,sd=sqrt(Sigma))
    }
    res=pnorm(alpha.new) #calculate prevalence using probit model
    coord.new[k,c('prev','lo','hi','missing')]=c(mean(res,na.rm=T),
                                                 quantile(res, c(0.025,0.975),na.rm=T),
                                                 sum(is.na(res)))
  }
  coord.new
}
