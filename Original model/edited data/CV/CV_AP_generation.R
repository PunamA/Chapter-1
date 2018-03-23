###### Subsetting data for 10-fold cross-validation
rm(list=ls())
set.seed(1)
#to round up
foo <- function(x) ceiling(max(x)/10)*10

for(l in 1:6){ 
  AP=paste("edited data/AP",l,".csv",sep="")
  #Bring in Data
  dat=read.csv(AP)
  summary(dat)
  
  comid.uniq=c(unique(dat$village), rep(NA, foo(max(dat$village))-max(dat$village)))
  comid.uniq1=sample(comid.uniq,size=length(comid.uniq), replace = FALSE) 
  tab=matrix(comid.uniq1,ncol=10) 
  #Perform 10 fold cross validation
  for(i in 1:10){
    CV_Data2=dat
    subset1=tab[,i]
    subset2=subset1[!is.na(subset1)]
    cond=dat$village%in%subset2
    #Segement your data by fold using the which() function 
    CV_Data2[cond, "Mal"] = NA
    #check number of missing added per dataset
    print(sum(is.na(CV_Data2$Mal)))
    #export data
    CV_Data2$test=is.na(CV_Data2$Mal)
    aggregate(test~village,data=CV_Data2,mean)
    cv = paste("edited data/CV/cv_AP",l, "data_",i,".csv",sep="")
    write.csv(x = CV_Data2, file = cv)
  }
} 
###NOTES: you should check the data to make sure it has truly sampled without replacement...