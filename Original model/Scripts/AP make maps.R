rm(list=ls(all=TRUE))
library(ggplot2)
library(RColorBrewer)
library(sp)
library(rgdal)
library(raster)
library(colorRamps)

set.seed(1)
nsurvey=6
cores=colorRampPalette(c('green','yellow','red'))(100)

#predicted prevalence rasters
for (i in 1:nsurvey){ 
  nome=paste('Original model/model predictions/pred',i,'.csv',sep='')
  dat=read.csv(nome,as.is=T)
  pts=read.csv(nome,as.is=T)[c("Longitude","Latitude", "prev")]
  coordinates(pts)=~Longitude+Latitude
  proj4string(pts)=CRS("+init=epsg:4326") # set it to lat-long
  pts = spTransform(pts,CRS("+init=epsg:4326"))
  gridded(pts) = TRUE
  r = raster(pts)
  projection(r) = CRS("+init=epsg:4326")
  plot(r, zlim=c(0,1), col=green2red(99))
  name=paste0("Original model/final rasters/","preds","prev",i,".tif")
  writeRaster(r,name, format="GTiff", overwrite=TRUE)
}

