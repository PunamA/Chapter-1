#########################################################
##Master file
##Author: Punam Amratia
##Comments:
#this script is to run each script in order to create the
#predictions for chapter 1.
#########################################################

#--------------------------------------------------------
#Original BYD single survey data model

#generate single survey datasets
source('Original model/edited data/AP generation.R')

#generate unsampled prediction dataset for each survey
source('Original model/edited data/Pred generation.R')

#run the main gibbs model (not for validation)
source('Original model/gibbs geospatial.R')

#create the inverse distance matrices for each survey
source('Original model/model results/inverses.R')

#create spatial prediction
source('Original model/spatial predictions generation.R')

#finally run the make AP maps to create the rasters
source('Original model/AP make maps.R')

#########################################################
#  From here move rasters to ArcGIS for figure creation

