
setwd("C:/Users/User/Desktop/2020_Corona_Fun/Coding/Rm1_ToadRDA/Code")

###### Common-Garden raised progeny #############

ToadData_D3 <- read.csv("RDA_MAIN_ALL_ProperNames.csv",stringsAsFactors=TRUE,sep=",")
str(ToadData_D3)

ToadVars_D3 <- ToadData_D3[c(9:15)] #morphology
ToadVars_D3 <- as.data.frame(lapply(ToadVars_D3, as.numeric))
str(ToadVars_D3)

PopCoords_D3 <- ToadData_D3[c(6,7)]
PopCoords_D3 <- as.data.frame(lapply(PopCoords_D3, as.numeric))
str(PopCoords_D3)

#####

#install.packages("raster")
#library(raster)
#climdata <- getData('worldclim',download=TRUE,var='bio',res=5)

points_D3 <- SpatialPoints(PopCoords_D3, proj4string=climdata@crs)
values_D3 <- extract(climdata,points_D3)
envdata_D3_NA <- cbind.data.frame(PopCoords_D3,values_D3)
str(envdata_D3_NA)
envdata_D3 <- envdata_D3_NA[rowSums(is.na(envdata_D3_NA)) == 0,]
str(envdata_D3)

colnames(envdata_D3) <- c("long","lat",
                          "AnnualMeanTemp","MeanDiurnalRange",
                          "Isothermality","TempSeasonality",
                          "MaxTempofWarmestMonth","MinTempofColdestMonth",
                          "TempAnnualRange","MeanTempofWettestQuarter",
                          "MeanTempofDriestQuarter","MeanTempofWarmestQuarter",
                          "MeanTempofColdestQuarter","AnnualPrecipitation",
                          "PrecipitationofWettestMonth","PrecipitationofDriestMonth",
                          "PrecipitationSeasonality","PrecipitationofWettestQuarter",
                          "PrecipitationofDriestQuarter","PrecipitationofWarmestQuarter",
                          "PrecipitationofColdestQuarter")
str(envdata_D3)

rowSums(is.na(envdata_D3_NA))

PopCoords_D3 <- PopCoords_D3[rowSums(is.na(envdata_D3_NA)) == 0,] ##REMOVE THE NA's that resulted from RASTER package

####

#install.packages("rgdal")

library(rgdal)

#extracting aridity index
raster_ai <- raster("global_ai_et0/ai_et0/ai_et0.tif")
points_ai_D3 <- SpatialPoints(PopCoords_D3, proj4string=raster_ai@crs)
identicalCRS(raster_ai,points_ai_D3)
aridity_index_D3 <- extract(raster_ai,points_ai_D3)
str(aridity_index_D3)
aridity_index_D3

#extracting PET
raster_pet <- raster("global-et0_annual.tif/et0_yr/et0_yr.tif")
points_pet_D3 <- SpatialPoints(PopCoords_D3, proj4string=raster_pet@crs)
identicalCRS(raster_pet,points_pet_D3)
potential_evap_D3 <- extract(raster_pet,points_pet_D3) #Potential Evapotranspiration
str(potential_evap_D3)

#calculating Daylength
#install.packages("geosphere")
library("geosphere")
day_length_D3 <- NULL
counter <- 1
for (val in PopCoords_D3$lat) {
  day_length_D3[[counter]] <-   print(
    mean(
      daylength(val, 1:365)
    )
  )
  names(day_length_D3[[counter]]) = names(val)
  counter <- counter + 1
}

str(day_length_D3)


#Combining all extra env variables
envdata_D3 <- cbind.data.frame(envdata_D3,aridity_index_D3,potential_evap_D3,day_length_D3)
envdata_D3 <- as.data.frame(lapply(envdata_D3, as.numeric))

str(envdata_D3)
envdata_D3_NoNA <- envdata_D3[rowSums(is.na(envdata_D3)) == 0,] #Remove NA's that resulted from 3 extra env. measures
str(envdata_D3_NoNA)
envdata_D3_NoNA <- envdata_D3_NoNA[,3:24]
str(envdata_D3_NoNA)

##Remove omitted NA rows from ToadVars list
ToadVars_D3 <- ToadVars_D3[rowSums(is.na(envdata_D3_NA)) == 0,] ##REMOVE THE NA's from original popcoords
ToadVars_D3 <- ToadVars_D3[rowSums(is.na(envdata_D3)) == 0,]
str(ToadVars_D3)

#RDA
#library(vegan)
# running with all env predictors

ToadVars_D3.rda <- rda(ToadVars_D3 ~. , data=envdata_D3_NoNA, scale=T)

RsquareAdj(ToadVars_D3.rda)
signif.full <- anova.cca(ToadVars_D3.rda, parallel=getOption("mc.cores"))
signif.full # test for model significance

