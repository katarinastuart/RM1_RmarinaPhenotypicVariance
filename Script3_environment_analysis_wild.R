
setwd("C:/Users/User/Desktop/2020_Corona_Fun/Coding/Rm1_ToadRDA/Code")

install.packages("DescTools")
library(ggplot2)

if (!require("remotes")) install.packages("remotes")
remotes::install_github("AndriSignorell/DescTools")
library(DescTools)

######Wild Caugth CANE TOAD#############

ToadData_D3 <- read.csv("RDA_MASTER_WA_ADULT.csv",stringsAsFactors=TRUE,sep=",")
str(ToadData_D3)

ToadVars_D3 <- ToadData_D3[c(5:14)]
ToadVars_D3 <- as.data.frame(lapply(ToadVars_D3, as.numeric))
str(ToadVars_D3)

PopCoords_D3 <- ToadData_D3[c(2,3)]
PopCoords_D3 <- as.data.frame(lapply(PopCoords_D3, as.numeric))
str(PopCoords_D3)

#####
#require(devtools)
#install_version("raster", version = "3.3-7", repos = "http://cran.us.r-project.org")

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

#library(rgdal)

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
for (val in PopCoords_D3$Latitude) {
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
str(envdata_D3_NoNA[3:24])

envdata_D3_NoNA <- envdata_D3_NoNA[3:24]

##Variable Predictions

#library(psych)
#pdf("pair_panels_1.pdf")
#pairs.panels(envdata_D3_NoNA[,3:24], scale=T)
#dev.off()

#VIFpredictions_D3 <- envdata_D3_NoNA[c(3,6,8,22)]
#str(VIFpredictions_D3)
#pairs.panels(VIFpredictions_D3, scale=T)

##Remove omitted NA rows from ToadVars list


ToadVars_D3 <- ToadVars_D3[rowSums(is.na(envdata_D3_NA)) == 0,] ##REMOVE THE NA's from original popcoords
ToadVars_D3 <- ToadVars_D3[rowSums(is.na(envdata_D3)) == 0,]
str(ToadVars_D3)

#RDA
#library(vegan)

#ToadVars_D3.rda <- rda(ToadVars_D3 ~ TempSeasonality +
#                         MeanTempofDriestQuarter + MeanTempofColdestQuarter +  AnnualPrecipitation + 
#                         PrecipitationSeasonality +  aridity_index_D3 + potential_evap_D3 + day_length_D3, data=envdata_D3_NoNA, scale=T)

ToadVars_D3.rda <- rda(ToadVars_D3 ~. , data=envdata_D3_NoNA, scale=T)

RsquareAdj(ToadVars_D3.rda)

signif.full <- anova.cca(ToadVars_D3.rda, parallel=getOption("mc.cores"))
signif.full # test for model significance

str(ToadVars_D3)
str(envdata_D3_NoNA)

##################
## Collected results from each RDA run below

state <- c("1-FG", "2-HI", "3-QLD", "4-NT", "5-WA","1-FG", "2-HI", "3-QLD", "4-NT", "5-WA","1-FG", "2-HI", "3-QLD", "4-NT", "5-WA")
sex <- c(rep("Female",5),rep("Male",5), rep("Adult",5))
r.squred <- c(0.07534771, 0.1677144, 0.05174366, 0.1050439, 0.0638219,0.349391, 0.1457928, 0.09462855, 0.1108208, 0.05929375, 0.3003761, 0.1090402, 0.04238388, 0.1064017, 0.04903424)

env.pred <- data.frame(state, sex, r.squred)
env.pred

PanelA <- ggplot(env.pred, aes(x=state, y=r.squred, colour=sex, group=sex)) + 
  geom_line(size=1.0) +
  geom_point(size=3, shape=21, fill="white") +
  xlab("Population") +
  ylab("Adjusted R-Squared") +
  scale_colour_manual("", values = c("Female" = "#000000", "Male" = "#A9A9A9", "Adult" = "#0000FF")) +
  theme_bw() +
  theme(legend.justification=c(1,0),
        legend.position=c(0.95,0.75))  +
  theme(legend.title = element_text(colour="black", size=12, 
                                    face="bold")) +
  theme(legend.text = element_text(colour="black", size=12, 
                                   face="bold")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

##################
################## ENv PCA ##################
##################

library(factoextra)
library(FD)
library(tibble)
library("ggpubr")

######Wild caught CANE TOAD data set #############

ToadData_D3 <- read.csv("Cane_Toad_Master_Sheet_SC_adult_NoNSW.CSV",stringsAsFactors=TRUE,sep=",")
str(ToadData_D3)

ToadVars_D3 <- ToadData_D3[c(10:19)]
ToadVars_D3 <- as.data.frame(lapply(ToadVars_D3, as.numeric))
str(ToadVars_D3)

PopCoords_D3 <- ToadData_D3[c(5,6)]
PopCoords_D3 <- as.data.frame(lapply(PopCoords_D3, as.numeric))
str(PopCoords_D3)

State <- ToadData_D3[c(3)]
State <- as.data.frame(lapply(State, as.factor))
str(State)

#####

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

#library(rgdal)

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
for (val in PopCoords_D3$Latitude) {
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
State_NoNA <- as.data.frame(State[rowSums(is.na(envdata_D3_NA)) == 0,])
str(State_NoNA)

envdata_D3 <- cbind.data.frame(State_NoNA,envdata_D3,aridity_index_D3,potential_evap_D3,day_length_D3)
str(envdata_D3)

envdata_D3_NoNA <- envdata_D3[rowSums(is.na(envdata_D3)) == 0,] #Remove NA's that resulted from 3 extra env. measures
str(envdata_D3_NoNA)

names(envdata_D3_NoNA)[names(envdata_D3_NoNA) == "State[rowSums(is.na(envdata_D3_NA)) == 0, ]"] <- "Sample_State"
str(envdata_D3_NoNA)

##PCA all env variables

wdbc.pr <- prcomp(envdata_D3_NoNA[,c(4:25)], center = TRUE, scale = TRUE)
wdbc.pr

summary(wdbc.pr)

#library("factoextra")
PanelB <- fviz_pca_ind(wdbc.pr, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             fill.ind = envdata_D3_NoNA$Sample_State, 
             col.ind = "black", 
             palette = c("#AA3377","#CCBB44", "#228833","#EE6677","#66CCEE"), 
             addEllipses = TRUE,
             label = "var",
             col.var = "black",
             repel = TRUE,
             legend.title = "Population") +
  theme(plot.title = element_text(color="white", hjust = 0.5))+
  xlab("Axis 1 (58.3%)") +
  ylab("Axis 1 (18.0%)") +
  theme(legend.justification=c(1,0),
        legend.position=c(1,0.80))  +
  theme(legend.title = element_text(colour="black", size=12, 
                                    face="bold")) +
  theme(legend.text = element_text(colour="black", size=12, 
                                   face="bold")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

#########
#########Combining Panel A and B into one fig
#########

library(ggplot2)
library(grid)


png(filename = paste("Figure2_EnvRelations.png", sep = ""), width = 1200, height = 600, units = "px", res = NA)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = c(600, 600))))

vplayout <- function(x, y)
  viewport(layout.pos.row = x, layout.pos.col = y)

print(PanelA, vp = vplayout(1,1))
print(PanelB, vp = vplayout(1,2))
dev.off()

