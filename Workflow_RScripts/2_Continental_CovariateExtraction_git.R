#Covariate Extraction for Movement Analyses
#Author: Michael Butler Brown
#Author Affiliations: Smithsonian Conservation Biology Institute and Giraffe Conservation Foundation

########################################
# Load and prepare tracking data
library(sf)
library(dplyr)
library(move)

#Set Working Directory
#setwd('C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data')
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Load Data
#Note: These data have already been filtered using the following criteria
##Start the day after the unit was deployed
##Stop the day before the unit went stationary or failed
##Stop on January 1, 2022
##All times are local time at the location the unit was deployed
data = read.csv(file="../Data/data_continental.csv", header = TRUE, sep =",")

#Format datetime to class POSIxCT
data$dt <- as.POSIXct(strptime(data$dt, "%Y-%m-%d %H:%M:%S", tz = "UTC")) #Turn date to class POSIXct

#Final filtering and formatting of master dataframe
#Note, this chunk of code allows you to break down the greater dataframe into more manageable chunks for subsequent ctmm analyses if so desired
data1= data #rename to maintain integrity
data = data1 %>%
  #dplyr::filter(ID %in% select)%>%
  #dplyr::filter(Park_location == "Majete")%>%
  #dplyr::filter(Park_location %in% c("MFNP_South","PianUpe"))%>%
  #dplyr::filter(Country == "Tanzania")%>%
  #dplyr::filter(Operation == "Namibia_2019")%>%
  dplyr::select(ID,Name,collarid,LATITUDE,LONGITUDE,dt,Operation,BATTERY,Park_location,Country,species,subspecies,sex,HACCURACY)#%>%
data$BATTERY = as.numeric(data$BATTERY)
data=droplevels(data)

#Define spatial coordinates
data$x = data$LONGITUDE
data$y=data$LATITUDE
coordinates(data) = c("x","y") 
proj4string(data)= CRS("+init=epsg:4326") # Convert to WGS 1984 in decimal degrees (4326 is the code for this reference system). To find your reference system, go to http://spatialreference.org/ref/epsg/

library(maptools)
dsn = ('C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates')
#pa = rgdal::readOGR("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Africa_ProtectedAreas.shp")
pa = rgdal::readOGR(dsn=dsn, layer="Africa_ProtectedAreas")
#shit = as.data.frame(pa)

#SpatialJoin Via Over command
#See for instructions: http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS2_MergingSpatialData_part1_Joins.html
library(sp)
library(dplyr)
library(tidyr)
giraffe.pa = over(data,pa)# Get region data for underlying points
giraffe_pa = spCbind(data,giraffe.pa)#recombine these data
giraffe_pa = as.data.frame(giraffe_pa)
#write.csv(giraffe_pa, file="giraffe_protected_area2.csv")

#Add some levels to subsequent sorting
levels(giraffe_pa$NAME) = c(levels(giraffe_pa$NAME),"None")
levels(giraffe_pa$DESIG) = c(levels(giraffe_pa$DESIG),"None")

df = giraffe_pa %>%
  tidyr::replace_na(list(NAME="None", DESIG="None"))%>%
  dplyr::select(ID,Name, LATITUDE,LONGITUDE,dt,Operation,Park_location,Country,species, subspecies, sex,TEMPERATURE,NAME,DESIG)

unique(df$NAME)  
unique(df$DESIG)  

#Look at Variation in Space Use
landuse = as.data.frame(df%>%
                          dplyr::group_by(ID,DESIG)%>%
                          dplyr::summarise(use =n()))
obs = as.data.frame(df %>%
                      dplyr::group_by(ID)%>%
                      dplyr::summarise(obs=n()))                    
landuse = left_join(landuse,obs,by="ID")
landuse$percent = landuse$use/landuse$obs

# #Creating a stacked barplot to show space use by protected area status
# library(ggplot2)
# ggplot(landuse, aes(fill=DESIG, y=percent, x=ID)) + 
#   geom_bar(position="stack", stat="identity")

################
##Extracting GlobCover
library(raster)
#Reading in GLOBVCOVER Landcover Classification
globcover = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/GLOBCOVER_L4_200901_200912_V2.3.tif")

#Reading in a shapefile to clip the raster layer
library(rgdal)
africa = rgdal::readOGR(dsn=dsn,layer="Africa")

#Crop the landcover classification to the extent of Africa
globcover_crop = crop(globcover,africa)
plot(globcover_crop, main = "Cropped GLOBCOVER to Africa")

#Extracting Raster associated with collar points
raster_extract = raster::extract(globcover_crop,data)
data$raster = raster_extract

#Translating raster Key with code
globcover_key = read.csv("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Globcover2009_Legend.csv", header = T)
data2 = as.data.frame(data)
globcover_2009 = left_join(data2,globcover_key, by="raster")
#write.csv(globcover_2009, file="globcover_2009.csv")

#Look at Variation in Globcover
globcover_summary = as.data.frame(globcover_2009%>%
                                    dplyr::group_by(ID,Label)%>%
                                    dplyr::summarise(use =n()))
obs = as.data.frame(globcover_2009 %>%
                      dplyr::group_by(ID)%>%
                      dplyr::summarise(obs=n()))                    
globcover_summary = left_join(globcover_summary,obs,by="ID")
globcover_summary$percent = globcover_summary$use/globcover_summary$obs

################
##Extracting Human Footprint
library(raster)
#Reading in Human Footprint Layer
HFP2009 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/HFP2009.tif")

# #Crop the landcover classification to the extent of Africa
# #hfp2009_crop = crop(HFP2009,africa)
# #plot(hfp2009_crop, main = "Cropped Human Footprint to Africa")
# #plot(HFP2009)

#Extracting Raster associated with collar points
raster_extract = raster::extract(HFP2009,data)
data$raster_HFP2009 = raster_extract
covariates = as.data.frame(data)
HFP_2009 = covariates %>%
  dplyr::select(ID, raster_HFP2009)
#write.csv(HFP_2009, file="hfp_2009.csv")

################
##Extracting Human Density
library(raster)
#Reading in Human Density Landcover Classification
human_density = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/gpw_v4_population_density_adjusted_to_2015_unwpp_country_totals_rev11_2020_2pt5_min.tif")

#Reading in a shapefile to clip the raster layer
library(rgdal)
#Crop the landcover classification to the extent of Africa
#hfp2009_crop = crop(HFP2009,africa)
#plot(hfp2009_crop, main = "Cropped Human Footprint to Africa")
#plot(human_density)

#Extracting Raster associated with collar points
raster_extract = raster::extract(human_density,data)
data$raster_human_density = raster_extract
covariates2 = as.data.frame(data)
human_density = covariates2 %>%
  dplyr::select(ID, raster_human_density)
#write.csv(human_density, file="human_density.csv")

# #Creating a stacked barplot to show space use by protected area status
# library(ggplot2)
# ggplot(globcover_summary, aes(fill=Label, y=percent, x=ID)) + 
#   geom_bar(position="stack", stat="identity")

##Extracting Above ground woody biomass
###Reading in the rasterfiles
##See Bouvet et al. 2018 for description of raster layers
woody_central_eastern = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/AGB_50m_Central_and_Eastern_Africa-003.tif")
woody_southern = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/AGB_50m_Southern_Africa.tif")
woody_western = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/AGB_50m_Western_Africa.tif")
woody_northern = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/AGB_50m_Northern_Africa.tif")
woody = list(woody_central_eastern,woody_southern,woody_western,woody_northern)

woody_merged <- mosaic(woody_central_eastern, woody_northern, fun=mean)

#Reading in a shapefile to clip the raster layer
library(rgdal)

#Extracting Raster of woody data associated with collar points
raster_east_central = raster::extract(woody_central_eastern,data)
raster_southern = raster::extract(woody_southern,data)
raster_western = raster::extract(woody_western,data)
raster_northern = raster::extract(woody_northern,data)

data$raster_woody_east_central = raster_east_central
data$raster_woody_southern = raster_southern
data$raster_woody_western = raster_western
data$raster_woody_northern = raster_northern
covariates2 = as.data.frame(data)

covariates2[is.na(covariates2)] <- 0
covariates2$raster_woody = covariates2$raster_woody_east_central + covariates2$raster_woody_southern + covariates2$raster_woody_northern + covariates2$raster_woody_western

woody = covariates2 %>%
  dplyr::select(ID, raster_woody)

#write.csv(woody, file="woody.csv")

# #Creating a stacked barplot to show space use by protected area status
# library(ggplot2)# ggplot(globcover_summary, aes(fill=Label, y=percent, x=ID)) + 
#   geom_bar(position="stack", stat="identity")


##Extracting Copernicus Data
###Reading in the rasterfiles
cop_tree_1 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Tree/1_Tree-CoverFraction-layer_EPSG-4326.tif")
cop_tree_2 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Tree/2_Tree-CoverFraction-layer_EPSG-4326.tif")
cop_tree_3 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Tree/3_Tree-CoverFraction-layer_EPSG-4326.tif")
cop_tree_4 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Tree/4_Tree-CoverFraction-layer_EPSG-4326.tif")
cop_tree_5 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Tree/5_Tree-CoverFraction-layer_EPSG-4326.tif")
cop_tree_6 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Tree/6_Tree-CoverFraction-layer_EPSG-4326.tif")

#Extract
library(rgdal)

#Extracting Raster of woody data associated with collar points
raster_cop_tree1 = raster::extract(cop_tree_1,data)
raster_cop_tree2 = raster::extract(cop_tree_2,data)
raster_cop_tree3 = raster::extract(cop_tree_3,data)
raster_cop_tree4 = raster::extract(cop_tree_4,data)
raster_cop_tree5 = raster::extract(cop_tree_5,data)
raster_cop_tree6 = raster::extract(cop_tree_6,data)

data$raster_cop_tree1 = raster_cop_tree1
data$raster_cop_tree2 = raster_cop_tree2
data$raster_cop_tree3 = raster_cop_tree3
data$raster_cop_tree4 = raster_cop_tree4
data$raster_cop_tree5 = raster_cop_tree5
data$raster_cop_tree6 = raster_cop_tree6

covariates3 = as.data.frame(data)

covariates3[is.na(covariates3)] <- 0
covariates3$raster_cop_tree = covariates3$raster_cop_tree1 + covariates3$raster_cop_tree2 + covariates3$raster_cop_tree3 + covariates3$raster_cop_tree4 + covariates3$raster_cop_tree5 + covariates3$raster_cop_tree6

cop_tree = covariates3 %>%
  dplyr::select(ID, raster_cop_tree)
write.csv(cop_tree, file="cop_tree.csv")

##Extracting Copernicus Data
###Reading in the rasterfiles
cop_shrub_1 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Shrub/1_shrub-CoverFraction-layer_EPSG-4326.tif")
cop_shrub_2 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Shrub/2_shrub-CoverFraction-layer_EPSG-4326.tif")
cop_shrub_3 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Shrub/3_shrub-CoverFraction-layer_EPSG-4326.tif")
cop_shrub_4 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Shrub/4_shrub-CoverFraction-layer_EPSG-4326.tif")
cop_shrub_5 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Shrub/5_shrub-CoverFraction-layer_EPSG-4326.tif")
cop_shrub_6 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_Shrub/6_shrub-CoverFraction-layer_EPSG-4326.tif")

#Extract
library(rgdal)

#Extracting Raster of woody data associated with collar points
raster_cop_shrub1 = raster::extract(cop_shrub_1,data)
raster_cop_shrub2 = raster::extract(cop_shrub_2,data)
raster_cop_shrub3 = raster::extract(cop_shrub_3,data)
raster_cop_shrub4 = raster::extract(cop_shrub_4,data)
raster_cop_shrub5 = raster::extract(cop_shrub_5,data)
raster_cop_shrub6 = raster::extract(cop_shrub_6,data)

data$raster_cop_shrub1 = raster_cop_shrub1
data$raster_cop_shrub2 = raster_cop_shrub2
data$raster_cop_shrub3 = raster_cop_shrub3
data$raster_cop_shrub4 = raster_cop_shrub4
data$raster_cop_shrub5 = raster_cop_shrub5
data$raster_cop_shrub6 = raster_cop_shrub6

covariates3 = as.data.frame(data)

covariates3[is.na(covariates3)] <- 0
covariates3$raster_cop_shrub = covariates3$raster_cop_shrub1 + covariates3$raster_cop_shrub2 + covariates3$raster_cop_shrub3 + covariates3$raster_cop_shrub4 + covariates3$raster_cop_shrub5 + covariates3$raster_cop_shrub6

cop_shrub = covariates3 %>%
  dplyr::select(ID, raster_cop_shrub)
write.csv(cop_shrub, file="cop_shrub.csv")


##Extracting Copernicus Data
###Reading in the rasterfiles
cop_bare_1 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_bare/1_bare-CoverFraction-layer_EPSG-4326.tif")
cop_bare_2 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_bare/2_bare-CoverFraction-layer_EPSG-4326.tif")
cop_bare_3 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_bare/3_bare-CoverFraction-layer_EPSG-4326.tif")
cop_bare_4 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_bare/4_bare-CoverFraction-layer_EPSG-4326.tif")
cop_bare_5 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_bare/5_bare-CoverFraction-layer_EPSG-4326.tif")
cop_bare_6 = raster("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data/Covariates/Copernicus_bare/6_bare-CoverFraction-layer_EPSG-4326.tif")

#Extract
library(rgdal)

#Extracting Raster of woody data associated with collar points
raster_cop_bare1 = raster::extract(cop_bare_1,data)
raster_cop_bare2 = raster::extract(cop_bare_2,data)
raster_cop_bare3 = raster::extract(cop_bare_3,data)
raster_cop_bare4 = raster::extract(cop_bare_4,data)
raster_cop_bare5 = raster::extract(cop_bare_5,data)
raster_cop_bare6 = raster::extract(cop_bare_6,data)

data$raster_cop_bare1 = raster_cop_bare1
data$raster_cop_bare2 = raster_cop_bare2
data$raster_cop_bare3 = raster_cop_bare3
data$raster_cop_bare4 = raster_cop_bare4
data$raster_cop_bare5 = raster_cop_bare5
data$raster_cop_bare6 = raster_cop_bare6

covariates3 = as.data.frame(data)

covariates3[is.na(covariates3)] <- 0
covariates3$raster_cop_bare = covariates3$raster_cop_bare1 + covariates3$raster_cop_bare2 + covariates3$raster_cop_bare3 + covariates3$raster_cop_bare4 + covariates3$raster_cop_bare5 + covariates3$raster_cop_bare6

cop_bare = covariates3 %>%
  dplyr::select(ID, raster_cop_bare)
write.csv(cop_bare, file="cop_bare.csv")






