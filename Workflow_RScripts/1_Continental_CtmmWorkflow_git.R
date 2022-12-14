##THE CONTINENTAL
##A Streamlined CTMM Workflow
#Author: Michael Butler Brown
#Author Affiliations: Smithsonian Conservation Biology Institute and Giraffe Conservation Foundation
#For Vignettes on Variograms and Model Selection, see: https://cran.r-project.org/web/packages/ctmm/vignettes/variogram.html

#Start Fresh
rm(list=ls())

#Load Relevant packages
library(plyr)
library(dplyr)
library(tidyverse)
library(raster)
library(rgeos)
library(rgbif)
library(viridis)
library(gridExtra)
library(rasterVis)
library(maps)
library(ggmap)
library(maptools)
library(Hmisc)
library(move)
library(lubridate)
library(mapview)
library(rgdal)
library(ctmm)

#Set Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Read in Movement Data
#Load Data
#Note: These data have already been filtered using the following criteria
##Start the day after the unit was deployed
##Stop the day before the unit went stationary or failed
##Stop on January 1, 2022
##All times are local time at the location the unit was deployed
data = read.csv(file="../Data/data_continental.csv", header = TRUE, sep =",")
#data2= read.csv(file="../Data/data_continental_add.csv", header = TRUE, sep =",")


data_export = data %>%
  #dplyr::filter()%>%
  dplyr::select(ID,LATITUDE,LONGITUDE,dt,Operation)
head(data_export)
library(tibble)
write.csv(data_export, file="TheContinental_GEE_Export.csv")# Write A csv #IF you want, you can write a .csv from the filtered dataset

#Names = unique(data$ID)
#write.csv(Names, file="names_doublecheck.csv")

#Format datetime to class POSIxCT
data$dt <- as.POSIXct(strptime(data$dt, "%Y-%m-%d %H:%M:%S", tz = "UTC")) #Turn date to class POSIXct

#From EXCEL
read.excel <- function(header=TRUE,...) {
  read.table("clipboard",sep="\t",header=header,...)
} 
select = read.excel()
select2 = as.list(select$ID)

select2 = c("IRI2016-3050","IRI2016-3052","IRI2016-3047","IRI2016-3046","IRI2016-3049","IRI2016-3042","IRI2016-3044","IRI2016-3055",   
            "IRI2016-3187","IRI2016-3269","IRI2016-3267")
select2 = c("IRI2016-3523","IRI2016-3530","IRI2016-3539","IRI2016-3544","IRI2016-3526",   
"IRI2016-3527","IRI2016-3529","IRI2016-3531","IRI2016-3528","SAT1631","SAT1630", "SAT91 GIRtwyfel", "SAT93 GIRtwyfel")

select3 = c("0-3601999", "0-3601919","0-3602234")


#Final filtering and formatting of master dataframe
#Note, this chunk of code allows you to break down the greater dataframe into more manageable chunks for subsequent ctmm analyses if so desired
data1= data #rename to maintain integrity
data = data1 %>%
  dplyr::filter(ID %in% select2)%>%
  #dplyr::filter(Park_location == "misc_1")%>%
  #dplyr::filter(Park_location %in% c("Zambezi","Otjonzondjupa","SouthKunene"))%>%
  #dplyr::filter(Park_location %in% c("Samburu","Melako","Sera","Leparua"))%>%
  #dplyr::filter(Country == "Tanzania")%>%
  #dplyr::filter(Unit_type == "CeresTag")%>%
  #dplyr::filter(Operation == "Namibia_2017")%>%
  dplyr::select(ID,Name,LATITUDE,LONGITUDE,dt,Operation,BATTERY,Park_location,Country,species,subspecies,HACCURACY)#%>%
data$BATTERY = as.numeric(data$BATTERY)
data=droplevels(data)

##Filter out Individuals with fewer than 400 fixes
data_filter_select = as.data.frame(data %>%
                                     dplyr::group_by(ID)%>%
                                     dplyr::summarise(total=n())%>%
                                     dplyr::filter(total > 400))
select_id = unique(data_filter_select$ID)
data = data %>%
  dplyr::filter(ID %in% select_id)

#Remove Duplicates in Datetime
nrow(data) # Check dimensions of original data
data<- data[!duplicated(data[c('ID','dt')]),] #Filter out duplicated timestamps
nrow(data) #Check dimensions of the filtered data

# Remove near timestamp repeats (fixes < 30min apart)
data = data %>%
  arrange(ID,dt)
timeinsec <- matrix(as.numeric(data$dt))
deltatime <- matrix(diff(timeinsec)) #create a new table of change in time measurements of data
colnames(deltatime, do.NULL = FALSE) #create collumn names
colnames(deltatime) <- "Delta_Time"
deltatime <- rbind(as.vector(3600, mode = "numeric"), deltatime) #add "3600" row to beginning to make length of deltatime = timeinsec
deltatime <- as_tibble(deltatime) #convert to a tibble to use filter, a dplyr function
data2 <- cbind(data, deltatime) #bind deltatime to clean_giraffe data
data2 <- data2 %>% filter(Delta_Time > 1800) #Only keep points with a time change of >30 min (1800 sec)
nrow(data) - nrow(data2) #number of datapoints removed
data <- data2 #set clean_giraffe equal to identical data frame without faulty readings
data=droplevels(data)

#Plot in space with dynamic map to find potential outlier points
#Define spatial coordinates
 data$x = data$LONGITUDE
 data$y=data$LATITUDE
# proj4string(data2)= CRS("+init=epsg:4326") # Convert to WGS 1984 in decimal degrees (4326 is the code for this reference system). To find your reference system, go to http://spatialreference.org/ref/epsg/
##Make a Dynamic Map
##If you have outlier points and want to specifically identify them, you can create a dynamic map an click on the point in question
##I've silenced this line because if the dataframes are large, it takes a long time to plot
#mapview(data2, zcol="ID")

##CTMM
library(ctmm)
library(sp)
library(rgdal)
library(pbapply) #required for progress bar only
library(gridBase)
library(ggplot2)

# Create output directory if needed
out_dir <- file.path("C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Output")
if (!dir.exists(out_dir))
  dir.create(out_dir)

akde_dir <- file.path(paste0(out_dir,"/AKDE/Plot"))
shp_dir <- file.path(paste0(out_dir,"/AKDE/Shapefiles"))
variogram_dir <- file.path(paste0(out_dir,"/Variogram"))

#Reformatting to create telemetry object
giraffe_df<-data.frame(
  individual.local.identifier = data$ID, 
  timestamp = data$dt, 
  location.long = data$x,
  location.lat = data$y)

giraffe_df2 = giraffe_df %>%
  arrange(individual.local.identifier,timestamp)%>%
  dplyr::select(individual.local.identifier, timestamp, location.long,location.lat) 
head(giraffe_df)

#Make some projection cheats to maintain north is up
#See: https://groups.google.com/forum/#!topic/ctmm-user/msBaXtIzHrg
datum <- "WGS84"
lon_0 <- stats::median(giraffe_df$location.long)
lat_0 <- stats::median(giraffe_df$location.lat)
proj <- paste0("+proj=aeqd +lon_0=",lon_0," +lat_0=",lat_0," +datum=",datum) # Manually change the porjection so that telemetry object maps out in appropriate direction

# Create the ctmm telemetry object from giraffe_df ----
giraffe_tel <- as.telemetry(giraffe_df2,timezone="UTC", projection = proj) # , timeformat="POSIXct")
giraffe_tel_backup = giraffe_tel
plot(giraffe_tel, col=rainbow(length(giraffe_tel)))

library(foreach)
foreach(i = 1:length(giraffe_tel)) %do% {
#Identify and remove outliers
OUT = outlie(giraffe_tel[[i]])
bad = OUT$speed > 1.05 # Change this parameter if you want to have a different filter speed
length(bad)
dim(giraffe_tel[[i]])
giraffe_tel[[i]] = giraffe_tel[[i]][!bad,]
}

##########################
filter
foreach(i = 1:length(giraffe_tel)) %do% {
  #Identify and remove outliers
  OUT = outlie(giraffe_tel[[i]])
  bad = OUT$speed > 1.05
  length(bad)
  dim(giraffe_tel[[i]])
  giraffe_tel[[i]] = giraffe_tel[[i]][!bad,]
} 
filter_summary = data.frame(file_name=vector(mode = "numeric", length = length(giraffe_tel)))

foreach(i = 1:length(giraffe_tel)) %do% {
  filter_summary$ID[i] = giraffe_tel[[i]]@info$identity
  OUT = outlie(giraffe_tel[[i]])
  bad = OUT$speed > 1.05
  filter_summary$outliers[i] = length(bad) - dim(giraffe_tel[[i]])[1]
}
filter_summary

#########################
#Making a Loop to generate, print, and save variograms
foreach(i = 1:length(giraffe_tel)) %do% {
  SVF = variogram(giraffe_tel[[i]])
  
  #Save plots of Variograms in the defined directory
  plot.variogram = file.path(paste0(variogram_dir,"/variograms_",giraffe_tel[[i]]@info$identity,".png"))
  png(file= paste0(plot.variogram),
      type="cairo",
      units = "in",
      width = 6.81, height = 5,
      pointsize = 10,
      res = 800) #
  plot(SVF,level =c(0.5,0.95),)
  title(main= paste('Variogram for', giraffe_tel[[i]]@info$identity))
  dev.off()
  }

# Making a loop to Fit movement models
FITS <- AKDES <- list()
for(i in 1:length(giraffe_tel))
{
  GUESS <- ctmm.guess(giraffe_tel[[i]],interactive=FALSE)
  FITS[[i]] <- ctmm.select(giraffe_tel[[i]],GUESS,trace=2)
}
save(FITS, file = "FITS_misc3.Rda")
#FITS = load("FITS_Namibia2016.Rda")
# calculate AKDES on a consistent grid
AKDES <- akde(giraffe_tel,FITS,trace=1)
summary(AKDES[[1]],level = 0.95, level.ud=0.50)

#Generate Summary for AKDES
library(dplyr)
library(foreach)

#Generate of Summary of AKDES
foreach(i = 1:length(AKDES)) %do% {
  summary(AKDES[[i]])
}
summary(AKDES[[1]])

#Extract and Summarize 
fits_summary_95 = data.frame(file_name=vector(mode = "numeric", length = length(FITS)))

foreach(i = 1:length(FITS)) %do% {
  #summary(FITS[[i]], level = 0.95, level.ud= 0.95)
  fits_summary_95$ID[i] = FITS[[i]]@info$identity
  fits_summary_95$model[i] = summary(FITS[[i]])$name
  fits_summary_95$akde_lclL_95[i] = summary(FITS[[i]])$CI[1,1]
  fits_summary_95$akde_ucl_95[i] = summary(FITS[[i]])$CI[1,3]
  fits_summary_95$akde_ml_95[i] = summary(FITS[[i]])$CI[1,2]
  fits_summary_95$tau_pos_lcl[i] = summary(FITS[[i]])$CI[2,1]
  fits_summary_95$tau_pos_ml[i] = summary(FITS[[i]])$CI[2,2]
  fits_summary_95$tau_pos_ucl[i] = summary(FITS[[i]])$CI[2,3]
  fits_summary_95$tau_vel_lcl[i] = summary(FITS[[i]])$CI[3,1]
  fits_summary_95$tau_vel_ucl[i] = summary(FITS[[i]])$CI[3,3]
  fits_summary_95$tau_vel_ml[i]= summary(FITS[[i]])$CI[3,2]
  fits_summary_95$speed_lcl[i] = summary(FITS[[i]])$CI[4,1]
  fits_summary_95$speed_ucl[i] = summary(FITS[[i]])$CI[4,3]
  fits_summary_95$speed_ml[i] = summary(FITS[[i]])$CI[4,2]
  fits_summary_95$DOF_mean[i] = summary(FITS[[i]])$DOF[1]
  fits_summary_95$DOF_area[i] = summary(FITS[[i]])$DOF[2]
  fits_summary_95$DOF_speed[i] = summary(FITS[[i]])$DOF[3]
  }
fits_summary_95
write.csv(fits_summary_95, file="fits_summary_misc3.csv")# Write A csv #IF you want, you can write a .csv from the filtered dataset

akde_summary = data.frame(file_name=vector(mode = "numeric", length = length(AKDES)))
foreach(i = 1:length(AKDES)) %do% {
akde_summary$ID[i] = AKDES[[i]]@info$identity
akde_summary$akde_95_lcl[i] = summary(AKDES[[i]], level=.95,level.ud =0.95)$CI[1,1]
akde_summary$akde_95_ML[i] = summary(AKDES[[i]], level=.95,level.ud =0.95)$CI[1,2]
akde_summary$akde_95_ucl[i] = summary(AKDES[[i]], level=.95,level.ud =0.95)$CI[1,3]
akde_summary$akde_95_DOF_area[i] = summary(AKDES[[i]], level=.95,level.ud =0.95)$DOF[1]
akde_summary$akde_95_DOF_bandwidth[i] = summary(AKDES[[i]], level=.95,level.ud =0.95)$DOF[2]
akde_summary$akde_50_lcl[i] = summary(AKDES[[i]], level=.95,level.ud =0.50)$CI[1,1]
akde_summary$akde_50_ML[i] = summary(AKDES[[i]], level=.95,level.ud =0.50)$CI[1,2]
akde_summary$akde_50_ucl[i] = summary(AKDES[[i]], level=.95,level.ud =0.50)$CI[1,3]
akde_summary$akde_50_DOF_area[i] = summary(AKDES[[i]], level=.95,level.ud =0.50)$DOF[1]
akde_summary$akde_50_DOF_bandwidth[i] = summary(AKDES[[i]], level=.95,level.ud =0.50)$DOF[2]
}
akde_summary
write.csv(akde_summary, file="akde_summary_misc3.csv")# Write A csv #IF you want, you can write a .csv from the filtered dataset
          
#Plot out the meta for this group
meta_95 = meta(AKDES, variable="area", level=0.95,level.UD=0.95)
meta_50 = meta(AKDES, variable="area", level=0.95,level.UD=0.50)

#plot AKDES
  #Save plots of Home Ranges in the defined directory
  plot.ctmm = file.path(paste0(akde_dir,"/ctmm_",giraffe_tel[[i]]@info$identity,".png"))
  png(file= paste0(plot.ctmm),
      type="cairo",
      units = "in",
      width = 6.81, height = 5,
      pointsize = 10,
      res = 800) #
  plot(AKDES[[i]])
  title(main= paste('Home Range Estimate for', giraffe_tel[[i]]@info$identity))
  plot(giraffe_tel[[i]],add=TRUE)
  dev.off()
}

# color to be spatially distinct
COL <- color(AKDES,by='individual')

#Calculate Extent
EXT <- extent(AKDES,level=0.95)

#Plotting all AKDES Over each other and trying to extract spatial polygon dataframe
plot(AKDES)

#Loop through AKDES to generate spatial polygons and export shapefiles
#Save plots of Home Ranges in the defined directory
library(rgdal)

foreach(i = 1:length(AKDES)) %do% {
  ud = AKDES[[i]]
  ud_95 = SpatialPolygonsDataFrame.UD(ud,level.UD=0.95)
  ud_50 = SpatialPolygonsDataFrame.UD(ud,level.UD=0.50)
  
  writeOGR(ud_95, shp_dir,paste0("ud_95_",AKDES[[i]]@info$identity) , driver="ESRI Shapefile")
  writeOGR(ud_50, shp_dir,paste0("ud_50_",AKDES[[i]]@info$identity) , driver="ESRI Shapefile")
}

#plot. = file.path(paste0(plot_dir,"/ctmm_",giraffe_tel[[i]]@info$identity,".png"))
png(file= paste0(plot.ctmm),
    type="cairo",
    units = "in",
    width = 6.81, height = 5,
    pointsize = 10,
    res = 800) #
plot(AKDES[[i]])
title(main= paste('Home Range Estimate for', giraffe_tel[[i]]@info$identity))
plot(giraffe_tel[[i]],add=TRUE)
dev.off()
