########################################
# Generating spatiotemporally explicit NDVI values for coordinate fixes in a telemetry dataset
#Author: Michael Butler Brown
#Author Affiliations: Smithsonian Conservation Biology Institute and Giraffe Conservation Foundation

#Load and prepare tracking data
#NOTE: This script takes a while to run and will probably crash if you run it locally without partitioning the data...I've added a block of code for separating the dataframe by country, site, or operation

#Start Fresh
rm(list=ls())

library(sf)
library(dplyr)
library(move)
library(geojsonio)

#Set Working Directory
#Change working Directory
setwd('C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data')

#Load Data
#Note: These data have already been filtered using the following criteria
##Start the day after the unit was deployed
##Stop the day before the unit went stationary or failed
##Stop on January 1, 2022
##All times are local time at the location the unit was deployed
data =read.csv(file="data_continental.csv", head=TRUE, sep=",") 
data =read.csv(file="ndvi_update.csv", head=TRUE, sep=",") 


#Format datetime to class POSIxCT
data$dt <- as.POSIXct(strptime(data$dt, "%Y-%m-%d %H:%M:%S", tz = "UTC")) #Turn date to class POSIXct

#Final filtering and formatting of master dataframe
data1= data #rename to maintain integrity
trackdata = data1 %>%
  #dplyr::filter(Park_location == "Northwest")%>%
  #dplyr::filter(Country == "Tanzania")%>%
  dplyr::filter(Operation == "Northwest_2021")%>%
  #dplyr::filter(Operation %in% c("Bots_2012","Etosha_2021","EtoshaHeights_2021","Northwest_2021","Phinda_2021"))%>%
  dplyr::select(ID,Name,LATITUDE,LONGITUDE,dt,Operation,BATTERY,Park_location,Country,species,subspecies,HACCURACY)#%>%
data$BATTERY = as.numeric(data$BATTERY)
data=droplevels(data)
dim(data)
dim(trackdata)

trackdata=droplevels(trackdata)

str(trackdata)
trackdata$Date <- as.Date(trackdata$dt)
trackdata$Date <- as.factor(trackdata$Date)
trackdata$index <- seq(1:nrow(trackdata))
head(trackdata)
datasf <- trackdata %>%
  dplyr::select(ID,index,Date,LATITUDE,LONGITUDE)
datasf <- st_as_sf(datasf, coords = c('LONGITUDE','LATITUDE'), crs = 4326)
head(datasf)

#######################################
# Initialize rgee
library(rgee)
#ee_Initialize(email='mbbrown62@gmail.com',drive=TRUE, gcs=TRUE)
ee_Initialize()
ee_check()

##############################
## NDVI extraction functions
#Function to add property with time in milliseconds
add_date<-function(feature) {
  date <- ee$Date(ee$String(feature$get("Date")))$millis()
  feature$set(list(date_millis=date))
}

#Join Image and Points based on a maxDifference Filter within a temporal window
tempwin <- 16 #set windows in days. This may depend on the remote sensing data used

maxDiffFilter<-ee$Filter$maxDifference(
  difference=tempwin*24*60*60*1000, #days * hr * min * sec * milliseconds
  leftField= "date_millis", #Date data was collected
  rightField="system:time_start" #Image date
)

# Define the join.
saveBestJoin<-ee$Join$saveBest(
  matchKey="bestImage",
  measureKey="timeDiff"
)

#Function to add property with raster pixel value from the matched image
add_value<-function(feature){
  #Get the best image and select band
  img1<-ee$Image(feature$get("bestImage")) #In this case we are using NDVI
  #Get ID
  id<-ee$Number(feature$get("ID"))
  #Extract geometry from the features
  point<-feature$geometry()
  #Get NDVI value for each point at the desired spatial resolution (argument scale)
  pixel_value<-img1$sample(region=point, scale=250, tileScale = 16, dropNulls = F) 
  #Return the data containing pixel value and ids. Here we are working with NDVI but the name can be changed 
  feature$setMulti(list(PixelVal = pixel_value$first()$get("NDVI"), ID=id))
}

# Function to remove image property from features
removeProperty<- function(feature) {
  #Get the properties of the data
  properties = feature$propertyNames()
  #Select all items except images
  selectProperties = properties$filter(ee$Filter$neq("item", "bestImage"))
  #Return selected features
  feature$select(selectProperties)
}

start<-"2000-01-01"
end<-"2022-01-01"
imagecoll<-ee$ImageCollection('MODIS/006/MOD13Q1')$filterDate(start,end)
band <- "NDVI" #Name of the band to use. You can change to EVI for instance.

datasf$uniq <- rep(1:1000, each=1000)[1:nrow(datasf)] #This is for up to 1 million points. To increase the max number of points, increase the value for max repetitions.

#datasf <- datasf[1:5000,]

start_time <- Sys.time()
dataoutput <- data.frame()
for(x in unique(datasf$uniq)){
  data1 <- datasf %>% filter(uniq == x)
  # Send sf to GEE
  data <- sf_as_ee(data1)
  # Transform day into milliseconds
  data<-data$map(add_date)
  # Apply the join.
  Data_match<-saveBestJoin$apply(data, imagecoll, maxDiffFilter)
  # Add pixel value to the data to the data
  DataFinal<-Data_match$map(add_value)
  # Remove image property from the data
  DataFinal<-DataFinal$map(removeProperty)
  # Transform GEE object in sf
  temp<- ee_as_sf(DataFinal)
  # append
  dataoutput <- rbind(dataoutput, temp)
}
end_time <- Sys.time()

head(dataoutput)
names(dataoutput)[4] <- band
dataoutput
dataoutput$NDVI2 = dataoutput$PixelVal/10000

#library(sfheaders)
#dataoutput.df = sf_to_df(dataoutput)
dataoutput.df = as.data.frame(dataoutput)
write.csv(dataoutput.df, file="Northwest_outputdf.csv")# Write A csv
dataoutput.df$ID_check = dataoutput.df$ID
dataoutput.df= dataoutput.df%>%
  dplyr::select(ID_check, PixelVal,geometry,NDVI2)
ndvi = cbind(trackdata,dataoutput.df)

# Build plot
library(cowplot)
library(ggmap)
library(ggplot2)
theme_set(theme_cowplot())
NDVI_plot <- ggplot(ndvi, aes(x=dt, y=NDVI2))+ 
  geom_point(aes(color=ID))+
  ylab("NDVI")+
  xlab("Date")+
  ggtitle("NDVI over Time")
#theme(axis.text.y = element_text(angle=90))#+ #NOTE: SILENCE THIS LINE OF CODE if you are looking at multiple animals 
#theme(text = element_text(size=20))
NDVI_plot 

#Write the Files to csv
write.csv(ndvi, file="NDVI_update_Nambia2021.csv")# Write A csv
