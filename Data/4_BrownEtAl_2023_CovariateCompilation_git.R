# Covariate Compilation for for Continental Scale Giraffe Movement Analyses
#Author: Michael Butler Brown
#Author Affiliations: Smithsonian Conservation Biology Institute and Giraffe Conservation Foundation

# Compiling datasets

rm(list=ls())
########################################
# Load and prepare tracking data
library(sf)
library(dplyr)
library(move)
library(ggplot2)
library(tidyr)
library(dplyr)

#trackdata <- read.csv("clean_giraffe.csv", header = T)
#Set Working Directory
setwd('C:/Users/mbbro/Documents/GitHub/TheContinental_Giraffe/Data')

#Load Data
collar_key = read.csv(file="ctmm_key2.csv",head=TRUE,sep=",")
collar_key$ID = collar_key$collarid

##Change the names of location for northern Kenya giraffe
collar_key$Park_location = replace(collar_key$Park_location, collar_key$Park_location =="Leparua","NorthernKenya")#For replacing characters with other characters
collar_key$Park_location = replace(collar_key$Park_location, collar_key$Park_location =="Loisaba","NorthernKenya")#For replacing characters with other characters
collar_key$Park_location = replace(collar_key$Park_location, collar_key$Park_location =="Mpala","NorthernKenya")#For replacing characters with other characters
collar_key$Park_location = replace(collar_key$Park_location, collar_key$Park_location =="WestGate","NorthernKenya")#For replacing characters with other characters
collar_key$Park_location = replace(collar_key$Park_location, collar_key$Park_location =="Sera","NorthernKenya")#For replacing characters with other characters
collar_key$Park_location = replace(collar_key$Park_location, collar_key$Park_location =="Samburu","NorthernKenya")#For replacing characters with other characters
collar_key$Park_location = replace(collar_key$Park_location, collar_key$Park_location =="Biliquo_Bulesa","NorthernKenya")#For replacing characters with other characters
collar_key$Park_location = replace(collar_key$Park_location, collar_key$Park_location =="Shaba","NorthernKenya")#For replacing characters with other characters

collar_key2 = collar_key %>%
  dplyr::filter(!Variogram_inspection =="bad")%>% #remove collars that visual inspection of variograms revealed to be bad
  dplyr::filter(!Number_of_fixes < 500) # remove collars with fewer than 500 fixes

#Make a table of Collar locations
library(sjPlot)
collar_key_summary = as.data.frame(collar_key2 %>%
  dplyr::group_by(Country, Park_location,taxa,sex)%>%
  dplyr::summarise(total = n()))

collar_summary_table = collar_key_summary %>%
  tidyr::spread(sex,total)

collar_summary_table[is.na(collar_summary_table)] = 0
collar_summary_table  

 tab_df(collar_summary_table,
        title = "Supplementary Table 1: Summary of Deployed Units", #always give your tables titles
        file = "Supp_Table1_deployed_units.doc",
        use.viewer = TRUE)
 
#####################################################
 
#Collate the model summaries calculated in step 1
fits_df1 = read.csv(file="fits_summary_Ethiopia.csv", head=TRUE, sep=",")
fits_df2 = read.csv(file="fits_summary_Kidepo.csv", head=TRUE, sep=",")
fits_df3 = read.csv(file="fits_summary_MFNP.csv", head=TRUE, sep=",")
fits_df4 = read.csv(file="fits_summary_Namibia_NW2019.csv", head=TRUE, sep=",")
fits_df5 = read.csv(file="fits_summary_NamibiaNotNw.csv", head=TRUE, sep=",")
fits_df6 = read.csv(file="fits_summary_Namibia_NW2018.csv", head=TRUE, sep=",")
fits_df7 = read.csv(file="fits_summary_Niger.csv", head=TRUE, sep=",")
fits_df8 = read.csv(file="fits_summary_NorthernKenya.csv", head=TRUE, sep=",")
fits_df9 = read.csv(file="fits_summary_Phinda.csv", head=TRUE, sep=",")
fits_df10 = read.csv(file="fits_summary_pian.csv", head=TRUE, sep=",")
fits_df11 = read.csv(file="fits_summary_SouthAfrica.csv", head=TRUE, sep=",")
fits_df12 = read.csv(file="fits_summary_zakouma.csv", head=TRUE, sep=",")
fits_df13 = read.csv(file="fits_summary_AWT_Ceres.csv", head=TRUE, sep=",")
fits_df14 = read.csv(file="fits_summary_Namibia_NW2021.csv", head=TRUE, sep=",")
fits_df15 = read.csv(file="fits_summary_Save_Amboseli.csv", head=TRUE, sep=",")
fits_df16 = read.csv(file="fits_summary_MFNP_South.csv", head=TRUE, sep=",")
fits_df17 = read.csv(file="fits_summary_CeresTag.csv", head=TRUE, sep=",")
fits_df18 = read.csv(file="fits_summary_Garamba2020.csv", head=TRUE, sep=",")
fits_df19 = read.csv(file="fits_summary_Serengeti.csv", head=TRUE, sep=",")
fits_df20 = read.csv(file="fits_summary_Namibia2017.csv", head=TRUE, sep=",")
fits_df21 = read.csv(file="fits_summary_Gambella_Tarangire.csv", head=TRUE, sep=",")
fits_df22 = read.csv(file="fits_summary_Zimbabwe.csv", head=TRUE, sep=",")
fits_df23 = read.csv(file="fits_summary_NorthernKenya3.csv", head=TRUE, sep=",")
fits_df24 = read.csv(file="fits_summary_southernAfrica1.csv", head=TRUE, sep=",")
fits_df25 = read.csv(file="fits_summary_Misc_Namibia.csv", head=TRUE, sep=",")
fits_df26 = read.csv(file="fits_summary_misc_12.csv", head=TRUE, sep=",")
fits_df27 = read.csv(file="fits_summary_misc3.csv", head=TRUE, sep=",")

##Remove tau_vel_ml if you want to match with the above
fits_all = rbind(fits_df1,fits_df2, fits_df3, fits_df4, fits_df5, fits_df6, fits_df7, fits_df8, fits_df9, fits_df10, fits_df11, fits_df12, fits_df13, fits_df15, fits_df16)

fits_all2 = rbind(fits_df14,fits_df17,fits_df18,fits_df19,fits_df20, fits_df21, fits_df22, fits_df23,fits_df24,fits_df25, fits_df26, fits_df27)
fits_all2 = dplyr::select(fits_all2, -tau_pos_ml)

fits_all = rbind(fits_all,fits_all2)
fits_all = fits_all %>%
  dplyr::select(ID,model,tau_pos_lcl,tau_pos_ucl, tau_vel_lcl, tau_vel_ucl, tau_vel_ml, speed_lcl, speed_ucl, speed_ml, DOF_mean, DOF_area, DOF_speed)
  
#Collate the akde summaries calculated in step 1
akde_df1 = read.csv(file="akde_summary_ethiopia.csv", head=TRUE, sep=",")
akde_df2 = read.csv(file="akde_summary_Kidepo.csv", head=TRUE, sep=",")
akde_df3 = read.csv(file="akde_summary_MFNP.csv", head=TRUE, sep=",")
akde_df4 = read.csv(file="akde_summary_Namibia_NW2019.csv", head=TRUE, sep=",")
akde_df5 = read.csv(file="akde_summary_NamibiaNotNw.csv", head=TRUE, sep=",")
akde_df6 = read.csv(file="akde_summary_NamibiaNW2018.csv", head=TRUE, sep=",")
akde_df7 = read.csv(file="akde_summary_Niger.csv", head=TRUE, sep=",")
akde_df8 = read.csv(file="akde_summary_NorthernKenya.csv", head=TRUE, sep=",")
akde_df9 = read.csv(file="akde_summary_Phinda.csv", head=TRUE, sep=",")
akde_df10 = read.csv(file="akde_summary_pian.csv", head=TRUE, sep=",")
akde_df11 = read.csv(file="akde_summary_save.csv", head=TRUE, sep=",")
akde_df12 = read.csv(file="akde_summary_zakouma.csv", head=TRUE, sep=",")
akde_df13 = read.csv(file="akde_summary_AWT_Ceres.csv", head=TRUE, sep=",")
akde_df14 = read.csv(file="akde_summary_Namibia_NW2021.csv", head=TRUE, sep=",")
akde_df15 = read.csv(file="akde_summary_Save_Amboseli.csv", head=TRUE, sep=",")
akde_df16 = read.csv(file="akde_summary_MFNP_South.csv", head=TRUE, sep=",")
akde_df17 = read.csv(file="akde_summary_CeresTag.csv", head=TRUE, sep=",")
akde_df18 = read.csv(file="akde_summary_Garamba2020.csv", head=TRUE, sep=",")
akde_df19 = read.csv(file="akde_summary_Serengeti.csv", head=TRUE, sep=",")
akde_df20 = read.csv(file="akde_summary_Namibia2017.csv", head=TRUE, sep=",")
akde_df21 = read.csv(file="akde_summary_Gambella_Tarangire.csv", head=TRUE, sep=",")
akde_df22 = read.csv(file="akde_summary_Zimbabwe.csv", head=TRUE, sep=",")
akde_df23 = read.csv(file="akde_summary_NorthernKenya3.csv", head=TRUE, sep=",")
akde_df24 = read.csv(file="akde_summary_southernAfrica1.csv", head=TRUE, sep=",")
akde_df25 = read.csv(file="akde_summary_Misc_Namibia.csv", head=TRUE, sep=",")
akde_df26 = read.csv(file="akde_summary_misc_12.csv", head=TRUE, sep=",")
akde_df27 = read.csv(file="akde_summary_misc3.csv", head=TRUE, sep=",")


akde_all = rbind(akde_df1,akde_df2, akde_df3, akde_df4, akde_df5, akde_df6, akde_df7, akde_df8, akde_df9, akde_df10, akde_df11, akde_df12, akde_df13, akde_df14, akde_df13,akde_df15, akde_df16,akde_df17, akde_df18,akde_df19,akde_df20, akde_df21, akde_df22, akde_df23,akde_df24, akde_df25, akde_df26, akde_df27)

#Combine Model fit, akde, and collar key to generate response variable dataframe
collar_key$ID = collar_key$collarid
response_df = left_join(collar_key, fits_all, by="ID")
response_df = left_join(response_df, akde_all, by="ID")

##Compile the predictor covariates
####NDVI
ndvi_PianMFNPSouth = read.csv("Pian_MFNPSouth_NDVI.csv",header=T)
ndvi_kvnp = read.csv("KVNP_NDVI.csv",header=T)
ndvi_namibia_noNW = read.csv("Namibia_NoNW_NDVI.csv",header=T)
ndvi_chad = read.csv("Chad_NDVI.csv",header=T)
ndvi_Tanzania = read.csv("Tanzania_NDVI.csv",header=T)
ndvi_NorthernKenya = read.csv("NorthernKenya_NDVI.csv",header=T)
ndvi_Amboseli = read.csv("Amboseli_NDVI.csv",header=T)
ndvi_Zimbabwe = read.csv("Zimbabwe_NDVI.csv",header=T)
ndvi_Botswana = read.csv("Botswana_NDVI.csv",header=T)
ndvi_Niger = read.csv("Niger_NDVI.csv",header=T)
ndvi_Ethiopia = read.csv("Ethiopia_NDVI.csv",header=T)
ndvi_DRC = read.csv("DRC_NDVI.csv",header=T)
ndvi_df2 = read.csv("dataframe2_update_NDVI.csv",header=T)
ndvi_update = read.csv("NDVI_update_NDVI.csv",header=T)
ndvi_update1 = read.csv("NDVI_update1.csv",header=T)
ndvi_MFNP_2018 = read.csv("MFNP_2018_NDVI_outputdf.csv",header=T)
ndvi_MFNP_2019 = read.csv("MFNP_2019_NDVI.csv",header=T)
ndvi_NamibNW = read.csv("Northwest_outputdf.csv",header=T)
ndvi_NamibNW2017 = read.csv("NDVI_update_NW2017.csv",header=T)
ndvi_NamibNW2018 = read.csv("NDVI_update_NW2018.csv",header=T)
ndvi_NamibNW2019 = read.csv("NDVI_update_NW2019.csv",header=T)

#Combine NDVI Data
ndvi_all = rbind(ndvi_PianMFNPSouth,ndvi_kvnp, ndvi_namibia_noNW,ndvi_chad, ndvi_Tanzania, ndvi_NorthernKenya,ndvi_Amboseli,ndvi_Zimbabwe,ndvi_Botswana,
                 ndvi_Niger, ndvi_Ethiopia,ndvi_DRC)
ndvi_all$Name_export = ndvi_all$ID.x
ndvi_all$dt = ndvi_all$LONGITUDE
ndvi_all = ndvi_all%>%
  dplyr::select(Name_export,dt,NDVI2)

ndvi_append = rbind(ndvi_df2,ndvi_update,ndvi_update1, ndvi_NamibNW2017, ndvi_NamibNW2018, ndvi_NamibNW2019)
ndvi_append$Name_export = ndvi_append$X
ndvi_append$dt = ndvi_append$LONGITUDE
ndvi_append = ndvi_append %>%
  dplyr::select(Name_export,dt,NDVI2)

ndvi_append2 = rbind(ndvi_MFNP_2018,ndvi_MFNP_2019,ndvi_NamibNW)
ndvi_append2 = as.data.frame(ndvi_append2)
ndvi_append2$Name_export = ndvi_append2$Date
ndvi_append2$dt = ndvi_append2$X
ndvi_append2$dt = as.POSIXct(strptime(ndvi_append2$dt,"%Y-%m-%d", tz = "UTC"))
ndvi_append2 = ndvi_append2 %>%
  dplyr::select(Name_export,dt,NDVI2)

ndvi_all = rbind(ndvi_all,ndvi_append,ndvi_append2)
ndvi_all$ID = ndvi_all$Name_export

# 
# trial_plot <- ggplot(crap, aes(x=dt, y=NDVI2))+ 
#   geom_point(aes(color=Name))+
#   ylab("NDVI")+
#   xlab("Date")+
#   ggtitle("NDVI over Time")
# #theme(axis.text.y = element_text(angle=90))#+ #NOTE: SILENCE THIS LINE OF CODE if you are looking at multiple animals 
# #theme(text = element_text(size=20))
# trial_plot 

#summarise NDVI results for covariate input
head(ndvi_all)
ndvi_all$NDVI2 = as.numeric(ndvi_all$NDVI2)
ndvi_summary = as.data.frame(ndvi_all%>%
                               dplyr::group_by(ID)%>%
                               dplyr::summarise(ndvi_mean = mean(NDVI2,na.rm=TRUE),
                                                ndvi_sd = sd(NDVI2,na.rm=TRUE)))

#summarise temperature results for covariate input
temp = read.csv("Temp_export.csv",header=TRUE)
temp_summary = as.data.frame(temp%>%
                               dplyr::group_by(ID)%>%
                               dplyr::summarise(temp_mean = mean(TEMPERATURE,na.rm=TRUE),
                                                temp_sd = sd(TEMPERATURE,na.rm=TRUE)))

#summarise human density results for covariate input
human_density = read.csv("human_density.csv",header=TRUE)
human_density_summary = as.data.frame(human_density %>%
                               dplyr::group_by(ID)%>%
                               dplyr::summarise(human_density_mean = mean(raster_human_density,na.rm=TRUE),
                                                human_density_sd = sd(raster_human_density,na.rm=TRUE)))

#summarise human footprint results for covariate input
hfp_2009= read.csv("hfp_2009.csv",header=TRUE)
hfp2009_summary = as.data.frame(hfp_2009%>%
                                        dplyr::group_by(ID)%>%
                                        dplyr::summarise(hfp2009_mean = mean(raster_HFP2009,na.rm=TRUE),
                                                         hfp2009_sd = sd(raster_HFP2009,na.rm=TRUE)))

#summarise woody vegetation results for covariate input
woody = read.csv("woody.csv",header=TRUE)
woody_summary = as.data.frame(woody%>%
                                  dplyr::group_by(ID)%>%
                                  dplyr::summarise(woody_mean = mean(raster_woody,na.rm=TRUE),
                                                   woody_sd = sd(raster_woody,na.rm=TRUE)))

#summarise protected area results for covariate input
pa = read.csv("giraffe_protected_area2.csv",header=TRUE)
pa$pa_type = ifelse(pa$DESIG =="None","None","ProtectedArea")

landuse = as.data.frame(pa%>%
                          dplyr::group_by(ID,pa_type)%>%
                          dplyr::summarise(use =n()))
obs = as.data.frame(pa %>%
                      dplyr::group_by(ID)%>%
                      dplyr::summarise(obs=n()))                    
landuse = left_join(landuse,obs,by="ID")
landuse$percent = landuse$use/landuse$obs

library()
pa_summary = landuse  %>% 
  spread (pa_type,percent)%>%
  tidyr::drop_na(ProtectedArea)%>%
  dplyr::select(ID,ProtectedArea)
  
pa_summary = left_join(collar_key,pa_summary,by="ID") %>%
  dplyr::select(ID,ProtectedArea)
pa_summary$ProtectedArea[is.na(pa_summary$ProtectedArea)] <- 0

#summarise landcover results for covariate input
globcover_2009 = read.csv("globcover_2009.csv",header=TRUE)
# globcover_key = read.csv("C:/Users/mbbro/Dropbox/Projects_Working/TwigaTracker_Layers/Globcover2009_Legend.csv", header = T)
# glob_key=globcover_key %>%
#   dplyr::select(raster,Label2)
# globcover_2009 =left_join(globcover_2009,glob_key,by="raster")

globcover_summary = as.data.frame(globcover_2009%>%
                                    dplyr::group_by(ID,raster,Label2)%>%
                                    dplyr::summarise(use =n()))
obs = as.data.frame(globcover_2009 %>%
                      dplyr::group_by(ID)%>%
                      dplyr::summarise(obs=n()))                    
globcover_summary = left_join(globcover_summary,obs,by="ID")
globcover_summary$percent = globcover_summary$use/globcover_summary$obs
globcover_summary$raster = as.character(globcover_summary$raster)
globcover_summary = globcover_summary %>%
  mutate(raster2 = paste0("raster_", raster))

globcover_summary = globcover_summary%>%
  dplyr::select(ID,Label2,percent)%>%
  tidyr::spread(Label2,percent,convert=TRUE)
globcover_summary[is.na(globcover_summary)] <- 0

#Find most Frequently used habitat type for each individual
globcover_summary2 = globcover_summary %>%
  dplyr::select(!ID)
habitat = colnames(globcover_summary2)[max.col(globcover_summary2, ties.method = "first")]  # Apply colnames & max.col functions
globcover_summary$habitat = habitat


#Combine all dataframes into a master dataframe for analyses
continental_df = response_df %>%
  dplyr::left_join(pa_summary,by="ID")%>%
  dplyr::left_join(temp_summary,by="ID")%>%
  dplyr::left_join(ndvi_summary,by="ID")%>%
  dplyr::left_join(globcover_summary,by="ID")%>%
  dplyr::left_join(hfp2009_summary, by = "ID")%>%
  dplyr::left_join(human_density_summary, by = "ID")%>%
  dplyr::left_join(woody_summary, by = "ID")%>%
  dplyr::distinct(ID,.keep_all=TRUE)

#Export the dataframe for subsequent modeling
#write.csv(continental_df, file="continental_df4.csv")

###################EX POST FACTO COVARIATE ADDITIONs
#Read in the Continental DF
continental_df4 = read.csv("continental_df4.csv",header=TRUE)
cop_bare = read.csv("cop_bare.csv", header = T)
cop_bare_summary = cop_bare %>%
  dplyr::group_by(ID)%>%
  dplyr::summarise(cop_bare_mean = mean(raster_cop_bare, na.rm=TRUE),
                   cop_bare_sd = sd(raster_cop_bare,na.rm = TRUE))
cop_bare_summary$collarid = cop_bare_summary$ID
cop_bare_summary = as.data.frame(cop_bare_summary %>%
  dplyr::select(collarid,cop_bare_mean, cop_bare_sd))


cop_tree = read.csv("cop_tree.csv", header = T)
cop_tree_summary = cop_tree %>%
  dplyr::group_by(ID)%>%
  dplyr::summarise(cop_tree_mean = mean(raster_cop_tree, na.rm=TRUE),
                   cop_tree_sd = sd(raster_cop_tree,na.rm = TRUE))
cop_tree_summary$collarid = cop_tree_summary$ID
cop_tree_summary = as.data.frame(cop_tree_summary %>%
                                   dplyr::select(collarid,cop_tree_mean, cop_tree_sd))


cop_shrub = read.csv("cop_shrub.csv", header = T)
cop_shrub_summary = cop_shrub %>%
  dplyr::group_by(ID)%>%
  dplyr::summarise(cop_shrub_mean = mean(raster_cop_shrub, na.rm=TRUE),
                   cop_shrub_sd = sd(raster_cop_shrub,na.rm = TRUE))
cop_shrub_summary$collarid = cop_shrub_summary$ID
cop_shrub_summary = as.data.frame(cop_shrub_summary %>%
                                   dplyr::select(collarid,cop_shrub_mean, cop_shrub_sd))

#Left Join to Continental Dataframe
continental_df5 = left_join(continental_df4,cop_bare_summary, by="collarid")
continental_df5 = left_join(continental_df5,cop_tree_summary, by="collarid")
continental_df5 = left_join(continental_df5,cop_shrub_summary, by="collarid")

#write.csv(continental_df5, file="continental_df5.csv")

