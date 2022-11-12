##THE CONTINENTAL
##Mixed Effects Models in the Meta Analysis Framework
#Author: Michael Butler Brown
#Author Affiliations: Smithsonian Conservation Biology Institute and Giraffe Conservation Foundation

#Start Fresh
rm(list=ls())

#Load Relevant packages
library(plyr)
library(dplyr)
library(tidyverse)
library(metafor)
library(forestplot)

#Set Working Directory
#Change working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in Movement Data
df = read.csv(file="../Data/continental_df.csv", header = TRUE, sep =",")
collar_key = read.csv(file="../Data/ctmm_key.csv",head=TRUE,sep=",")

#Changing some sites to northern Kenya
df$Park_location = replace(df$Park_location, df$Park_location =="Biliqo_Bulesa","NorthernKenya")#For replacing characters with other characters
df$Park_location = replace(df$Park_location, df$Park_location =="Melako","NorthernKenya")#For replacing characters with other characters

#Transforming some of the variables for gaussian modeling
df$log_akde_95 = log(df$akde_95_ML) #log transformation fo r
df$log_akde_50 = log(df$akde_50_ML)
df$var_log_area_95 = log(df$akde_95_ML)^2/df$akde_95_DOF_area
df$var_log_area_50 = log(df$akde_50_ML)^2/df$akde_50_DOF_area

#Some Simple filtering to exclude giraffe with bad models
df2 = df %>%
  dplyr::filter(!Variogram_inspection =="bad")%>% #If variograms didn't converge on home range behaviour as determined by subjective visual inspection
  dplyr::filter(!Number_of_fixes < 500)%>% # If there are fewer than 500 fixes...I didn't run ctmm models on individuals with fewer than 500 fixes, so there shouldn't be too many of these
  dplyr::filter(!DOF_area < 4)%>% #If DOF area is less than 4
  dplyr::filter(!akde_95_ML > 20000)%>% #Just to get rid of that one weird one up in Niger
  drop_na(model) #Get rid of all collars without sufficient model data

#Generating a quick plot to visualize the HR size by site
ggplot(data=df2,aes(x=akde_95_ML, y= ID))+
  geom_point()+
  geom_errorbar(aes(xmin=akde_95_lcl,xmax=akde_95_ucl))+
  xlab("AKDE 95% Size")+
  ylab(" Collar ID")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size")

#Generating a forest Plot to look at HR Variation Across Sites
summary_site = df2 %>%
  dplyr::group_by(Park_location)%>%
  dplyr::summarise(mean_akde = mean(akde_95_ML,na.rm=TRUE),
            sd.akde =sd(akde_95_ML, na.rm=TRUE),
            n.akde = n())%>%
  dplyr::mutate(se.akde = sd.akde/sqrt(n.akde),
         lower.ci.akde = mean_akde - qt(1-(0.05/2),n.akde-1)*se.akde,
         upper.ci.akde = mean_akde + qt(1-(0.05/2),n.akde-1)*se.akde)
         
#Generating a quick plot to visualize the HR size by site
ggplot(data=summary_site,aes(x=mean_akde, y= Park_location))+
  geom_point()+
  geom_errorbar(aes(xmin=lower.ci.akde,xmax=upper.ci.akde))+
  xlab("AKDE 95% Size")+
  ylab("Study Site")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size by Site")

#Generating a forest Plot to look at HR Variation Across Countries
summary_country = df2 %>%
  dplyr::group_by(Country)%>%
  dplyr::summarise(mean_akde = mean(akde_95_ML,na.rm=TRUE),
                   sd.akde =sd(akde_95_ML, na.rm=TRUE),
                   n.akde = n())%>%
  dplyr::mutate(se.akde = sd.akde/sqrt(n.akde),
                lower.ci.akde = mean_akde - qt(1-(0.05/2),n.akde-1)*se.akde,
                upper.ci.akde = mean_akde + qt(1-(0.05/2),n.akde-1)*se.akde)

summary_country = summary_country %>%
  drop_na(sd.akde)

#Generating a quick plot to visualize the HR size by Country
ggplot(data=summary_country,aes(x=mean_akde, y= Country))+
  geom_point()+
  geom_errorbar(aes(xmin=lower.ci.akde,xmax=upper.ci.akde))+
  xlab("AKDE 95% Size")+
  ylab("Country")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("Country")

summary_species = df2 %>%
  dplyr::group_by(species)%>%
  dplyr::summarise(mean_akde = mean(akde_95_ML,na.rm=TRUE),
                   sd.akde =sd(akde_95_ML, na.rm=TRUE),
                   n.akde = n())%>%
  dplyr::mutate(se.akde = sd.akde/sqrt(n.akde),
                lower.ci.akde = mean_akde - qt(1-(0.05/2),n.akde-1)*se.akde,
                upper.ci.akde = mean_akde + qt(1-(0.05/2),n.akde-1)*se.akde)

summary_species = summary_species %>%
  drop_na(sd.akde)

#Generating a quick plot to visualize the HR size by Species
ggplot(data=summary_species,aes(x=mean_akde, y= species))+
  geom_point()+
  geom_errorbar(aes(xmin=lower.ci.akde,xmax=upper.ci.akde))+
  xlab("AKDE 95% Size")+
  ylab("Species")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE Variation by Species")

#Plotting some simple linear relationships
#Plot out some relationships
# #Variation by Taxa
ggplot(data=df2,aes(x= ID, y= akde_95_ML , color = taxa))+
  geom_point()+
  geom_errorbar(aes(ymin=akde_95_lcl,ymax= akde_95_ucl))+
  xlab("GiraffeID")+
  ylab("AKDE area in km^2")+
  facet_grid(col=vars(taxa),scales="free_x")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size by Taxa")

# #Variation by Taxa 2
# akde_species = as.data.frame(df %>%
#                                dplyr::group_by(taxa) %>%
#                                dplyr::summarise (akde_mean=mean(ML_95),
#                                                  akde_sd =sd(ML_95)))
# ggplot(data=akde_species,aes(x=taxa, y= akde_mean))+
#   geom_point()+
#   geom_errorbar(aes(ymin=akde_mean-akde_sd,ymax=akde_mean+akde_sd))+
#   #geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
#   xlab("Taxa")+
#   ylab("AKDE area in km^2")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle("AKDE_95 Size by Taxa")
# 
# #Variation by percent protected area
# ggplot(data=df,aes(x= ProtectedArea, y= ML_95))+
#   geom_point()+
#   geom_errorbar(aes(ymin=lcl_95,ymax=ucl_95))+
#   geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
#   xlab("Percent of Time in Protected Area")+
#   ylab("AKDE area in km^2")+
#   #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle("AKDE_95 Size by Time in Protected Area")
# 
#Variation by Mean NDVI
ggplot(data=df2,aes(x= ndvi_mean, y= log_akde_95))+
  geom_point()+
  geom_errorbar(aes(ymin=(log_akde_95-var_log_area_95),ymax=(log_akde_95+var_log_area_95)))+
  geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
  xlab("NDVI")+
  ylab("Log AKDE area in km^2")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size by NDVI")
# 
# #Variation by SD NDVI
 ggplot(data=df2,aes(x= ndvi_sd, y= log_akde_95))+
  geom_point()+
  geom_errorbar(aes(ymin=(log_akde_95-var_log_area_95),ymax=(log_akde_95+var_log_area_95)))+
  geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
  xlab("NDVI SD")+
  ylab("Log AKDE area in km^2")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size by NDVI SD")
# 
#Variation by Mean NDVI by Country
ggplot(data=df2,aes(x= ndvi_mean, y= log_akde_95, color = Country))+
  geom_point()+
  geom_errorbar(aes(ymin=(log_akde_95-var_log_area_95),ymax=(log_akde_95+var_log_area_95)))+
  xlab("Mean NDVI")+
  ylab("Log AKDE area")+
  facet_grid(col=vars(Country),scales="free_x")+
  geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
  ggtitle("AKDE_95 Size by Country")

#Variation by Mean NDVI by Habitat Type
#summarize the number of measurements per dominant habitat type
habitat_summary = as.data.frame(df2 %>%
                                  drop_na(ndvi_mean)%>% #Get rid of all collars without sufficient ndvi data
                                  dplyr::group_by(habitat)%>%
                                  dplyr::summarise(total = n())%>%
                                  dplyr::filter(total>3))
habitat_filter = as.vector(habitat_summary$habitat) #create a vecotr to only filter habitats with >3 measurements
df3 = df2 %>%
  dplyr::filter(habitat %in% habitat_filter)

ggplot(data=df3,aes(x= ndvi_mean, y= log_akde_95, color = habitat))+
  geom_point()+
  geom_errorbar(aes(ymin=(log_akde_95-var_log_area_95),ymax=(log_akde_95+var_log_area_95)))+
  xlab("Mean NDVI")+
  ylab("Log 95% AKDE area")+
  facet_grid(col=vars(habitat),scales="free_x")+
  geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
  ggtitle("AKDE_95 Size by Habitat Type")
# 
# #Variation by STDEV NDVI
ggplot(data=df2,aes(x=ndvi_sd, y= log_akde_95))+
  geom_point()+
  geom_errorbar(aes(ymin=(log_akde_95-var_log_area_95),ymax=(log_akde_95+var_log_area_95)))+
  geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
  xlab("NDVI Standard Deviation")+
  ylab("Log 95% AKDE area")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size by NDVI STDEV")
# 

# #Variation by Protected Area
ggplot(data=df2,aes(x=ProtectedArea, y= log_akde_95))+
  geom_point()+
  geom_errorbar(aes(ymin=(log_akde_95-var_log_area_95),ymax=(log_akde_95+var_log_area_95)))+
  geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
  xlab("NDVI Standard Deviation")+
  ylab("Log 95 % AKDE area")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size by % Time in Protected Area")
# 

# #Variation by Density
df3=df2 %>%
  dplyr::filter(!human_density_mean >200)
ggplot(data=df3,aes(x=human_density_mean, y= log_akde_95))+
  geom_point()+
  geom_errorbar(aes(ymin=(log_akde_95-var_log_area_95),ymax=(log_akde_95+var_log_area_95)))+
  geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
  xlab("Mean Human Density")+
  ylab("Log 95 % AKDE area")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size by Mean Human Density")
# 

# #Variation by Density
ggplot(data=df3,aes(x=hfp2009_mean, y= log_akde_95))+
  geom_point()+
  geom_errorbar(aes(ymin=(log_akde_95-var_log_area_95),ymax=(log_akde_95+var_log_area_95)))+
  geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
  xlab("Mean Human Footprint Index")+
  ylab("Log 95 % AKDE area")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size by Mean Human Footprint Index")

# #Variation by Taxa 
# akde_species = as.data.frame(df %>%
#                                dplyr::group_by(taxa) %>%
#                                dplyr::summarise (akde_mean=mean(ML_95),
#                                                   akde_sd =sd(ML_95)))
# ggplot(data=akde_species,aes(x=taxa, y= akde_mean))+
#   geom_point()+
#   geom_errorbar(aes(ymin=akde_mean-akde_sd,ymax=akde_mean+akde_sd))+
#   #geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
#   xlab("Taxa")+
#   ylab("AKDE area in km^2")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle("AKDE_95 Size by Taxa")
# 
# #Variation by Habitat Type 
# akde_habitat = as.data.frame(df %>%
#                                dplyr::group_by(habitat) %>%
#                                dplyr::summarise (akde_mean=mean(log_ML_95),
#                                                  akde_sd =sd(log_ML_95)))
# 
# ggplot(data=akde_habitat,aes(x=habitat, y= akde_mean))+
#   geom_point()+
#   geom_errorbar(aes(ymin=akde_mean-akde_sd,ymax=akde_mean+akde_sd))+
#   #geom_smooth(method="lm", se=TRUE, formula= y ~ x)+
#   xlab("Taxa")+
#   ylab("AKDE area in km^2")+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   ggtitle("AKDE_95 Size by Habitat Type")

############################################
############################################
############################################
#Looking at some quick forest plots per species to compare HR estimates
head(df2)
dat = df2
mean(dat$akde_95_ML)
dat$akde_95_se = (dat$akde_95_ucl-dat$akde_95_lcl)/3.92

dat2 = dat %>%
  dplyr::select(collarid, unit_ID, species, subspecies,taxa,Country,Park_location,sex,DOF_area, akde_95_lcl, akde_95_ML,akde_95_ucl,akde_95_DOF_area,akde_95_se, log_akde_95,var_log_area_95)
dat_reticulata = dat2 %>%
  dplyr::filter(species == "reticulata")
dat_reticulata = dat2 %>%
  dplyr::filter(species =="reticulata")
dat_camelopardalis = dat2 %>%
  dplyr::filter(species =="camelopardalis")
dat_tippelskirchi = dat2 %>%
  dplyr::filter(species =="tippelskirchi")

#reticulata
#Plot out the raw data
ForestPlot_reticulata=
  ggplot(data=dat_reticulata, aes(x=akde_95_ML,y=collarid))+
  geom_point(aes(color=collarid))+
  geom_errorbar(aes(xmin=(akde_95_lcl), xmax=(akde_95_ucl), color=collarid))+
  theme(legend.position="none")+
  geom_vline(xintercept=247.89)+
  geom_vline(xintercept=271.6,col="red")
ForestPlot_reticulata

#Trying log transformed estimates
res_reticulata <- rma(yi=log_akde_95, vi=var_log_area_95, data=dat_reticulata)
res_reticulata
forest(res_reticulata)
predict(res_reticulata,transf=exp,digits=2)

#Trying with raw estimates
res2_reticulata <- rma(yi=akde_95_ML, sei = akde_95_se, data=dat_reticulata)
res2_reticulata
forest(res2_reticulata)

#Calculating the mean without any meta analyses
mean(dat_reticulata$akde_95_ML)
dat$akde_95_se = (dat$akde_95_ucl-dat$akde_95_lcl)/3.92

#giraffa
#Plot out the raw data
ForestPlot_giraffa=
  ggplot(data=dat_giraffa, aes(x=akde_95_ML,y=collarid))+
  geom_point(aes(color=collarid))+
  geom_errorbar(aes(xmin=(akde_95_lcl), xmax=(akde_95_ucl), color=collarid))+
  theme(legend.position="none")+
  geom_vline(xintercept=372.63)+
  geom_vline(xintercept=692,col="red")
ForestPlot_giraffa

#Trying log transformed estimates
res_giraffa <- rma(yi=log_akde_95, vi=var_log_area_95, data=dat_giraffa)
res_giraffa
forest(res_giraffa)
predict(res_giraffa,transf=exp,digits=2)

#Trying with raw estimates
res2_giraffa <- rma(yi=akde_95_ML, sei = akde_95_se, data=dat_giraffa)
res2_giraffa
forest(res2_giraffa)

#Calculating the mean without any meta analyses
mean(dat_giraffa$akde_95_ML)
dat$akde_95_se = (dat$akde_95_ucl-dat$akde_95_lcl)/3.92

#tippelskirchi
#Plot out the raw data
ForestPlot_tippelskirchi=
  ggplot(data=dat_tippelskirchi, aes(x=akde_95_ML,y=collarid))+
  geom_point(aes(color=collarid))+
  geom_errorbar(aes(xmin=(akde_95_lcl), xmax=(akde_95_ucl), color=collarid))+
  theme(legend.position="none")+
  geom_vline(xintercept= 77.61)+
  geom_vline(xintercept=174.99,col="red")
ForestPlot_tippelskirchi

#Trying log transformed estimates
res_tippelskirchi <- rma(yi=log_akde_95, vi=var_log_area_95, data=dat_tippelskirchi)
res_tippelskirchi
forest(res_tippelskirchi)
predict(res_tippelskirchi,transf=exp,digits=2)

#Trying with raw estimates
res2_tippelskirchi <- rma(yi=akde_95_ML, sei = akde_95_se, data=dat_tippelskirchi)
res2_tippelskirchi
forest(res2_tippelskirchi)

#Calculating the mean without any meta analyses
mean(dat_tippelskirchi$akde_95_ML)
dat$akde_95_se = (dat$akde_95_ucl-dat$akde_95_lcl)/3.92


#camelopardalis
#Plot out the raw data
ForestPlot_camelopardalis=
  ggplot(data=dat_camelopardalis, aes(x=akde_95_ML,y=collarid))+
  geom_point(aes(color=collarid))+
  geom_errorbar(aes(xmin=(akde_95_lcl), xmax=(akde_95_ucl), color=collarid))+
  theme(legend.position="none")+
  geom_vline(xintercept= 414.68)+
  geom_vline(xintercept=585.11,col="red")
ForestPlot_camelopardalis

#Trying log transformed estimates
res_camelopardalis <- rma(yi=log_akde_95, vi=var_log_area_95, data=dat_camelopardalis)
res_camelopardalis
forest(res_camelopardalis)
predict(res_camelopardalis,transf=exp,digits=2)

#Trying with raw estimates
res2_camelopardalis <- rma(yi=akde_95_ML, sei = akde_95_se, data=dat_camelopardalis)
res2_camelopardalis
forest(res2_camelopardalis)

#Calculating the mean without any meta analyses
mean(dat_camelopardalis$akde_95_ML)
dat$akde_95_se = (dat$akde_95_ucl-dat$akde_95_lcl)/3.92


#COMBINING ALL OF THE SUMMARIES FROM RAW META ANALYSES 
species = c("giraffa","reticulata","tippelskirchi","camelopardalis")
est = c(372.63, 247.89,77.61,414.68)
se = c(37.79,23.12,13.76,33.88)
lcl=c(298.55,202.58,50.65,348.29)
ucl =c(446.71,293.19,104.57,481.08)
summary = as.data.frame(cbind(species,est,se,lcl,ucl))
summary$est = as.numeric(summary$est)
summary$se = as.numeric(summary$se)
summary$ucl = as.numeric(summary$ucl)
summary$lcl = as.numeric(summary$lcl)

#Generating a quick plot to visualize the HR size by site
ggplot(data=summary,aes(x=est, y= species))+
  geom_point()+
  geom_errorbar(aes(xmin=lcl,xmax=ucl))+
  xlab("AKDE 95% Size")+
  ylab("Species")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE Variation by Species")

#Generating a quick plot to visualize the HR size by site
ggplot(data=summary,aes(x=species, y=est))+
  geom_point()+
  geom_errorbar(aes(ymin=lcl,ymax=ucl))+
  ylab("AKDE 95% Size")+
  xlab("Species")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE Variation by Species")


#Calculating Mean using meta::metafor()
### fit fixed model
library(metafor)
res <- rma(yi=log_95, vi=var_log_95, data=dat)
res
#predited effects 
predict(res,transf=exp,digits=2)

res2 = rma(yi = akde_95_ML, sei=akde_95_se, data=dat)
res2

forest(res2)


############################################
############################################
############################################
library(lme4)
require(AICcmodavg)
require(kimisc)

#Generating some models
#Simple linear regressions
m1 = glm(df2$log_akde_95 ~ df2$ndvi_mean) 
m2 = glm(df2$log_akde_95 ~ df2$ndvi_sd) 
m3 = glm(df2$log_akde_95 ~ df2$ProtectedArea) 
m4 = glm(df2$log_akde_95 ~ df2$habitat) 
m5 = glm(df2$log_akde_95 ~ df2$ndvi_mean) 
m6 = glm(df2$log_akde_95 ~ df2$temp_mean) 
m7 = glm(df2$log_akde_95 ~ df2$Park_location) 
m9 = glm(df2$log_akde_95 ~ df2$sex) 
m10 = glm(df2$log_akde_95 ~ df2$species) 
m11= glm(df2$log_akde_95 ~ df2$taxa) 
m12 = glm(df2$log_akde_95 ~ df2$subspecies) 
m13 = glm(df2$log_akde_95 ~ df2$Country ) 
m14 = glm(df2$log_akde_95 ~ df2$human_density_mean)
m16 = glm(df2$log_akde_95 ~ df2$human_density_sd)
m17 = glm(df2$log_akde_95 ~ df2$hfp2009_mean)
m18 = glm(df2$log_akde_95 ~ df2$hfp2009_sd)
m19 = glm(df2$log_akde_95 ~ df2$human_density_sd)

#With Interactions
m20 = glm()

models_lm = nlist(m1,m2,m3,m4,m5,m6,m7,m9,m10,m11,m12,m13,m14,m16,m17,m18,m19)
aic_table= aictab(models_lm)

#Trying Global Models
require(MuMIn)
df3 = df2 %>%
  dplyr::select(log_akde_95,var_log_area_95,ndvi_mean,ndvi_sd,ProtectedArea,temp_mean, temp_sd, habitat, hfp2009_mean, hfp2009_sd, human_density_mean, human_density_sd,species,Park_location) %>%
  dplyr::filter(!is.na(log_akde_95))%>%
  dplyr::filter(!is.na(temp_mean))%>%
  dplyr::filter(!is.na(ProtectedArea))%>%
  dplyr::filter(!is.na(ndvi_mean))%>%
  dplyr::filter(!is.na(ndvi_sd))%>%
  dplyr::filter(!is.na(habitat))%>%
  dplyr::filter(!is.na(hfp2009_mean))%>%
  dplyr::filter(!is.na(human_density_mean))

#Find rows with na so you can remove them if neccesary
new_DF <- df3[rowSums(is.na(df3)) > 0,]
new_DF

#Create a Global model and run all combinations of covariates
globalmodel1 = lm(log_akde_95 ~ ndvi_mean + ndvi_sd + ProtectedArea + temp_mean + temp_sd + habitat + hfp2009_mean + hfp2009_sd + human_density_mean + human_density_sd, data=df3, na.action = na.fail) 
combinations = dredge(globalmodel1)
print(combinations)

globalmodel2 = lm(log_akde_95 ~ ndvi_mean + ndvi_sd + ProtectedArea + temp_mean + temp_sd + habitat + hfp2009_mean + hfp2009_sd + human_density_mean + human_density_sd + ndvi_mean * habitat + ndvi_mean * human_density_mean + ndvi_sd * habitat + ndvi_sd * human_density_mean + ndvi_mean * ProtectedArea + ndvi_sd * ProtectedArea, data=df3, na.action = na.fail) 
combinations2 = dredge(globalmodel2)
print(combinations2)

#'Best' model
summary(get.models(combinations2, 1)[[1]])

#Look at 50
globalmodel3 = lm(log_akde_50 ~ ndvi_mean + ndvi_sd + ProtectedArea + temp_mean + temp_sd + habitat + hfp2009_mean + hfp2009_sd + human_density_mean + human_density_sd + ndvi_mean * habitat + ndvi_mean * human_density_mean + ndvi_sd * habitat + ndvi_sd * human_density_mean + ndvi_mean * ProtectedArea + ndvi_sd * ProtectedArea, data=df, na.action = na.fail) 
combinations3 = dredge(globalmodel3)
print(combinations3)

#globalmodel3 = lm(log_akde_95 ~ species + ndvi_mean + ndvi_sd + ProtectedArea + temp_mean + temp_sd + habitat + hfp2009_mean + hfp2009_sd + human_density_mean + human_density_sd + ndvi_mean * habitat + ndvi_mean * human_density_mean + ndvi_sd * habitat + ndvi_sd * human_density_mean + ndvi_mean * ProtectedArea + ndvi_sd * ProtectedArea + species * ndvi_mean + species * ndvi_sd + species * hfp2009_mean + species * ProtectedArea, data=df3, na.action = na.fail) 
#combinations3 = dredge(globalmodel3)
#print(combinations3)


###Looking at Mixed Effects Models (without the meta framework)
#Creating a mixed effects model to look at the effect of time on lesion size
library(lme4)
glmm_1 = lmer(log_akde_95 ~ habitat +  hfp2009_mean + hfp2009_sd + 
                human_density_mean + ndvi_mean + ProtectedArea + temp_mean + 
                habitat:ndvi_mean + ndvi_mean:ProtectedArea + (1|Park_location),data = df3, REML = FALSE, na.action = 'na.exclude')

#Report the model summary
library(report)
glmm_results <- report(glmm_1)
glmm_results
summary(glmm_1)

#Report the model summary
library(report)
growth_model_results <- report(growth_model1)
growth_model_results

#Attempt at a mixed effects meta analysis
or.rma.mixed = rma(log_akde_95,var_log_area_95, mods = ~ habitat +  hfp2009_mean + hfp2009_sd + 
                     human_density_mean + ndvi_mean + ProtectedArea + temp_mean + 
                     habitat:ndvi_mean + ndvi_mean:ProtectedArea,
                     method= "REML",
                   data=df3)