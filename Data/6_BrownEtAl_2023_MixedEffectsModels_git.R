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
df = read.csv(file="../Data/continental_df5.csv", header = TRUE, sep =",")
collar_key = read.csv(file="../Data/ctmm_key.csv",head=TRUE,sep=",")

#Changing some sites to northern Kenya
df$Park_location = replace(df$Park_location, df$Park_location =="Biliqo_Bulesa","NorthernKenya")#For replacing characters with other characters
df$Park_location = replace(df$Park_location, df$Park_location =="Melako","NorthernKenya")#For replacing characters with other characters

#Transforming some of the variables for gaussian modeling
df$log_akde_95 = log(df$akde_95_ML) #log transformation for the 95%akde 
df$log_akde_50 = log(df$akde_50_ML) # log transformation for the 50% akde
df$log_speed = log(df$speed_ml)
df$var_log_area_95 = log(df$akde_95_ML)^2/df$akde_95_DOF_area #Proper variance for subsequent metafor analyses
df$var_log_area_50 = log(df$akde_50_ML)^2/df$akde_50_DOF_area #Proper variance for subsewuent metafor analyses
df$var_log_speed = log(df$speed_ml)^2/df$DOF_speed #Proper variance for subsewuent metafor analyses


#Some Simple filtering to exclude giraffe with bad models
df2 = df %>%
  dplyr::filter(!Number_of_fixes < 500)%>% # If there are fewer than 500 fixes...I didn't run ctmm models on individuals with fewer than 500 fixes, so there shouldn't be too many of these
  tidyr::drop_na(model)%>%
  drop_na(model)%>% #Get rid of all collars without sufficient model data
  dplyr::filter(!collarid =="IRI2016-3241")%>% #Just to get rid of that one weird one up in Niger
  dplyr::filter(!DOF_area < 4)%>% #If DOF area is less than 4
  dplyr::filter(!Variogram_inspection =="bad") #If variograms didn't converge on home range behaviour as determined by subjective visual inspection

#How many giraffe had home range residency
dim(df)
dim(df2)

#Generating a quick plot to visualize the HR size by Species
ggplot(data=df2,aes(x=akde_95_ML, y= ID))+
  geom_point(aes(color=species))+
  geom_errorbar(aes(xmin=akde_95_lcl,xmax=akde_95_ucl,color=species))+
  xlab("AKDE 95% Size")+
  ylab(" Collar ID")+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size")

#Creating estimate for tau_pos_ml
df2$tau_pos_ml = ((df2$tau_pos_lcl+ df2$tau_pos_ucl)/2)

#Calculating Overall summaries
summary_overall = df2 %>%
  dplyr::summarise(akde_mean = mean(akde_95_ML, na.rm=TRUE),
                   akde_sd = sd(akde_95_ML,na.rm=TRUE),
                   tau_pos_mean = mean(tau_pos_ml, na.rm=TRUE),
                   tau_pos_sd = sd(tau_pos_ml, na.rm=TRUE),
                   speed_mean = mean(speed_ml, na.rm=TRUE),
                   speed_sd = sd(speed_ml, na.rm=TRUE),
                   tau_vel_mean = mean(tau_vel_ml,na.rm=TRUE),
                   tau_vel_sd = sd(tau_vel_ml,na.rm=TRUE))
summary_overall

summary_species = as.data.frame(df2 %>%
  dplyr::group_by(species)%>%
  dplyr::summarise(total=n(),
    akde_mean = mean(akde_95_ML, na.rm=TRUE),
                   akde_sd = sd(akde_95_ML,na.rm=TRUE),
                   tau_pos_mean = mean(tau_pos_ml, na.rm=TRUE),
                   tau_pos_sd = sd(tau_pos_ml, na.rm=TRUE),
                   speed_mean = mean(speed_ml, na.rm=TRUE),
                   speed_sd = sd(speed_ml, na.rm=TRUE),
                   tau_vel_mean = mean(tau_vel_ml,na.rm=TRUE),
                   tau_vel_sd = sd(tau_vel_ml,na.rm=TRUE)))
summary_species

############################################
#Looking at some quick forest plots per species to compare HR estimates
head(df2)
dat = df2
mean(dat$akde_95_ML, na.rm=TRUE)
sd(dat$akde_95_ML, na.rm=TRUE)

dat$akde_95_se = (dat$akde_95_ucl-dat$akde_95_lcl)/3.92
dat$speed_se = (dat$speed_ucl-dat$speed_lcl)/3.92

# dat2 = dat %>%
#   dplyr::select(collarid, unit_ID, species, subspecies,taxa,Country,Park_location,sex,DOF_area, akde_95_lcl, akde_95_ML,akde_95_ucl,akde_95_DOF_area,akde_95_se, log_akde_95,var_log_area_95)
dat2=dat

dat_reticulata = dat2 %>%
  dplyr::filter(species == "reticulata")
dat_giraffa = dat2 %>%
  dplyr::filter(species =="giraffa")
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
# res_reticulata <- rma(yi=log_akde_95, vi=var_log_area_95, data=dat_reticulata)
# res_reticulata
# forest(res_reticulata)
# predict(res_reticulata,transf=exp,digits=2)

#Trying with raw estimates
res2_reticulata <- rma(yi=akde_95_ML, sei = akde_95_se, data=dat_reticulata)
res2_reticulata
forest(res2_reticulata)

#Trying with raw estimates
res_speed_reticulata <- rma(yi=speed_ml, sei = speed_se, data=dat_reticulata)
res_speed_reticulata
forest(res_speed_reticulata)

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

#Trying with raw estimates
res2_reticulata <- rma(yi=akde_95_ML, sei = akde_95_se, data=dat_reticulata)
res2_reticulata
forest(res2_reticulata)

#Trying with raw estimates
res_speed_reticulata <- rma(yi=speed_ml, sei = speed_se, data=dat_reticulata)
res_speed_reticulata
forest(res_speed_reticulata)

#Trying with raw estimates
res2_giraffa <- rma(yi=akde_95_ML, sei = akde_95_se, data=dat_giraffa)
res2_giraffa
forest(res2_giraffa)

#Trying with raw estimates
res_speed_giraffa <- rma(yi=speed_ml, sei = speed_se, data=dat_giraffa)
res_speed_giraffa
forest(res_speed_giraffa)

#Trying with raw estimates
res2_tippelskirchi <- rma(yi=akde_95_ML, sei = akde_95_se, data=dat_tippelskirchi)
res2_tippelskirchi
forest(res2_tippelskirchi)

#Trying with raw estimates
res_speed_tippelskirchi <- rma(yi=speed_ml, sei = speed_se, data=dat_tippelskirchi)
res_speed_tippelskirchi
forest(res_speed_tippelskirchi)

#Trying with raw estimates
res2_camelopardalis <- rma(yi=akde_95_ML, sei = akde_95_se, data=dat_camelopardalis)
res2_camelopardalis
forest(res2_camelopardalis)

#Trying with raw estimates
res_speed_camelopardalis <- rma(yi=speed_ml, sei = speed_se, data=dat_camelopardalis)
res_speed_camelopardalis
forest(res_speed_camelopardalis)


#Calculating Mean using meta::metafor()
### fit fixed model
# library(metafor)
# res <- rma(yi=log_akde_95, vi=var_log_area_95, mods = ~species, data=dat)
# res

#predicted effects 
#predict(res,transf=exp,digits=2)
res2 = rma(yi = akde_95_ML, sei=akde_95_se, mods = ~ species, data=dat2)
res2
forest(res2)

res_speed = rma(yi = speed_ml, sei=speed_se, mods = ~ species, data=dat2)
res_speed
forest(res_speed)

#Overall
res_overall = rma(yi = akde_95_ML, sei=akde_95_se, data=dat2)
res_overall
forest(res_overall)

res_speed_overall = rma(yi = speed_ml, sei=speed_se, data=dat2)
res_speed_overall
forest(res_speed_overall)


# 
# res3 = rma(log_akde_95,var_log_area_95,mods = ~ species, data=dat)
# res3
# forest(res3)
# results = as.data.frame(predict(res3,transf=exp,digits=2))
# results = results %>%
#   dplyr::distinct(pred, .keep_all = TRUE)

#COMBINING ALL OF THE SUMMARIES FROM RAW META ANALYSES 
species = c("G. reticulata","G. giraffa","G. tippelskirchi","G. camelopardalis")
est = c(262.1,360.8,112.89,380.01)
se = c(21.95,37.93,19.14,30.96)
lcl=c(219.1,286.5,75.37,319.32)
ucl =c(305.1,435.2,150.43,440.70)
# est = c(377.5272, 333.3975,288.5072,196.3472)
# se = c(30.18,42.5296,55.2618,64.446)
# lcl=c(318.3755,250.0409,180.1955,69.6408)
# ucl =c(436.679,416.754,396.8177,323.043)
summary = as.data.frame(cbind(species,est,se,lcl,ucl))
summary$est = as.numeric(summary$est)
summary$se = as.numeric(summary$se)
summary$ucl = as.numeric(summary$ucl)
summary$lcl = as.numeric(summary$lcl)

#Generating a quick plot to visualize the HR size by site
ggplot(data=summary,aes(x=species, y=est))+
  geom_errorbar(aes(ymin=lcl,ymax=ucl),size=.8,width=.5, colour="gray")+
  geom_point(aes(),size=5)+
  ylab(bquote('AKDE 95% Size'~km^2))+
  xlab("Species")+
  theme_bw()+
  theme(text = element_text(size=rel(5)),
        strip.text.x = element_text(size=rel(5)),
        strip.text.y = element_text(size=rel(5)))
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  #ggtitle("AKDE Variation by Species")


#COMBINING ALL OF THE SUMMARIES FROM RAW META ANALYSES 
species = c("G. reticulata","G. giraffa","G. tippelskirchi","G. camelopardalis")
est = c(13.65,13.55,11.90,16.01)
se = c(.27,.42,.57,41)
lcl=c(13.10,12.71,10.77,15.2)
ucl =c(14.2,14.39,13.02,16.82)
# est = c(377.5272, 333.3975,288.5072,196.3472)
# se = c(30.18,42.5296,55.2618,64.446)
# lcl=c(318.3755,250.0409,180.1955,69.6408)
# ucl =c(436.679,416.754,396.8177,323.043)
summary_speed = as.data.frame(cbind(species,est,se,lcl,ucl))
summary_speed$est = as.numeric(summary_speed$est)
summary_speed$se = as.numeric(summary_speed$se)
summary_speed$ucl = as.numeric(summary_speed$ucl)
summary_speed$lcl = as.numeric(summary_speed$lcl)

#Generating a quick plot to visualize the HR size by site
ggplot(data=summary_speed,aes(x=species, y=est))+
  geom_errorbar(aes(ymin=lcl,ymax=ucl),size=.8,width=.5, colour="gray")+
  geom_point(aes(),size=5)+
  ylab('Mean Daily Distance km')+
  xlab("Species")+
  theme_bw()+
  theme(text = element_text(size=rel(5)),
        strip.text.x = element_text(size=rel(5)),
        strip.text.y = element_text(size=rel(5)))
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#ggtitle("AKDE Variation by Species")


############################################
############################################
############################################
library(lme4)
require(AICcmodavg)
require(kimisc)

#Generating some models
#Simple linear regressions
#m1 = glm(df2$log_akde_95 ~ df2$ndvi_mean) 
#m2 = glm(df2$log_akde_95 ~ df2$ndvi_sd) 
m3 = glm(df2$log_akde_95 ~ df2$ProtectedArea) 
m4 = glm(df2$log_akde_95 ~ df2$habitat) 
#m5 = glm(df2$log_akde_95 ~ df2$ndvi_mean) 
#m6 = glm(df2$log_akde_95 ~ df2$temp_mean) 
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

models_lm = nlist(m3,m4,m7,m9,m10,m11,m12,m13,m14,m16,m17,m18,m19)
models_lm = nlist(m1,m2,m3,m4,m5,m6,m7,m9,m10,m11,m12,m13,m14,m16,m17,m18,m19)
aic_table= aictab(models_lm)
aic_table

#Trying Global Models
require(MuMIn)
df3 = df2 %>%
  dplyr::select(log_akde_95,log_speed, var_log_speed,var_log_area_95,ndvi_mean,ndvi_sd,ProtectedArea,temp_mean, temp_sd, habitat, hfp2009_mean, hfp2009_sd, human_density_mean, human_density_sd,species,Park_location) %>%
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
# globalmodel1 = lm(log_akde_95 ~ Park_location,ndvi_mean + ndvi_sd + ProtectedArea + temp_mean + temp_sd + habitat + hfp2009_mean + hfp2009_sd + human_density_mean + human_density_sd, data=df3, na.action = na.fail) 
# combinations = dredge(globalmodel1)
# print(combinations)

globalmodel2 = lm(log_akde_95 ~ ndvi_mean + ndvi_sd + ProtectedArea + temp_mean + temp_sd + habitat + hfp2009_mean + hfp2009_sd + human_density_mean + human_density_sd + ndvi_mean * habitat + ndvi_mean * human_density_mean + ndvi_sd * habitat + ndvi_sd * human_density_mean + ndvi_mean * ProtectedArea + ndvi_sd * ProtectedArea, data=df3, na.action = na.fail) 
combinations2 = dredge(globalmodel2)
print(combinations2)

#'Best' model
summary(get.models(combinations2, 1)[[1]])

###Looking at Mixed Effects Models (without the meta framework)
#Creating a mixed effects model to look at the effect of time on lesion size
library(lme4)
# glmm_full = lmer(log_akde_95 ~ ndvi_mean + 
#                 ndvi_sd + 
#                 ProtectedArea +
#                 temp_mean +
#                 temp_sd +
#                 habitat +
#                 hfp2009_mean +
#                 hfp2009_sd +
#                 human_density_mean +
#                 human_density_sd +
#                 ndvi_mean * habitat +
#                 ndvi_mean * human_density_mean +
#                 ndvi_sd * habitat +
#                 ndvi_sd * human_density_mean +
#                 ndvi_mean * ProtectedArea + 
#                 ndvi_sd * ProtectedArea + 
#                 (1|Park_location),
#               data = df3, REML = FALSE, na.action = 'na.exclude')

glmm_full = lmer(log_akde_95 ~ ndvi_mean + 
                   ndvi_sd + 
                   ProtectedArea +
                   temp_mean +
                   temp_sd +
                   habitat +
                   hfp2009_mean +
                   hfp2009_sd +
                   human_density_mean +
                   human_density_sd +
                   ndvi_mean * habitat +
                   ndvi_mean * human_density_mean +
                   ndvi_sd * habitat +
                   ndvi_sd * human_density_mean +
                   ndvi_mean * ProtectedArea + 
                   ndvi_sd * ProtectedArea + 
                   (1|Park_location),na.action= "na.fail",
                 data = df3, REML = FALSE)
dredge_glmm = dredge(glmm_full)
print(dredge_glmm)

#'Best' model
summary(get.models(dredge_glmm, 1)[[1]])

#Report the model summary
library(report)
#glmm_results <- report(glmm_1)
glmm_results <- report(get.models(dredge_glmm, 1)[[1]])
glmm_results


#Attempt at a mixed effects meta1_ analysis
meta1_1 = rma(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean, 
                
               method= "REML",
               data=df3)

meta1_2 = rma(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_sd, 
                
               method= "REML",
               data=df3)

meta1_3 = rma(log_akde_95,var_log_area_95, mods = ~ 
                 ProtectedArea,
                
               method= "REML",
               data=df3)

meta1_4 = rma(log_akde_95,var_log_area_95, mods = ~ 
                 temp_mean,
                
               method= "REML",
               data=df3)

meta1_5 = rma(log_akde_95,var_log_area_95, mods = ~
                 temp_sd, 
                
               method= "REML",
               data=df3)

meta1_6 = rma(log_akde_95,var_log_area_95, mods = ~ 
                 habitat, 
                
               method= "REML",
               data=df3)

meta1_7 = rma(log_akde_95,var_log_area_95, mods = ~ 
                 hfp2009_mean,
                
               method= "REML",
               data=df3)

meta1_8 = rma(log_akde_95,var_log_area_95, mods = ~ 
                 hfp2009_sd,
                
               method= "REML",
               data=df3)

meta1_9 = rma(log_akde_95,var_log_area_95, mods = ~
                 human_density_mean,
                
               method= "REML",
               data=df3)

meta1_10 = rma(log_akde_95,var_log_area_95, mods = ~ 
                  human_density_sd,
                 
                method= "REML",
                data=df3)

meta1_11 = rma(log_akde_95,var_log_area_95, mods = ~ 
                  ndvi_mean * habitat,
                 
                method= "REML",
                data=df3)

meta1_12 = rma(log_akde_95,var_log_area_95, mods = ~ 
                  ndvi_mean +
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd +
                  ndvi_mean * habitat +
                  ndvi_mean * human_density_mean +
                  ndvi_sd * habitat +
                  ndvi_sd * human_density_mean +
                  ndvi_mean * ProtectedArea + 
                  ndvi_sd * ProtectedArea, 
                 
                method= "REML",
                data=df3)

meta1_13 = rma(log_akde_95,var_log_area_95, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd +
                  ndvi_mean * habitat +
                  ndvi_mean * human_density_mean +
                  ndvi_sd * habitat +
                  ndvi_sd * human_density_mean +
                  ndvi_mean * ProtectedArea + 
                  ndvi_sd * ProtectedArea, 
                 
                method= "REML",
                data=df3)

meta1_14 = rma(log_akde_95,var_log_area_95, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd +
                  ndvi_mean * habitat +
                  ndvi_mean * human_density_mean +
                  ndvi_sd * habitat +
                  ndvi_sd * human_density_mean +
                  ndvi_mean * ProtectedArea, 
                 
                method= "REML",
                data=df3)

meta1_15 = rma(log_akde_95,var_log_area_95, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd +
                  ndvi_mean * habitat +
                  ndvi_mean * human_density_mean +
                  ndvi_sd * habitat +
                  ndvi_sd * human_density_mean, 
                 
                method= "REML",
                data=df3)

meta1_16 = rma(log_akde_95,var_log_area_95, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd, 
                 
                method= "REML",
                data=df3)

meta1_17 = rma(log_akde_95,var_log_area_95, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean, 
                 
                method= "REML",
                data=df3)

models_meta1__summary = nlist(meta1_1,meta1_2,meta1_3,meta1_4,meta1_5,meta1_6,meta1_7,meta1_8,meta1_9,meta1_10,meta1_11,meta1_12,meta1_13,meta1_14,meta1_15,meta1_16,meta1_17)
aic_table2 = aictab(models_meta1__summary)
aic_table2

AIC(meta1_1)
AIC(meta1_2)
AIC(meta1_3)
AIC(meta1_4)
AIC(meta1_5)
AIC(meta1_6)
AIC(meta1_7)
AIC(meta1_8)
AIC(meta1_9)
AIC(meta1_10)
AIC(meta1_11)
AIC(meta1_12)
AIC(meta1_13)
AIC(meta1_14)
AIC(meta1_15)
AIC(meta1_16)
AIC(meta1_17)




#Attempt at a mixed effects meta analysis
meta1 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta2 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_sd, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta3 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ProtectedArea,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta4 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 temp_mean,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta5 = rma.mv(log_akde_95,var_log_area_95, mods = ~
                 temp_sd, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta6 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 habitat, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta7 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 hfp2009_mean,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta8 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 hfp2009_sd,
               random= ~1|Park_location,
                     method= "REML",
                   data=df3)

meta9 = rma.mv(log_akde_95,var_log_area_95, mods = ~
                 human_density_mean,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta10 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 human_density_sd,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta11 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean * habitat,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta12 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean +
                 ProtectedArea +
                 temp_mean +
                 temp_sd +
                 habitat +
                 hfp2009_mean +
                 hfp2009_sd +
                 human_density_mean +
                 human_density_sd +
                 ndvi_mean * habitat +
                 ndvi_mean * human_density_mean +
                 ndvi_sd * habitat +
                 ndvi_sd * human_density_mean +
                 ndvi_mean * ProtectedArea + 
                 ndvi_sd * ProtectedArea, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta13 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean +
                 ndvi_sd + 
                 ProtectedArea +
                 temp_mean +
                 habitat +
                 hfp2009_mean +
                 hfp2009_sd +
                 human_density_mean +
                 human_density_sd +
                 ndvi_mean * habitat +
                 ndvi_mean * human_density_mean +
                 ndvi_sd * habitat +
                 ndvi_sd * human_density_mean +
                 ndvi_mean * ProtectedArea + 
                 ndvi_sd * ProtectedArea, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta14 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean +
                 ndvi_sd + 
                 ProtectedArea +
                 temp_mean +
                 temp_sd +
                 habitat +
                 hfp2009_mean +
                 hfp2009_sd +
                 human_density_mean +
                 human_density_sd +
                 ndvi_mean * habitat +
                 ndvi_mean * human_density_mean +
                 ndvi_sd * habitat +
                 ndvi_sd * human_density_mean +
                 ndvi_mean * ProtectedArea, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta15 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean +
                 ndvi_sd + 
                 ProtectedArea +
                 temp_mean +
                 temp_sd +
                 habitat +
                 hfp2009_mean +
                 hfp2009_sd +
                 human_density_mean +
                 human_density_sd +
                 ndvi_mean * habitat +
                 ndvi_mean * human_density_mean +
                 ndvi_sd * habitat +
                 ndvi_sd * human_density_mean, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta16 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean +
                 ndvi_sd + 
                 ProtectedArea +
                 temp_mean +
                 temp_sd +
                 habitat +
                 hfp2009_mean +
                 hfp2009_sd +
                 human_density_mean +
                 human_density_sd, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta17 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean +
                 ndvi_sd + 
                 ProtectedArea +
                 temp_mean +
                 temp_sd +
                 habitat +
                 hfp2009_mean, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

models_meta_summary = nlist(meta1,meta2,meta3,meta4,meta5,meta6,meta7,meta8,meta9,meta10,meta11,meta12,meta13,meta14,meta15,meta16,meta17)
aic_table2 = aictab(models_meta_summary)
aic_table2

AIC(meta1)
AIC(meta2)
AIC(meta3)
AIC(meta4)
AIC(meta5)
AIC(meta6)
AIC(meta7)
AIC(meta8)
AIC(meta9)
AIC(meta10)
AIC(meta11)
AIC(meta12)
AIC(meta13)
AIC(meta14)
AIC(meta15)
AIC(meta16)
AIC(meta17)

glmm_meta_results = report(meta1)

meta1.df = as.data.frame(meta1[1])
meta1.df <- cbind(rownames(meta1.df), data.frame(meta1.df, row.names=NULL))
meta1.df['lcl'] = meta1[6]
meta1.df['lucl'] = meta1[7]

###############################################
###############################################
###############################################

#Attempt at a mixed effects meta_speed analysis
meta_speed1 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 ndvi_mean, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta_speed2 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 ndvi_sd, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta_speed3 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 ProtectedArea,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta_speed4 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 temp_mean,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta_speed5 = rma.mv(log_speed,var_log_speed, mods = ~
                 temp_sd, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta_speed6 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 habitat, 
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta_speed7 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 hfp2009_mean,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta_speed8 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 hfp2009_sd,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta_speed9 = rma.mv(log_speed,var_log_speed, mods = ~
                 human_density_mean,
               random= ~1|Park_location,
               method= "REML",
               data=df3)

meta_speed10 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  human_density_sd,
                random= ~1|Park_location,
                method= "REML",
                data=df3)

meta_speed11 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean * habitat,
                random= ~1|Park_location,
                method= "REML",
                data=df3)

meta_speed12 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean +
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd + ndvi_mean * habitat + ndvi_mean * human_density_mean +ndvi_sd * habitat +ndvi_sd * human_density_mean +ndvi_mean * ProtectedArea + ndvi_sd * ProtectedArea, random= ~1|Park_location,
                method= "REML",
                data=df3)

meta_speed13 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd +
                  ndvi_mean * habitat +
                  ndvi_mean * human_density_mean +
                  ndvi_sd * habitat +
                  ndvi_sd * human_density_mean +
                  ndvi_mean * ProtectedArea + 
                  ndvi_sd * ProtectedArea, 
                random= ~1|Park_location,
                method= "REML",
                data=df3)

meta_speed14 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd +
                  ndvi_mean * habitat +
                  ndvi_mean * human_density_mean +
                  ndvi_sd * habitat +
                  ndvi_sd * human_density_mean +
                  ndvi_mean * ProtectedArea, 
                random= ~1|Park_location,
                method= "REML",
                data=df3)

meta_speed15 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd +
                  ndvi_mean * habitat +
                  ndvi_mean * human_density_mean +
                  ndvi_sd * habitat +
                  ndvi_sd * human_density_mean, 
                random= ~1|Park_location,
                method= "REML",
                data=df3)

meta_speed16 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean +
                  hfp2009_sd +
                  human_density_mean +
                  human_density_sd, 
                random= ~1|Park_location,
                method= "REML",
                data=df3)

meta_speed17 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean +
                  ndvi_sd + 
                  ProtectedArea +
                  temp_mean +
                  temp_sd +
                  habitat +
                  hfp2009_mean, 
                random= ~1|Park_location,
                method= "REML",
                data=df3)

models_meta_speed_summary = nlist(meta_speed1,meta_speed2,meta_speed3,meta_speed4,meta_speed5,meta_speed6,meta_speed7,meta_speed8,meta_speed9,meta_speed10,meta_speed11,meta_speed12,meta_speed13,meta_speed14,meta_speed15,meta_speed16,meta_speed17)
aic_table2 = aictab(models_meta_speed_summary)
aic_table2

AIC(meta_speed1)
AIC(meta_speed2)
AIC(meta_speed3)
AIC(meta_speed4)
AIC(meta_speed5)
AIC(meta_speed6)
AIC(meta_speed7)
AIC(meta_speed8)
AIC(meta_speed9)
AIC(meta_speed10)
AIC(meta_speed11)
AIC(meta_speed12)
AIC(meta_speed13)
AIC(meta_speed14)
AIC(meta_speed15)
AIC(meta_speed16)
AIC(meta_speed17)


#plotting out Effect size in forest plot
top_speed_model = meta_speed14
speed_sum_df = as.data.frame(top_speed_model$beta)
library(tibble)
speed_sum_df <- tibble::rownames_to_column(speed_sum_df, "Effect")
speed_sum_df$Beta = speed_sum_df$V1
speed_sum_df$LCI = top_speed_model$ci.lb
speed_sum_df$UCI = top_speed_model$ci.ub

speed_sum_df = speed_sum_df %>%
  dplyr::filter(!Effect =="intrcpt")

#Generating a quick plot to visualize the HR size by Species
ggplot(data=speed_sum_df,aes(x=Beta, y= Effect))+
  geom_point()+
  geom_errorbar(aes(xmin=LCI,xmax=UCI))+
  xlab("Effect Size")+
  ylab("Moderator")+
  geom_vline(xintercept = 0, colour="red")+
  theme_minimal()+
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size")



