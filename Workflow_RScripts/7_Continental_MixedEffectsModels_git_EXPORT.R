##THE CONTINENTAL
##Mixed Effects Models in the Meta Analysis Framework
#Author: Michael Butler Brown
#Author Affiliations: Smithsonian Conservation Biology Institute and Giraffe Conservation Foundation

#See Here for a good tutorial: https://www.youtube.com/watch?v=E7-EI13FGKc
#See also: https://www.youtube.com/watch?v=IkduL5iRdqo&t=1013s

#Start Fresh
rm(list=ls())

#Load Relevant packages
library(plyr)
library(dplyr)
library(tidyverse)
library(metafor)
library(forestplot)
library(forcats)
library(ggplot2)

#Set Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read in Movement Data
df1 = read.csv(file="../Data/continental_df5.csv", header = TRUE, sep =",")
df=df1
collar_key = read.csv(file="../Data/ctmm_key2.csv",head=TRUE,sep=",")
count= length(dplyr::filter(df,df$Variogram_inspection == "bad"))
#Changing some site names in northern Kenya
df$Park_location = replace(df$Park_location, df$Park_location =="Biliqo_Bulesa","NorthernKenya")#For replacing characters with other characters
df$Park_location = replace(df$Park_location, df$Park_location =="Melako","NorthernKenya")#For replacing characters with other characters

#Transforming some of the response variables for gaussian modeling
df$log_akde_95 = log(df$akde_95_ML) #log transformation for the 95% akde 
df$log_akde_50 = log(df$akde_50_ML) # log transformation for the 50% akde
df$log_speed = log(df$speed_ml)# log transformation for the speed
df$var_log_area_95 = 1/df$akde_95_DOF_area #Proper variance for subsequent metafor analyses
df$var_log_area_50 = 1/df$akde_50_DOF_area #Proper variance for subsequent metafor analyses
df$var_log_speed = (1/4)/df$DOF_speed

# df$var_log_area_95 = log(df$df$akde_95_DOF_area)^2/df$akde_50_DOF_area #Proper variance for subsequent metafor analyses
# df$var_log_area_50 = log(df$akde_50_ML)^2/df$akde_50_DOF_area #Proper variance for subsequent metafor analyses
# df$var_log_speed = log(df$speed_ml)^2/df$DOF_speed #Proper variance for subsequent metafor analyses

#Some Simple filtering to exclude giraffe with bad models
df = df %>%
  dplyr::filter(!Duration < 152)%>% # If there are fewer than 500 fixes...I didn't run many ctmm models on individuals with fewer than 500 fixes, so there shouldn't be too many of these
  #dplyr::filter(!Number_of_fixes < 500)%>% # If there are fewer than 500 fixes...I didn't run many ctmm models on individuals with fewer than 500 fixes, so there shouldn't be too many of these
  dplyr::filter(!DOF_area < 4)%>% #If DOF area is less than 4
  dplyr::filter(!duplicate == "yes")%>% #If DOF area is less than 4
  #dplyr::filter(Country=="Niger")%>%
  dplyr::filter(!akde_95_ML > 20000)%>% #Just to get rid of that one weird one up in Niger
  dplyr::filter(!Variogram_inspection =="bad")%>% #If variograms didn't converge on home range behaviour as determined by subjective visual inspection
  drop_na(model) #Get rid of all collars without sufficient model data

species_summary = df %>%
  group_by(species)%>%
  dplyr::summarise(total=n())
species_summary

movement_summary = df %>%
  group_by(model)%>%
  dplyr::summarise(total=n())
movement_summary

sex_summary = df %>%
  group_by(sex)%>%
  dplyr::summarise(total=n())
sex_summary

head(df)


#Generating a quick plot to visualize the HR size by Species
ggplot(data=df,aes(x=akde_95_ML, y= ID))+
  geom_point(aes(color= Park_location))+
  geom_errorbar(aes(xmin=akde_95_lcl,xmax=akde_95_ucl,color=Park_location))+
  xlab("AKDE 95% Size")+
  ylab(" Collar ID")+
  facet_wrap(~species, scale = "free")+
  xlim(0,600) #excluding point estimates above 600km2 in the graph just so we can visually inspect the smaller sizes
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size")
  
  
library(ggpmisc)
#Generating a quick plot to visualize the HR size by Species
  ggplot(data=df,aes(x=ndvi_mean, y= log_akde_95))+
    geom_point(aes(color= species))+
        #geom_errorbar(aes(ymin=akde_95_lcl,ymax=akde_95_ucl,color=Park_location))+
    xlab("Mean NDVI")+
    ylab("AKDE 95% Size")+
    stat_smooth(aes(fill = species, color = species), method = "lm", formula = y ~ x) +
    #stat_regline_equation(aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~"))) +
    #geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
    facet_wrap(~species, scale = "free")
    #xlim(0,600) #excluding point estimates above 600km2 in the graph just so we can visually inspect the smaller sizes
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggtitle("AKDE_95 Size and NDVI by Species")
  
df_trial = df %>%
  dplyr::filter(species =="tippelskirchi")
test = lm(df_trial$log_akde_95 ~ df_trial$ndvi_mean)  
plot(df_trial$log_akde_95 ~ df_trial$ndvi_mean)    
abline(test)
summary(test)
    
  ### copy BCG vaccine data into 'dat'
  dat <- dat.bcg
  
  ### calculate log risk ratios and corresponding sampling variances
  dat <- escalc(measure="RR", ai=tpos, bi=tneg, ci=cpos, di=cneg, data=dat,
                slab=paste(author, ", ", year, sep=""))
  
  ### fit random-effects model
  res <- rma(yi, vi, data=dat)
  
  # \dontrun{
  ### generate report
 # reporter(res)
  #> 
  #> Directory for generating the report is: /tmp/RtmpAWjSbW
  #> Copying references.bib and apa.csl to report directory ...
  #> Saving model object to report_res.rdata ...
  #> Creating report_res.rmd file ...
  #> Rendering report_res.rmd file ...
  #> Generated /tmp/RtmpAWjSbW/report_res.html ...
  #> Opening report ...
  # }
  
  
#########################################
#Looking at some quick forest plots per species to compare HR estimates
df = df %>%
  dplyr::filter(!is.na(log_akde_95))%>%
  #dplyr::filter(!is.na(temp_mean))%>%
  dplyr::filter(!is.na(ProtectedArea))%>%
  dplyr::filter(!is.na(ndvi_mean))%>%
  dplyr::filter(!is.na(ndvi_sd))%>%
  dplyr::filter(!is.na(habitat))%>%
  dplyr::filter(!is.na(hfp2009_mean))%>%
  dplyr::filter(!is.na(human_density_mean))

  
#Scale the covariates  
df$human_density_mean_scale = scale(df$human_density_mean)
df$human_density_sd_scale = scale(df$human_density_sd)
df$hfp2009_mean_scale = scale(df$hfp2009_mean)
df$hfp2009_sd_scale = scale(df$hfp2009_sd)
df$ndvi_mean_scale = scale(df$ndvi_mean)
df$ndvi_sd_scale = scale(df$ndvi_sd)
df$ProtectedArea_scale = scale(df$ProtectedArea)
df$woody_mean_scale = scale(df$woody_mean)
df$woody_sd_scale = scale(df$woody_sd)
df$cop_tree_mean_scale = scale(df$cop_tree_mean)
df$cop_tree_sd_scale = scale(df$cop_tree_sd)
df$cop_shrub_mean_scale = scale(df$cop_shrub_mean)
df$cop_shrub_sd_scale = scale(df$cop_shrub_sd)

species_summary = df%>%
  dplyr::group_by(species)%>%
  dplyr::summarise(total=n())
species_summary
  
ctmm_summary = df%>%
  dplyr::group_by(model)%>%
  dplyr::summarise(total=n())
ctmm_summary

#Calculating mean AKDE from raw ML point estimates...for comparison purposes
mean(df$akde_95_ML) #Looking at mean value 
summary= df %>%
  dplyr::summarise(mean = mean(akde_95_ML),
                   total = n(),
                   sd = sd(akde_95_ML),
                   se = sd/sqrt(total))
alpha = 0.05
degrees.freedom = summary$total - 1
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
print(t.score)
summary$ucl = summary$mean + (t.score * summary$se)
summary$lcl = summary$mean + (t.score - summary$se)

#summary
#ci_model = lm(df$akde_95_ML ~ 1)
#confint(ci_model, level= 0.95)

#Summarise by Species
summary_species= df %>%
  dplyr::group_by(species)%>%
  dplyr::summarise(mean = mean(akde_95_ML),
                   total = n(),
                   sd = sd(akde_95_ML),
                   se = sd/sqrt(total))
alpha = 0.05
degrees.freedom = summary_species$total - 1
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
print(t.score)
summary_species$ucl = summary_species$mean + (t.score * summary_species$se)
summary_species$lcl = summary_species$mean - (t.score * summary_species$se)
summary_species

#Calculating Global Mean AKDE using metafor::rma()
library(metafor)
res1_akde <- rma(yi=log_akde_95, vi=var_log_area_95, data=df)
res1_akde
mean_global = predict(res1_akde, transf=exp,digits=2)
mean_global
forest(res1_akde)

#Calculating Species Mean AKDE using metafor::rma()
### Fit Fixed Model with species as moderator
library(metafor)
res2_akde <- rma(yi=log_akde_95, vi=var_log_area_95,mods = ~species, data=df)
res2_akde
forest(res2_akde)
#predicted effects 
df_species = as.data.frame(df %>% dplyr::select(species))
mean_akde_species = cbind(df_species,as.data.frame(predict(res2_akde,transf=exp,digits=2)))%>%
  dplyr::distinct(species,.keep_all = TRUE)
mean_akde_species
mean_akde_species$species = replace(mean_akde_species$species, mean_akde_species$species =="camelopardalis","G. camelopardalis")#For replacing characters with other characters
mean_akde_species$species = replace(mean_akde_species$species, mean_akde_species$species =="giraffa","G. giraffa")#For replacing characters with other characters
mean_akde_species$species = replace(mean_akde_species$species, mean_akde_species$species =="reticulata","G. reticulata")#For replacing characters with other characters
mean_akde_species$species = replace(mean_akde_species$species, mean_akde_species$species =="tippelskirchi","G. tippelskirchi")#For replacing characters with other characters

speed_summary$Effect = replace(speed_summary$Effect, speed_summary$Effect =="woody_sd_scale","Woody Biomass (SD)")#For replacing characters with other characters

#Generating a quick plot to visualize the HR size by site
ggplot(data=mean_akde_species,aes(x=species, y=pred))+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub),size=.8,width=.5, colour="black")+
  geom_point(aes(),size=5)+
  ylab(bquote('AKDE 95% Size'~km^2))+
  xlab("Species")+
  geom_abline(intercept=371.1, slope=0, size=0.75, color="red")+
  geom_abline(intercept=438.4, slope=0,size=0.2, linetype="dashed",color="red")+
  geom_abline(intercept=314.2, slope=0,size=0.2, linetype="dashed",color="red")+
  theme_bw()+
  theme(text = element_text(size=rel(5.5)),
        strip.text.x = element_text(size=rel(5)),
        strip.text.y = element_text(size=rel(5)))
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#ggtitle("AKDE Variation by Species")

theme(text = element_text(size=rel(3.5)),
      strip.text.x = element_text(size=rel(3.5)),
      strip.text.y = element_text(size=rel(3.5)))

#######################
##SPEED##
#######################
#Calculating mean AKDE from raw ML point estimates...for comparison purposes
mean(df$speed_ml) #Looking at mean value 

#Calculating Global Mean AKDE using metafor::rma()
library(metafor)
res1_speed <- rma(yi=log_speed, vi=var_log_speed, data=df)
res1_speed
mean_global = predict(res1_speed, transf=exp,digits=2)
mean_global
forest(res1_speed)

#Calculating Species Mean AKDE using metafor::rma()
### Fit Fixed Model with species as moderator
library(metafor)
res2_speed <- rma(yi=log_speed, vi=var_log_speed,mods = ~species, data=df)
res2_speed
forest(res2_speed)

reporter(res2_speed)
#predicted effects 
df_species = as.data.frame(df %>% dplyr::select(species))
mean_speed_species = cbind(df_species,as.data.frame(predict(res2_speed,transf=exp,digits=2)))%>%
  dplyr::distinct(species,.keep_all = TRUE)
mean_speed_species

mean_speed_species$species = replace(mean_speed_species$species, mean_speed_species$species =="camelopardalis","G. camelopardalis")#For replacing characters with other characters
mean_speed_species$species = replace(mean_speed_species$species, mean_speed_species$species =="giraffa","G. giraffa")#For replacing characters with other characters
mean_speed_species$species = replace(mean_speed_species$species, mean_speed_species$species =="reticulata","G. reticulata")#For replacing characters with other characters
mean_speed_species$species = replace(mean_speed_species$species, mean_speed_species$species =="tippelskirchi","G. tippelskirchi")#For replacing characters with other characters


#Calculating mean AKDE from raw ML point estimates...for comparison purposes
mean(df$speed_ml) #Looking at mean value 
summary_speed= df %>%
  dplyr::group_by(species)%>%
  dplyr::summarise(mean = mean(speed_ml),
                   total = n(),
                   sd = sd(speed_ml),
                   se = sd/sqrt(total))
alpha = 0.05
degrees.freedom = summary$total - 1
t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
print(t.score)
summary_speed$ucl = summary_speed$mean + (t.score * summary_speed$se)
summary_speed$lcl = summary_speed$mean - (t.score + summary_speed$se)
summary_speed

#Generating a quick plot to visualize the HR size by site
ggplot(data=mean_speed_species,aes(x=species, y=pred))+
  geom_errorbar(aes(ymin=ci.lb,ymax=ci.ub),size=.8,width=.5, colour="black")+
  geom_point(aes(),size=5)+
  ylab(bquote('Daily Speed'~km))+
  xlab("Species")+
  geom_abline(intercept=14.2, slope=0, size=0.75, color="red")+
  geom_abline(intercept=14.7, slope=0,size=0.2, linetype="dashed",color="red")+
  geom_abline(intercept=13.7, slope=0,size=0.2, linetype="dashed",color="red")+
  theme_bw()+
  theme(text = element_text(size=rel(5.5)),
        strip.text.x = element_text(size=rel(5)),
        strip.text.y = element_text(size=rel(5)))
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#ggtitle("AKDE Variation by Species")

################################################
################################################

library(lme4)
require(kimisc)
#Attempt at a mixed effects meta analysis with site as a random effect
meta1 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean_scale, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta2 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_sd_scale, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

library(lme4)
lmm1 <- lmer(log_akde_95 ~ ndvi_sd_scale + (1|Park_location), data=df)
lmm1

meta3 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ProtectedArea,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta5 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 woody_mean_scale, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta6 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 hfp2009_mean_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta7 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 hfp2009_sd_scale,
               random= ~1|Park_location,
                     method= "REML",
                   data=df)

meta8 = rma.mv(log_akde_95,var_log_area_95, mods = ~
                 human_density_mean_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta9 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 human_density_sd_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta10 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
               woody_mean_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta11 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
               woody_sd_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta12 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean_scale * woody_mean_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta13 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean_scale +
                 ndvi_sd_scale + 
                 ProtectedArea +
                 woody_mean_scale +
                 woody_sd_scale +
                 #habitat +
                 hfp2009_mean_scale +
                 hfp2009_sd_scale +
                 #human_density_mean_scale +
                 #human_density_sd_scale +
                 ndvi_mean_scale * woody_mean_scale + 
                 #ndvi_mean_scale * woody_sd_scale +
                 #ndvi_sd_scale * woody_mean_scale + 
                 #ndvi_sd_scale * woody_sd_scale+
                 #ndvi_mean_scale * habitat +
                 ndvi_mean_scale * hfp2009_mean_scale +
                 #ndvi_sd_scale * habitat +
                 #ndvi_sd_scale * hfp2009_mean_scale +
                 #ndvi_sd_scale * ProtectedArea + 
                 ndvi_mean_scale * ProtectedArea, 
               random= ~1|Park_location,
               method= "REML",
               data=df)


trial = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                  ndvi_mean_scale +
                  ndvi_sd_scale + 
                  ProtectedArea +
                  woody_mean_scale +
                  woody_sd_scale +
                  #habitat +
                  hfp2009_mean_scale +
                  hfp2009_sd_scale +
                  #human_density_mean_scale +
                  #human_density_sd_scale +
                  ndvi_mean_scale * woody_mean_scale + 
                  #ndvi_mean_scale * woody_sd_scale +
                  #ndvi_sd_scale * woody_mean_scale + 
                  #ndvi_sd_scale * woody_sd_scale+
                  #ndvi_mean_scale * habitat +
                  ndvi_mean_scale * hfp2009_mean_scale +
                  #ndvi_sd_scale * habitat +
                  #ndvi_sd_scale * hfp2009_mean_scale +
                  #ndvi_sd_scale * ProtectedArea + 
                  ndvi_mean_scale * ProtectedArea, 
                #random= ~1|Park_location,
                method= "REML",
                data=df)


meta14 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean_scale +
                 ndvi_sd_scale + 
                 ProtectedArea +
                 #habitat +
                 woody_mean_scale + 
                 woody_sd_scale + 
                 hfp2009_mean_scale +
                 hfp2009_sd_scale +
                 ndvi_mean_scale * woody_mean_scale +
                 ndvi_mean_scale * hfp2009_mean_scale +
                 ndvi_sd_scale * woody_mean_scale +
                 ndvi_sd_scale * hfp2009_mean_scale +
                 ndvi_mean_scale * ProtectedArea + 
                 ndvi_sd_scale * ProtectedArea, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta14.2 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                  ndvi_mean_scale +
                  ndvi_sd_scale + 
                  ProtectedArea +
                  #habitat +
                  cop_tree_mean_scale + 
                  cop_tree_sd_scale + 
                  hfp2009_mean_scale +
                  hfp2009_sd_scale +
                  ndvi_mean_scale * woody_mean_scale +
                  ndvi_mean_scale * hfp2009_mean_scale +
                  ndvi_sd_scale * woody_mean_scale +
                  ndvi_sd_scale * hfp2009_mean_scale +
                  ndvi_mean_scale * ProtectedArea + 
                  ndvi_sd_scale * ProtectedArea, 
                random= ~1|Park_location,
                method= "REML",
                data=df)
AIC(meta14)
AIC(meta14.2)

meta15 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean_scale +
                 ndvi_sd_scale + 
                 ProtectedArea +
                 #habitat +
                 woody_mean_scale +
                 woody_sd_scale + 
                 hfp2009_mean_scale +
                 hfp2009_sd_scale +
                 ndvi_mean_scale * woody_mean_scale +
                 ndvi_mean_scale * hfp2009_mean_scale +
                 ndvi_sd_scale * woody_mean_scale +
                 ndvi_sd_scale * hfp2009_mean_scale +
                 ndvi_mean_scale * ProtectedArea, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta16 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean_scale +
                 ndvi_sd_scale + 
                 ProtectedArea +
                 #habitat +
                 woody_mean_scale + 
                 woody_sd_scale +
                 hfp2009_mean_scale +
                 hfp2009_sd_scale +
                 ndvi_mean_scale * woody_mean_scale +
                 ndvi_mean_scale * hfp2009_mean_scale +
                 ndvi_sd_scale * woody_mean_scale +
                 ndvi_sd_scale * hfp2009_mean_scale, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta17 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean_scale +
                 ndvi_sd_scale + 
                 ProtectedArea +
                 woody_mean_scale +
                 hfp2009_mean_scale +
                 hfp2009_sd_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta18 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                 ndvi_mean_scale +
                 ndvi_sd_scale + 
                 ProtectedArea +
                 woody_mean_scale +
                 hfp2009_mean_scale, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta19 = rma.mv(log_akde_95,var_log_area_95, mods = ~ 
                   ndvi_mean_scale +
                   ndvi_sd_scale + 
                   ProtectedArea +
                   #habitat +
                   hfp2009_mean_scale +
                   woody_mean_scale+
                   woody_sd_scale, 
                 random= ~1|Park_location,
                 method= "REML",
                 data=df)


###############################################
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
AIC(meta18)
AIC(meta19)

#Model Summary
top_model1 = meta14
akde_summary = as.data.frame(top_model1[2])
akde_summary <- as.data.frame(tibble::rownames_to_column(akde_summary, "Effect"))
akde_summary$se = top_model1[3]
akde_summary$se = top_model1$se
akde_summary$CI_lower = top_model1$ci.lb 
akde_summary$CI_upper = top_model1$ci.ub

#Remove the intercept from the model
akde_summary = akde_summary %>%
  dplyr::filter(!Effect =="intrcpt")
akde_summary = akde_summary[-c(1),]
akde_summary=as.data.frame(akde_summary)
# akde_summary2 = akde_summary%>%
#   dplyr::select(!Effect=="intrcpt")

#Generating a quick plot to visualize effect sizes as a forest plot
ggplot(data=akde_summary,aes(x=beta, y= Effect))+
  geom_point()+
  geom_vline(xintercept = 0, linetype="dotted", color="red", size=1.5)+
  geom_errorbar(aes(xmin = CI_lower, xmax =CI_upper))+
  xlab("Effect Size on 95% AKDE")+
  ylab("Covariate")

##Developing Variance Components for Meta13
#SEE THIS HELP PAGE: https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
#How much Variance is due to heterogeneity?
W <- diag(1/meta13$vi)
X <- model.matrix(meta13)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(meta13$sigma2) / (sum(meta13$sigma2) + (meta13$k-meta13$p)/sum(diag(P)))

#How much variance is due to between cluster heterogeneity
100 * meta13$sigma2 / (sum(meta13$sigma2) + (meta13$k-meta13$p)/sum(diag(P)))

#Attempt at a mixed effects meta_speed analysis with site as a random effect
meta_speed1 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 ndvi_mean_scale, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta_speed2 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 ndvi_sd_scale, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta_speed3 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 ProtectedArea,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta_speed5 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 woody_mean_scale, 
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta_speed6 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 hfp2009_mean_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta_speed7 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 hfp2009_sd_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta_speed8 = rma.mv(log_speed,var_log_speed, mods = ~
                 human_density_mean_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta_speed9 = rma.mv(log_speed,var_log_speed, mods = ~ 
                 human_density_sd_scale,
               random= ~1|Park_location,
               method= "REML",
               data=df)

meta_speed10 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  woody_mean_scale,
                random= ~1|Park_location,
                method= "REML",
                data=df)

meta_speed11 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  woody_sd_scale,
                random= ~1|Park_location,
                method= "REML",
                data=df)

meta_speed12 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean_scale * woody_mean_scale,
                random= ~1|Park_location,
                method= "REML",
                data=df)

meta_speed13.2 = rma.mv(log_speed,var_log_speed, mods = ~ 
                          ndvi_mean_scale +
                          ndvi_sd_scale + 
                          ProtectedArea +
                          woody_mean_scale +
                          woody_sd_scale +
                          #habitat +
                          hfp2009_mean_scale +
                          hfp2009_sd_scale +
                          #human_density_mean_scale +
                          #human_density_sd_scale +
                          ndvi_mean_scale * woody_mean_scale + 
                          #ndvi_mean_scale * woody_sd_scale +
                          #ndvi_sd_scale * woody_mean_scale + 
                          #ndvi_sd_scale * woody_sd_scale+
                          #ndvi_mean_scale * habitat +
                          ndvi_mean_scale * hfp2009_mean_scale +
                          #ndvi_sd_scale * habitat +
                          #ndvi_sd_scale * hfp2009_mean_scale +
                          #ndvi_sd_scale * ProtectedArea + 
                          ndvi_mean_scale * ProtectedArea,  
                random= ~1|Park_location,
                method= "REML",
                data=df)
meta_speed13.2

meta_speed14 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean_scale +
                  ndvi_sd_scale + 
                  ProtectedArea +
                  woody_mean_scale + 
                  woody_sd_scale + 
                  hfp2009_mean_scale +
                  hfp2009_sd_scale +
                  ndvi_mean_scale * woody_mean_scale +
                  ndvi_mean_scale * hfp2009_mean_scale +
                  ndvi_sd_scale * woody_mean_scale +
                  ndvi_sd_scale * hfp2009_mean_scale +
                  ndvi_mean_scale * ProtectedArea + 
                  ndvi_sd_scale * ProtectedArea, 
                random= ~1|Park_location,
                method= "REML",
                data=df)

meta_speed15 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean_scale +
                  ndvi_sd_scale + 
                  ProtectedArea +
                  #habitat +
                  woody_mean_scale +
                  woody_sd_scale + 
                  hfp2009_mean_scale +
                  hfp2009_sd_scale +
                  ndvi_mean_scale * woody_mean_scale +
                  ndvi_mean_scale * hfp2009_mean_scale +
                  ndvi_sd_scale * woody_mean_scale +
                  ndvi_sd_scale * hfp2009_mean_scale +
                  ndvi_mean_scale * ProtectedArea, 
                random= ~1|Park_location,
                method= "REML",
                data=df)

meta_speed16 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean_scale +
                  ndvi_sd_scale + 
                  ProtectedArea +
                  #habitat +
                  woody_mean_scale + 
                  woody_sd_scale +
                  hfp2009_mean_scale +
                  hfp2009_sd_scale +
                  ndvi_mean_scale * woody_mean_scale +
                  ndvi_mean_scale * hfp2009_mean_scale +
                  ndvi_sd_scale * woody_mean_scale +
                  ndvi_sd_scale * hfp2009_mean_scale, 
                random= ~1|Park_location,
                method= "REML",
                data=df)

meta_speed17 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean_scale +
                  ndvi_sd_scale + 
                  ProtectedArea +
                  woody_mean_scale +
                  hfp2009_mean_scale +
                  hfp2009_sd_scale,
                random= ~1|Park_location,
                method= "REML",
                data=df)

meta_speed18 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean_scale +
                  ndvi_sd_scale + 
                  ProtectedArea +
                  woody_mean_scale +
                  hfp2009_mean_scale, 
                random= ~1|Park_location,
                method= "REML",
                data=df)

meta_speed19 = rma.mv(log_speed,var_log_speed, mods = ~ 
                  ndvi_mean_scale +
                  ndvi_sd_scale + 
                  ProtectedArea +
                  #habitat +
                  hfp2009_mean_scale +
                  woody_mean_scale+
                  woody_sd_scale, 
                random= ~1|Park_location,
                method= "REML",
                data=df)

###############################################
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
AIC(meta_speed18)
AIC(meta_speed19)

#Model Summary
top_model = meta_speed14
speed_summary = as.data.frame(top_model[2])
speed_summary <- as.data.frame(tibble::rownames_to_column(speed_summary, "Effect"))
speed_summary$se = top_model[3]
speed_summary$se = top_model$se
speed_summary$CI_lower = top_model$ci.lb 
speed_summary$CI_upper = top_model$ci.ub

#Remove the intercept from the model
speed_summary = speed_summary %>%
  dplyr::filter(!Effect=="intrcpt")
#speed_summary = speed_summary[-c(1),]
speed_summary=as.data.frame(speed_summary)

#Renaming Covariates for Plotting
speed_summary$Effect = replace(speed_summary$Effect, speed_summary$Effect =="ndvi_mean_scale","NDVI (mean)")#For replacing characters with other characters
speed_summary$Effect = replace(speed_summary$Effect, speed_summary$Effect =="woody_sd_scale","Woody Biomass (SD)")#For replacing characters with other characters
speed_summary$Effect = replace(speed_summary$Effect, speed_summary$Effect =="woody_mean_scale","Woody Biomass (mean)")#For replacing characters with other characters
speed_summary$Effect = replace(speed_summary$Effect, speed_summary$Effect =="ProtectedArea","Protected Area Overlap")#For replacing characters with other characters
speed_summary$Effect = replace(speed_summary$Effect, speed_summary$Effect =="ndvi_sd_scale","NDVI (SD)")#For replacing characters with other characters
speed_summary$Effect = replace(speed_summary$Effect, speed_summary$Effect =="hfp2009_sd_scale","HFI (SD)")#For replacing characters with other characters
speed_summary$Effect = replace(speed_summary$Effect, speed_summary$Effect =="hfp2009_mean_scale","HFI (Mean)")#For replacing characters with other characters

#Generating a quick plot to visualize effect sizes as a forest plot
ggplot(data=speed_summary,aes(x=beta, y= Effect))+
  geom_point()+
  geom_vline(xintercept = 0, linetype="dotted", color="red", size=1.5)+
  geom_errorbar(aes(xmin = CI_lower, xmax =CI_upper))+
  xlab("Effect Size on Daily Speed")+
  ylab("Covariate")+
  theme_minimal()+
  #xlim(0,600) #excluding point estimates above 600km2 in the graph just so we can visually inspect the smaller sizes
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
ggtitle("Movement Rate")

##Developing Variance Componenets for Meta13.2
#SEE THIS HELP PAGE: https://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
#How much Variance is due to heterogeneity?
W <- diag(1/meta_speed13.2$vi)
X <- model.matrix(meta_speed13.2)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
100 * sum(meta_speed13.2$sigma2) / (sum(meta_speed13.2$sigma2) + (meta_speed13.2$k-meta_speed13.2$p)/sum(diag(P)))

#How much variacen is du to between cluster heterogeneity
100 * meta_speed13.2$sigma2 / (sum(meta_speed13.2$sigma2) + (meta_speed13.2$k-meta_speed13.2$p)/sum(diag(P)))

