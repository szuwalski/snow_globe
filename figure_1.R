#==loking at NBS by size
#++combine the NBS and EBS data from 2017-2019
#==have a multipanel plot with columns of different slices of the population

library(maps)
library("rnaturalearth")
library(interp)
library(RColorBrewer)
library(reshape2) # for melt
library(mgcv)  
library(PBSmapping)
library(mapdata)    #some additional hires data
library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(ggplot2)
library(patchwork)
library(dplyr)
#==============================
# EBS data
#==============================
survDAT<-read.csv("C:/Users/cody.szuwalski/Work/snow_2023_9/data/survey/EBSCrab_Haul/EBSCrab_Haul.csv",header=T,skip=5)
survDAT<-filter(survDAT,HAUL_TYPE==3&nchar(GIS_STATION)<5)
drvYear<-as.numeric(substr(survDAT$CRUISE,1,4))
SurvYR<-unique(drvYear)
AllStation<-unique(survDAT$GIS_STATION)

#==plot GIS stations
AllStnLoc<-matrix(ncol=2,nrow=length(AllStation))
for(w in 1:length(AllStation))
{
  temp<-survDAT[survDAT$GIS_STATION==AllStation[w],]
  AllStnLoc[w,1]<-temp$MID_LATITUDE[1]
  AllStnLoc[w,2]<-temp$MID_LONGITUDE[1]
}

# plot(AllStnLoc[,1]~AllStnLoc[,2],cex=.3)
# text(AllStation,x=AllStnLoc[,2],y=AllStnLoc[,1])
# nbs_bound<-data.frame(lat=c(62.1,62.1,61.1,61.1,60.5,60.5),
#                       lon=c(-175,-172.5,-172.5,-171,-171,-165))
#========================
#==find survey densities (total)
#========================
nmiSurv<-140350
DensityM_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityF_45_55<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_55_65<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_65_75<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM_78_100<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityM101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityMge45<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
DensityF_mat<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

num_M_45_85<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
num_M_101<-matrix(nrow=length(SurvYR),ncol=length(AllStation))

StationYr<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
station_bot_temp<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
station_depth<-matrix(nrow=length(SurvYR),ncol=length(AllStation))
for(y in 1:length(SurvYR))
{
  yrDAT<-survDAT[drvYear==SurvYR[y],]
  fileyr<-SurvYR[y]
  stationsUNQ<-(unique(yrDAT$GIS_STATION))
  #==density at station
  for(j in 1:length(stationsUNQ))
  {
    stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j],]
    StationYr[y,j]<-as.character(stationsUNQ[j])
    Hauls<-(unique(stationALL$HAUL))
    female_45_55<-stationALL[stationALL$SEX==2 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    female_mat<-stationALL[stationALL$SEX==2 &  stationALL$EGG_COLOR>0,]
    male_45_55<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_55_65<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>55 & stationALL$WIDTH<66,]
    male_65_75<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>65 & stationALL$WIDTH<76,]
    male_45_85<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<86,]
    male_78_100<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>78 & stationALL$WIDTH<101,]
    Males101<-stationALL[stationALL$SEX==1& stationALL$WIDTH>101,]
    male_ge45<-stationALL[stationALL$SEX==1& stationALL$WIDTH>44,]
    #==densities across hauls in crabs per km^2
    tempDensM45<-NULL
    tempDensF45<-NULL
    tempDensM55<-NULL
    tempDensM65<-NULL
    tempDensM85<-NULL
    tempDensM75101<-NULL
    tempDensM101<-NULL
    tempDensF<-NULL
    tempDensFmat<-NULL
    temp_num_45<-0
    temp_num_101<-0
    tempDensMge45<-NULL
    for(k in 1:length(Hauls))
    {
      
      SampFactM<-female_mat$SAMPLING_FACTOR[which(female_mat$HAUL==Hauls[k])[1]] 
      AreaSweptM<-female_mat$AREA_SWEPT[which(female_mat$HAUL==Hauls[k])[1]]
      tempDensFmat<-length(female_mat$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-female_45_55$SAMPLING_FACTOR[which(female_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-female_45_55$AREA_SWEPT[which(female_45_55$HAUL==Hauls[k])[1]]
      tempDensF45<-length(female_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_55$SAMPLING_FACTOR[which(male_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_55$AREA_SWEPT[which(male_45_55$HAUL==Hauls[k])[1]]
      tempDensM45<-length(male_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_55_65$SAMPLING_FACTOR[which(male_55_65$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_55_65$AREA_SWEPT[which(male_55_65$HAUL==Hauls[k])[1]]
      tempDensM55<-length(male_55_65$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_65_75$SAMPLING_FACTOR[which(male_65_75$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_65_75$AREA_SWEPT[which(male_65_75$HAUL==Hauls[k])[1]]
      tempDensM65<-length(male_65_75$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_85$SAMPLING_FACTOR[which(male_45_85$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_85$AREA_SWEPT[which(male_45_85$HAUL==Hauls[k])[1]]
      tempDensM85<-length(male_45_85$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_45<-temp_num_45+length(Males101$HAUL==Hauls[k])*SampFactM
      
      SampFactM<-male_78_100$SAMPLING_FACTOR[which(male_78_100$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_78_100$AREA_SWEPT[which(male_78_100$HAUL==Hauls[k])[1]]
      tempDensM78101<-length(male_78_100$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_ge45$SAMPLING_FACTOR[which(male_ge45$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_ge45$AREA_SWEPT[which(male_ge45$HAUL==Hauls[k])[1]]
      tempDensMge45<-length(male_ge45$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-Males101$SAMPLING_FACTOR[which(Males101$HAUL==Hauls[k])[1]] 
      AreaSweptM<-Males101$AREA_SWEPT[which(Males101$HAUL==Hauls[k])[1]]
      tempDensM101<-length(Males101$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_101<-temp_num_101+length(Males101$HAUL==Hauls[k])*SampFactM
    }
    DensityF_45_55[y,j]<-mean(tempDensF45)
    DensityM_45_55[y,j]<-mean(tempDensM45)
    DensityM_55_65[y,j]<-mean(tempDensM55)
    DensityM_65_75[y,j]<-mean(tempDensM65)
    DensityM_45_85[y,j]<-mean(tempDensM85)
    DensityM_78_100[y,j]<-mean(tempDensM78101)
    DensityMge45[y,j]<-mean(tempDensMge45)
    DensityM101[y,j]<-mean(tempDensM101)
    num_M_45_85[y,j]<-temp_num_45    
    num_M_101[y,j]<-temp_num_101
    DensityF_mat[y,j]<-mean(tempDensFmat)
    station_bot_temp[y,j]<-stationALL$GEAR_TEMPERATURE[1]
    station_depth[y,j]<-stationALL$BOTTOM_DEPTH[1]
  }
}
plot(apply(num_M_101,1,sum,na.rm=T)[-1]~SurvYR[-1],type='l',ylab='total number of observed crab >101mm')

grr<-DensityM101
grr[32,95]<-NA
plot(apply(grr,1,sum,na.rm=T)[-1]~SurvYR[-1],type='l',ylab='total density of observed crab >101mm')


#================================
# Northern Bering Sae
#================================
#nbs_survDAT<-read.csv("C:/data/NBS crab data/nbs_opilio_crabhaul.csv",header=T)
#nbs_survDAT<-read.csv("data/CRABHAUL_OPILIO_NBS.csv",header=T)
#survDAT<-read.csv("C:/Users/cody.szuwalski/Work/snow_2023_9/data/survey/EBSCrab_Haul/EBSCrab_Haul.csv",header=T,skip=5)
nbs_survDAT<-read.csv("C:/data/NBS crab data/CRABHAUL_OPILIO_NBS.csv",header=T)
nbs_drvYear<-as.numeric(substr(nbs_survDAT$CRUISE,1,4))
nbs_survDAT$AKFIN_SURVEY_YEAR <-as.numeric(substr(nbs_survDAT$CRUISE,1,4))
nbs_SurvYR<-unique(nbs_drvYear)
nbs_AllStation<-unique(nbs_survDAT$GIS_STATION)

colnames(nbs_survDAT)

#==plot GIS stations
nbs_AllStnLoc<-matrix(ncol=2,nrow=length(nbs_AllStation))
for(w in 1:length(nbs_AllStation))
{
  temp<-nbs_survDAT[nbs_survDAT$GIS_STATION==nbs_AllStation[w],]
  nbs_AllStnLoc[w,1]<-temp$MID_LATITUDE[1]
  nbs_AllStnLoc[w,2]<-temp$MID_LONGITUDE[1]
}

nbs_num_M_45_85<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_num_M_101<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))

nbs_DensityF_45_55<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityF_mat<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_45_55<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_55_65<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_65_75<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_45_85<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_78_100<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM101<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_StationYr<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_station_bot_temp<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_station_depth<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
nbs_DensityM_ge45<-matrix(nrow=length(nbs_SurvYR),ncol=length(nbs_AllStation))
for(y in 1:length(nbs_SurvYR))
{
  yrDAT<-nbs_survDAT[nbs_drvYear==nbs_SurvYR[y],]
  fileyr<-nbs_SurvYR[y]
  stationsUNQ<-(unique(yrDAT$GIS_STATION))
  #==density at station
  for(j in 1:length(stationsUNQ))
  {
    stationALL<-yrDAT[yrDAT$GIS_STATION==stationsUNQ[j],]
    nbs_StationYr[y,j]<-as.character(stationsUNQ[j])
    Hauls<-(unique(stationALL$HAUL))
    female_45_55<-stationALL[stationALL$SEX==2 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    female_mat<-stationALL[stationALL$SEX==2 &  stationALL$EGG_COLOR>0,]
    male_45_55<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<56,]
    male_55_65<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>55 & stationALL$WIDTH<66,]
    male_65_75<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>65 & stationALL$WIDTH<76,]
    male_45_85<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>44 & stationALL$WIDTH<86,]
    male_78_100<-stationALL[stationALL$SEX==1 &  stationALL$WIDTH>78 & stationALL$WIDTH<101,]
    Males101<-stationALL[stationALL$SEX==1& stationALL$WIDTH>101,]
    male_ge45<-stationALL[stationALL$SEX==1& stationALL$WIDTH>45,]
    
    #==densities across hauls in crabs per km^2
    tempDensF45<-NULL
    tempDensM45<-NULL
    tempDensM55<-NULL
    tempDensM65<-NULL
    tempDensM85<-NULL
    tempDensM75101<-NULL
    tempDensM101<-NULL
    tempDensF<-NULL
    tempDensFmat<-NULL
    temp_num_45<-0
    temp_num_101<-0
    tempDensMge45<-NULL
    
    for(k in 1:length(Hauls))
    {
      SampFactM<-female_mat$SAMPLING_FACTOR[which(female_mat$HAUL==Hauls[k])[1]] 
      AreaSweptM<-female_mat$AREA_SWEPT[which(female_mat$HAUL==Hauls[k])[1]]
      tempDensFmat<-length(female_mat$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-female_45_55$SAMPLING_FACTOR[which(female_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-female_45_55$AREA_SWEPT[which(female_45_55$HAUL==Hauls[k])[1]]
      tempDensF45<-length(female_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_55$SAMPLING_FACTOR[which(male_45_55$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_55$AREA_SWEPT[which(male_45_55$HAUL==Hauls[k])[1]]
      tempDensM45<-length(male_45_55$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_55_65$SAMPLING_FACTOR[which(male_55_65$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_55_65$AREA_SWEPT[which(male_55_65$HAUL==Hauls[k])[1]]
      tempDensM55<-length(male_55_65$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_65_75$SAMPLING_FACTOR[which(male_65_75$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_65_75$AREA_SWEPT[which(male_65_75$HAUL==Hauls[k])[1]]
      tempDensM65<-length(male_65_75$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_45_85$SAMPLING_FACTOR[which(male_45_85$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_45_85$AREA_SWEPT[which(male_45_85$HAUL==Hauls[k])[1]]
      tempDensM85<-length(male_45_85$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_45<-temp_num_45+length(Males101$HAUL==Hauls[k])*SampFactM
      
      SampFactM<-male_78_100$SAMPLING_FACTOR[which(male_78_100$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_78_100$AREA_SWEPT[which(male_78_100$HAUL==Hauls[k])[1]]
      tempDensM78101<-length(male_78_100$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-male_ge45$SAMPLING_FACTOR[which(male_ge45$HAUL==Hauls[k])[1]] 
      AreaSweptM<-male_ge45$AREA_SWEPT[which(male_ge45$HAUL==Hauls[k])[1]]
      tempDensMge45<-length(male_ge45$HAUL==Hauls[k])*SampFactM/AreaSweptM
      
      SampFactM<-Males101$SAMPLING_FACTOR[which(Males101$HAUL==Hauls[k])[1]] 
      AreaSweptM<-Males101$AREA_SWEPT[which(Males101$HAUL==Hauls[k])[1]]
      tempDensM101<-length(Males101$HAUL==Hauls[k])*SampFactM/AreaSweptM
      temp_num_101<-temp_num_101+length(Males101$HAUL==Hauls[k])*SampFactM
    }
    nbs_DensityF_mat[y,j]<-mean(tempDensFmat)
    nbs_DensityF_45_55[y,j]<-mean(tempDensF45)
    nbs_DensityM_45_55[y,j]<-mean(tempDensM45)
    nbs_DensityM_55_65[y,j]<-mean(tempDensM55)
    nbs_DensityM_65_75[y,j]<-mean(tempDensM65)
    nbs_DensityM_45_85[y,j]<-mean(tempDensM85)
    nbs_DensityM_78_100[y,j]<-mean(tempDensM78101)
    nbs_DensityM101[y,j]<-mean(tempDensM101)
    nbs_station_bot_temp[y,j]<-stationALL$GEAR_TEMPERATURE[1]
    nbs_station_depth[y,j]<-stationALL$BOTTOM_DEPTH[1]
    nbs_num_M_45_85[y,j]<-temp_num_45    
    nbs_num_M_101[y,j]<-temp_num_101
    nbs_DensityM_ge45[y,j]<-mean(tempDensMge45)
    
  }
}

#===============================
# mappity map map data
#===============================
ebs_dat<-data.frame(log_abund_101=log(c(t(DensityM101[-1,]))),
                    fish_recruits=log(c(t(DensityM_78_100[-1,]))),
                    bot_temp=c(t(station_bot_temp[-nrow(station_bot_temp),])),
                    depth=c(t(station_depth[-nrow(station_depth),])),
                    station=c(t(StationYr[-1,])),
                    lat=rep(AllStnLoc[,1],(nrow(DensityM_78_100)-1)),
                    lon=rep(AllStnLoc[,2],(nrow(DensityM_78_100)-1)),
                    year=rep(SurvYR[-1],each=ncol(DensityM101)),
                    log_abund_101=log(c(t(DensityM101[-1,]))),
                    log_abund_45=log(c(t(DensityM_45_55[-1,]))),
                    log_abund_45_f=log(c(t(DensityF_45_55[-1,]))),
                    log_abund_55=log(c(t(DensityM_55_65[-1,]))),
                    log_abund_65=log(c(t(DensityM_65_75[-1,]))),
                    log_abund_85=log(c(t(DensityM_45_85[-1,]))),
                    log_abund_ge45=log(c(t(DensityMge45[-1,]))),
                    loc="EBS",
                    obs_num_101=c(t(num_M_101[-1,])),
                    obs_num_45_85=c(t(num_M_45_85[-1,])),
                    log_abund_45_mat=log(c(t(DensityF_mat[-1,]))))

for(x in 1:length(AllStation))
{
  ebs_dat$lat[ebs_dat$station==AllStation[x]]<-AllStnLoc[x,1]
  ebs_dat$lon[ebs_dat$station==AllStation[x]]<-AllStnLoc[x,2]
}

nbs_dat<-data.frame(log_abund_101=log(c(t(nbs_DensityM101))),
                    fish_recruits=log(c(t(nbs_DensityM_78_100))),
                    bot_temp=c(t(nbs_station_bot_temp)),
                    depth=c(t(nbs_station_depth)),
                    station=c(t(nbs_StationYr)),
                    lat=rep(nbs_AllStnLoc[,1],(nrow(nbs_DensityM_78_100))),
                    lon=rep(nbs_AllStnLoc[,2],(nrow(nbs_DensityM_78_100))),
                    year=rep(nbs_SurvYR,each=ncol(nbs_DensityM101)),
                    log_abund_101=log(c(t(nbs_DensityM101))),
                    log_abund_45=log(c(t(nbs_DensityM_45_55))),
                    log_abund_45_f=log(c(t(nbs_DensityF_45_55))),
                    log_abund_55=log(c(t(nbs_DensityM_55_65))),
                    log_abund_65=log(c(t(nbs_DensityM_65_75))),
                    log_abund_85=log(c(t(nbs_DensityM_45_85))),
                    log_abund_ge45=log(c(t(nbs_DensityM_ge45))),
                    loc='NBS',
                    obs_num_101=c(t(nbs_num_M_101)),
                    obs_num_45_85=c(t(nbs_num_M_45_85)),
                    log_abund_45_mat=log(c(t(nbs_DensityF_mat))))

for(x in 1:length(nbs_AllStation))
{
  nbs_dat$lat[nbs_dat$station==nbs_AllStation[x]]<-nbs_AllStnLoc[x,1]
  nbs_dat$lon[nbs_dat$station==nbs_AllStation[x]]<-nbs_AllStnLoc[x,2]
}

ebs_dat$surv<-"EBS"
nbs_dat$surv<-"NBS"
in_dat<-rbind(ebs_dat,nbs_dat)
world <- ne_countries(scale = "medium", returnclass = "sf")
lon_1<-min(in_dat$lon,na.rm=T)
lon_2<-max(in_dat$lon,na.rm=T)*.99
lat_1<-min(in_dat$lat,na.rm=T)
lat_2<-max(in_dat$lat,na.rm=T)

lon_1<- -45
lon_2<- 140
lat_1<-42
lat_2<-67

plotted<-ggplot() + 
  geom_tile(data=filter(in_dat,year==2019&!is.na(log_abund_ge45)), aes(x = lon, y = lat, fill = -1*log_abund_ge45 ),width=.5,height=.25) +
  scale_fill_distiller(palette="Blues",direction= -1, na.value="white") +
  geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA),legend.position='none')+labs(fill="Depth (m)")

#==need to specify the color for region, scale the density by region,
#==then use the scaled color as alpha

lg_ml<-read.csv("data/comp_region_lg_males.csv")
reg_col<-c("#F8766D","#619CFF")
lg_ml$incol[lg_ml$region=="Newfoundland"]<-reg_col[1]
lg_ml$incol[lg_ml$region=="Bering Sea"]<-reg_col[2]

NL_bot_tmp<-read.csv("data/btemps_long.csv")
use_nl_tmp<-filter(NL_bot_tmp,year==2020)%>%
  group_by(lat,long)%>%
  summarize(value=mean(btemp))

#=bottom temps
nl_tmp<-NL_bot_tmp%>%
  group_by(year)%>%
  summarize(tmp=mean(btemp))
nl_tmp$region<-"Newfoundland"
nl_tmp$tmp<-(nl_tmp$tmp)
ak_dat<-read.csv("data/AK_m_gam_dat.csv")

mg_ak<-ak_dat[,c(2,13)]
colnames(mg_ak)<-c("year","tmp")
mg_ak$tmp<-(mg_ak$tmp)
mg_ak$region<-"Bering Sea"

plot_tmp<-rbind(mg_ak,nl_tmp)

bot_temps<-ggplot(plot_tmp)+
  geom_point(aes(x=year,y=tmp,col=region,group=region),size=3)+
  geom_line(aes(x=year,y=tmp,col=region,group=region),lwd=1.2)+
  theme_bw()+ylab("")+xlab("")+
  theme(legend.position='none',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background=element_blank(),
        legend.background=element_rect(fill=rgb(1,1,1,0)),
        legend.key=element_rect(fill=rgb(1,1,1,0)),
        panel.border = element_blank(),
        plot.background=element_blank())+ggtitle("Bottom temperature")


ebs_dat<-filter(in_dat,year==2019&!is.na(log_abund_ge45))[,c(6,7,15)]
ebs_dat$log_abund_ge45<-ebs_dat$log_abund_ge45/max(ebs_dat$log_abund_ge45)
colnames(ebs_dat)<-c("lat","long","value")
ebs_dat$col<-reg_col[2]
use_nl_tmp$value<-use_nl_tmp$value/max(use_nl_tmp$value)
use_nl_tmp$col<-reg_col[1]
new_dat<-rbind(use_nl_tmp,ebs_dat)


lon_1<-min(in_dat$lon,na.rm=T)
lon_2<- -47
lat_1<-42
lat_2<-67

plotted<-ggplot() + 
  geom_tile(data=new_dat, aes(x = long, y = lat, fill = col,alpha=value ),width=.5,height=.25) +
  scale_color_identity()+geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+xlab("")+ylab("")+
  theme(axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, size = 10),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA),legend.position='none')+labs(fill="Depth (m)")


scl_lg_dat<-ggplot(lg_ml)+
  geom_point(aes(x=years,y=values,col=region,group=region),size=3)+
  geom_line(aes(x=years,y=values,col=region,group=region),lwd=1.2)+
  theme_bw()+ylab("")+xlab("")+
  theme(legend.position=c(.65,.75),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background=element_blank(),
        legend.background=element_rect(fill=rgb(1,1,1,0)),
        legend.key=element_rect(fill=rgb(1,1,1,0)),
        panel.border = element_blank(),
        plot.background=element_blank())+ggtitle("Large male abundance")





library(grid)
library(patchwork)
sm_grob<-ggplotGrob(scl_lg_dat)
sm_grob2<-ggplotGrob(bot_temps)
figure_1<-plotted + annotation_custom(grob=sm_grob,
                            xmin=-127,
                            xmax=-65,
                            ymin=42,
                            ymax=67)+
  annotation_custom(grob=sm_grob2,
                    xmin=-178,
                    xmax=-125,
                    ymin=42,
                    ymax=55.5)

png("plots/figure_1.png",res=300,units='in',height=4,width=12)
print(figure_1)
dev.off()

design<-"23
         11
         11"
png("plots/figure_1_alt.png",res=300,units='in',height=7,width=12)
plotted + scl_lg_dat + bot_temps + plot_layout(design=design)
dev.off()



lon_1<--180
lon_2<- 180
lat_1<-30
lat_2<-85

plotted2<-ggplot() + 
  geom_tile(data=new_dat, aes(x = long, y = lat, fill = col,alpha=value ),width=.5,height=.25) +
  scale_color_identity()+geom_sf(data=world) +
  coord_sf(xlim = c(lon_1,lon_2), ylim = c(lat_1,lat_2), expand = FALSE) +
  theme_bw()+xlab("")+ylab("")+
  theme(axis.text.x = element_blank(), #element_text(angle = 90, vjust = 0.5, size = 10),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(margin= margin(1,0,1,0)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(color="white",fill="white"))+
  theme(legend.background = element_rect(fill='transparent',color=NA),
        legend.box.background = element_rect(fill='transparent',color=NA),legend.position='none')+labs(fill="Depth (m)")
png("plots/global.png",res=300,units='in',height=4,width=12)
print(plotted2)
dev.off()
