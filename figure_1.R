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
#library(maptools)   #useful tools such as reading shapefiles
library(mapproj)
library(ggplot2)
library(patchwork)
library(dplyr)
#devtools::install_github("AFSC-Shellfish-Assessment-Program/crabpack")
library(crabpack)
channel <- "API"
specimen_data <- crabpack::get_specimen_data(species = "SNOW",
                                             region = "EBS",
                                             years = c(1975:2025),
                                             channel = channel)

specimen_data_nbs <- crabpack::get_specimen_data(species = "SNOW",
                                             region = "NBS",
                                             years = c(1975:2025),
                                             channel = channel)

#==make a spatial foot print file
cpue <- crabpack::calc_cpue(crab_data = specimen_data,
                            species = "SNOW",
                            region = "EBS",
                            district = "ALL",
                            years = 2019,
                            crab_category = "all_categories")

cpue_NBS <- crabpack::calc_cpue(crab_data = specimen_data_nbs,
                            species = "SNOW",
                            region = "NBS",
                            district = "ALL",
                            years = 2019,
                            crab_category = "all_categories")



tot_sp_pop<-cpue%>%
  group_by(YEAR,LATITUDE,LONGITUDE)%>%
  summarize(tot_cpue=sum(CPUE))

tot_sp_pop_nbs<-cpue_NBS%>%
  group_by(YEAR,LATITUDE,LONGITUDE)%>%
  summarize(tot_cpue=sum(CPUE))
  
#==colors for different regions
reg_col<-c("#F8766D","#619CFF","#00BA39")
reg_col<-c("#9900ff","#9900ff","#9900ff")
no_data_col<-"#9900ff"
no_dat_val<-1
in_dat<-rbind(tot_sp_pop,tot_sp_pop_nbs)
in_dat$col<-reg_col[2]
colnames(in_dat)<-c('year','lat','long','value','col')
write.csv(in_dat,"data/AK_spatial_footprint.csv")


#==global map
#==add survey stations with no data
#==(and then go persuade folks to collaborate!)
#==okhost
okh_loc<-read.csv('data/spatial/okhost_loc.csv')
okh_loc<-okh_loc[,c(2,1)]
okh_loc$value<-no_dat_val
okh_loc$col<-no_data_col
colnames(okh_loc)<-c('lat','long','value','col')
new_dat2<-rbind(in_dat,okh_loc)

#==Barents
lats<-seq(68,82,0.66)
lons<-seq(10,55,10/5)

fake_dat<-expand.grid(lats,lons)
colnames(fake_dat)<-c("lat","long")
fake_dat$value<-no_dat_val
fake_dat$col<-no_data_col


lon_1<--180
lon_2<- 180
lat_1<-30
lat_2<-85

for(x in 1:nrow(fake_dat))
{
  if(fake_dat[x,1]>79 & fake_dat[x,2]>40 |fake_dat[x,1]>72 & fake_dat[x,2]>54 |
     fake_dat[x,1]<76 & fake_dat[x,2]<16 )
  {
    fake_dat[x,3]<-NA
    fake_dat[x,4]<-0
  }
}


colnames(fake_dat)<-c('lat','long','value','col')
new_dat2<-rbind(new_dat2,fake_dat)


#==japan West
jp_loc<-read.csv('data/spatial/JP_west_2020.csv')
jp_loc<-jp_loc[,c(2,1)]
jp_loc$value<-no_dat_val
jp_loc$col<-no_data_col
colnames(jp_loc)<-c('lat','long','value','col')
new_dat2<-rbind(new_dat2,jp_loc)


#==japan east
jp_loc<-read.csv('data/spatial/JP_east_2020.csv')
jp_loc<-jp_loc[,c(2,1)]
jp_loc$value<-no_dat_val
jp_loc$col<-no_data_col
colnames(jp_loc)<-c('lat','long','value','col')
new_dat2<-rbind(new_dat2,jp_loc)

#==korea
kor_loc<-read.csv('data/spatial/korea_fish.csv')
kor_loc<-kor_loc[,c(2,1)]
kor_loc$value<-no_dat_val
kor_loc$col<-no_data_col
colnames(kor_loc)<-c('lat','long','value','col')
new_dat2<-rbind(new_dat2,kor_loc)


#==western Bering Sea
wbs_loc<-read.csv('data/spatial/west_bering_1990.csv')
wbs_loc<-wbs_loc[,c(2,1)]
wbs_loc$value<-no_dat_val
wbs_loc$col<-no_data_col
wbs_loc[wbs_loc[,2]<180,2]<- -1*wbs_loc[wbs_loc[,2]<180,2]
wbs_loc[wbs_loc[,2]>180,2]<- 180-(wbs_loc[wbs_loc[,2]>180,2]-180)
colnames(wbs_loc)<-c('lat','long','value','col')
new_dat2<-rbind(new_dat2,wbs_loc)

#==add GSL
gsl_loc<-as.data.frame(read.csv("data/spatial/gsl_avg_dens_spatial.csv",check.names=FALSE))
colnames(gsl_loc)[1]<-"long"
gsl_dat<-melt(gsl_loc,id.vars=c("long"))
gsl_dat<-gsl_dat[,c(2,1,3)]
gsl_dat$value<-gsl_dat$value/max(gsl_dat$value,na.rm=T)
gsl_dat$col<-reg_col[3]
gsl_dat<-data.frame(lat=as.numeric(as.character(unlist(gsl_dat[,1]))),
                    long=as.numeric(unlist(gsl_dat[,2])),
                    value=unlist(gsl_dat[,3]),
                    col=unlist(gsl_dat[,4]))

ugh<-unique(gsl_dat$long)
tmp_gsl<-filter(gsl_dat,lat%in%seq(45,49,by=0.2))
ugh<-unique(tmp_gsl$long)
takers<-seq(1,length(ugh),length.out=20)
tmp_gsl<-filter(tmp_gsl,long%in%ugh[takers])
tmp_gsl<-tmp_gsl[!is.na(tmp_gsl$value),]
new_dat2<-rbind(new_dat2,tmp_gsl)

#==add NFL
#==FIX THIS, CURRENTLY BASED ON TEMPERATURE
NL_bot_tmp<-read.csv("data/spatial/NFL_btemps_long.csv")
use_nl_tmp<-filter(NL_bot_tmp,year==2020)%>%
  group_by(lat,long)%>%
  summarize(value=mean(btemp))
use_nl_tmp$value<-use_nl_tmp$value/max(use_nl_tmp$value)
use_nl_tmp$col<-reg_col[1]
new_dat2<-rbind(use_nl_tmp,new_dat2)

lon_1<--180
lon_2<- 180
lat_1<-30
lat_2<-85
world <- ne_countries(scale = "medium", returnclass = "sf")
plotted2<-ggplot() + 
  geom_tile(data=new_dat2, aes(x = long, y = lat, fill = col,alpha=log(value+0.0001) ),width=.5,height=.25) +
  scale_fill_identity()+geom_sf(data=world) +
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
#png("plots/global.png",res=300,units='in',height=4,width=12)

layout <- "
A##
BBB
BBB
"
png("plots/global_2.png",res=300,units='in',height=4,width=12)
plotted2 + inset_element(fao_catch, left = 0.37, bottom = 0.0, right = .83, top = .6)
dev.off()
