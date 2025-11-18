library(ggplot2)
library(dplyr)
library(reshape2)

fao_dat<-read.csv("C:\\Users\\Cody\\OneDrive\\Desktop\\Capture_2025.1.0\\Capture_Quantity.csv")
cnt_dat<-read.csv("C:\\Users\\Cody\\OneDrive\\Desktop\\Capture_2025.1.0\\CL_FI_COUNTRY_GROUPS.csv")
spp_dat<-read.csv("C:\\Users\\Cody\\OneDrive\\Desktop\\Capture_2025.1.0\\CL_FI_SPECIES_GROUPS.csv")
chiono_dat<-filter(fao_dat,SPECIES.ALPHA_3_CODE%in%c("CRQ","CWJ","CVB","PCR"))

#==chionocetes spp = CRQ (queen crab--opilio), CWJ (red snow crab--japonicus, CVB (Tanner --bairdi), PCR (Chionoecetes spp))
opie_ct<-unique(chiono_dat$COUNTRY.UN_CODE)
cnt_dat[match(opie_ct,cnt_dat$UN_Code),5]

chiono_dat$country<-cnt_dat[match(chiono_dat$COUNTRY.UN_CODE,cnt_dat$UN_Code),5]

ggplot(chiono_dat, aes(x = PERIOD, y = VALUE, fill = as.factor(country))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Stacked Bar Plot Example",
       x = "Main Category",
       y = "Value",
       fill = "Sub-Category") +
  theme_bw()+
  facet_wrap(~SPECIES.ALPHA_3_CODE)


all_chio<-filter(chiono_dat,!country%in%c("Ireland","Latvia","Lithuania","Saint Pierre and Miquelon","Spain"))%>%
  group_by(country,PERIOD)%>%
  summarize(tot_catch=sum(VALUE))

fao_catch<-ggplot(all_chio, aes(x = PERIOD, y = tot_catch/1000, fill = as.factor(country))) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_brewer(palette="Spectral")+
  labs( x = "Year",
       y = "Catch (1,000 t)",
       fill = "Country") +
  theme_bw()+theme(plot.background = element_rect(fill = "#FFFFFF11", color = NA),
                   legend.background = element_rect(fill = "#FFFFFF66", color = NA) # Set fill color and alpha
  )



#==============================
# king crab
#========================


yrg<-filter(spp_dat,ISSCAAP_Group_En=="King crabs, squat-lobsters")$X3A_Code 
king_dat<-filter(fao_dat,SPECIES.ALPHA_3_CODE%in%yrg)

#==chionocetes spp = CRQ (queen crab--opilio), CWJ (red snow crab--japonicus, CVB (Tanner --bairdi), PCR (Chionoecetes spp))
king_ct<-unique(king_dat$COUNTRY.UN_CODE)
cnt_dat[match(king_ct,cnt_dat$UN_Code),5]

king_dat$country<-cnt_dat[match(king_dat$COUNTRY.UN_CODE,cnt_dat$UN_Code),5]

all_king<-filter(king_dat)%>%
  group_by(country,PERIOD)%>%
  summarize(tot_catch=sum(VALUE))

top_10<-all_king%>%
  group_by(country)%>%
  summarize(mean_cat=mean(tot_catch))%>%
  arrange(desc(mean_cat))

fao_catch_king<-ggplot(filter(all_king,country%in%top_10$country[1:10]), aes(x = PERIOD, y = tot_catch/1000, fill = as.factor(country))) +
  geom_bar(stat = "identity", position = "stack") +
  labs( x = "Year",
        y = "Catch (1,000 t)",
        fill = "Country") +
  theme_bw()+theme(plot.background = element_rect(fill = "#FFFFFF11", color = NA),
                   legend.background = element_rect(fill = "#FFFFFF66", color = NA) # Set fill color and alpha
  )


king_spp_dat<-filter(king_dat,country%in%top_10$country[1:10])
ggplot(king_spp_dat, aes(x = PERIOD, y = VALUE, fill = as.factor(country))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Stacked Bar Plot Example",
       x = "Main Category",
       y = "Value",
       fill = "Sub-Category") +
  theme_bw()+
  facet_wrap(~SPECIES.ALPHA_3_CODE)
