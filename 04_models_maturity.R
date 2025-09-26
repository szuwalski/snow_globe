library(reshape2)
library(ggplot2)
library(mgcv)  
library(dplyr)
library(ggridges)
library(png)
library(PBSmodelling)
library(patchwork)
library(mgcv)
library(png)
library(grid)
library(DHARMa)
library(itsadug)
library(mgcViz)
library(directlabels)

pmat_all<-NULL
use_stocks<-c("Bering Sea","Newfoundland","Gulf of St Lawrence")
for(x in 1:length(outs_in))
{
  tmp<-outs_in[[x]]$prob_term_molt  
  rownames(tmp)<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  colnames(tmp)<-outs_in[[x]]$sizes
  tmps2<-melt(tmp)
  tmps2$region<-use_stocks[x]
  colnames(tmps2)<-c("Year","Size","probability","Region")
  
  #==get densities in  
  use_mat<-outs_in[[x]]$`pred_mat_pop_num`
  use_imm<-outs_in[[x]]$`pred_imm_pop_num`
  
  rownames(use_mat)<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  rownames(use_imm)<-seq(outs_in[[x]]$styr,outs_in[[x]]$endyr)
  colnames(use_mat)<-outs_in[[x]]$sizes
  colnames(use_imm)<-outs_in[[x]]$sizes

  all_crab<-use_mat+use_imm
  large_crab<-apply(all_crab[,as.numeric(colnames(all_crab))>95],1,sum)
  small_crab<-apply(all_crab[,as.numeric(colnames(all_crab))<95],1,sum) 
  
  sm_df<-data.frame(Year=names(small_crab),
                    sm_abund=small_crab/100000000)
  lg_df<-data.frame(Year=names(large_crab),
                    lg_abund=large_crab/10000000)  

  tmp_big<-merge(tmps2,sm_df,by="Year")
  tmp_big<-merge(tmp_big,lg_df,by="Year")  
  
  pmat_all<-rbind(pmat_all,tmp_big)
}

med_mat<-pmat_all%>%
  group_by(Size,Region)%>%
  summarize(med_mat=median(probability))

alldat<-ggplot()+
  geom_line(data=pmat_all,aes(x=Size,y=probability,group=interaction(Year,Region),col=Region),alpha=0.1,lwd=1.1)+
  geom_line(data=med_mat,aes(x=Size,y=med_mat,group=(Region),col=Region),lwd=1.3)+
  theme_bw()+
  theme(legend.position=c(.3,.85))+
  ylab("p(maturing)")+
  xlab("Carapace width (mm)")

save_dat<-NULL
for(x in 1:length(use_stocks))
{
  gam_dat<-filter(pmat_all,Region==use_stocks[x])
  if(x==1)
  {
   ref<-filter(gam_dat,Year==1982,Size==67.5)
   comp<-filter(gam_dat,Size==67.5)
   removers<-comp$Year[which(!is.na(match(comp$probability,ref$probability)))]
   gam_dat<-filter(gam_dat,!Year%in%removers)
  }
  mod<-gam(probability~s(Size)+s(lg_abund,sm_abund,k=9),data = gam_dat,family =betar(link='logit'))
  summary(mod)
  #plot(mod,pages=1,too.far=2,scheme=2,cex=3)

#==model selected
term_molt_mod<-mod

new_data_1 <- expand.grid(lg_abund = seq(min(gam_dat$lg_abund), max(gam_dat$lg_abund), length.out = 30),
                          sm_abund = seq(min(gam_dat$sm_abund), max(gam_dat$sm_abund), length.out = 30),
                          Size=72.5)
pred_1 <- predict(term_molt_mod, newdata = new_data_1, type = "response", se.fit = TRUE)
df_pred_1 <- data.frame(new_data_1, fit = pred_1$fit, se = pred_1$se.fit, model = use_stocks[x])
save_dat <- rbind(df_pred_1, save_dat)
#==model consistency checks
# summary(term_molt_mod)
# concurvity(term_molt_mod)
# gam.check(term_molt_mod)
# simout<-simulateResiduals(term_molt_mod,n=200)
# png("plots/dharma_mat.png",height=6,width=8,res=350,units='in') 
# plot(simout)
# dev.off()

#==plotting
# gam_plot_term<-plot(term_molt_mod,pages=1)
# b<-getViz(term_molt_mod)
# if(x==1)
# ebs_molt<-plot(sm(b, 2),too.far=2) +
#   ylab("Small mature males")+
#   xlab("Large mature males")+
#   ggtitle("Bering Sea")+ l_fitRaster() +
#   l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") +
#   l_points()
# if(x==2)
#   nfl_molt<-plot(sm(b, 2),too.far=2) +
#   ylab("Small mature males")+
#   xlab("Large mature males")+
#   ggtitle("Newfoundland")+ l_fitRaster() +
#   l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") +
#   l_points()
# if(x==3)
#   gsl_molt<-plot(sm(b, 2),too.far=2) +
#   ylab("Small mature males")+
#   xlab("Large mature males")+
#   ggtitle("Gulf of St Lawrence")+ l_fitRaster() +
#   l_fitContour(mapping = aes(z = tz, colour = ..level..), colour = "black") +
#   l_points()
}


yoink<-ggplot(save_dat, aes(x = lg_abund, y = sm_abund, z = fit)) +
  geom_raster(aes(fill = fit), interpolate = TRUE) + # Or geom_contour()
  scale_fill_viridis_c() + # Or other color scales
  facet_wrap(~ model,scales='free',ncol=1) + # For separate panels per model
  labs(title = "p(maturing) at 75 mm",
       x = ">95mm males",
       y = "<95mm males",
       fill = "Probability") +
  theme_bw()+
  geom_contour()

png("plots/all_maturity.png",height=6,width=9,res=350,units='in') 
alldat + yoink
dev.off()
