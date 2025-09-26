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

outs<-read.csv("data/all_output.csv")
outs<-filter(outs,Year<2023)
unc_mort<-read.csv("data/uncertainty_mort.csv")
other_met<-read.csv("data/coldpool_ice.csv")
gsl_env<-read.csv("data/gsl_gam_dat.csv")
nfl_env<-read.csv("data/NFL_gam_dat.csv")
ak_env<-read.csv("data/AK_gam_dat_2.csv")

#============================
# chionoecetes species
#============================
use_stocks<-c("Bering Sea","Newfoundland","Gulf of St Lawrence")
rec_term<-NULL
mort_term<-NULL
out_plot_r<-NULL
out_plot_m<-NULL
dev_expl_m<-NULL
keep_AIC_r<-NULL
keep_AIC_m<-NULL
keep_AIC_m_big<-NULL
keep_conc<-list(list())
conc_cnt<-1

for(y in 1:length(use_stocks))
{
  #==this annoying...try to find a way to standardize these in the data munging stage
  if(y==1) # Bering Sea
  {
    set1<-filter(outs,species==use_stocks[y])[,-c(1,6)]
    tmp_unc<-unc_mort[grep(use_stocks[y],unc_mort$stock),]
    imm_sd<-tmp_unc[grep('imm',tmp_unc$stock),c(6,7)]
    colnames(imm_sd)[1]<-'imm_sd'
    mat_sd<-tmp_unc[grep('mat',tmp_unc$stock),c(6,7)]
    colnames(mat_sd)[1]<-'mat_sd'
    
    casted<-dcast(set1,Year~process,value.var="values")
    mod_dat<-merge(casted,ak_env[,c(2,8,9,10,11,12,13,19,20,21)],by="Year",all=TRUE)
    mod_dat<-merge(mod_dat,imm_sd,by="Year",all=TRUE)
    mod_dat<-merge(mod_dat,mat_sd,by="Year",all=TRUE)
    mod_dat$p_mort_mat<-1-exp(-mod_dat$M_mat)
    mod_dat$p_mort_imm<-1-exp(-mod_dat$M_imm) 
    colnames(mod_dat)[c(9,15)]<-c("temp","large_scale")
  }
  if(y==2) # Newf
  {
    set1<-filter(outs,species==use_stocks[y])[,-c(1,6)]
    tmp_unc<-unc_mort[grep(use_stocks[y],unc_mort$stock),]
    imm_sd<-tmp_unc[grep('imm',tmp_unc$stock),c(6,7)]
    colnames(imm_sd)[1]<-'imm_sd'
    mat_sd<-tmp_unc[grep('mat',tmp_unc$stock),c(6,7)]
    colnames(mat_sd)[1]<-'mat_sd'
    
    casted<-dcast(set1,Year~process,value.var="values")
    mod_dat<-merge(casted,nfl_env[,c(2,7,8)],by="Year",all=TRUE)
    mod_dat<-merge(mod_dat,imm_sd,by="Year",all=TRUE)
    mod_dat<-merge(mod_dat,mat_sd,by="Year",all=TRUE)
    mod_dat$p_mort_mat<-1-exp(-mod_dat$M_mat)
    mod_dat$p_mort_imm<-1-exp(-mod_dat$M_imm) 
    colnames(mod_dat)[c(9,10)]<-c("temp","large_scale")
  }  
  if(y==3) # gsl
  {
    set1<-filter(outs,species==use_stocks[y])[,-c(1,6)]
    tmp_unc<-unc_mort[grep(use_stocks[y],unc_mort$stock),]
    imm_sd<-tmp_unc[grep('imm',tmp_unc$stock),c(6,7)]
    colnames(imm_sd)[1]<-'imm_sd'
    mat_sd<-tmp_unc[grep('mat',tmp_unc$stock),c(6,7)]
    colnames(mat_sd)[1]<-'mat_sd'
    
    casted<-dcast(set1,Year~process,value.var="values")
    mod_dat<-merge(casted,gsl_env[,c(2,7,8,9)],by="Year",all=TRUE)
    mod_dat<-merge(mod_dat,imm_sd,by="Year",all=TRUE)
    mod_dat<-merge(mod_dat,mat_sd,by="Year",all=TRUE)
    mod_dat$p_mort_mat<-1-exp(-mod_dat$M_mat)
    mod_dat$p_mort_imm<-1-exp(-mod_dat$M_imm) 
    colnames(mod_dat)[c(11)]<-c("large_scale")
  }  
  
  
  #==p(mort) mat
  # mod_dat_base<-mod_dat[,-which(colnames(mod_dat)=="Recruitment")]
  # mod_dat_base<-mod_dat_base[complete.cases(mod_dat_base),]
  # base_mod_m<-gam(data=mod_dat_base,p_mort_mat~1,family = betar(link = "logit"),weights=1/mat_sd)
  if(y==2)
    mod_mat_m<-gam(data=mod_dat,p_mort_mat~s(N_mat,k=4)+s(temp,k=3)+s(avg_size,k=3)+s(large_scale,k=3),family = betar(link = "logit"),weights=1/mat_sd)
  if(y!=2)
    mod_mat_m<-gam(data=mod_dat,p_mort_mat~s(N_mat,k=4)+s(temp_mat,k=3)+s(avg_size,k=3)+s(large_scale,k=3),family = betar(link = "logit"),weights=1/mat_sd)
  
  # mod_mat_a<-gam(data=mod_dat_base,p_mort_mat~s(N_mat,k=4),family = betar(link = "logit"),weights=1/mat_sd)
  # mod_mat_t<-gam(data=mod_dat_base,p_mort_mat~s(temp,k=4),family = betar(link = "logit"),weights=1/mat_sd)
  # mod_mat_s<-gam(data=mod_dat_base,p_mort_mat~s(avg_size,k=3),family = betar(link = "logit"),weights=1/mat_sd)
  # mod_mat_i<-gam(data=mod_dat_base,p_mort_mat~s(large_scale,k=3),family = betar(link = "logit"),weights=1/mat_sd)
  summary(mod_mat_m)
  # simout<-simulateResiduals(mod_mat_m,n=250)
  # 
  # png(paste("plots/dharma_m",use_stocks[y],"_mat.png",sep=''),height=6,width=8,res=350,units='in') 
  # plot(simout)
  # dev.off()
  # 
  # keep_conc[[conc_cnt]]<- concurvity(mod_mat_m,full=FALSE)$observed[-1,-1]
  # conc_cnt<-conc_cnt+1
  dev_expl_m<-c(dev_expl_m,round(summary(mod_mat_m)$dev,2))
  plotted<-plot(mod_mat_m,pages=1)
  
  # keep_AIC_m<-rbind(keep_AIC_m,c(AIC(base_mod_m),AIC(mod_mat_m)))
  # keep_AIC_m_big<-rbind(keep_AIC_m_big,c(AIC(base_mod_m),AIC(mod_mat_m),AIC(mod_mat_a),AIC(mod_mat_t),AIC(mod_mat_s),AIC(mod_mat_i)))
  # 
  preds<-predict.gam(mod_mat_m,type='response',se.fit=TRUE)
  
  plo_gam_m<-data.frame(obs=mod_dat$p_mort_mat[as.numeric(names(preds$fit))],
                        preds=preds$fit,
                        y_up=preds$fit+2*preds$se,
                        y_dn=preds$fit-2*preds$se,
                        year=mod_dat$Year[as.numeric(names(preds$fit))],
                        stock=paste(use_stocks[y],"_mat",sep=''))
  
  out_plot_m<-rbind(out_plot_m,plo_gam_m)
  
  for(x in 1:(length(plotted)))
  {
    temp<-data.frame(x=plotted[[x]]$x,
                     y=plotted[[x]]$fit,
                     y_up=plotted[[x]]$fit+2*plotted[[x]]$se,
                     y_dn=plotted[[x]]$fit-2*plotted[[x]]$se,
                     covar=plotted[[x]]$xlab,
                     stock=paste(use_stocks[y],"_mat",sep=''))
    mort_term<-rbind(mort_term,temp)
  }
  
  #==p(mort) imm
  # base_mod_m<-gam(data=mod_dat_base,p_mort_imm~1,family = betar(link = "logit"),weights=1/imm_sd)
  # mod_imm_m<-gam(data=mod_dat,p_mort_imm~s(Abundance,k=4)+s(Temperature,k=4)+s(Size,k=3)+s(ice,k=3),family = betar(link = "logit"),weights=1/imm_sd)
  # mod_imm_a<-gam(data=mod_dat_base,p_mort_imm~s(Abundance,k=4),family = betar(link = "logit"),weights=1/imm_sd)
  # mod_imm_t<-gam(data=mod_dat_base,p_mort_imm~s(Temperature,k=4),family = betar(link = "logit"),weights=1/imm_sd)
  # mod_imm_s<-gam(data=mod_dat_base,p_mort_imm~s(Size,k=3),family = betar(link = "logit"),weights=1/imm_sd)
  # mod_imm_i<-gam(data=mod_dat_base,p_mort_imm~s(ice,k=3),family = betar(link = "logit"),weights=1/imm_sd)
  # 
  if(y==2)
    mod_imm_m<-gam(data=mod_dat,p_mort_imm~s(N_imm,k=4)+s(temp,k=3)+s(avg_size,k=3)+s(large_scale,k=3),family = betar(link = "logit"),weights=1/mat_sd)
  if(y!=2)
    mod_imm_m<-gam(data=mod_dat,p_mort_imm~s(N_imm,k=4)+s(temp_imm,k=3)+s(avg_size,k=3)+s(large_scale,k=3),family = betar(link = "logit"),weights=1/mat_sd)
  
  # keep_conc[[conc_cnt]]<- concurvity(mod_imm_m,full=FALSE)$observed[-1,-1]
  # conc_cnt<-conc_cnt+1
  # simout<-simulateResiduals(mod_imm_m,n=250)
  # 
  # png(paste("plots/dharma_m",use_stocks[y],"_imm.png",sep=''),height=6,width=8,res=350,units='in') 
  # plot(simout)
  # dev.off()
  
  dev_expl_m<-c(dev_expl_m,round(summary(mod_imm_m)$dev,2))
  plotted<-plot(mod_imm_m,pages=1)
  
  # keep_AIC_m<-rbind(keep_AIC_m,c(AIC(base_mod_m),AIC(mod_imm_m)))
  # keep_AIC_m_big<-rbind(keep_AIC_m_big,c(AIC(base_mod_m),AIC(mod_imm_m),AIC(mod_imm_a),AIC(mod_imm_t),AIC(mod_imm_s),AIC(mod_imm_i)))
  # 
  preds<-predict.gam(mod_imm_m,type='response',se.fit=TRUE)
  
  plo_gam_m<-data.frame(obs=mod_dat$M_imm[as.numeric(names(preds$fit))],
                        preds=preds$fit,
                        y_up=preds$fit+2*preds$se,
                        y_dn=preds$fit-2*preds$se,
                        year=mod_dat$Year[as.numeric(names(preds$fit))],
                        stock=paste(use_stocks[y],"_imm",sep=''))
  
  out_plot_m<-rbind(out_plot_m,plo_gam_m)
  
  for(x in 1:(length(plotted)))
  {
    temp<-data.frame(x=plotted[[x]]$x,
                     y=plotted[[x]]$fit,
                     y_up=plotted[[x]]$fit+2*plotted[[x]]$se,
                     y_dn=plotted[[x]]$fit-2*plotted[[x]]$se,
                     covar=plotted[[x]]$xlab,
                     stock=paste(use_stocks[y],"_imm",sep=''))
    mort_term<-rbind(mort_term,temp)
  }
  
 
  
}

colnames(keep_AIC_m_big)<-c("No covars","All covars","Density","Temp","Size","Ice")
rownames(keep_AIC_m_big)<-c("BBRKC","PIRKC","SMBKC","PIBKC","Snow (mat)","Snow (imm)","Tanner (mat)","Tanner (imm)")

names(dev_expl_m)<-c("ESB (mat)","EBS (imm)","NFL (mat)","NFL (imm)","GSL (mat)","GSL (imm)")

# plot_conc<-NULL
# for(x in 1:length(keep_conc))
# {
#  tmp<-melt(keep_conc[[x]]) 
#  #levels(tmp$Var1)<-c("s(Abundance)","s(Size)","s(Temperature)","s(Ice)")
#  #levels(tmp$Var2)<-c("s(Abundance)","s(Size)","s(Temperature)","s(Ice)")
#  tmp$stock<-names(dev_expl_m)[x]
#  plot_conc<-rbind(plot_conc,tmp)
# }
# 
# plot_conc[plot_conc==1]<-NA
# png("plots/model_concurvity.png",height=8,width=8,res=350,units='in') 
# ggplot(plot_conc)+
#   geom_tile(aes(x=Var1,y=Var2,fill=value))+
#   geom_text(aes(x=Var1,y=Var2,label=round(value,2)))+
#   scale_fill_gradient(low='white',high='red',na.value='whitesmoke')+
#   facet_wrap(~stock)+theme_bw()+
#   theme(legend.position=c(.85,.15))+
#   ggtitle("Estimated concurvity among variables by model")+
#   ylab("")+xlab("")+guides(fill=guide_legend(title="Concurvity"))
# dev.off()
# 
# 
# alt_AIC<-keep_AIC_m_big
# for(x in 1:nrow(alt_AIC))
#   for(y in 2:ncol(alt_AIC))
#     alt_AIC[x,y]<--1*(alt_AIC[x,1]-alt_AIC[x,y])
# 
# in_dat_t1<-as.data.frame(round(alt_AIC,2))
# 
# library(tidyr)
# library(scales)
# # Specify the columns to be colored
# columns_to_color <- colnames(in_dat_t1[2:6])
# 
# # Reshape the data for ggplot2
# data_long <- in_dat_t1 %>%
#   mutate(RowID = rownames(in_dat_t1)) %>%
#   pivot_longer(cols = all_of(columns_to_color), names_to = "Variable", values_to = "Value")
# 
# # Create the heatmap-like graphic
# aic_table<-ggplot(data_long, aes(x = Variable, y = as.factor(RowID), fill = Value)) +
#   geom_tile(color = "black") + # Add borders to cells
#   scale_fill_gradient2(
#     low = "forestgreen",
#     mid = "white",
#     high = "tomato3",
#     midpoint = 0,
#     space = "Lab",
#     limits=c(-10,10),
#     oob=squish
#   ) +
#   geom_text(aes(label = Value), color = "black") + # Add cell values as text
#   labs(title = "Change in AIC from null model", x = "",y='') +
#   theme_minimal()+
#   theme(legend.position="none")
# 
# use_dat<-melt(dev_expl_m)
# colnames(use_dat)[1]<-"Deviance_explained"
# use_dat$stock<-rownames(use_dat)
# 
# my_values <- dev_expl_m
# 
# # Convert the vector to a data frame, necessary for ggplot2
# df <- data.frame(
#   Value = my_values,
#   Row = seq_along(my_values),  # Assign a row number for each value
#   Row_Name = names(my_values)         # Add the row names to the data frame
# )
# 
# # Create the ggplot
# dev_ex_plot<-ggplot(df, aes(x = 1, y = as.factor(Row_Name), fill = Value)) + # Map x to a single column (1), y to Row, and fill to Value
#   geom_tile(color = "black") +  # Create tiles, add a black border for clarity
#   scale_fill_gradient(low = "white", high = "cornflowerblue") +  # Apply the gradient
#   geom_text(aes(label = Value), color = "black") + # Add the value as text inside each tile
#  # geom_text(aes(x = 0.5, label = Row_Name), hjust = 1, color = "black") + # Add row names to the left of the tiles
#   labs(title = "Deviance explained", x = NULL, y = NULL) +  # Add title and remove axis labels
#   theme_minimal() + 
#   xlim(0,2)+
#   theme(axis.text.x = element_blank(), # Remove x-axis text
#         axis.ticks.x = element_blank(), # Remove x-axis ticks
#         axis.text.y = element_blank(), # Remove x-axis text
#         axis.ticks.y = element_blank(), # Remove x-axis ticks
#         panel.grid = element_blank(),  # Remove grid lines
#         legend.position = "none")  # Place the legend on the right
# 
# png("plots/model_select.png",height=6,width=7,res=350,units='in') 
# aic_table + dev_ex_plot +plot_layout(width=c(2,1))
# dev.off()

in_col2<-c("#0034c3","#0034c388","#3da550","#3da55088","#ff8f38","#ff8f3888")
mort_term$covar[mort_term$covar=="N_imm"]<-"Abundance"
mort_term$covar[mort_term$covar=="N_mat"]<-"Abundance"
mort_term$covar[mort_term$covar=="temp_imm"]<-"temp"
mort_term$covar[mort_term$covar=="temp_mat"]<-"temp"

ughh<-out_plot_m[,c(1,5,6)]
colnames(ughh)<-c("est_m","year","stock'")

mod_mort<-unc_mort
mod_mort$up_m[mod_mort$up_m>max(mod_mort$est_m)]<-1.5*max(mod_mort$est_m)

mod_mort_cap <- mod_mort %>%
  group_by(stock) %>%
  mutate(up_m = pmin(up_m, 1.5 * max(est_m, na.rm = TRUE))) %>% # The `pmin` function is used here to take the minimum of the two values element-wise, effectively capping the up_m values.
  ungroup()

mort_plot_alt_trans<-ggplot()+
  geom_point(data=filter(unc_mort,Year<2023&Year!=2020),aes(x=Year,y=1-exp(-est_m)),col='darkgrey')+
  geom_errorbar(data=filter(unc_mort,Year<2023&Year!=2020),aes(x=Year,ymin=1-exp(-dn_m),ymax=1-exp(-up_m)),col='grey')+
  geom_ribbon(data=out_plot_m,aes(x=year,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.3,lwd=2)+
  geom_line(data=out_plot_m,aes(x=year,y=preds,col=stock),lwd=1.15)+  theme_bw()+ylab("p(mortality)")+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)+
  facet_wrap(~stock,ncol=1,scales='free_y')+
  theme(legend.position='none')+ylim(0,1)


mort_plot_trm<-ggplot(mort_term)+
  geom_line(aes(x=x,y=y,col=stock),lwd=2)+
  geom_ribbon(aes(x=x,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.3,lwd=2)+
  theme_bw()+
  theme(legend.position='none',
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.key.size=unit(.7, 'lines'))+
  geom_hline(yintercept=0,lty=2)+xlab("Observed value")+ylab("Smooth")+
  facet_grid(rows=vars(stock),cols=vars(covar),scales='free_x')+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)

mort_plot_1<-ggplot(mort_term)+
  geom_line(aes(x=x,y=y,col=stock),lwd=1.25)+
  geom_ribbon(aes(x=x,ymin=y_dn,ymax=y_up,fill=stock),alpha=0.2,lwd=2)+
  theme_bw()+
  theme(legend.position='none',
        legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank(),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.key.size=unit(.7, 'lines'))+
  geom_hline(yintercept=0,lty=2)+xlab("Observed value")+ylab("Smooth")+
  facet_wrap(~covar,scales='free',ncol=1)+
  scale_color_manual(values=in_col2)+
  scale_fill_manual(values=in_col2)

# mortality <- density + temperature + competitor + size
library(patchwork)
design<-"112
         112
         112"
png("plots/all_mort_all3.png",height=9,width=6,res=350,units='in') 
mort_plot_alt_trans +mort_plot_1 + plot_layout(nrow=2, design=design)
dev.off()



#==================================
# among pop cors
#========================
unique(outs$process)
library(GGally)
library(forecast)
#==by prGGally#==by process
#==put ccf in upper triangle
use_gg<-filter(outs,process=='Recruitment')
casted<-dcast(use_gg,Year~species,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

# Correlation matrix plot
p2 <- ggcorr(casted, label = TRUE, label_round = 2)
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p<-ncol(casted)
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}
p1 = ggpairs(casted, lower = list(continuous = my_fn))
png("plots/all_rec_corr.png",height=10,width=10,res=350,units='in') 
print(p1)
dev.off()

#=========================
# other mortality

use_gg<-filter(outs,process%in%c('Other mortality',"M_mat"))
casted<-dcast(use_gg,Year~species,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

# Correlation matrix plot
p2 <- ggcorr(casted, label = TRUE, label_round = 2)
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p<-ncol(casted)
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))
png("plots/all_mort_corr.png",height=10,width=10,res=350,units='in') 
print(p1)
dev.off()

# other mortality

use_gg<-filter(outs,process%in%c('Other mortality',"M_imm"))
casted<-dcast(use_gg,Year~species,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

# Correlation matrix plot
p2 <- ggcorr(casted, label = TRUE, label_round = 2)
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p<-ncol(casted)
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))
png("plots/all_mort_imm_corr.png",height=10,width=10,res=350,units='in') 
print(p1)
dev.off()

#=========================
# abundance

use_gg<-filter(outs,process%in%c('Abundance'))
casted<-dcast(use_gg,Year~species,value.var="values")[,-1]
my_fn <- function(data, mapping, ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=loess, fill="red", color="red", ...) +
    geom_smooth(method=lm, fill="blue", color="blue", ...)
  p+theme_bw()
}

# Correlation matrix plot
p2 <- ggcorr(casted, label = TRUE, label_round = 2)
g2 <- ggplotGrob(p2)
colors <- g2$grobs[[6]]$children[[3]]$gp$fill

# Change background color to tiles in the upper triangular matrix of plots 
idx <- 1
p<-ncol(casted)
for (k1 in 1:(p-1)) {
  for (k2 in (k1+1):p) {
    plt <- getPlot(p1,k1,k2) +
      theme(panel.background = element_rect(fill = colors[idx], color="white"),
            panel.grid.major = element_line(color=colors[idx]))
    p1 <- putPlot(p1,plt,k1,k2)
    idx <- idx+1
  }
}

p1 = ggpairs(casted, lower = list(continuous = my_fn))
png("plots/all_abund_corr.png",height=10,width=10,res=350,units='in') 
print(p1)
dev.off()





