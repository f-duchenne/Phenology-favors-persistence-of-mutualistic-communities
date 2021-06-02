library(MASS)
require(sp)
require(rgdal)
require(maps)
require(doBy)
require(ggplot2)
require(gridExtra)
require(lubridate)
library(chron)
library(dplyr)
library(chron)
library(geosphere)
library(reshape2)
library(car)
library(data.table)
library(jtools)
library(plot3D)
library(cowplot)
library(viridis)
library(RColorBrewer)
library(igraph)
library(semPlot)
library(lavaan)
require(scales)
library(expm)
library(cowplot)
library(lme4)
library(nlme)

memory.limit(size=250000)
setwd(dir="C:/Users/duche/Documents/5eme - papier modelo/75x75")
lili=list.files()
lili=lili[grep("75x75",lili,fixed=T)]
lili=lili[grep("samedistrib",lili,fixed=T,invert=T)]
tabnet=NULL
tabsp=NULL
for(i in 1:length(lili)){
bidon=fread(lili[i],sep="\t",header=T)
if(i %in% grep("netcar_",lili)){tabnet=rbind(tabnet,bidon)}else{tabsp=rbind(tabsp,bidon)}
}

tabsp$surv=0
tabsp$surv[tabsp$N_eq1>0]=1
tabsp$mfd[tabsp$mfd>365]=tabsp$mfd[tabsp$mfd>365]-365
tabsp$mfd[tabsp$mfd<1]=365+tabsp$mfd[tabsp$mfd<1]
#model=glm(surv~sd*beta+trait_mu*alpha+alpha*beta+interf+type,family=binomial,data=tabsp)
tabsp$mfd2=abs(tabsp$mfd-190)
tabsp$emi_indfp=tabsp$emi_totfp-tabsp$emi_directfp
tabsp$emi_indpp=tabsp$emi_totpp-tabsp$emi_directpp
tabsp$recu_indfp=tabsp$recu_totfp-tabsp$recu_directfp
tabsp$recu_indpp=tabsp$recu_totpp-tabsp$recu_directpp
tabsp=tabsp %>% dplyr::group_by(alpha,beta,interf,random) %>% dplyr::mutate(
I_mean3=I_mean/max(I_mean))
tabsp$I_mean2=plyr::round_any(tabsp$I_meaneq1,0.1,f=ceiling)
tabsp$I_mean3=plyr::round_any(tabsp$I_mean3,0.1)
tabsp$P_mean2=plyr::round_any(tabsp$P_mean,0.1)
tabsp$A_mean2=plyr::round_any(tabsp$A_mean,0.1)
tabsp$sd2=plyr::round_any(tabsp$sd,2)
tabsp$trait_sd2=plyr::round_any(tabsp$trait_sd,0.1)
tabsp$type[tabsp$type=="flow"]="plants"

datt=subset(tabsp,beta==0 & alpha>0) %>% dplyr::group_by(interf,trait_sd2,sd2,type) %>%
dplyr::summarise(N_eq1_sd=sd(surv,na.rm=T)/sqrt(length(surv)),N_eq1=mean(surv,na.rm=T))
datp=subset(tabsp,alpha==0 & beta>0) %>% dplyr::group_by(interf,trait_sd2,sd2,type) %>%
dplyr::summarise(N_eq1_sd=sd(surv,na.rm=T)/sqrt(length(surv)),N_eq1=mean(surv,na.rm=T))
dattp=subset(tabsp,alpha>0 & beta>0) %>% dplyr::group_by(interf,trait_sd2,sd2,type) %>%
dplyr::summarise(N_eq1_sd=sd(surv,na.rm=T)/sqrt(length(surv)),N_eq1=mean(surv,na.rm=T))

pl1=ggplot(data=datt,aes(x=trait_sd2,y=sd2,fill=N_eq1))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
strip.text=element_text(size=12),plot.subtitle=element_text(size=14))+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),na.value="white")+
labs(fill='')+
coord_cartesian(expand=F)+xlab("Morphological niche width")+ylab("Flowering/flight period duration")+
facet_wrap(~paste0("c = ",interf)+type,ncol=2)+ggtitle("a")

pl2=ggplot(data=datp,aes(x=trait_sd2,y=sd2,fill=N_eq1))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
strip.text=element_text(size=12),plot.subtitle=element_text(size=14))+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),na.value="white")+
labs(fill='')+
coord_cartesian(expand=F)+xlab("Morphological niche width")+ylab("Flowering/flight period duration")+
facet_wrap(~paste0("c = ",interf)+type,ncol=2)+ggtitle("b")


pl3=ggplot(data=dattp,aes(x=trait_sd2,y=sd2,fill=N_eq1))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
strip.text=element_text(size=12),plot.subtitle=element_text(size=14))+
scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),na.value="white")+
labs(fill='')+
coord_cartesian(expand=F)+xlab("Morphological niche width")+ylab("Flowering/flight period duration")+
facet_wrap(~paste0("c = ",interf)+type,ncol=2)+ggtitle("c")

png("figS1.png",width=1400,height=1000,res=140)
grid.arrange(pl1,pl2,pl3,ncol=3)
dev.off();


#EFFETS indirect propagés par généralistes
datt=subset(tabsp,alpha!=0 & beta!=0) %>% group_by(interf,I_mean2,type) %>%
summarise(emi_indfp_sd=sd(emi_indfp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),emi_indpp_sd=sd(emi_indpp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),
emi_directfp_sd=sd(emi_directfp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),emi_directpp_sd=sd(emi_directpp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),
emi_indfp=mean(emi_indfp,na.rm=T),emi_indpp=mean(emi_indpp,na.rm=T),surv=mean(surv),
emi_directfp=mean(emi_directfp,na.rm=T),emi_directpp=mean(emi_directpp,na.rm=T),surv=mean(surv))
datt$emi_directpp[is.na(datt$emi_directpp)]=0
pl1=ggplot(datt,aes(shape=type))+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point(data=datt,aes(x=I_mean2+0.02,y=emi_indpp),fill="red")+
geom_errorbar(data=datt,aes(x=I_mean2+0.02,ymin=emi_indpp-1.96*emi_indpp_sd,
ymax=emi_indpp+1.96*emi_indpp_sd),colour="red",width=0)+
geom_point(data=datt,aes(x=I_mean2+0.02,y=emi_directpp),fill="black")+
geom_errorbar(data=datt,aes(x=I_mean2+0.02,ymin=emi_directpp-1.96*emi_directpp_sd,
ymax=emi_directpp+1.96*emi_directpp_sd),colour="black",width=0)+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Generated effects within guilds")+
xlab("Generalism level at the equilibrium")+
scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+
theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.text=element_text(size=12),axis.title=element_text(size=12),legend.position="none")+ggtitle("")+
scale_shape_manual(values=c(21,22))
# pl2=ggplot(datt,aes(shape=type))+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
# geom_point(data=datt,aes(x=I_mean2-0.02,y=emi_indfp),fill="red")+
# geom_errorbar(data=datt,aes(x=I_mean2-0.02,ymin=emi_indfp-1.96*emi_indfp_sd,
# ymax=emi_indfp+1.96*emi_indfp_sd),colour="red",width=0)+
# geom_point(data=datt,aes(x=I_mean2-0.02,y=emi_directfp),fill="black")+
# geom_errorbar(data=datt,aes(x=I_mean2-0.02,ymin=emi_directfp-1.96*emi_directfp_sd,
# ymax=emi_directfp+1.96*emi_directfp_sd),colour="black",width=0)+
# facet_wrap(~paste0("c = ",interf),ncol=4)+
# ylab("Propagated effects among guilds")+
# xlab("Degree at the equilibrium")+
# scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+
# theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
# strip.text=element_text(size=12),axis.title=element_text(size=12),legend.position="none")+ggtitle("b")+
# scale_shape_manual(values=c(21,22))
# grid.arrange(pl1,pl2)
png("figS6.png",width=1000,height=700,res=140)
pl1
dev.off();


####PROPAGATIOn tot
datt=subset(tabsp,beta>0 & alpha>0) %>% group_by(interf,I_mean2,type) %>%
summarise(emi_totfp_sd=sd(emi_totfp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),emi_totpp_sd=sd(emi_totpp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),
emi_totfp=mean(emi_totfp,na.rm=T),emi_totpp=mean(emi_totpp,na.rm=T),surv=mean(surv))
datt$emi_directpp[is.na(datt$emi_directpp)]=0
datt$type[datt$type=="flow"]="Plants"
datt$type[datt$type=="poll"]="Pollinators"

pl1=ggplot(datt,aes(shape=type))+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_totpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point(data=datt,aes(x=I_mean2,y=emi_totpp),fill="black",alpha=0.7)+
geom_errorbar(data=datt,aes(x=I_mean2,ymin=emi_totpp-1.96*emi_totpp_sd,
ymax=emi_totpp+1.96*emi_totpp_sd),colour="black",width=0,)+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Propagated effects within guilds")+
xlab("Degree at equilibrium")+
scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+
theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.text=element_text(size=12),axis.title=element_text(size=12),legend.position="none")+ggtitle("b")+
scale_shape_manual(values=c(21,22))

# #MAIS EFFETS intotal émis plus positifs pour spécialistes
# datt=subset(tabsp,beta==0 & alpha>0) %>% group_by(interf,I_mean2,type) %>%
# summarise(prop_pp=mean(emi_totpp-recu_totpp,na.rm=T),prop_pp_sd=sd(emi_totpp-recu_totpp,na.rm=T)/sqrt(length(na.omit(recu_totpp))),
# prop_fp=mean(emi_totfp-recu_totfp,na.rm=T),prop_fp_sd=sd(emi_totfp-recu_totfp,na.rm=T)/sqrt(length(na.omit(recu_totpp))))
# datp=subset(tabsp,alpha==0 & beta>0) %>% group_by(interf,I_mean2,type) %>%
# summarise(prop_pp=mean(emi_totpp-recu_totpp,na.rm=T),prop_pp_sd=sd(emi_totpp-recu_totpp,na.rm=T)/sqrt(length(na.omit(recu_totpp))),
# prop_fp=mean(emi_totfp-recu_totfp,na.rm=T),prop_fp_sd=sd(emi_totfp-recu_totfp,na.rm=T)/sqrt(length(na.omit(recu_totpp))),
# recu_totpp_moy=mean(recu_totpp,na.rm=T),emi_totpp_moy=mean(emi_totpp,na.rm=T))
# # datt=subset(datt,I_mean3!=0 & I_mean3!=1)
# ggplot(data=datp)+
# geom_point(data=datp,aes(x=I_mean2,y=recu_totpp_moy),colour="gold3")+
# geom_point(data=datp,aes(x=I_mean2,y=emi_totpp_moy),colour="dodgerblue3")+
# facet_wrap(~paste0("c = ",interf),ncol=4)

# pl2=ggplot(datt,aes(shape=type))+
# geom_point(data=datt,aes(x=I_mean2+0.02,y=prop_pp),position = position_dodge(width = 0.02),fill="dodgerblue3",alpha=0.7)+
# geom_errorbar(data=datt,aes(x=I_mean2+0.02,ymin=prop_pp-1.96*prop_pp_sd,
# ymax=prop_pp+1.96*prop_pp_sd),position = position_dodge(width = 0.02),colour="dodgerblue3",width=0)+
# geom_point(data=datp,aes(x=I_mean2-0.02,y=prop_pp),position = position_dodge(width = 0.02),fill="gold3",alpha=0.7)+
# geom_errorbar(data=datp,aes(x=I_mean2-0.02,ymin=prop_pp-1.96*prop_pp_sd,
# ymax=prop_pp+1.96*prop_pp_sd),position = position_dodge(width = 0.02),colour="gold3",width=0)+
# facet_wrap(~paste0("c = ",interf),ncol=4)+
# ylab("Propagated - Received effects")+
# xlab("Degree at the equilibrium")+
# scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+
# theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
# strip.text=element_text(size=12),axis.title=element_text(size=12),legend.position="none")+ggtitle("a")+
# scale_shape_manual(values=c(21,22))
# #Among guilds
# pl2b=ggplot(datt,aes(shape=type))+
# geom_point(data=datt,aes(x=I_mean2+0.02,y=prop_fp),colour="dodgerblue3")+
# geom_errorbar(data=datt,aes(x=I_mean2+0.02,ymin=prop_fp-prop_fp_sd,
# ymax=prop_fp+prop_fp_sd),colour="dodgerblue3",width=0)+
# geom_point(data=datp,aes(x=I_mean2-0.02,y=prop_fp),colour="gold3")+
# geom_errorbar(data=datp,aes(x=I_mean2-0.02,ymin=prop_fp-prop_fp_sd,
# ymax=prop_fp+prop_fp_sd),colour="gold3",width=0)+
# facet_wrap(~paste0("c = ",interf),ncol=4)+
# ylab("Propagated - Received effects")+
# xlab("Degree at the equilibrium")+
# scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+
# theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
# strip.text=element_text(size=12),axis.title=element_text(size=12))+ggtitle("a")


#SPECIALISTES TRAITS MEURS PLUS
datt=subset(tabsp,beta==0 & alpha>0) %>% dplyr::group_by(interf,I_mean3,type) %>%
dplyr::summarise(N_eq1_sd=sd(surv,na.rm=T)/sqrt(length(surv)),N_eq1=mean(surv,na.rm=T))
datp=subset(tabsp,alpha==0 & beta>0) %>% dplyr::group_by(interf,I_mean3,type) %>%
dplyr::summarise(N_eq1_sd=sd(surv,na.rm=T)/sqrt(length(surv)),N_eq1=mean(surv,na.rm=T))
pl1b=ggplot(datt,aes(shape=type))+
geom_point(data=datt,aes(x=I_mean3+0.02,y=N_eq1),position = position_dodge(width = 0.02),fill="dodgerblue3",alpha=0.7)+
geom_errorbar(data=datt,aes(x=I_mean3+0.02,ymin=N_eq1-1.96*N_eq1_sd,
ymax=N_eq1+1.96*N_eq1_sd),position = position_dodge(width = 0.02),colour="dodgerblue3",width=0)+
geom_point(data=datp,aes(x=I_mean3-0.02,y=N_eq1),position = position_dodge(width = 0.02),fill="gold3",alpha=0.7)+
geom_errorbar(data=datp,aes(x=I_mean3-0.02,ymin=N_eq1-1.96*N_eq1_sd,
ymax=N_eq1+1.96*N_eq1_sd),position = position_dodge(width = 0.02),colour="gold3",width=0)+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Average persistence probability")+
xlab("Initial Degree")+
scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+scale_y_continuous(breaks=c(0,0.5,1))+
theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.text=element_text(size=12),axis.title=element_text(size=12),legend.position="none")+ggtitle("a")+
scale_shape_manual(values=c(21,22))

grid.arrange(pl1,pl1b)

png("fig5.png",width=1100,height=1000,res=140)
grid.arrange(pl1b,pl1)
dev.off();


#DECORTIQUER TOTAL DES SPECIALISTES 
obj=tabsp %>% dplyr::group_by(alpha,beta,interf,random,type) %>% dplyr::mutate(
contr_indpp=emi_indpp/(abs(emi_indpp)+abs(emi_directpp)),
contr_indfp=emi_indfp/(abs(emi_indfp)+abs(emi_directfp)))

datt=subset(obj,alpha!=0 | beta!=0) %>% group_by(interf,I_mean2) %>%
summarise(contr_indpp=mean(contr_indpp,na.rm=T),contr_indfp=mean(contr_indfp,na.rm=T))
# datt=subset(datt,I_mean2!=0 & I_mean2!=1)
pl1=ggplot(datt)+
# geom_point(aes(x=sd,y=recu_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=recu_totpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point(data=datt,aes(x=I_mean2,y=contr_indfp),colour="black")+
geom_point(data=datt,aes(x=I_mean2,y=contr_indpp),colour="red")+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Standardized strenght of totirect effects \n on other species (inter-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))

#Si on prend que les effets inds pour spécialistes
datt=subset(tabsp,beta==0 & alpha>0) %>% group_by(interf,I_mean2,type) %>%
summarise(prop_pp=mean(emi_directpp-recu_directpp,na.rm=T),prop_pp_sd=sd(emi_directpp-recu_directpp,na.rm=T),
prop_pp_ind=mean(emi_indpp-recu_indpp,na.rm=T),prop_pp_sd_ind=sd(emi_indpp-recu_indpp,na.rm=T))
datp=subset(tabsp,alpha==0 & beta>0) %>% group_by(interf,I_mean2,type) %>%
summarise(prop_pp=mean(emi_directpp-recu_directpp,na.rm=T),prop_pp_sd=sd(emi_directpp-recu_directpp,na.rm=T),
prop_pp_ind=mean(emi_indpp-recu_indpp,na.rm=T),prop_pp_sd_ind=sd(emi_indpp-recu_indpp,na.rm=T))
pl2=ggplot(datt,aes(shape=type))+
geom_point(data=datt,aes(x=I_mean2+0.02,y=prop_pp),fill="dodgerblue3",alpha=0.7)+
geom_errorbar(data=datt,aes(x=I_mean2+0.02,ymin=prop_pp-prop_pp_sd,
ymax=prop_pp+prop_pp_sd),colour="dodgerblue3",width=0)+
geom_point(data=datp,aes(x=I_mean2-0.02,y=prop_pp),fill="gold3",alpha=0.7)+
geom_errorbar(data=datp,aes(x=I_mean2-0.02,ymin=prop_pp-prop_pp_sd,
ymax=prop_pp+prop_pp_sd),colour="gold3",width=0)+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Propagated - Received\ndirect effects")+
xlab("Degree at the equilibrium")+
scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+
theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.text=element_text(size=12),axis.title=element_text(size=12),legend.position="none")+ggtitle("a")+
scale_shape_manual(values=c(21,22))

pl1=ggplot(datt,aes(shape=type))+
geom_point(data=datt,aes(x=I_mean2+0.02,y=prop_pp_ind),fill="dodgerblue3",alpha=0.7)+
geom_errorbar(data=datt,aes(x=I_mean2+0.02,ymin=prop_pp_ind-prop_pp_sd_ind,
ymax=prop_pp_ind+prop_pp_sd_ind),colour="dodgerblue3",width=0)+
geom_point(data=datp,aes(x=I_mean2-0.02,y=prop_pp_ind),fill="gold3",alpha=0.7)+
geom_errorbar(data=datp,aes(x=I_mean2-0.02,ymin=prop_pp_ind-prop_pp_sd_ind,
ymax=prop_pp_ind+prop_pp_sd_ind),colour="gold3",width=0)+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Propagated - Received\nindirect effects")+
xlab("Degree at the equilibrium")+
scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+
theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.text=element_text(size=12),axis.title=element_text(size=12),legend.position="none")+ggtitle("b")+
scale_shape_manual(values=c(21,22))


grid.arrange(pl2,pl1)


png("figS4.png",width=1100,height=700,res=140)
grid.arrange(pl2,pl1)
dev.off();
















pl1=ggplot(data=subset(tabsp,alpha>0 & beta>0))+
# geom_point(aes(x=sd,y=emi_totpp-recu_totpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
facet_wrap(~interf,scale="free_y",ncol=4)+
stat_smooth(aes(x=sd,y=emi_totpp-recu_totpp,col="type"))
pl2=ggplot(data=subset(tabsp,alpha==0.5 & beta==0.5 & type=="poll"))+
# geom_point(aes(x=sd,y=emi_totfp-recu_totfp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
stat_smooth(aes(x=sd,y=emi_totfp-recu_totfp))+
facet_wrap(~interf,scale="free_y",ncol=4)
grid.arrange(pl1,pl2)


#EFFETS indirect propagés par généralistes
obj=tabsp %>% dplyr::group_by(alpha,beta,interf,random,type) %>% dplyr::mutate(
recu_indpp=scale(abs(recu_indpp)),
recu_indfp=scale(abs(recu_indfp)),emi_indpp=scale(abs(emi_indpp)),emi_indfp=scale(abs(emi_indfp)),
recu_totpp=scale(abs(recu_totpp)),recu_totfp=scale(abs(recu_totfp)),emi_totpp=scale(abs(emi_totpp)),
emi_totfp=scale(abs(emi_totfp)))

obj$I_mean=plyr::round_any(obj$I_mean,0.1)
obj$P_mean=plyr::round_any(obj$P_mean,0.1)
obj$A_mean=plyr::round_any(obj$A_mean,0.1)

datt=subset(obj,alpha>0 & beta>0) %>% group_by(interf,alpha,beta,A_mean) %>%
summarise(emi_indfp_sd=sd(emi_indfp,na.rm=T),emi_indpp_sd=sd(emi_indpp,na.rm=T),
emi_indfp=mean(emi_indfp,na.rm=T),emi_indpp=mean(emi_indpp,na.rm=T),surv=mean(surv))
datp=subset(obj,alpha>0 & beta>0) %>% group_by(interf,alpha,beta,P_mean) %>%
summarise(emi_indfp_sd=sd(emi_indfp,na.rm=T),emi_indpp_sd=sd(emi_indpp,na.rm=T),
emi_indfp=mean(emi_indfp,na.rm=T),emi_indpp=mean(emi_indpp,na.rm=T),surv=mean(surv))

pl1=ggplot(datt)+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
stat_smooth(data=datt,aes(x=A_mean,y=surv,col=as.factor(alpha),fill=as.factor(alpha)),linetype="dashed")+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Standardized strenght of indirect effects \n on other species (inter-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))
pl2=ggplot(datp)+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
stat_smooth(data=datp,aes(x=P_mean,y=surv,col=as.factor(beta),fill=as.factor(beta)))+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Standardized strenght of indirect effects \n on other species (inter-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))
grid.arrange(pl1,pl2)

pl1=ggplot(datt,aes(x=I_mean,y=surv,col=as.factor(beta),fill=as.factor(beta)))+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point()+stat_smooth()+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Standardized strenght of indirect effects \n on other species (inter-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))





pl2=ggplot(data=datt)+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point(aes(x=I_mean,y=emi_indpp),colour="black")+
geom_errorbar(aes(x=I_mean,ymin=emi_indpp-emi_indpp_sd,
ymax=emi_indpp+emi_indpp_sd),colour="black",width=0)+
facet_wrap(~paste0("c = ",interf),scale="free_y",ncol=4)+
ylab("Standardized strenght of indirect effects \n on other species (intra-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))
grid.arrange(pl1,pl2)
###########################

#EFFETS indirect propagés par généralistes
obj=tabsp %>% dplyr::group_by(alpha,beta,interf,random,type) %>% dplyr::mutate(
recu_indpp=scale(abs(recu_indpp)),
recu_indfp=scale(abs(recu_indfp)),emi_indpp=scale(abs(emi_indpp)),emi_indfp=scale(abs(emi_indfp)),
recu_totpp=scale(abs(recu_totpp)),recu_totfp=scale(abs(recu_totfp)),emi_totpp=scale(abs(emi_totpp)),
emi_totfp=scale(abs(emi_totfp)))

obj$sd=plyr::round_any(obj$sd,2)
obj$trait_sd=plyr::round_any(obj$trait_sd,0.05)


datt=subset(obj,alpha>0 & beta>0) %>% group_by(interf,trait_sd) %>%
summarise(emi_indfp_sd=sd(emi_indfp,na.rm=T),emi_indpp_sd=sd(emi_indpp,na.rm=T),
emi_indfp=mean(emi_indfp,na.rm=T),emi_indpp=mean(emi_indpp,na.rm=T))
datp=subset(obj,alpha>0 & beta>0) %>% group_by(interf,sd) %>%
summarise(emi_indfp_sd=sd(emi_indfp,na.rm=T),emi_indpp_sd=sd(emi_indpp,na.rm=T),
emi_indfp=mean(emi_indfp,na.rm=T),emi_indpp=mean(emi_indpp,na.rm=T))
pl1=ggplot(datt)+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point(data=datt,aes(x=(trait_sd-0.1)/0.8,y=emi_indfp),colour="black")+
geom_errorbar(data=datt,aes(x=(trait_sd-0.1)/0.8,ymin=emi_indfp-emi_indfp_sd,
ymax=emi_indfp+emi_indfp_sd),colour="black",width=0)+
geom_point(data=datp,aes(x=(sd-5)/35,y=emi_indfp),colour="red")+
geom_errorbar(data=datp,aes(x=(sd-5)/35,ymin=emi_indfp-emi_indfp_sd,
ymax=emi_indfp+emi_indfp_sd),colour="red",width=0)+
facet_wrap(~paste0("c = ",interf),scale="free_y",ncol=4)+
ylab("Standardized strenght of indirect effects \n on other species (inter-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))
pl2=ggplot(data=obj)+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_indpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point(data=datt,aes(x=(trait_sd-0.1)/0.8,y=emi_indpp),colour="black")+
geom_errorbar(data=datt,aes(x=(trait_sd-0.1)/0.8,ymin=emi_indpp-emi_indpp_sd,
ymax=emi_indpp+emi_indpp_sd),colour="black",width=0)+
geom_point(data=datp,aes(x=(sd-5)/35,y=emi_indpp),colour="red")+
geom_errorbar(data=datp,aes(x=(sd-5)/35,ymin=emi_indpp-emi_indpp_sd,
ymax=emi_indpp+emi_indpp_sd),colour="red",width=0)+
facet_wrap(~paste0("c = ",interf),scale="free_y",ncol=4)+
ylab("Standardized strenght of indirect effects \n on other species (intra-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))
grid.arrange(pl1,pl2)
###########################


obj=tabsp %>% dplyr::group_by(alpha,beta,interf,random,type) %>% dplyr::mutate(
recu_indpp=scale(recu_indpp),
recu_indfp=scale(recu_indfp),emi_indpp=scale(emi_indpp),emi_indfp=scale(emi_indfp),
recu_totpp=scale(recu_totpp),recu_totfp=scale(recu_totfp),emi_totpp=scale(emi_totpp),
emi_totfp=scale(emi_totfp))

obj$sd=plyr::round_any(obj$sd,2)
obj$trait_sd=plyr::round_any(obj$trait_sd,0.05)
#Total effects
datt=subset(obj,alpha>0 & beta>0) %>% dplyr::group_by(interf,trait_sd) %>%
summarise(emi_totfp_sd=sd(emi_totfp,na.rm=T),emi_totpp_sd=sd(emi_totpp,na.rm=T),
emi_totfp=mean(emi_totfp,na.rm=T),emi_totpp=mean(emi_totpp,na.rm=T))

datp=subset(obj,alpha>0 & beta>0) %>% dplyr::group_by(interf,sd) %>%
summarise(emi_totfp_sd=sd(emi_totfp,na.rm=T),emi_totpp_sd=sd(emi_totpp,na.rm=T),
emi_totfp=mean(emi_totfp,na.rm=T),emi_totpp=mean(emi_totpp,na.rm=T))
pl1=ggplot(datt)+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_totpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point(data=datt,aes(x=(trait_sd-0.1)/0.8,y=emi_totfp),colour="black")+
geom_errorbar(data=datt,aes(x=(trait_sd-0.1)/0.8,ymin=emi_totfp-emi_totfp_sd,
ymax=emi_totfp+emi_totfp_sd),colour="black",width=0)+
geom_point(data=datp,aes(x=(sd-5)/35,y=emi_totfp),colour="red")+
geom_errorbar(data=datp,aes(x=(sd-5)/35,ymin=emi_totfp-emi_totfp_sd,
ymax=emi_totfp+emi_totfp_sd),colour="red",width=0)+
facet_wrap(~paste0("c = ",interf),scale="free_y",ncol=4)+
ylab("Standardized sum of total effects \n on other species (inter-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))
pl2=ggplot(data=obj)+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_totpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point(data=datt,aes(x=(trait_sd-0.1)/0.8,y=emi_totpp),colour="black")+
geom_errorbar(data=datt,aes(x=(trait_sd-0.1)/0.8,ymin=emi_totpp-emi_totpp_sd,ymax=emi_totpp+emi_totpp_sd),colour="black",width=0)+
geom_point(data=datp,aes(x=(sd-5)/35,y=emi_totpp),colour="red")+
geom_errorbar(data=datp,aes(x=(sd-5)/35,ymin=emi_totpp-emi_totpp_sd,
ymax=emi_totpp+emi_totpp_sd),colour="red",width=0)+
facet_wrap(~paste0("c = ",interf),scale="free_y",ncol=4)+
ylab("Standardized sum of total effects \n on other species (intra-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))
grid.arrange(pl1,pl2)


#Direct effects
datt=subset(obj,alpha>0 & beta>0) %>% dplyr::group_by(interf,trait_sd) %>%
summarise(emi_directfp_sd=sd(emi_directfp,na.rm=T),emi_directpp_sd=sd(emi_directpp,na.rm=T),
emi_directfp=mean(emi_directfp,na.rm=T),emi_directpp=mean(emi_directpp,na.rm=T))

datp=subset(obj,alpha>0 & beta>0) %>% dplyr::group_by(interf,sd) %>%
summarise(emi_directfp_sd=sd(emi_directfp,na.rm=T),emi_directpp_sd=sd(emi_directpp,na.rm=T),
emi_directfp=mean(emi_directfp,na.rm=T),emi_directpp=mean(emi_directpp,na.rm=T))
pl1=ggplot(datt)+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_directpp),size=0.6,alpha=0.4)+
geom_point(data=datt,aes(x=(trait_sd-0.1)/0.8,y=emi_directfp),colour="black")+
geom_errorbar(data=datt,aes(x=(trait_sd-0.1)/0.8,ymin=emi_directfp-emi_directfp_sd,
ymax=emi_directfp+emi_directfp_sd),colour="black",width=0)+
geom_point(data=datp,aes(x=(sd-5)/35,y=emi_directfp),colour="red")+
geom_errorbar(data=datp,aes(x=(sd-5)/35,ymin=emi_directfp-emi_directfp_sd,
ymax=emi_directfp+emi_directfp_sd),colour="red",width=0)+
facet_wrap(~paste0("c = ",interf),scale="free_y",ncol=4)+
ylab("Standardized sum of directal effects \n on other species (inter-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))
pl2=ggplot(data=obj)+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_directpp),size=0.6,alpha=0.4)+
geom_point(data=datt,aes(x=(trait_sd-0.1)/0.8,y=emi_directpp),colour="black")+
geom_errorbar(data=datt,aes(x=(trait_sd-0.1)/0.8,ymin=emi_directpp-emi_directpp_sd,ymax=emi_directpp+emi_directpp_sd),colour="black",width=0)+
geom_point(data=datp,aes(x=(sd-5)/35,y=emi_directpp),colour="red")+
geom_errorbar(data=datp,aes(x=(sd-5)/35,ymin=emi_directpp-emi_directpp_sd,
ymax=emi_directpp+emi_directpp_sd),colour="red",width=0)+
facet_wrap(~paste0("c = ",interf),scale="free_y",ncol=4)+
ylab("Standardized sum of directal effects \n on other species (intra-guild)")+
xlab("Phenology/trait generalism")+
scale_x_continuous(breaks=c(0,0.5,1))
grid.arrange(pl1,pl2)





library(Hmisc)
obj=tabsp %>% group_by(interf,alpha,beta,random,type) %>% summarise(
moy_sd=(wtd.mean(P_mean,N_eq1,na.rm=TRUE)-wtd.mean(P_mean,N,na.rm=TRUE))/wtd.mean(P_mean,N,na.rm=TRUE),
moy_tr=(wtd.mean(A_mean,N_eq1,na.rm=TRUE)-wtd.mean(A_mean,N,na.rm=TRUE))/wtd.mean(A_mean,N,na.rm=TRUE))
b=na.omit(obj) %>% group_by(interf,alpha,beta) %>% summarise(
moy_sd=mean(moy_sd,na.rm=T),
moy_tr=mean(moy_tr,na.rm=T))
#fwrite(b,"gene_par_net.txt",sep="\t",row.names=F)

# b=subset(tabsp,interf<=0.75) %>% group_by(interf,alpha,beta) %>% summarise(
# moy_sd=(mean(sd[surv>0])-mean(sd))/mean(sd),
# moy=(wtd.mean(mfd2,N_eq1)-wtd.mean(mfd2,N)),
# moy_tr=(mean(trait_sd[surv>0])-mean(trait_sd))/mean(trait_sd),nb=sum(surv),
# sd_moy_sd=sd(sd[surv>0]),
# sd_trait_sd=sd(trait_sd[surv>0]))

pl1=ggplot(data=b, aes(x=alpha, y=beta, fill=moy_sd))+
geom_raster(interpolate=F,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=14))+
scale_fill_gradientn(colours=c("firebrick4","pink","white","lightblue1","dodgerblue2","darkblue","darkblue"),
values=c(0,0.01,0.02,0.03,0.5,0.8,1),na.value="grey")+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+
labs(fill='')+
facet_wrap(~interf,ncol=4)
pl2=ggplot(data=b, aes(x=alpha, y=beta, fill=moy_tr))+
geom_raster(interpolate=F,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=14))+
scale_fill_gradientn(colours=c("firebrick4","pink","white","lightblue1","dodgerblue2","darkblue","darkblue"),
values=c(0,0.01,0.02,0.03,0.5,0.8,1),na.value="grey")+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+
labs(fill='')+
facet_wrap(~interf,ncol=4)

grid.arrange(pl1,pl2)

tabsp$trait_sd=(tabsp$trait_sd-0.1)/0.8
tabsp$sd=(tabsp$sd-5)/35
tabsp$alpha=as.factor(tabsp$alpha)
tabsp$beta=as.factor(tabsp$beta)
tabsp$mfd2=abs(tabsp$mfd-190)/190
tabsp$mfd2=plyr::round_any(tabsp$mfd,5)
obj=tabsp %>% dplyr::group_by(type) %>% dplyr::mutate(N_eq1=scale(N_eq1))


ggplot(data=subset(tabsp,alpha==0 & surv==0))+geom_density(aes(x=mfd2))+
stat_function(fun = dnorm, n = 101, args = list(mean = 190, sd = 70),col="red")+
facet_wrap(~interf+beta,ncol=5,scales="free")

model=lm(N_eq1~as.factor(interf)*(trait_sd*alpha+sd*beta+alpha*beta+mfd2*beta),
data=subset(obj,random<100))
newdat=data.frame(interf=as.factor(c(0,0.25,0.5,0.75)),trait_sd=rep(seq(0,1,0.1),each=4),
sd=mean(tabsp$sd),
beta=as.factor(0.5),alpha=as.factor(rep(c(0,0.25,0.5,0.75,1),each=11*4)),mfd2=0)
newdat$fit=predict(model,newdata=newdat,type="response")
ggplot(data=newdat,aes(x=trait_sd,y=fit,col=as.factor(alpha)))+geom_line()+
facet_wrap(~interf,scales="free")

newdat=data.frame(interf=as.factor(c(0,0.25,0.5,0.75)),sd=rep(seq(0,1,0.1),each=4),
trait_sd=mean(tabsp$trait_sd),
alpha=as.factor(0),beta=as.factor(rep(c(0,0.25,0.5,0.75,1),each=11*4)),mfd2=0)
newdat$fit=predict(model,newdata=newdat,type="response")
ggplot(data=newdat,aes(x=sd,y=fit,col=as.factor(beta)))+geom_line()+
facet_wrap(~interf,scales="free")


newdat=data.frame(interf=as.factor(c(0,0.25,0.5,0.75)),sd=0.5,
trait_sd=mean(tabsp$trait_sd),
alpha=as.factor(0),beta=as.factor(rep(c(0,0.25,0.5,0.75,1),each=11*4)),mfd2=rep(seq(0,1,0.1),each=4))
newdat$fit=predict(model,newdata=newdat,type="response")
ggplot(data=newdat,aes(x=mfd2,y=fit,col=as.factor(beta)))+geom_line()+
facet_wrap(~interf,scales="free")

grid.arrange(pl1,pl2)

ggplot(data = b, aes(x=alpha, y=beta, fill=moy_tr))+
geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=14))+
scale_fill_viridis(na.value="white")+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+
labs(fill='')+
facet_wrap(~interf,ncol=4)

b=subset(tabnet,stade=="eq1" & interff>0) %>% group_by(interff,alpha,beta) %>% summarise(
moy=mean(mod_I,na.rm=T),nb=mean(ntot),moy2=mean(indirect_inter,na.rm=T))
ggplot(data = b, aes(x=alpha, y=beta, fill=moy2))+
geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=14))+
scale_fill_viridis(na.value="white")+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+
labs(fill='')+
facet_wrap(~interff,ncol=4)


tabsp$sdr=plyr::round_any(tabsp$sd,1)
tabsp$mfdr=plyr::round_any(tabsp$mfd,10)
tabsp$I_meanr=plyr::round_any(tabsp$I_mean,0.02)
tabsp$I_sdr=plyr::round_any(tabsp$I_sd,0.02)
tabsp$trait_sdr=plyr::round_any(tabsp$trait_sd,0.05)
b=tabsp %>% group_by(interf,beta,alpha,sdr) %>% summarise(moy=mean(surv))
ggplot(b,aes(x=sdr,y=moy,col=as.factor(beta),fill=as.factor(beta)))+
geom_point()+facet_wrap(~alpha+interf,ncol=4,scales="free")+
stat_smooth()

b=tabsp %>% group_by(interf,beta,alpha,I_meanr,type) %>% summarise(moy=mean(surv))
ggplot(b,aes(x=I_meanr,y=moy,col=as.factor(beta),shape=type,linetype=type,fill=as.factor(beta)))+
geom_point()+facet_wrap(~alpha+interf,ncol=4,scales="free")+
stat_smooth()
b=subset(tabsp,interf==0.5) %>% group_by(beta,alpha,I_sdr,I_meanr,type) %>% summarise(moy=mean(surv))
ggplot(data = b, aes(x=I_sdr, y=I_meanr, fill=moy))+
geom_raster(interpolate=F,hjust=1,vjust=1)+theme_bw()+facet_wrap(~alpha+beta,ncol=5)+
scale_fill_viridis(na.value="white")

ggplot(b,aes(x=I_sdr,y=moy,col=as.factor(beta),shape=type,linetype=type,fill=as.factor(beta)))+
facet_wrap(~alpha+interf,ncol=4,scales="free")+
stat_smooth()

b=tabsp %>% group_by(interf,alpha,beta,trait_sdr,type) %>% summarise(moy=mean(surv))
ggplot(b,aes(x=trait_sdr,y=moy,col=as.factor(beta),shape=type,linetype=type,fill=as.factor(beta)))+
facet_wrap(~interf,ncol=4,scales="free")+
stat_smooth(alpha=0.1)

b=subset(tabsp,N_eq1>0 & interf==0.25 & alpha>0 & beta>0 & random<100) %>%
group_by(interf,alpha,beta,mfd,type,sp,random) %>%
do(data.frame(abond=dnorm(1:365,.$mfd,.$sd)*.$N_eq1,jour=1:365,sp=.$sp,type=.$type,
alpha=.$alpha,beta=.$beta,interf=.$interf,random=.$random))
b2=b %>% group_by(interf,alpha,beta,jour,type,random) %>% summarise(moy=sum(abond)/75)
b3=b2 %>% group_by(interf,alpha,beta,jour,type) %>% summarise(moy=mean(moy))
ggplot(b3,aes(x=jour,y=moy,col=as.factor(beta),shape=type,linetype=type))+
facet_wrap(~type+interf+alpha,ncol=4,scales="free")+
geom_point()

b=tabsp %>% group_by(interf,alpha,trait_sdr) %>% summarise(moy=mean(surv))
ggplot(b,aes(x=trait_sdr,y=moy,col=as.factor(alpha)))+geom_point()+facet_wrap(~interf,ncol=4)+
stat_smooth()



b=subset(tabsp,interf>0) %>% group_by(interf,beta,sdr) %>% summarise(moy=mean(surv))
pl1=ggplot(data = b, aes(x=beta, y=sdr, fill=moy))+
geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=14))+
scale_fill_viridis(na.value="white")+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+
labs(fill='')+
facet_wrap(~interf,ncol=4)

b=subset(tabsp,interf>0) %>% group_by(interf,alpha,trait_sdr) %>% summarise(moy=mean(surv))
pl1=ggplot(data = b, aes(x=alpha, y=trait_sdr, fill=moy))+
geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=14))+
scale_fill_viridis(na.value="white")+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+
labs(fill='')+
facet_wrap(~interf,ncol=4)

