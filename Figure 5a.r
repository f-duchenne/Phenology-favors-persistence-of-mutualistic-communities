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
library(piecewiseSEM)
library(qgraph)


setwd(dir="C:/Users/Duchenne/Documents/5eme - papier modelo/75x75")
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
tabnet$interff=tabnet$interfp
tabnet$interfp=paste0("c = ",tabnet$interfp)
tabsp$specialist=0
tabsp$specialist[tabsp$I_meaneq1<0.2]=1
# tabsp2=tabsp %>% group_by(interf,alpha,beta,random) %>% summarise(spec_moy=mean(I_meaneq1,na.rm=T),
# spec_somme=sum(specialist))
# names(tabsp2)[1]="interff"
# tabnet=merge(tabnet,tabsp2,by=c("interff","alpha","beta","random"),all.x=T,all.y=F)


####PROPAGATIOn tot
datt=subset(tabsp,beta==0 & alpha>0) %>% dplyr::group_by(interf,I_mean2,type) %>%
dplyr::summarise(emi_totfp_sd=sd(emi_totfp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),emi_totpp_sd=sd(emi_totpp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),
emi_totfp=mean(emi_totfp,na.rm=T),emi_totpp=mean(emi_totpp,na.rm=T),surv=mean(surv))
datt$colo="blue"
datp=subset(tabsp,alpha==0 & beta>0) %>% dplyr::group_by(interf,I_mean2,type) %>%
dplyr::summarise(emi_totfp_sd=sd(emi_totfp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),emi_totpp_sd=sd(emi_totpp,na.rm=T)/sqrt(length(na.omit(emi_directfp))),
emi_totfp=mean(emi_totfp,na.rm=T),emi_totpp=mean(emi_totpp,na.rm=T),surv=mean(surv))
datp$colo="yellow"
datt=rbind(datt,datp)
datt$emi_directpp[is.na(datt$emi_directpp)]=0
datt$type[datt$type=="flow"]="Plants"
datt$type[datt$type=="poll"]="Pollinators"

pl1=ggplot(datt,aes(shape=type,col=colo,fill=colo))+
# geom_point(aes(x=sd,y=emi_directpp),size=0.6,alpha=0.4,col="red")+
# geom_point(aes(x=sd,y=emi_totpp),size=0.6,alpha=0.4)+
# geom_point(aes(x=trait_sd,y=recu_totpp),size=0.6,alpha=0.4)+
geom_point(aes(x=I_mean2,y=emi_totpp),alpha=0.7)+
geom_errorbar(aes(x=I_mean2,ymin=emi_totpp-1.96*emi_totpp_sd,
ymax=emi_totpp+1.96*emi_totpp_sd),width=0)+
facet_wrap(~paste0("c = ",interf),ncol=4)+
ylab("Propagated effects within guilds")+
xlab("Degree at equilibrium")+
scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+
theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.text=element_text(size=12),axis.title=element_text(size=12),legend.position="none")+ggtitle("b")+
scale_shape_manual(values=c(21,22))+scale_color_manual(values=c("dodgerblue3","gold3"))+
scale_fill_manual(values=c("dodgerblue3","gold3"))

pdf("fig5b.pdf",width=12,height=3.5)
pl1
dev.off();


tabnet2=subset(tabnet,stade=="eq1")
biche=subset(tabnet2,interff>0 & !is.na(ai_contrib_ff))
biche$NODF_5=biche$NODF_5/100
biche$ntot=biche$ntot-75
biche$total=biche$alpha*biche$beta
biche=as.data.frame(biche)

# ggplot(data=biche,aes(x=ai_indirect_pp,y=ai_direct_pp,col=ai_contrib_pp))+geom_point()+scale_color_viridis()+
# facet_wrap(~interff,scales="free")

 model1=lme(ai_contrib_pp~alpha+beta+total+ntot+NODF_5,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 model2=lme(NODF_5~alpha+beta+total+ntot,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 model4=lme(ntot~alpha+beta+total,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 obj=piecewiseSEM::psem(model1,model2,model4,data=biche)
 objb=multigroup(obj, group = "interff",data=biche)
 

pdf("fig5.pdf",width=12,height=3)
par(mfrow=c(1,3))
for(i in c(0.25,0.5,0.75)){
center=10
left=0
right=20
haut=20
bas=0
l=objb$group.coefs[[paste(i)]][,c("Predictor","Response","Std.Estimate")]
l$colo="black"
l$colo[which(l$Std.Estimate<0)]="red"
l$Std.Estimate[l$Std.Estimate<0.01 & l$Std.Estimate>0]=l$Std.Estimate[l$Std.Estimate<0.01 & l$Std.Estimate>0]+0.001
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
g= g %>% set_vertex_attr("name", value =c("MS","PS","MS:PS","Diversity","Nestedness",
"\nIE contribution to total effects\namong pollinators"))
coord=data.frame(label=vertex_attr(g, "name"),
lab2=c("MS","PS","MS:PS","Diversity","Nestedness",
"\nIE contribution to total effects\namong pollinators"),
x=c(left,right,center,left+2,right-2,center),y=c(haut,haut,haut,bas+10,bas+10,bas))
EL=as_edgelist(g)
EL=cbind(EL,l[,3])
asi=abs(l[,3])/0.05
asi[asi<5]=5
qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=2.2,label.scale=F,
edge.label.cex = 2.2,edge.label.position=0.2,vsize2=6,vsize=49,title=paste0("c = ",i),
title.cex=2.2,shape="ellipse",edge.labels=T,fade=F,esize=max(abs(l[,3]))/0.05,asize=asi,
mar=c(5,3,3,3))
bidon=data.frame(total=c(l[1,3]*l[6,3]+l[5,3]*l[3,3],l[2,3]*l[6,3]+l[5,3]*l[4,3]),
type="interguild",interff=i,eff=c("TS","PS"))

# l=objb2$group.coefs[[paste(i)]][,c("Predictor","Response","Std.Estimate")]
# l$colo="black"
# l$colo[which(l$Std.Estimate<0)]="red"
# l[3,3]=l[3,3]+0.001
# g <- graph.data.frame(l, directed=T)
# g= g %>% set_edge_attr("color", value =l$colo)
# g= g %>% set_vertex_attr("name", value =c("MS","PS","MS:PS","Diversity","Nestedness",
# "Contribution of indirect effects\namong guilds"))
# coord=data.frame(label=vertex_attr(g, "name"),
# lab2=c("MS","PS","MS:PS","Diversity","Nestedness",
# "Contribution of indirect effects\namong guilds"),
# x=c(left,right,center,left+2,right-2,center),y=c(haut,haut,haut,bas+10,bas+10,bas))
# EL=as_edgelist(g)
# EL=cbind(EL,l[,3])
# qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
# border.color="white",label.cex=2.2,label.scale=F,
# edge.label.cex = 2.2,edge.label.position=0.2,vsize2=9,vsize=32,title=paste0("c = ",i),
# title.cex=2.2,shape="ellipse",edge.labels=T,fade=F,esize=max(abs(l[,3]))/0.1,asize=8)
# bidon=data.frame(total=c(l[1,3]*l[6,3]+l[5,3]*l[3,3],l[2,3]*l[6,3]+l[5,3]*l[4,3]),
# type="interguild",interff=i,eff=c("TS","PS"))

}
dev.off();


tabnet2=subset(tabnet,stade=="eq1")
biche=subset(tabnet2,interff>0 & !is.na(ai_contrib_ff))
biche$NODF_5=biche$NODF_5
biche$ntot=biche$ntot/150
biche$total=biche$alpha*biche$beta
library(piecewiseSEM)
 model1=lme(ai_direct_pp~alpha+beta+total+ntot+NODF_5,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 model2=lme(NODF_5~alpha+beta+total+ntot,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 model4=lme(ntot~alpha+beta+total,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 model5=lme(ai_indirect_pp~alpha+beta+total+ntot+NODF_5,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 obj=piecewiseSEM::psem(model1,model2,model4,model5)
 objb=multigroup(obj, group = "interff")
 
pdf("figS4.pdf",width=12,height=7)
par(mfrow=c(2,3))
for(i in c(0.25,0.5,0.75)){
center=10
left=0
right=20
haut=20
bas=0
l=objb$group.coefs[[paste(i)]][,c("Predictor","Response","Std.Estimate")]
l=l[grep("ai_indirect_pp",l$Response,invert=T),]
l$colo="black"
l$colo[which(l$Std.Estimate<0)]="red"
l$Std.Estimate[l$Std.Estimate<0.01 & l$Std.Estimate>0]=l$Std.Estimate[l$Std.Estimate<0.01 & l$Std.Estimate>0]+0.001
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
g= g %>% set_vertex_attr("name", value =c("MS","PS","MS:PS","Diversity","Nestedness",
"Direct effects among pollinators"))
coord=data.frame(label=vertex_attr(g, "name"),
lab2=c("MS","PS","MS:PS","Diversity","Nestedness",
"Direct effects among pollinators"),
x=c(left,right,center,left+2,right-2,center),y=c(haut,haut,haut,bas+10,bas+10,bas))
EL=as_edgelist(g)
EL=cbind(EL,l[,3])
asi=abs(l[,3])/0.05
asi[asi<5]=5
qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=2.2,label.scale=F,
edge.label.cex = 2.2,edge.label.position=0.2,vsize2=6,vsize=49,title=paste0("c = ",i),
title.cex=2.2,shape="ellipse",edge.labels=T,fade=F,esize=max(abs(l[,3]))/0.05,asize=asi)
bidon=data.frame(total=c(l[1,3]*l[6,3]+l[5,3]*l[3,3],l[2,3]*l[6,3]+l[5,3]*l[4,3]),
type="interguild",interff=i,eff=c("TS","PS"))
}
for(i in c(0.25,0.5,0.75)){
# l=objb2$group.coefs[[paste(i)]][,c("Predictor","Response","Std.Estimate")]
# l$colo="black"
# l$colo[which(l$Std.Estimate<0)]="red"
# l[3,3]=l[3,3]+0.001
# g <- graph.data.frame(l, directed=T)
# g= g %>% set_edge_attr("color", value =l$colo)
# g= g %>% set_vertex_attr("name", value =c("MS","PS","MS:PS","Diversity","Nestedness",
# "Contribution of indirect effects\namong guilds"))
# coord=data.frame(label=vertex_attr(g, "name"),
# lab2=c("MS","PS","MS:PS","Diversity","Nestedness",
# "Contribution of indirect effects\namong guilds"),
# x=c(left,right,center,left+2,right-2,center),y=c(haut,haut,haut,bas+10,bas+10,bas))
# EL=as_edgelist(g)
# EL=cbind(EL,l[,3])
# qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
# border.color="white",label.cex=2.2,label.scale=F,
# edge.label.cex = 2.2,edge.label.position=0.2,vsize2=9,vsize=32,title=paste0("c = ",i),
# title.cex=2.2,shape="ellipse",edge.labels=T,fade=F,esize=max(abs(l[,3]))/0.1,asize=8)
# bidon=data.frame(total=c(l[1,3]*l[6,3]+l[5,3]*l[3,3],l[2,3]*l[6,3]+l[5,3]*l[4,3]),
# type="interguild",interff=i,eff=c("TS","PS"))
center=10
left=0
right=20
haut=20
bas=0
l=objb$group.coefs[[paste(i)]][,c("Predictor","Response","Std.Estimate")]
l=l[grep("ai_direct_pp",l$Response,invert=T),]
l$colo="black"
l$colo[which(l$Std.Estimate<0)]="red"
l$Std.Estimate[l$Std.Estimate<0.01 & l$Std.Estimate>0]=l$Std.Estimate[l$Std.Estimate<0.01 & l$Std.Estimate>0]+0.001
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
g= g %>% set_vertex_attr("name", value =c("MS","PS","MS:PS","Diversity","Nestedness",
"Indirect effects among pollinators"))
coord=data.frame(label=vertex_attr(g, "name"),
lab2=c("MS","PS","MS:PS","Diversity","Nestedness",
"Indirect effects among pollinators"),
x=c(left,right,center,left+2,right-2,center),y=c(haut,haut,haut,bas+10,bas+10,bas))
EL=as_edgelist(g)
EL=cbind(EL,l[,3])
asi=abs(l[,3])/0.05
asi[asi<5]=5
asi[asi>12]=12
qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=2.2,label.scale=F,
edge.label.cex = 2.2,edge.label.position=0.2,vsize2=6,vsize=49,title=paste0("c = ",i),
title.cex=2.2,shape="ellipse",edge.labels=T,fade=F,esize=max(abs(l[,3]))/0.05,asize=asi)
bidon=data.frame(total=c(l[1,3]*l[6,3]+l[5,3]*l[3,3],l[2,3]*l[6,3]+l[5,3]*l[4,3]),
type="interguild",interff=i,eff=c("TS","PS"))
}
dev.off();






 model1=lme(ai_contrib_ff~alpha+beta+total+ntot+NODF_5,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 model2=lme(NODF_5~alpha+beta+total+ntot,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 model4=lme(ntot~alpha+beta+total,data=biche,random=~1|random,
 control = lmeControl(opt = "optim"))
 obj=piecewiseSEM::psem(model1,model2,model4)
 objb=multigroup(obj, group = "interff")
 
pdf("figS5.pdf",width=12,height=3)
par(mfrow=c(1,3))
for(i in c(0.25,0.5,0.75)){
center=10
left=0
right=20
haut=20
bas=0
l=objb$group.coefs[[paste(i)]][,c("Predictor","Response","Std.Estimate")]
l$colo="black"
l$colo[which(l$Std.Estimate<0)]="red"
l$Std.Estimate[l$Std.Estimate<0.01 & l$Std.Estimate>0]=l$Std.Estimate[l$Std.Estimate<0.01 & l$Std.Estimate>0]+0.001
g <- graph.data.frame(l, directed=T)
g= g %>% set_edge_attr("color", value =l$colo)
g= g %>% set_vertex_attr("name", value =c("MS","PS","MS:PS","Diversity","Nestedness",
"\nIE contribution to total effects\namong plants"))
coord=data.frame(label=vertex_attr(g, "name"),
lab2=c("MS","PS","MS:PS","Diversity","Nestedness",
"\nIE contribution to total effects\namong plants"),
x=c(left,right,center,left+2,right-2,center),y=c(haut,haut,haut,bas+10,bas+10,bas))
EL=as_edgelist(g)
EL=cbind(EL,l[,3])
asi=abs(l[,3])/0.05
asi[asi<5]=5
qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
border.color="white",label.cex=2.2,label.scale=F,
edge.label.cex = 2.2,edge.label.position=0.2,vsize2=6,vsize=49,title=paste0("c = ",i),
title.cex=2.2,shape="ellipse",edge.labels=T,fade=F,esize=max(abs(l[,3]))/0.05,asize=asi,
mar=c(5,3,3,3))
bidon=data.frame(total=c(l[1,3]*l[6,3]+l[5,3]*l[3,3],l[2,3]*l[6,3]+l[5,3]*l[4,3]),
type="interguild",interff=i,eff=c("TS","PS"))

# l=objb2$group.coefs[[paste(i)]][,c("Predictor","Response","Std.Estimate")]
# l$colo="black"
# l$colo[which(l$Std.Estimate<0)]="red"
# l[3,3]=l[3,3]+0.001
# g <- graph.data.frame(l, directed=T)
# g= g %>% set_edge_attr("color", value =l$colo)
# g= g %>% set_vertex_attr("name", value =c("MS","PS","MS:PS","Diversity","Nestedness",
# "Contribution of indirect effects\namong guilds"))
# coord=data.frame(label=vertex_attr(g, "name"),
# lab2=c("MS","PS","MS:PS","Diversity","Nestedness",
# "Contribution of indirect effects\namong guilds"),
# x=c(left,right,center,left+2,right-2,center),y=c(haut,haut,haut,bas+10,bas+10,bas))
# EL=as_edgelist(g)
# EL=cbind(EL,l[,3])
# qgraph(EL,layout=as.matrix(coord[,c("x","y")]),edge.color=l$colo,
# border.color="white",label.cex=2.2,label.scale=F,
# edge.label.cex = 2.2,edge.label.position=0.2,vsize2=9,vsize=32,title=paste0("c = ",i),
# title.cex=2.2,shape="ellipse",edge.labels=T,fade=F,esize=max(abs(l[,3]))/0.1,asize=8)
# bidon=data.frame(total=c(l[1,3]*l[6,3]+l[5,3]*l[3,3],l[2,3]*l[6,3]+l[5,3]*l[4,3]),
# type="interguild",interff=i,eff=c("TS","PS"))

}
dev.off();






tabnet2$spec_moy2=plyr::round_any(tabnet2$spec_moy,0.05)
tabnet2$ntot2=plyr::round_any(tabnet2$ntot,5)
b=tabnet2 %>% group_by(ntot2,spec_moy2,alpha,beta) %>% summarise(NODF_5=mean(NODF_5,na.rm=T))
ggplot(data=b,aes(x=spec_moy2,y=ntot2,fill=NODF_5))+geom_raster()+facet_wrap(~paste("MS = ",alpha)+
paste("PS = ",beta))+scale_color_viridis()


model=lme(contrib_intra~ntot*as.factor(interff),data=biche,random=~1|random)
biche$contrib_intra_c=residuals(model)

model=lme(contrib_inter~ntot*as.factor(interff),data=biche,random=~1|random)
biche$contrib_inter_c=residuals(model)

b=biche %>% group_by(interfp,alpha,beta) %>% summarise(moy=mean(contrib_intra_c),contrib_intra=mean(contrib_intra),
moy2=mean(contrib_inter_c),contrib_inter=mean(contrib_inter),NODF_5_r=mean(NODF_5_r),NODF_5=mean(NODF_5))
b=b %>% group_by(interfp) %>% mutate(average=mean(contrib_intra),average2=mean(contrib_inter))


pl1=ggplot(data = b, aes(x=alpha, y=beta, fill=NODF_5))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=14),axis.title.x=element_blank())+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
breaks=round(quantile(b$NODF_5_r,probs=c(0.05,0.5,0.95),names=F),digits=2),name="")+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+ggtitle("a",subtitle ="Intra-guild")+
coord_fixed(ratio = 1)+xlab("Morphological structure strength")+ylab("Phenological structure strength")

pl2=ggplot(data = b, aes(x=alpha, y=beta, fill=moy2+average2))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
axis.title=element_text(size=14,hjust=1),
plot.subtitle=element_text(size=14),axis.title.y=element_blank())+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
breaks=round(quantile(b$moy2+b$average2,probs=c(0.05,0.5,0.95),names=F),digits=2))+scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='')+ggtitle("b",subtitle ="Inter-guild")+
coord_fixed(ratio = 1)+xlab("Morphological structure strength")+ylab("Phenological structure strength")


plot_grid(pl1,pl2,align = "h",ncol=2)


