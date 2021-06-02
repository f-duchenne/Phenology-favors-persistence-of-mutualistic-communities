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
setwd(dir="C:/Users/duche/Documents/5eme - papier modelo")
jacob=as.matrix(read.table("jacob.txt",sep="\t",header=T))
inv=-1*solve(jacob)
inv2=inv
for(i in 1:nrow(jacob)){
for(j in 1:ncol(jacob)){
inv2[i,j]=inv[i,j]/(inv[i,i]*inv[j,j]-inv[i,j]*inv[j,i])
}}
indi=inv2-jacob
diag(jacob)=NA
diag(inv2)=NA
diag(indi)=NA
jacob=as.data.frame(jacob)
inv2=as.data.frame(inv2)
indi=as.data.frame(indi)
jacob$var1=names(jacob)
jacob=melt(jacob,id.vars="var1")
jacob$type="direct"
inv2$var1=names(inv2)
inv2=melt(inv2,id.vars="var1")
inv2$type="net"
indi$var1=names(indi)
indi=melt(indi,id.vars="var1")
indi$type="indirect"
jacob=rbind(jacob,inv2,indi)
setwd(dir="C:/Users/duche/Documents/5eme - papier modelo")
dyna=read.table("dyna.txt",sep="\t",header=T)
dyna=melt(dyna,id.vars="Time")
dyna$var2=as.numeric(gsub("X","",dyna$variable))
dyna$type="poll."
dyna$type[dyna$var2>75]="flow."
np=nrow(subset(dyna,Time==max(dyna$Time) & type=="poll." & value>0))+1
nf=nrow(subset(dyna,Time==max(dyna$Time) & type=="flow." & value>0))
jacob$var1=as.numeric(gsub("V","",jacob$var1))
jacob$variable=as.numeric(gsub("V","",jacob$variable))
milieu=min(jacob$value,na.rm=T)/(min(jacob$value,na.rm=T)-abs(max(jacob$value,na.rm=T)))


jaco=ggplot(data =jacob, aes(x=variable, y=var1, fill=value))+geom_raster(hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.text=element_blank(),axis.title=element_blank(),axis.text=element_blank(),
axis.ticks=element_blank(),legend.position="bottom",plot.title=element_text(size=18,face="bold",hjust = 0))+
facet_wrap(~type,ncol=1)+scale_fill_gradientn(colours=c("darkred","firebrick1","lightcoral","pink","white",
"darkturquoise","deepskyblue","dodgerblue4","darkblue"),na.value = "black",
values=c(0,milieu-0.2,milieu-0.01,milieu-0.0001,milieu,milieu+0.001,milieu+0.08,0.70,1))+labs(fill="")+
coord_fixed(ratio = 1,expand=F)+ggtitle("c")+
#annotate("rect", ymin=np,ymax=np+nf,xmin=0,xmax=np,col=NA,fill="white",alpha=0.8)+
#annotate("rect", ymin=0,ymax=np,xmin=np,xmax=np+nf,col=NA,fill="white",alpha=0.8)+
geom_vline(xintercept=np,size=1.2)+geom_hline(yintercept=np,size=1.2)

pdf("fig_methods.pdf",width=5,height=8)
jaco
dev.off();
#######################################################################
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
tabnet$interfp=paste0("c = ",tabnet$interfp)

biche=subset(tabnet,stade=="eq1")

bidon=melt(biche,id.vars=c("alpha","beta","interff"),measure.vars=c("ai_net_pp","ai_net_fp","ai_net_ff","ai_net_pf"))
bidon$cate="Between guilds"
bidon$cate[grep("fp",bidon$variable)]="Between guilds"
bidon$cate[grep("pp",bidon$variable)]="Within guild"
bidon$cate[grep("ff",bidon$variable)]="Within guild"
bidon$cate2="Received by pollinators"
bidon$cate2[grep("ff",bidon$variable)]="Received by plants"
bidon$cate2[grep("_pf",bidon$variable)]="Received by plants"
coucou=c("magenta2",brewer.pal(5,"Set2")[5],"red","dodgerblue4")
coucou=brewer.pal(9,"PuRd")[c(3,5,7,9)]
bidon=bidon %>% group_by(cate2,cate) %>% mutate(mini=boxplot.stats(value)$stats[1],
maxi=boxplot.stats(value)$stats[5])
bidon=subset(bidon,value>=mini & value<=maxi)

pl1=ggplot(data=bidon,aes(x=cate,y=value,col=as.factor(paste0("c = ",interff))))+geom_hline(yintercept=0)+
scale_color_manual(values=coucou)+
geom_boxplot(outlier.shape = NA)+facet_wrap(~cate2,scales="free")+
labs(color="")+xlab("")+ylab("Average total effects\n")+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=18,face="bold",hjust=0),strip.text=element_text(size=16),
axis.title=element_text(size=16),legend.title=element_text(size=14),axis.text.y=element_text(size=16),
axis.text.x=element_blank(),
plot.subtitle=element_text(size=14))+ggtitle("a")+labs(color="Competition strength")+
guides(colour =guide_legend(title.position="top", title.hjust = 0.5))
# Extract the legend. Returns a gtable
leg <- get_legend(pl1)
pl1=pl1+theme(legend.position="none")

biche$contrib_dir_pp=biche$ai_direct_pp/(abs(biche$ai_indirect_pp)+abs(biche$ai_direct_pp))
biche$contrib_dir_fp=biche$ai_direct_fp/(abs(biche$ai_indirect_fp)+abs(biche$ai_direct_fp))
biche$contrib_dir_pf=biche$ai_direct_pf/(abs(biche$ai_indirect_pf)+abs(biche$ai_direct_pf))
biche$contrib_dir_ff=biche$ai_direct_ff/(abs(biche$ai_indirect_ff)+abs(biche$ai_direct_ff))

bidon=melt(biche,id.vars=c("alpha","beta","interff"),measure.vars=c("ai_contrib_pp","ai_contrib_fp","ai_contrib_ff","ai_contrib_pf",
"contrib_dir_pp","contrib_dir_fp","contrib_dir_pf","contrib_dir_ff"))
bidon$cate="Between guilds"
bidon$cate[grep("fp",bidon$variable)]="Between guilds"
bidon$cate[grep("pp",bidon$variable)]="Within guild"
bidon$cate[grep("ff",bidon$variable)]="Within guild"
bidon$cate2="Received by pollinators"
bidon$cate2[grep("ff",bidon$variable)]="Received by plants"
bidon$cate2[grep("_pf",bidon$variable)]="Received by plants"
bidon$cate3="indirect"
bidon$cate3[grep("dir_",bidon$variable)]="direct"
coucou=c("magenta2",brewer.pal(5,"Set2")[5],"red","dodgerblue4")
coucou=brewer.pal(9,"PuRd")[c(3,5,7,9)]
#bidon=bidon %>% group_by(cate2,cate) %>% mutate(mini=boxplot.stats(value)$stats[1],
#maxi=boxplot.stats(value)$stats[5])
#bidon=subset(bidon,value>=mini & value<=maxi)
bidon=bidon %>% group_by(cate,cate2,cate3,interff,variable) %>%
summarise(moy=mean(value,na.rm=T),conf.low=mean(value,na.rm=T)-sd(value,na.rm=T),conf.high=mean(value,na.rm=T)+sd(value,na.rm=T))
bidon$interff=as.factor(paste0("c = ",bidon$interff))

pl2=ggplot(data=bidon,aes(x=cate,y=moy,fill=interff,col=interff,
group=interff))+
geom_hline(yintercept=0)+
scale_color_manual(values=coucou)+
scale_fill_manual(values=coucou)+
geom_bar(data=subset(bidon,cate3=="direct"),position=position_dodge(width=0.6),stat = "identity",width=0.5)+
facet_wrap(~cate2,scales="free")+
geom_errorbar(data=subset(bidon,cate3=="direct"),aes(ymin=conf.low,ymax=conf.high),width=0.1,position=position_dodge(width=0.6),col="black",alpha=1)+
geom_bar(data=subset(bidon,cate3=="indirect"),position=position_dodge(width=0.6),stat = "identity",width=0.5,fill="white")+
facet_wrap(~cate2,scales="free")+
geom_errorbar(data=subset(bidon,cate3=="indirect"),aes(ymin=conf.low,ymax=conf.high),width=0.1,position=position_dodge(width=0.6),col="black",alpha=1)+
labs(color="")+xlab("")+ylab("Average contributions of\nindirect effects to total effect")+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="none",
plot.title=element_text(size=18,face="bold",hjust=0),strip.text=element_text(size=16),
axis.title=element_text(size=16),legend.title=element_text(size=10),axis.text.y=element_text(size=16),
axis.text.x=element_text(size=16,angle=45,hjust=1),
plot.subtitle=element_text(size=14))+ggtitle("b")+scale_y_continuous(labels=percent)#+
#guides(fill="none", colour = "none",alpha=guide_legend(override.aes = list(color="dodgerblue4",fill=c("dodgerblue4","white"))))+
labs(alpha="Effect:")

pfi=plot_grid(pl1,pl2,ncol=1,rel_heights=c(3,5),align="v")

grid.arrange(pfi,ggpubr::as_ggplot(leg),ncol=1,heights=c(9,0.6))

pdf("fig4.pdf",width=7,height=9)
grid.arrange(pfi,ggpubr::as_ggplot(leg),ncol=1,heights=c(9,0.6))
dev.off();



##########FIGURE3:
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
library(ggplotify)
memory.limit(size = 1e6)
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
tabsp$I_mean2=plyr::round_any(tabsp$I_meaneq1,0.1,f=ceiling)-0.05
tabsp$I_mean3=plyr::round_any(tabsp$I_mean3,0.1,f=ceiling)-0.05
tabsp$P_mean2=plyr::round_any(tabsp$P_mean,0.1)
tabsp$A_mean2=plyr::round_any(tabsp$A_mean,0.1)
tabsp$sd2=plyr::round_any(tabsp$sd,2)
tabsp$trait_sd2=plyr::round_any(tabsp$trait_sd,0.1)
tabnet$interfp=paste0("c = ",tabnet$interfp)

subset(tabsp,N_eq1>0 & beta==1) %>% group_by(interf,type) %>% summarise(mfd_moy=mean(mfd),mfd_sd=sd(mfd),sd_moy=mean(sd),sd_sd=sd(sd)) 

setwd(dir="C:/Users/Duchenne/Documents/5eme - papier modelo")
dyna=read.table("dyna.txt",sep="\t",header=T)
dyna=melt(dyna,id.vars="Time")
dyna$var2=as.numeric(gsub("X","",dyna$variable))
dyna$type="poll."
dyna$type[dyna$var2>75]="flow."
pld=ggplot(data=dyna,aes(x=Time,y=value,col=as.factor(type),group=variable))+geom_line(alpha=0.5,size=0.7)+ scale_x_continuous(trans='log2')+
facet_wrap(~type,scales="free_y",ncol=1)+theme_bw()+ylab("Abundance")+theme(strip.text=element_blank(),legend.position="none",
panel.border = element_blank(),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),axis.text.x=element_text(size=12))+
scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x)))+
annotation_logticks(sides = "b")+scale_colour_manual(values=c("forestgreen","dodgerblue3"))+ggtitle("a")

library(network)
library(sna)
library(ggplot2)
library(ggnet)
res=fread("reseaux_figure.txt",sep="\t",header=T)
size=3
par(mfrow=c(2,2))
par(mar=c(0,0,0,0)+.1)
b=b[order(b$value),]
b=subset(res,alpha==0 & beta==1)
b=b[b$value>1,]
b=b[order(b$value),]
g=graph_from_data_frame(b, directed =F)
vec=vertex_attr(g,"name")
vec[grep("p",vec)]="dodgerblue3"
vec[grep("f",vec)]="forestgreen"
vertex_attr(g,"color")=vec
edge_attr(g,"weight")=b$value
E(g)$width <- E(g)$weight/10
vec=vertex_attr(g,"name")
vec[grep("p",vec)]=FALSE
vec[grep("f",vec)]=TRUE
vertex_attr(g,"type")=vec
edge_attr(g,"color")=adjustcolor(paste0("gray",round(101-100*E(g)$weight/max(E(g)$weight))))
lay=layout_as_bipartite(g,maxiter = 1,hgap=10000,vgap=1)
gr1=plot_grid(base2grob(~plot(g,layout=lay,vertex.label=NA,vertex.size=size,label.cex=0,size2=size,curved=T)))+
theme(plot.margin = unit(c(0,0,0,0), "cm"))+coord_cartesian(expand=F)

b=subset(res,alpha==1 & beta==1)
b=b[b$value>1,]
b=b[order(b$value),]
g=graph_from_data_frame(b, directed =F)
vec=vertex_attr(g,"name")
vec[grep("p",vec)]="dodgerblue3"
vec[grep("f",vec)]="forestgreen"
vertex_attr(g,"color")=vec
edge_attr(g,"weight")=b$value
E(g)$width <- E(g)$weight/5
vec=vertex_attr(g,"name")
vec[grep("p",vec)]=FALSE
vec[grep("f",vec)]=TRUE
vertex_attr(g,"type")=vec
#vertex_attr(g,"name")=""
edge_attr(g,"color")=adjustcolor(paste0("gray",round(101-100*E(g)$weight/max(E(g)$weight))))
lay=layout_as_bipartite(g,maxiter = 1000,hgap=10000,vgap=1)
gr2=plot_grid(base2grob(~plot(g,layout=lay,vertex.label=NA,vertex.size=size,label.cex=0,size2=size,curved=T)))+
theme(plot.margin = unit(c(0,0,0,0), "cm"))+coord_cartesian(expand=F)

b=subset(res,alpha==0 & beta==0)
b=b[b$value>1,]
b=b[order(b$value),]
g=graph_from_data_frame(b, directed =F)
vec=vertex_attr(g,"name")
vec[grep("p",vec)]="dodgerblue3"
vec[grep("f",vec)]="forestgreen"
vertex_attr(g,"color")=vec
edge_attr(g,"weight")=b$value
E(g)$width <- E(g)$weight/5
vec=vertex_attr(g,"name")
vec[grep("p",vec)]=FALSE
vec[grep("f",vec)]=TRUE
vertex_attr(g,"type")=vec
#vertex_attr(g,"name")=""
edge_attr(g,"color")=adjustcolor(paste0("gray",round(101-100*E(g)$weight/max(E(g)$weight))))
lay=layout_as_bipartite(g,maxiter = 1000,hgap=10000,vgap=1)
gr3=plot_grid(base2grob(~plot(g,layout=lay,vertex.label=NA,vertex.size=size,label.cex=0,size2=size,curved=T)))+
theme(plot.margin = unit(c(0,0,0,0), "cm"))+coord_cartesian(expand=F)


b=subset(res,alpha==1 & beta==0)
b=b[b$value>1,]
b=b[order(b$value),]
g=graph_from_data_frame(b, directed =F)
vec=vertex_attr(g,"name")
vec[grep("p",vec)]="dodgerblue3"
vec[grep("f",vec)]="forestgreen"
vertex_attr(g,"color")=vec
edge_attr(g,"weight")=b$value
E(g)$width <- E(g)$weight/5
vec=vertex_attr(g,"name")
vec[grep("p",vec)]=FALSE
vec[grep("f",vec)]=TRUE
vertex_attr(g,"type")=vec
#vertex_attr(g,"name")=""
edge_attr(g,"color")=adjustcolor(paste0("gray",round(101-100*E(g)$weight/max(E(g)$weight))))
lay=layout_as_bipartite(g,maxiter = 1000,hgap=10000,vgap=1)
gr4=plot_grid(base2grob(~plot(g,layout=lay,vertex.label=NA,vertex.size=size,label.cex=0,size2=size,curved=T)))+
theme(plot.margin = unit(c(0,0,0,0), "cm"))+coord_cartesian(expand=F)

bla=ggplot() + theme_void()
gr=plot_grid(bla,bla,bla,gr1,bla,gr2,bla,bla,bla,gr3,bla,gr4,bla,bla,bla,ncol=3,rel_heights=c(0.05,1,-0.3,1,-0.3),align = "hv",
rel_widths = c(1, -0.3, 1),labels=c("","","","PF = 1 & MF = 0","","PF = 1 & MF = 1","","","","PF = 0 & MF = 0","",
"PF = 0 & MF = 1"),vjust=c(rep(0.5,8),rep(2,8)),label_size=10)+
ggtitle("b")+theme(plot.title=element_text(size=14,face="bold",hjust = 0))
gr



b=subset(res,alpha==0 & beta==0)
mat=as.data.frame(dcast(b,Var1~Var2,value.var="value"))
rownames(mat)=as.character(mat[,1])
obj=nestedrank(as.matrix(mat[,-1]), method = "NODF", weighted=TRUE, normalise=TRUE, return.matrix=T)
b$Var1=factor(b$Var1,rownames(obj$nested.matrix))
b$Var2=factor(b$Var2,rev(colnames(obj$nested.matrix)))
pl1=ggplot(data=b,aes(x=Var1,y=Var2,fill=sqrt(value)))+geom_tile()+
scale_fill_viridis()+theme(axis.text=element_blank(),axis.ticks=element_blank())


coucou=c("magenta2",brewer.pal(5,"Set2")[5],"red","dodgerblue4")
b=subset(tabnet,stade=="eq1")
pl_flow=ggplot(data=b, aes(x=ntot, y=nflow, color=interfp))+geom_jitter(alpha=0.5,size=0.2)+
scale_color_manual(values=coucou)+geom_abline(intercept=0,slope=0.5)+theme_bw()+
theme(strip.text=element_blank(),legend.position="right",
panel.border = element_blank(),axis.line = element_line(colour = "black"),panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),axis.text=element_text(size=12))+
labs(color="")+ylab("Number of plant species at equilibrium")+xlab("Number of species at equilibrium")+
ggtitle("b")+guides(color = guide_legend(override.aes = list(size=2)))

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
facet_wrap(~paste0("c = ",interf),ncol=2)+
ylab("Average persistence probability")+
xlab("Initial generalism level")+
scale_x_continuous(breaks=c(0,0.5,1))+theme_bw()+scale_y_continuous(breaks=c(0,0.5,1))+
theme(strip.background=element_blank(),plot.title=element_text(size=14,face="bold",hjust = 0),
strip.text=element_text(size=12),axis.title=element_text(size=12),legend.position="none")+ggtitle("c")+
scale_shape_manual(values=c(21,22))


b=subset(tabnet,stade=="eq1") %>% dplyr::group_by(interfp,alpha,beta) %>% dplyr::summarise(moy=mean(ntot)/150)
pl1=ggplot(data = b, aes(x=alpha, y=beta, fill=moy))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title.x=element_text(size=14),
plot.subtitle=element_text(size=14),axis.title.y=element_blank())+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
labels = percent(c(0.5*0:2)),breaks = 0.5*0:2,limits=c(0,1))+scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='')+ggtitle("",subtitle ="Persistence")+
coord_fixed(ratio = 1,expand = FALSE)+xlab("Morphological forcing")

b=subset(tabnet,stade=="eq1") %>% dplyr::group_by(interfp,alpha,beta) %>% dplyr::summarise(moy=length(ntot[ntot>0])/length(ntot))
pl2=ggplot(data = b, aes(x=alpha, y=beta, fill=moy))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="bottom",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),axis.title=element_text(size=14),
plot.subtitle=element_text(size=14),axis.title.x=element_blank())+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
labels = percent(c(0.5*0:2)),breaks = 0.5*0:2,limits=c(0,1))+scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='')+ggtitle("a",subtitle ="Viability")+
coord_fixed(ratio = 1,expand = FALSE)+xlab("Morphological forcing")+ylab("Phenological forcing")


virgo=c(viridis(100)[1:50],rep(viridis(100)[50],200),viridis(100)[51:100])
b=subset(tabnet,stade=="eq1") %>% dplyr::group_by(interfp,alpha,beta) %>% dplyr::summarise(moy=mean(NODF_5,na.rm=T))
pl3=ggplot(data=b, aes(x=alpha, y=beta, fill=moy))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title.y=element_blank(),legend.position="bottom",
strip.text=element_text(size=12),axis.title.x=element_blank(),plot.subtitle=element_text(size=14))+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),na.value="white")+ggtitle("",subtitle ="Nestedness")+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='')+
coord_fixed(ratio = 1,expand = FALSE)+xlab("Morphological forcing")


pfi=plot_grid(pl2,pl1,pl3,ncol=3,align = "h")
grid.arrange(pfi,gr,pl1b,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(1.6,1.6))

pdf("fig3.pdf",width=12,height=9)
grid.arrange(pfi,gr,pl1b,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(1.6,1.6))
dev.off();

png("fig3.png",width=1500,height=1100,res=140)
grid.arrange(pfi,gr,pl1b,layout_matrix = rbind(c(1,2),c(1,3)),widths=c(1,1),heights=c(1,1.3))
dev.off();


###FIGURE S1:
b=subset(tabnet,tabnet$stade=="eq1" & interff==0) %>% group_by(interfp,alpha,beta) %>% 
summarise(moy=mean(ntot)/150)
pl1=ggplot(data = b, aes(x=alpha, y=beta, fill=moy))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=14),axis.title=element_blank())+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),
labels = percent)+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='')+ggtitle("",subtitle ="Persistence")+
coord_fixed(ratio = 1)+xlab("Trait structure strength")

b=subset(tabnet,tabnet$stade=="eq1" & interff==0) %>% group_by(interfp,alpha,beta) %>% summarise(moy=length(ntot[ntot>0])/length(ntot))
pl2=ggplot(data = b, aes(x=alpha, y=beta, fill=moy))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),legend.position="right",
plot.title=element_text(size=14,face="bold",hjust = 0),strip.text=element_text(size=12),
plot.subtitle=element_text(size=14),axis.title=element_blank())+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours="#440154FF",labels=percent(1))+scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='')+ggtitle("",subtitle ="Coexistence")+
coord_fixed(ratio = 1)


virgo=c(viridis(100)[1:50],rep(viridis(100)[50],200),viridis(100)[51:100])
b=subset(tabnet,tabnet$stade=="eq1" & interff==0) %>% group_by(interfp,alpha,beta) %>% summarise(moy=mean(valprop,na.rm=T))
pl3=ggplot(data=b, aes(x=alpha, y=beta, fill=log(abs(moy))))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title.y=element_blank(),legend.position="right",
strip.text=element_text(size=12),,axis.title=element_blank(),plot.subtitle=element_text(size=14))+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))),breaks=c(-5.85,-5.9),na.value="white")+
ggtitle("",subtitle ="Resilience")+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='')+
coord_fixed(ratio = 1)+xlab("Trait structure strength")


grid.arrange(pl2,pl1,pl3,ncol=3,left="Phenological structure strength",
bottom="Trait structure strength")

pdf("fig_s1.pdf",width=8,height=9)
pfi
dev.off();



b=subset(tabnet,tabnet$stade=="eq1") %>% group_by(interfp,alpha,beta) %>% summarise(moy=mean(connectance))
pl2=ggplot(data = b, aes(x=alpha, y=beta, fill=moy))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title=element_blank())+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))))+scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='NODF')+
coord_fixed(ratio = 1)



b=subset(tabnet,tabnet$stade=="eq1") %>% group_by(interfp,alpha,beta) %>% summarise(moy=mean(indirect_poll))
pl4=ggplot(data = b, aes(x=alpha, y=beta, fill=sqrt(moy)))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title=element_blank())+
facet_wrap(~interfp,ncol=1)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))))+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='Intra-guild NE')+
coord_fixed(ratio = 1)

b=subset(tabnet,tabnet$stade=="eq1") %>% group_by(interfp,alpha,beta) %>% summarise(moy=mean(indirect_flow))
pl5=ggplot(data = b, aes(x=alpha, y=beta, fill=sqrt(moy)))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+
theme(panel.grid=element_blank(),strip.background=element_blank(),axis.title=element_blank())+
facet_wrap(~interfp,ncol=4)+scale_fill_gradientn(colours=rev(c(viridis(100),
colorRampPalette(c("yellow","orange","red","darkred"))(66))))+
scale_y_continuous(breaks=c(0,1.25),labels=c("0","1"))+
scale_x_continuous(breaks=c(0,1.25),labels=c("0","1"))+labs(fill='LPE strength')+
coord_fixed(ratio = 1)

pfi=plot_grid(pl1,pl3,ncol=2,align = "h")
#create common x and y labels
y.grob <- textGrob("Common Y",gp=gpar(fontface="bold", col="blue", fontsize=15), rot=90)
x.grob <- textGrob("Common X",gp=gpar(fontface="bold", col="blue", fontsize=15))
#add to plot
grid.arrange(arrangeGrob(plot, left = y.grob, bottom = x.grob))

#SPECIES:
ggplot()+geom_histogram(data=subset(tabsp,interf==0.75 & alpha==1 & beta==0.5),aes(x=sd,col=type,fill=type),linetype="dashed")+theme_bw()+
geom_histogram(data=subset(tabsp,interf==0.75 & alpha==1 & beta==0.5 & N_eq1>0),aes(x=sd,col=type),alpha=1,fill="white")+theme(panel.grid=element_blank())+
xlab("Mean flight/flowering date")+scale_color_brewer(palette="Set1")

tabsp$surv=0
tabsp$surv[tabsp$N_eq1>0]=1
tabsp$mfd_cen=abs(tabsp$mfd-190)
tabsp$trait_mu=abs(tabsp$trait_mu-0)
tabsp$trait_sd=plyr::round_any(tabsp$trait_sd,0.1)


newdat=tabsp %>% group_by(trait_sd,alpha,interf,type) %>% summarise(moy=mean(surv),sdev=sd(surv),moy_N=mean(N_eq1))
ggplot(data =subset(newdat,interf==0.5), aes(x=trait_sd, y=moy,col=as.factor(alpha),shape=type))+geom_line()+theme_bw()+theme(panel.grid=element_blank())+
facet_wrap(~type+beta,ncol=5,scales="free_y")

newdat=tabsp %>% group_by(trait_sd,alpha,interf,type) %>% summarise(moy=mean(surv),sdev=sd(surv),moy_N=mean(N_eq1))
ggplot(data =subset(newdat,interf==0.5), aes(x=alpha, y=trait_sd, fill=moy))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+theme(panel.grid=element_blank())+
facet_wrap(~type)+scale_fill_viridis()+scale_x_continuous(breaks=c(0,1.2),labels=c("0","1"))+xlab("Trait structure strength")+
ylab("Generalism")

newdat=tabsp %>% group_by(trait_sd,beta,interf,type) %>% summarise(moy=mean(surv),sdev=sd(surv),moy_N=mean(N_eq1))
ggplot(data =subset(newdat,interf==0.5), aes(x=beta, y=trait_sd, fill=moy))+geom_raster(interpolate=T,hjust=1,vjust=1)+theme_bw()+theme(panel.grid=element_blank())+
facet_wrap(~type)+scale_fill_viridis()+scale_x_continuous(breaks=c(0,1.2),labels=c("0","1"))+xlab("Pheno. structure strength")+
ylab("Generalism")

model=glm(surv~((scale(mfd_cen)+scale(sd))*beta+(scale(trait_sd)+scale(trait_mu))*alpha+beta*alpha)*interf,data=tabsp,family=binomial)
model=step(model,direction = c("backward"))


newdat=data.frame(alpha=c(0,0.25,0.5,0.75,1),beta=rep(c(0,0.25,0.5,0.75,1),each=5),interf=rep(c(0,0.25,0.5,0.75),each=25),
sd=30,mfd_cen=190,trait_sd=rep(seq(0.1,1,0.1),each=25*4),trait_mu=rep(0:60,25*4*10))
newdat$fit=predict(model,newdata=newdat,type="response")
newdat$interff=paste0("c=",newdat$interf)
newdat$colo=paste0(newdat$beta," - ",newdat$alpha)
newdat=newdat[order(newdat$beta,newdat$alpha),]
newdat$colo=factor(newdat$colo,unique(newdat$colo))
library(RColorBrewer)
ggplot(data=subset(newdat,trait_mu==0),aes(x=trait_sd,y=fit,col=colo))+geom_line(size=1.2)+facet_wrap(~interff)+
theme_bw()+theme(strip.background =element_rect(fill="white"))+scale_color_manual(values=c(brewer.pal(9,"YlOrRd")[3:7],brewer.pal(9,"Greens")[3:7],
brewer.pal(9,"Blues")[3:7],brewer.pal(9,"BuPu")[3:7],brewer.pal(9,"Greys")[3:7]),name="PS - TS values")+ylab("Persistence probability")







newdat=data.frame(mfd=1:365,sd=mean(tab$sd),alpha=rep(c(0,0.25,0.5,0.75,1),each=365),trait_sd=mean(tab$trait_sd),
trait_mu=mean(tab$trait_mu),beta=0.25,compet=rep(c(0,0.25),each=365*5))
newdat$fit=predict(model,newdata=newdat,type="response")
ggplot(data=newdat,aes(x=mfd,y=fit,col=as.factor(alpha)))+geom_line()+facet_wrap(~compet)


newdat=data.frame(mfd=mean(tab$mfd),sd=seq(min(tab$sd),max(tab$sd),length.out=365),alpha=rep(c(0,0.25,0.5,0.75,1),each=365),trait_sd=mean(tab$trait_sd),
trait_mu=mean(tab$trait_mu),beta=0.25,compet=rep(c(0,0.25),each=365*5))
newdat$fit=predict(model,newdata=newdat,type="response")
ggplot(data=newdat,aes(x=sd,y=fit,col=as.factor(alpha)))+geom_line()+facet_wrap(~compet)

























