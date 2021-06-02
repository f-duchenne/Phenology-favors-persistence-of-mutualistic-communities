.libPaths(c("/home/duchenne/R/x86_64-pc-linux-gnu-library/3.4/",.libPaths()))
library(truncnorm)
library(bipartite)
library(Rmpfr)
library(circular)
library(CircStats)
library(doBy)
library(gridExtra)
library(odeintr)
library(data.table)
library(doParallel)
library(foreach)
library(parallel)
library(doMPI)
library(igraph)
cl<-startMPIcluster()
registerDoMPI(cl)


comb <- function(x, ...) {  
      mapply(rbind,x,...,SIMPLIFY=FALSE)
}


interfp=0.75
interff=0.75


alpha_vec=c(0,0.25,0.5,0.75,1)
beta_vec=alpha_vec
moys=expand.grid(beta_vec,alpha_vec)
fini=1000
resultat=foreach(jj=1:fini,.combine=comb)%dopar%{
library(truncnorm)
library(bipartite)
library(plot3D)
library(Rmpfr)
library(circular)
library(CircStats)
library(doBy)
library(phangorn)
library(gridExtra)
library(deSolve)
library(rootSolve)
library(expm)
library(igraph)
deb=1
list_P=list()
setwd(dir="/home/duchenne/part_1/initial_networks_75x75")
final=read.table(paste("pops_eq_i_",jj,".txt",sep=""),sep="\t",header=T)
precision=0.1 #precision for integrals
popf_p=subset(final,type=="poll" & random==jj)$N
popf_f=subset(final,type=="flow" & random==jj)$N
mu_p=subset(final,type=="poll" & random==jj)$mfd
mu_f=subset(final,type=="flow" & random==jj)$mfd
sd_p=subset(final,type=="poll" & random==jj)$sd
sd_f=subset(final,type=="flow" & random==jj)$sd
genmu_p=subset(final,type=="poll" & random==jj)$trait_mu
genmu_f=subset(final,type=="flow" & random==jj)$trait_mu
gensd_p=subset(final,type=="poll" & random==jj)$trait_sd
gensd_f=subset(final,type=="flow" & random==jj)$trait_sd

nbsp_p=length(popf_p)
nbsp_f=length(popf_f)
Nini=c(popf_p,popf_f)


P=matrix(0,nbsp_f,nbsp_p)
colnames(P)=1:nbsp_p
rownames(P)=1:nbsp_f
for(i in 1:nbsp_p){
for(i2 in 1:nbsp_f){
func <- function(x) {
  f1 <- dwrpnorm(circular(rad(x/(365/360))),circular(rad(mu_p[i]/(365/360))),sd=circular(rad(sd_p[i]/(365/360))))
  f2 <- dwrpnorm(circular(rad(x/(365/360))),circular(rad(mu_f[i2]/(365/360))),sd=circular(rad(sd_f[i2]/(365/360))))
  pmin(f1, f2)
}
#P[i2,i]=integrateR(func,0,365,abs.tol=precision)[[1]]/(365/(2*pi))
P[i2,i]=sum(func(seq(0,365,precision)))*precision/(365/(2*pi))
}
}
P2=P

A=P
for(i in 1:nbsp_p){
for(i2 in 1:nbsp_f){
func <- function(x) {
 f1 <- dnorm(x,genmu_p[i],gensd_p[i])
 f2 <- dnorm(x,genmu_f[i2],gensd_f[i2])
  return(list(inte=pmin(f1, f2),poll=f1,flow=f2))
}
inte=func(seq(-10,10,0.01))
A[i2,i]=sum(inte$inte)*0.01#*((1-sqrt(gensd_p[i]*gensd_f[i2]))^0.3)
}
}
#A[,]=sort(P)[rank(A)]
A2=A


#Competition matrix:
C_p_phe=matrix(0,nbsp_p,nbsp_p)
colnames(C_p_phe)=1:nbsp_p
rownames(C_p_phe)=1:nbsp_p
for(i in 1:nbsp_p){
for(i2 in i:nbsp_p){
func <- function(x) {
  f1 <- dwrpnorm(circular(rad(x/(365/360))),circular(rad(mu_p[i]/(365/360))),sd=circular(rad(sd_p[i]/(365/360))))
  f2 <- dwrpnorm(circular(rad(x/(365/360))),circular(rad(mu_p[i2]/(365/360))),sd=circular(rad(sd_p[i2]/(365/360))))
  pmin(f1, f2)
}
#pm=integrateR(func,0,365,abs.tol=precision)[[1]]/(365/(2*pi))
pm=sum(func(seq(1,365,precision)))*precision/(365/(2*pi))
C_p_phe[i2,i]=pm
}
}
#C_p=round(C_p,digits=digi)
# image2D(C_p)
# image2D(t(apply(C_p,2,rev)),col=rev(heat.colors(n=150)))
# image2D(t(apply(C_p_phe,2,rev)))
C_p_phe=C_p_phe+t(C_p_phe)
diag(C_p_phe)=1
#C_p_phe=round(C_p_phe,digits=digi)


C_f_phe=matrix(0,nbsp_f,nbsp_f)
colnames(C_f_phe)=1:nbsp_f
rownames(C_f_phe)=1:nbsp_f
for(i in 1:nbsp_f){
for(i2 in i:nbsp_f){
func <- function(x) {
  f1 <- dwrpnorm(circular(rad(x/(365/360))),circular(rad(mu_f[i]/(365/360))),sd=circular(rad(sd_f[i]/(365/360))))
  f2 <- dwrpnorm(circular(rad(x/(365/360))),circular(rad(mu_f[i2]/(365/360))),sd=circular(rad(sd_f[i2]/(365/360))))
  pmin(f1, f2)
}
#pm=integrateR(func,0,365,abs.tol=precision)[[1]]/(365/(2*pi))
pm=sum(func(seq(1,365,precision)))*precision/(365/(2*pi))
C_f_phe[i2,i]=pm
}
}
#C_f=round(C_f,digits=digi)
# image2D(C_f)
C_f_phe=C_f_phe+t(C_f_phe)
diag(C_f_phe)=1

C_p_phe2=C_p_phe
C_f_phe2=C_f_phe

for(advi in deb:nrow(moys)){
setwd(dir="/home/duchenne/part_1/initial_networks_75x75")
final=read.table(paste("pops_eq_i_",jj,".txt",sep=""),sep="\t",header=T)

popf_p=subset(final,type=="poll" & random==jj)$N
popf_f=subset(final,type=="flow" & random==jj)$N

mu_p=subset(final,type=="poll" & random==jj)$mfd
mu_f=subset(final,type=="flow" & random==jj)$mfd
sd_p=subset(final,type=="poll" & random==jj)$sd
sd_f=subset(final,type=="flow" & random==jj)$sd
genmu_p=subset(final,type=="poll" & random==jj)$trait_mu
genmu_f=subset(final,type=="flow" & random==jj)$trait_mu
gensd_p=subset(final,type=="poll" & random==jj)$trait_sd
gensd_f=subset(final,type=="flow" & random==jj)$trait_sd

seuil=1e-5 #extinction cutoff
precision=0.1 #precision for integrals
conv=1e-9 #convergence criterion

nbsp_p=length(popf_p)
nbsp_f=length(popf_f)

rmax <- runif(nbsp_p,1,1); rmax2<-runif(nbsp_f,1,1)
m=subset(final,type=="poll" & random==jj)$m
m2=subset(final,type=="flow" & random==jj)$m
K=subset(final,type=="poll" & random==jj)$K
K2=subset(final,type=="flow" & random==jj)$K
cs=subset(final,type=="poll" & random==jj)$cs
cs2=subset(final,type=="flow" & random==jj)$cs
beta=moys$Var1[advi]
alpha=moys$Var2[advi]
hangling=subset(final,random==jj)$hangling
efficience=subset(final,random==jj)$efficience

P=P2^beta
A=A2^alpha
I=A*P
C_p_phe=C_p_phe2^beta
C_f_phe=C_f_phe2^beta


#fonctions
derivs <-function(t, y,parms){
dy=rep(0,1)
for(i in 1:length(Nini)){
if(i<=nbsp_p){
I2=I[,i]*I*y[(1+nbsp_p):(nbsp_p+nbsp_f)]
vec=(C_p_phe[,i]*(apply(I2,2,sum))/sum(y[(1+nbsp_p):(nbsp_p+nbsp_f)]*I[,i]))
vec[which(is.na(vec))]=0
eq1=((-y[i]/K[i])+ #densite dependance incluant competition pour la place
efficience[i]*sum(I[,i]*y[(1+nbsp_p):(nbsp_p+nbsp_f)])/(1+hangling[i]*sum(I[,i]*y[(1+nbsp_p):(nbsp_p+nbsp_f)])+
interfp*sum(vec*y[1:nbsp_p]))-  #De Angelis-Beddington functional response
m[i])*  #mortalitee
y[i] #croissance logistique incluant compÃƒÂ©tition pour les sites de nidif (incluant phenologie et signal phylo)

myenv=new.env(i)
environment(eq1)=myenv
myenv$i=i
myenv$Nini=Nini
myenv$I=I
myenv$I2=I2
myenv$vec=vec
myenv$K=K
myenv$rmax=rmax
myenv$nbsp_p=nbsp_p
myenv$nbsp_f=nbsp_f
myenv$m=m
myenv$cs=cs
myenv$C_f_phe=C_f_phe
myenv$C_p_phe=C_p_phe
myenv$interfp=interfp
}else{
I2=t(I)[,(i-nbsp_p)]*t(I)*y[1:nbsp_p]
vec=C_f_phe[,(i-nbsp_p)]*(apply(I2,2,sum))/sum(y[1:nbsp_p]*I[(i-nbsp_p),])
vec[which(is.na(vec))]=0
eq1=(-y[i]/K2[(i-nbsp_p)]+ #densite dependance incluant competition pour la place
efficience[i]*sum(I[(i-nbsp_p),]*y[1:nbsp_p])/(1+hangling[i]*sum(I[(i-nbsp_p),]*y[1:nbsp_p])+
interff*sum(vec*y[(1+nbsp_p):(nbsp_p+nbsp_f)]))- #De Angelis-Beddington functional response
m2[(i-nbsp_p)])*  #mortalite
y[i]

myenv=new.env(i)
environment(eq1)=myenv
myenv$i=i
myenv$I=I
myenv$Nini=Nini
myenv$I2=I2
myenv$vec=vec
myenv$K2=K2
myenv$rmax2=rmax2
myenv$nbsp_p=nbsp_p
myenv$nbsp_f=nbsp_f
myenv$m2=m2
myenv$cs2=cs2
myenv$C_f_phe=C_f_phe
myenv$C_p_phe=C_p_phe
myenv$interff=interff
}
dy[i]=eq1}
return(list(dy))}

a=1
b=0
Nini=c(popf_p,popf_f)


####NETWORK LEVEL INFORMATIONS TO EXPORT
dig_num=16

### AGREGER:
mati=popf_f%*%t(popf_p)
resf=mati*I
nbsp_p=length(popf_p)
nbsp_f=length(popf_f)
Nini=c(popf_p,popf_f)




jacob=rootSolve::jacobian.full(Nini,derivs,pert = 1e-6)
valp=max(abs(eigen(jacob,only.values=T)$values))
eig=eigen(jacob,only.values=T)$values
#inflow to abundance
inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(i in 1:nrow(jacob)){
for(j in 1:ncol(jacob)){
jacob2[i,j]=jacob[i,j]/(jacob[i,i]*jacob[j,j]-jacob[i,j]*jacob[j,i])
}}
diag(jacob2)=NA
diag(inv2)=NA
ia_direct_pp=mean(jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ia_direct_ff=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ia_direct_fp=mean(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ia_direct_pf=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ia_net_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ia_net_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ia_net_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ia_net_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
#abundance to abundance
inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(i in 1:nrow(jacob)){
for(j in 1:ncol(jacob)){
jacob2[i,j]=-jacob[i,j]/jacob[i,i]
inv2[i,j]=inv[i,j]/inv[j,j]
}}
diag(jacob2)=NA
diag(inv2)=NA
aa_direct_pp=mean(jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
aa_direct_ff=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
aa_direct_fp=mean(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
aa_direct_pf=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
aa_net_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
aa_net_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
aa_net_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
aa_net_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
#inflow to inflow
inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(i in 1:nrow(jacob)){
for(j in 1:ncol(jacob)){
jacob2[i,j]=-jacob[i,j]/jacob[j,j]
inv2[i,j]=inv[i,j]/inv[i,i]
}}
diag(jacob2)=NA
diag(inv2)=NA
ii_direct_pp=mean(jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ii_direct_ff=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ii_direct_fp=mean(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ii_direct_pf=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ii_net_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ii_net_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ii_net_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ii_net_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
#abundance to inflow
inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(i in 1:nrow(jacob)){
for(j in 1:ncol(jacob)){
inv2[i,j]=inv[i,j]/(inv[i,i]*inv[j,j]-inv[i,j]*inv[j,i])
}}
diag(jacob2)=NA
diag(inv2)=NA
ai_direct_pp=mean(jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ai_direct_ff=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_direct_fp=mean(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_direct_pf=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ai_net_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ai_net_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ai_net_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_net_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_indirect_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_indirect_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ai_indirect_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)]-jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ai_indirect_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_contrib_fp=mean((inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)])/(abs(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)])+abs(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)])),na.rm=T)
ai_contrib_pf=mean((inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)])/(abs(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)])+abs(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)])),na.rm=T)
ai_contrib_pp=mean((inv2[(1:nbsp_p),(1:nbsp_p)]-jacob2[(1:nbsp_p),(1:nbsp_p)])/(abs(inv2[(1:nbsp_p),(1:nbsp_p)]-jacob2[(1:nbsp_p),(1:nbsp_p)])+abs(jacob2[(1:nbsp_p),(1:nbsp_p)])),na.rm=T)
ai_contrib_ff=mean((inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)])/(abs(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)])+abs(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)])),na.rm=T)

inv2ff=inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]
inv2pp=inv2[(1:nbsp_p),(1:nbsp_p)]
indirect_intra=(length(inv2ff[inv2ff>0 & !is.na(inv2ff)])+length(inv2pp[inv2pp>0 & !is.na(inv2pp)]))/
(nbsp_f*nbsp_f+nbsp_p*nbsp_p-nbsp_f-nbsp_p)
inv2pf=inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)]
inv2fp=inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)]
indirect_inter=(length(inv2fp[inv2fp>0 & !is.na(inv2fp)])+length(inv2pf[inv2pf>0 & !is.na(inv2pf)]))/
(nbsp_f*nbsp_p*2)

netcareq1=data.frame(valprop=max(Re(eig)),stab=Re(eig)[which.max(abs(Re(eig)))],valp=valp,
connectance=networklevel(round(resf,digits=dig_num),index="weighted connectance")[[1]],NODF=networklevel(round(resf,digits=dig_num),index="weighted NODF",weighted=TRUE)[[1]],
connectance_I=networklevel(round(I,digits=dig_num),index="weighted connectance",
weighted=TRUE)[[1]],NODF_I=networklevel(round(I,digits=dig_num),index="weighted NODF",weighted=TRUE)[[1]],
NODF_5=networklevel(round(resf,digits=5),index="weighted NODF",weighted=TRUE)[[1]],
npoll=nbsp_p,nflow=nbsp_f,ntot=nbsp_f+nbsp_p,ovtot=mean(P),
ovpoll=mean(C_p_phe[upper.tri(C_p_phe)],na.rm=T),ovflow=mean(C_f_phe[upper.tri(C_f_phe)],na.rm=T),
denspoll=sum(popf_p),densflow=sum(popf_f),
divpoll=vegan::diversity(popf_p,index="shannon"),divflow=vegan::diversity(popf_f,index="shannon"),
vulnerability=networklevel(resf,index="vulnerability",weighted=TRUE)[[1]],generability=networklevel(resf,index="generality",weighted=TRUE)[[1]],
nicheovpollobs=networklevel(round(resf,digits=dig_num),index="niche overlap",weighted=TRUE)[[1]],
nicheovflowobs=networklevel(round(resf,digits=dig_num),index="niche overlap",weighted=TRUE)[[2]],
nicheovpollth=networklevel(I,index="niche overlap",weighted=TRUE)[[1]],
nicheovflowth=networklevel(I,index="niche overlap",weighted=TRUE)[[2]],
interfp=interfp,interff=interff,
indirect_intra=indirect_intra,
indirect_inter=indirect_inter,
modularity=computeModules(resf, method="Beckett")@likelihood,random=jj,mod_I=computeModules(I, method="Beckett")@likelihood,
alpha=alpha,beta=beta,stade="initial",A_mean=mean(A),A_sd=sd(A),P_mean=mean(P),P_sd=sd(P),I_mean=mean(I),I_sd=sd(I),
ia_direct_pp=ia_direct_pp,ia_direct_ff=ia_direct_ff,ia_direct_fp=ia_direct_fp,
ia_direct_pf=ia_direct_pf,ia_net_pf=ia_net_pf,ia_net_pp=ia_net_pp,ia_net_ff=ia_net_ff,ia_net_fp=ia_net_fp,
aa_direct_pp=aa_direct_pp,aa_direct_ff=aa_direct_ff,aa_direct_fp=aa_direct_fp,
aa_direct_pf=aa_direct_pf,aa_net_pf=aa_net_pf,aa_net_pp=aa_net_pp,aa_net_ff=aa_net_ff,aa_net_fp=aa_net_fp,
ii_direct_pp=ii_direct_pp,ii_direct_ff=ii_direct_ff,ii_direct_fp=ii_direct_fp,
ii_direct_pf=ii_direct_pf,ii_net_pf=ii_net_pf,ii_net_pp=ii_net_pp,ii_net_ff=ii_net_ff,ii_net_fp=ii_net_fp,
ai_direct_pp=ai_direct_pp,ai_direct_ff=ai_direct_ff,ai_direct_fp=ai_direct_fp,
ai_direct_pf=ai_direct_pf,ai_net_pf=ai_net_pf,ai_net_pp=ai_net_pp,ai_net_ff=ai_net_ff,ai_net_fp=ai_net_fp,
ai_indirect_pp=ai_indirect_pp,ai_indirect_ff=ai_indirect_ff,ai_indirect_fp=ai_indirect_fp,ai_indirect_pf=ai_indirect_pf,
ai_contrib_pf=ai_contrib_pf,ai_contrib_pp=ai_contrib_pp,ai_contrib_ff=ai_contrib_ff,ai_contrib_fp=ai_contrib_fp)
if(advi==deb){netcar=netcareq1}else{netcar=rbind(netcar,netcareq1)}



while(a>conv){
if(a==1){out <- ode(y=Nini, times=seq(0,20,0.05), func = derivs,parms = NULL,method="lsoda")}else{
out <- ode(y=Nini, times=seq(0,20,1), func = derivs,parms = NULL,method="lsoda")}
colnames(out)[1]="Time"
# if(a==1){out=integrate_sys(derivs2,Nini,20,0.05)}else{out=integrate_sys(derivs2,Nini,20,1)}
Nini=t(out)[2:(1+nbsp_p+nbsp_f),nrow(out)]
Nini[which(Nini<seuil)]=0
Nini[which(is.na(Nini))]=0
if(b==0){dyna=out}else{
out[,"Time"]=out[,"Time"]+b
dyna=rbind(dyna,out[-1,])}
if(nrow(dyna)>10){a=max(apply(t(dyna[(nrow(dyna)-10):nrow(dyna),2:(1+nbsp_p+nbsp_f)]),1,var))}
b=b+20
}



popf_f=Nini[(nbsp_p+1):(nbsp_p+nbsp_f)]
popf_p=Nini[1:nbsp_p]
index=c(which(Nini>=seuil))

final$N_eq1=c(popf_p,popf_f)
final$interf=c(rep(interfp,nbsp_p),rep(interff,nbsp_f))
final$random=jj
final$beta=beta
final$alpha=alpha
final$I_mean=c(apply(I,2,mean),apply(I,1,mean))
final$A_mean=c(apply(A,2,mean),apply(A,1,mean))
final$P_mean=c(apply(P,2,mean),apply(P,1,mean))
final$C_mean=c(apply(C_p_phe,1,mean),apply(C_f_phe,1,mean))
final$recu_directfp=NA
final$emi_directfp=NA
final$recu_totfp=NA
final$emi_totfp=NA
final$recu_directpp=NA
final$emi_directpp=NA
final$recu_totpp=NA
final$emi_totpp=NA
#final$closeness_j=NA
#final$betweenness_j=NA
#final$closeness_t=NA
#final$betweenness_t=NA
final$I_meaneq1=NA
final$A_meaneq1=NA
final$P_meaneq1=NA



if(sum(popf_f)>0 & sum(popf_p)>0){
####NETWORK LEVEL INFORMATIONS TO EXPORT
popf_f=c(t(dyna[nrow(dyna),(nbsp_p+2):(1+nbsp_p+nbsp_f)]))
popf_p=c(t(dyna[nrow(dyna),2:(nbsp_p+1)]))

### AGREGER:
I=I[which(popf_f>=seuil),which(popf_p>=seuil)]
C_p_phe=C_p_phe[which(popf_p>=seuil),which(popf_p>=seuil)]
C_f_phe=C_f_phe[which(popf_f>=seuil),which(popf_f>=seuil)]
A=A[which(popf_f>=seuil),which(popf_p>=seuil)]
P=P[which(popf_f>=seuil),which(popf_p>=seuil)]
m=m[which(popf_p>=seuil)]
m2=m2[which(popf_f>=seuil)]
K=K[which(popf_p>=seuil)]
K2=K2[which(popf_f>=seuil)]
hangling=hangling[Nini>seuil]
efficience=efficience[Nini>seuil]
popf_f=popf_f[which(popf_f>=seuil)]
popf_p=popf_p[which(popf_p>=seuil)]
mati=sqrt(popf_f%*%t(popf_p))
resf=mati*I
nbsp_p=length(popf_p)
nbsp_f=length(popf_f)
Nini=c(popf_p,popf_f)

if(class(I)!="matrix"){
I=matrix(I,nbsp_f,nbsp_p)
C_p_phe=matrix(C_p_phe,nbsp_p,nbsp_p)
C_f_phe=matrix(C_f_phe,nbsp_f,nbsp_f)
A=matrix(A,nbsp_f,nbsp_p)
P=matrix(P,nbsp_f,nbsp_p)
}


jacob=rootSolve::jacobian.full(Nini,derivs,pert = 1e-6)
valp=max(abs(eigen(jacob,only.values=T)$values))
eig=eigen(jacob,only.values=T)$values
#inflow to abundance
inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(i in 1:nrow(jacob)){
for(j in 1:ncol(jacob)){
jacob2[i,j]=jacob[i,j]/(jacob[i,i]*jacob[j,j]-jacob[i,j]*jacob[j,i])
}}
diag(jacob2)=NA
diag(inv2)=NA
ia_direct_pp=mean(jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ia_direct_ff=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ia_direct_fp=mean(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ia_direct_pf=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ia_net_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ia_net_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ia_net_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ia_net_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
#abundance to abundance
inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(i in 1:nrow(jacob)){
for(j in 1:ncol(jacob)){
jacob2[i,j]=-jacob[i,j]/jacob[i,i]
inv2[i,j]=inv[i,j]/inv[j,j]
}}
diag(jacob2)=NA
diag(inv2)=NA
aa_direct_pp=mean(jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
aa_direct_ff=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
aa_direct_fp=mean(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
aa_direct_pf=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
aa_net_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
aa_net_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
aa_net_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
aa_net_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
#inflow to inflow
inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(i in 1:nrow(jacob)){
for(j in 1:ncol(jacob)){
jacob2[i,j]=-jacob[i,j]/jacob[j,j]
inv2[i,j]=inv[i,j]/inv[i,i]
}}
diag(jacob2)=NA
diag(inv2)=NA
ii_direct_pp=mean(jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ii_direct_ff=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ii_direct_fp=mean(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ii_direct_pf=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ii_net_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ii_net_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ii_net_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ii_net_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
#abundance to inflow
inv=-1*solve(jacob)
jacob2=jacob
inv2=inv
for(i in 1:nrow(jacob)){
for(j in 1:ncol(jacob)){
inv2[i,j]=inv[i,j]/(inv[i,i]*inv[j,j]-inv[i,j]*inv[j,i])
}}
diag(jacob2)=NA
diag(inv2)=NA
ai_direct_pp=mean(jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ai_direct_ff=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_direct_fp=mean(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_direct_pf=mean(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ai_net_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ai_net_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ai_net_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_net_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_indirect_fp=mean(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_indirect_pf=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],na.rm=T)
ai_indirect_pp=mean(inv2[(1:nbsp_p),(1:nbsp_p)]-jacob2[(1:nbsp_p),(1:nbsp_p)],na.rm=T)
ai_indirect_ff=mean(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],na.rm=T)
ai_contrib_fp=mean((inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)])/(abs(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)])+abs(jacob2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)])),na.rm=T)
ai_contrib_pf=mean((inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)])/(abs(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)])+abs(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)])),na.rm=T)
ai_contrib_pp=mean((inv2[(1:nbsp_p),(1:nbsp_p)]-jacob2[(1:nbsp_p),(1:nbsp_p)])/(abs(inv2[(1:nbsp_p),(1:nbsp_p)]-jacob2[(1:nbsp_p),(1:nbsp_p)])+abs(jacob2[(1:nbsp_p),(1:nbsp_p)])),na.rm=T)
ai_contrib_ff=mean((inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)])/(abs(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]-jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)])+abs(jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)])),na.rm=T)

inv2ff=inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]
inv2pp=inv2[(1:nbsp_p),(1:nbsp_p)]
indirect_intra=(length(inv2ff[inv2ff>0 & !is.na(inv2ff)])+length(inv2pp[inv2pp>0 & !is.na(inv2pp)]))/
(nbsp_f*nbsp_f+nbsp_p*nbsp_p-nbsp_f-nbsp_p)
inv2pf=inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)]
inv2fp=inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)]
indirect_inter=(length(inv2fp[inv2fp>0 & !is.na(inv2fp)])+length(inv2pf[inv2pf>0 & !is.na(inv2pf)]))/
(nbsp_f*nbsp_p*2)

if(nbsp_p>1 & nbsp_f>1){
diag(jacob)=NA
jacob_norm=jacob/matrix(Nini,ncol=length(Nini),nrow=length(Nini),byrow=F)
inv2=inv2/matrix(Nini,ncol=length(Nini),nrow=length(Nini),byrow=F)
#effet emis et recu entre poll et plantes, directs et totaux par sp
final$recu_directfp[index]=c(apply(jacob[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],1,sum),
apply(jacob[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],1,sum))
final$recu_totfp[index]=c(apply(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],1,sum),
apply(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],1,sum))
final$emi_totfp[index]=c(apply(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],2,sum),
apply(inv2[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],2,sum))
final$emi_directfp[index]=c(apply(jacob[(nbsp_p+1):(nbsp_p+nbsp_f),(1:nbsp_p)],2,sum),
apply(jacob[(1:nbsp_p),(nbsp_p+1):(nbsp_p+nbsp_f)],2,sum))

#effet emis et recu intraguilds, directs et totaux par sp
final$recu_directpp[index]=c(apply(jacob[(1:nbsp_p),(1:nbsp_p)],1,sum,na.rm=T),
apply(jacob[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],1,sum,na.rm=T))
final$emi_directpp[index]=c(apply(jacob[(1:nbsp_p),(1:nbsp_p)],2,sum,na.rm=T),
apply(jacob[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],2,sum,na.rm=T))
final$recu_totpp[index]=c(apply(inv2[(1:nbsp_p),(1:nbsp_p)],1,sum,na.rm=T),
apply(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],1,sum,na.rm=T))
final$emi_totpp[index]=c(apply(inv2[(1:nbsp_p),(1:nbsp_p)],2,sum,na.rm=T),
apply(inv2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)],2,sum,na.rm=T))


jacob2=abs(jacob)
jacob2[(1:nbsp_p),(1:nbsp_p)]=0
jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]=0
g=graph_from_adjacency_matrix(jacob2,mode="directed",weighted=T,diag =F)
#final$betweenness_j[index]=igraph::betweenness(g,v=V(g),directed=T)
#final$closeness_j[index]=igraph::closeness(g,v=V(g),mode="out")

jacob2=abs(inv2)
jacob2[(1:nbsp_p),(1:nbsp_p)]=0
jacob2[(nbsp_p+1):(nbsp_p+nbsp_f),(nbsp_p+1):(nbsp_p+nbsp_f)]=0
g=graph_from_adjacency_matrix(jacob2,mode="directed",weighted=T,diag =F)
#final$betweenness_t[index]=igraph::betweenness(g,v=V(g),directed=T)
#final$closeness_t[index]=igraph::closeness(g,v=V(g),mode="out")
final$I_meaneq1[index]=c(apply(I,2,mean),apply(I,1,mean))

netcareq1=data.frame(valprop=max(Re(eig)),stab=Re(eig)[which.max(abs(Re(eig)))],valp=valp,
connectance=networklevel(round(resf,digits=dig_num),index="weighted connectance")[[1]],NODF=networklevel(round(resf,digits=dig_num),index="weighted NODF",weighted=TRUE)[[1]],
connectance_I=networklevel(round(I,digits=dig_num),index="weighted connectance",
weighted=TRUE)[[1]],NODF_I=networklevel(round(I,digits=dig_num),index="weighted NODF",weighted=TRUE)[[1]],
NODF_5=networklevel(round(resf,digits=5),index="weighted NODF",weighted=TRUE)[[1]],
npoll=nbsp_p,nflow=nbsp_f,ntot=nbsp_f+nbsp_p,ovtot=mean(P),
ovpoll=mean(C_p_phe[upper.tri(C_p_phe)],na.rm=T),ovflow=mean(C_f_phe[upper.tri(C_f_phe)],na.rm=T),
denspoll=sum(popf_p),densflow=sum(popf_f),
divpoll=vegan::diversity(popf_p,index="shannon"),divflow=vegan::diversity(popf_f,index="shannon"),
vulnerability=networklevel(resf,index="vulnerability",weighted=TRUE)[[1]],generability=networklevel(resf,index="generality",weighted=TRUE)[[1]],
nicheovpollobs=networklevel(round(resf,digits=dig_num),index="niche overlap",weighted=TRUE)[[1]],
nicheovflowobs=networklevel(round(resf,digits=dig_num),index="niche overlap",weighted=TRUE)[[2]],
nicheovpollth=networklevel(I,index="niche overlap",weighted=TRUE)[[1]],
nicheovflowth=networklevel(I,index="niche overlap",weighted=TRUE)[[2]],
interfp=interfp,interff=interff,
indirect_intra=indirect_intra,
indirect_inter=indirect_inter,
modularity=computeModules(resf, method="Beckett")@likelihood,random=jj,mod_I=computeModules(I, method="Beckett")@likelihood,
alpha=alpha,beta=beta,stade="eq1",A_mean=mean(A),A_sd=sd(A),P_mean=mean(P),P_sd=sd(P),I_mean=mean(I),I_sd=sd(I),
ia_direct_pp=ia_direct_pp,ia_direct_ff=ia_direct_ff,ia_direct_fp=ia_direct_fp,
ia_direct_pf=ia_direct_pf,ia_net_pf=ia_net_pf,ia_net_pp=ia_net_pp,ia_net_ff=ia_net_ff,ia_net_fp=ia_net_fp,
aa_direct_pp=aa_direct_pp,aa_direct_ff=aa_direct_ff,aa_direct_fp=aa_direct_fp,
aa_direct_pf=aa_direct_pf,aa_net_pf=aa_net_pf,aa_net_pp=aa_net_pp,aa_net_ff=aa_net_ff,aa_net_fp=aa_net_fp,
ii_direct_pp=ii_direct_pp,ii_direct_ff=ii_direct_ff,ii_direct_fp=ii_direct_fp,
ii_direct_pf=ii_direct_pf,ii_net_pf=ii_net_pf,ii_net_pp=ii_net_pp,ii_net_ff=ii_net_ff,ii_net_fp=ii_net_fp,
ai_direct_pp=ai_direct_pp,ai_direct_ff=ai_direct_ff,ai_direct_fp=ai_direct_fp,
ai_direct_pf=ai_direct_pf,ai_net_pf=ai_net_pf,ai_net_pp=ai_net_pp,ai_net_ff=ai_net_ff,ai_net_fp=ai_net_fp,
ai_indirect_pp=ai_indirect_pp,ai_indirect_ff=ai_indirect_ff,ai_indirect_fp=ai_indirect_fp,
ai_indirect_pf=ai_indirect_pf,ai_contrib_pf=ai_contrib_pf,ai_contrib_pp=ai_contrib_pp,ai_contrib_ff=ai_contrib_ff,ai_contrib_fp=ai_contrib_fp)
}else{
netcareq1=data.frame(valprop=max(Re(eig)),stab=Re(eig)[which.max(abs(Re(eig)))],valp=valp,
connectance=NA,NODF=NA,connectance_I=NA,NODF_I=NA,NODF_5=NA,
npoll=nbsp_p,nflow=nbsp_f,ntot=nbsp_f+nbsp_p,ovtot=mean(P),
ovpoll=NA,ovflow=NA,
denspoll=sum(popf_p),densflow=sum(popf_f),
divpoll=NA,divflow=NA,
vulnerability=NA,generability=NA,
nicheovpollobs=NA,
nicheovflowobs=NA,
nicheovpollth=NA,
nicheovflowth=NA,
interfp=interfp,interff=interff,
indirect_intra=indirect_intra,
indirect_inter=indirect_inter,
modularity=NA,random=jj,mod_I=NA,
alpha=alpha,beta=beta,stade="eq1",A_mean=mean(A),A_sd=sd(A),P_mean=mean(P),P_sd=sd(P),I_mean=mean(I),I_sd=sd(I),
ia_direct_pp=ia_direct_pp,ia_direct_ff=ia_direct_ff,ia_direct_fp=ia_direct_fp,
ia_direct_pf=ia_direct_pf,ia_net_pf=ia_net_pf,ia_net_pp=ia_net_pp,ia_net_ff=ia_net_ff,ia_net_fp=ia_net_fp,
aa_direct_pp=aa_direct_pp,aa_direct_ff=aa_direct_ff,aa_direct_fp=aa_direct_fp,
aa_direct_pf=aa_direct_pf,aa_net_pf=aa_net_pf,aa_net_pp=aa_net_pp,aa_net_ff=aa_net_ff,aa_net_fp=aa_net_fp,
ii_direct_pp=ii_direct_pp,ii_direct_ff=ii_direct_ff,ii_direct_fp=ii_direct_fp,
ii_direct_pf=ii_direct_pf,ii_net_pf=ii_net_pf,ii_net_pp=ii_net_pp,ii_net_ff=ii_net_ff,ii_net_fp=ii_net_fp,
ai_direct_pp=ai_direct_pp,ai_direct_ff=ai_direct_ff,ai_direct_fp=ai_direct_fp,
ai_direct_pf=ai_direct_pf,ai_net_pf=ai_net_pf,ai_net_pp=ai_net_pp,ai_net_ff=ai_net_ff,ai_net_fp=ai_net_fp,
ai_indirect_pp=ai_indirect_pp,ai_indirect_ff=ai_indirect_ff,ai_indirect_fp=ai_indirect_fp,
ai_indirect_pf=ai_indirect_pf,
ai_contrib_pf=ai_contrib_pf,ai_contrib_pp=ai_contrib_pp,ai_contrib_ff=ai_contrib_ff,ai_contrib_fp=ai_contrib_fp)
}
}else{
netcareq1=data.frame(valprop=NA,stab=NA,valp=NA,
connectance=NA,NODF=NA,connectance_I=NA,NODF_I=NA,NODF_5=NA,
npoll=0,nflow=0,ntot=0,ovtot=NA,
ovpoll=NA,ovflow=NA,
denspoll=0,densflow=0,
divpoll=NA,divflow=NA,
vulnerability=NA,generability=NA,
nicheovpollobs=NA,
nicheovflowobs=NA,
nicheovpollth=NA,
nicheovflowth=NA,
interfp=interfp,interff=interff,
indirect_intra=NA,
indirect_inter=NA,
modularity=NA,random=jj,mod_I=NA,
alpha=alpha,beta=beta,stade="eq1",A_mean=NA,A_sd=NA,P_mean=NA,P_sd=NA,I_mean=NA,I_sd=NA,
ia_direct_pp=NA,ia_direct_ff=NA,ia_direct_fp=NA,
ia_direct_pf=NA,ia_net_pf=NA,ia_net_pp=NA,ia_net_ff=NA,ia_net_fp=NA,
aa_direct_pp=NA,aa_direct_ff=NA,aa_direct_fp=NA,
aa_direct_pf=NA,aa_net_pf=NA,aa_net_pp=NA,aa_net_ff=NA,aa_net_fp=NA,
ii_direct_pp=NA,ii_direct_ff=NA,ii_direct_fp=NA,
ii_direct_pf=NA,ii_net_pf=NA,ii_net_pp=NA,ii_net_ff=NA,ii_net_fp=NA,
ai_direct_pp=NA,ai_direct_ff=NA,ai_direct_fp=NA,
ai_direct_pf=NA,ai_net_pf=NA,ai_net_pp=NA,ai_net_ff=NA,ai_net_fp=NA,
ai_indirect_pp=NA,ai_indirect_ff=NA,ai_indirect_fp=NA,
ai_indirect_pf=NA,ai_contrib_pf=NA,ai_contrib_pp=NA,ai_contrib_ff=NA,ai_contrib_fp=NA)
}
netcar=rbind(netcar,netcareq1)
if(advi==deb){resspecies=final}else{resspecies=rbind(resspecies,final)}
}
return(list(netcar,resspecies))
}

setwd(dir="/home/duchenne/part_1/75x75")
fwrite(resultat[[1]],paste0("netcar_75x75_interfp_",interfp,"_interff_",interff,"_",fini,".txt"),row.names=F,sep="\t")
fwrite(resultat[[2]],paste0("resspecies_75x75_interfp_",interfp,"_interff_",interff,"_",fini,".txt"),row.names=F,sep="\t")

closeCluster(cl)
mpi.quit()

