library(truncnorm)
library(bipartite)
library(plot3D)
library(Rmpfr)
library(circular)
library(CircStats)
library(phylolm)
library(ggplot2)
library(phytools)
library(ape)
library(picante)
library(doBy)
library(phytools)
library(phylotools)
library(phangorn)
library(gridExtra)
library(doParallel)
library(foreach)
# cl<-makeCluster(6) 
# registerDoParallel(cl)


foreach(jj=1:1000,.combine=rbind)%dopar%{
library(truncnorm)
library(bipartite)
library(plot3D)
library(Rmpfr)
library(circular)
library(CircStats)
library(phylolm)
library(ggplot2)
library(phytools)
library(ape)
library(picante)
library(doBy)
library(phytools)
library(phylotools)
library(phangorn)
library(gridExtra)
seuil=1e-5 #extinction cutoff
precision=0.1 #precision for integrals
digi=5 #round interaction matrix
conv=1e-11 #convergence criterion
dig_num=16 #round values for calcul
nbsp_p=50
nbsp_f=50
interfp=0
interff=0

vec_p=c()
vec_f=c()
indi=0
duree=70
indi=indi+1
mu_p=sort(rnorm(nbsp_p,190,duree))
sd_p=runif(nbsp_p,5,40)
mu_f=sort(rnorm(nbsp_f,190,duree))
sd_f=runif(nbsp_f,5,40)
N_p=runif(nbsp_p,0.1,4)
N_f=runif(nbsp_f,5,50)
genmu_p=runif(nbsp_p,-1.5,1.5)
genmu_f=runif(nbsp_f,-1.5,1.5)
gensd_f=runif(nbsp_f,0.1,0.9)
gensd_p=runif(nbsp_p,0.1,0.9)
names(genmu_p)=paste0("p",1:nbsp_p)
names(genmu_f)=paste("f",1:nbsp_f)


#vecteurs parametres
rmax <- runif(nbsp_p,1,1); rmax2<-runif(nbsp_f,1,1)
K <- runif(nbsp_p,1,60);K2<-runif(nbsp_f,10,600) 
Nini <- c(N_p,N_f); 
m <- runif(nbsp_p,0.8,1);m2<-runif(nbsp_f,0.2,0.4)
cs <- rep(0.3,nbsp_p);cs2<-rep(0.7,nbsp_f)
hangling=runif(nbsp_f+nbsp_p,0.9,0.9)
efficience=runif(nbsp_f+nbsp_p,0.8,1)

#####SPECIES LEVEL INFORMATIONS TO EXPORT

final=data.frame(sp=c(names(genmu_p),names(genmu_f)),type=c(rep("poll",nbsp_p),rep("flow",nbsp_f)),mfd=c(mu_p,mu_f),
sd=c(sd_p,sd_f),trait_sd=c(gensd_p,gensd_f),trait_mu=c(genmu_p,genmu_f),N=c(N_p,N_f),K=c(K,K2),
m=c(m,m2),cs=c(cs,cs2),
interf=c(rep(interfp,nbsp_p),rep(interff,nbsp_f)),hangling=hangling,efficience=efficience)
final$stade="initial"
final$random=jj


setwd(dir=paste0("C:/Users/Francois/Documents/modélisation part1/initial_networks_",nbsp_p,"x",nbsp_f))
write.table(final,paste("pops_eq_i_",jj,".txt",sep=""),sep="\t",row.names=F)
}


stopCluster(cl)


