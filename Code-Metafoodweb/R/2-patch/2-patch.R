require(gtools)
library(ggplot2)
library(reshape2)
library(plyr)
library(sfsmisc)
library(gridExtra)
library(grid)
library(EMD)
require(doParallel)
######

rm(list=ls())
cl<-makeCluster(5)    ##### WARNING! Number of cores may depend on your machine!
registerDoParallel(cl)


replicates<-seq(1,5,1)
nutrient<-seq(0,20,1)
loss<-seq(0.001,1,0.05)

data.end<-foreach (i=replicates,.combine = "rbind") %dopar% {

  bdd <- data.frame(
    patch = numeric(),
    rep = numeric(),
    nut = numeric(),
    loss = numeric(),
    div = numeric()
  )

  for (j in nutrient)
  {
    for (k in loss)
    {
#######
      
dat<-read.table(paste("../../Simulations/Heatmap/timeseries_REP_",i,"_NUT_",j,"_LOSS_",k,".out",sep=""),header=T,sep="")

### fluxes
subset1 <- subset(dat,dat$time>20000)



div=sum(subset1[,2:7]>0)
div1=sum(subset1[,2:4]>0)
div2=sum(subset1[,5:7]>0)
  
bddx.p1 <- data.frame(
    patch = 1,
    rep = i,
    nut = 20-j,
    loss = k,
    div = div1
  )

bddx.p2 <- data.frame(
  patch = 2,
  rep = i,
  nut = 20-j,
  loss = k,
  div = div2
 )

bddx <- data.frame(
  patch = 3,
  rep = i,
  nut = 20-j,
  loss = k,
  div = div
)

bdd <- rbind(bdd,bddx)
bdd<- rbind(bdd,bddx.p1)
bdd <- rbind(bdd,bddx.p2)

    }
  }
return(bdd)  
}

write.table(data.end,"Diversity.txt")



######









replicates<-seq(1,5,1)
nutrient<-seq(0,20,1)
loss<-seq(0.001,1.000,0.05)

data.end2<-foreach (i=replicates,.combine = "rbind",.packages = c("EMD")) %dopar% {

  bddamp <- data.frame(
    patch = numeric(),
    rep = numeric(),
    nut = numeric(),
    loss = numeric(),
    ampl = numeric()
  )
  
  for (j in nutrient)
  {
    for (k in loss)
    {
      #######
      dat<-read.table(paste("../../Simulations/Heatmap/timeseries_REP_",i,"_NUT_",j,"_LOSS_",k,".out",sep=""),header=T,sep="")
      
      ### fluxes
     
      ifelse(
        extrema(dat$B2P0)$nextreme<2,
        b1.min <- rep(min(dat$B2P0),25),
        b1.min <- dat$B2P0[extrema(dat$B2P0)$minindex[1:25]]
      )
      ## b3max
      ifelse(
        extrema(dat$B2P0)$nextreme<2,
        b1.max <- rep(min(dat$B2P0),25),
        b1.max <- dat$B2P0[extrema(dat$B2P0)$maxindex[1:25]]
      ) 
      ifelse(
      b1.max[1]>0,
      ampl1<-mean(I(b1.max-b1.min), na.rm=T),
      ampl1<-"NA"
      )
      
      ifelse(
        extrema(dat$B2P1)$nextreme<2,
        b2.min <- rep(min(dat$B2P1),25),
        b2.min <- dat$B2P1[extrema(dat$B2P1)$minindex[1:25]]
      )
      ## b3max
      ifelse(
        extrema(dat$B2P1)$nextreme<2,
        b2.max <- rep(min(dat$B2P1),25),
        b2.max <- dat$B2P1[extrema(dat$B2P1)$maxindex[1:25]]
      )
      ifelse(
      b2.max[1]>0,
      ampl2<-mean(I(b2.max-b2.min), na.rm=T),
      ampl2<-"NA"
      )
      
      
      bddx.p1 <- data.frame(
        patch=1,
        rep = i,
        nut = 20-j,
        loss = k,
        ampl = ampl1
      )
      bddx.p2 <- data.frame(
        patch=2,
        rep = i,
        nut = 20-j,
        loss = k,
        ampl = ampl2
      )
      
  
      bddamp <- rbind(bddamp,bddx.p1)
      bddamp <- rbind(bddamp,bddx.p2)
      
    }
    
  }
  return(bddamp)
}

write.table(data.end2,"Ampl.txt")



bdd<-read.table("Diversity.txt")
bdd.ls<-subset(bdd, bdd$patch==3)
bdd.p1<-subset(bdd, bdd$patch==1)
bdd.p2<-subset(bdd, bdd$patch==2)


g<-ggplot(bdd.ls , aes(x=nut,y=loss)) +
  scale_fill_viridis_c(name="# Populations", breaks=c(0,1,2,3), limits=c(0,3), option = "C") +
  geom_raster(aes(fill = div/2), interpolate=F) +
  labs(x="Heterogeneity",y="Dispersal Loss") +
  theme_classic() + 
  xlim(0,20) +
  theme(axis.text=element_text(size=15),text=element_text(size=15, face="bold"),
        axis.title=element_text(size=15,face="bold"))

g1<-ggplot(bdd.p1 , aes(x=nut,y=loss)) +
  scale_fill_gradientn(colours = mycol, name=NULL, breaks=c(0,1,2,3), limits=c(0,3)) +
  geom_raster(aes(fill = div), interpolate=T) +
  labs(x="Heterogeneity",y="Hostility", title="Variable Patch") +
  theme_classic() + 
  xlim(0,20)+
  theme(axis.text=element_text(size=22),text=element_text(size=22, face="bold"),
        axis.title=element_text(size=22,face="bold"))

g2<-ggplot(bdd.p2 , aes(x=nut,y=loss)) +
  scale_fill_gradientn(colours = mycol, name=NULL, breaks=c(0,1,2,3), limits=c(0,3)) +
  geom_raster(aes(fill = div), interpolate=T) +
  labs(x="Heterogeneity",y="Hostility", title="Eutrophic Patch") +
  theme_classic() + 
  xlim(0,20)+
  theme(axis.text=element_text(size=22),text=element_text(size=22, face="bold"),
                       axis.title=element_text(size=22,face="bold"))


pdf("Heatmap.pdf",width=8,height=6)
grid.arrange(g)
dev.off()

data<-read.table("Ampl.txt")
p1<-subset(data,data$patch==1)
p2<-subset(data,data$patch==2)

mycol <- c("blue","green", "yellow", "red")

g1<-ggplot(p1 , aes(x=nut,y=loss)) +
  scale_fill_viridis_c(name="Amplitude", option = "C") +
  geom_raster(aes(fill = ampl), interpolate=T) +
  labs(x="Heterogeneity",y="Hostility", title="Variable Patch") +
  theme_classic() + 
  xlim(0,20)+
  theme(axis.text=element_text(size=22),text=element_text(size=22, face="bold"),
        axis.title=element_text(size=22,face="bold"))

g2<-ggplot(p2 , aes(x=nut,y=loss)) +
  scale_fill_viridis_c(name="Amplitude", option = "C") +
  geom_raster(aes(fill = ampl), interpolate=T) +
  labs(x="Heterogeneity",y="Hostility", title="Eutrophic Patch") +
  theme_classic() + 
  xlim(0,20)+
  theme(axis.text=element_text(size=15),text=element_text(size=15, face="bold"),
        axis.title=element_text(size=15,face="bold"))


pdf("Amplitude.pdf",width=8,height=6)
grid.arrange(g2)
dev.off()
