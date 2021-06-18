## Analysis of Full Foodwebs in different nutrient landscapes
library(plyr)
library(ggplot2)
library(igraph)
library(raster)
library(matrixStats)
library("viridis")

require(doParallel)
######


cl<-makeCluster(16)    ##### WARNING! Number of cores may depend on your machine!
registerDoParallel(cl)



data.end<-foreach (i=c(1:15,30:44),.combine = "rbind") %dopar% {
  
  persistance<-data.frame(Scenario=numeric(0),Nutrient.factor=numeric(0),Nutrient=numeric(0),Web=numeric(0),Landscape=numeric(0),Patch=numeric(0), Species=numeric(0), Persist=numeric(0), Growth1=numeric(0),Growth2=numeric(0), Plant=numeric(0))
  
  local<-read.table("../../Simulations/FullWeb/SummaryMass.out", sep=",", header=T)

  subset.nutrient<-subset(local,local$nut==i)
Nutrients<- read.table(paste("../../Simulations/FullWeb/Nutrients/Nutrient_", i, ".out",sep=""),header=F)

  for (w in 1)
  {
    subset.web<-subset(subset.nutrient,subset.nutrient$web==w)
    
    for (l in 1:5)
    {
      subset.landscape<-subset(subset.web,subset.web$landscape==l)
    
  for (p in 0:49)
  {
      subset.patch<- subset(subset.landscape, subset.landscape$patch==p)
      
  for (k in 0:39)
  {
    subset.sp<-subset(subset.patch,subset.patch$species==k)
    
    if (subset.sp$biomass.tend.10k==0 && subset.sp$biomass.tend.20k==0 && subset.sp$biomass.tend!=0) 
    {
      persist<-NA
      mean.net.growth1<-NA
      mean.net.growth2<-NA
      
    } else {
      persist<-(sum(subset.sp$biomass.tend)>0)*1
      mean.net.growth1<-subset.sp$Mean.net.growth1
      mean.net.growth2<-subset.sp$Mean.net.growth2
    }
    if.basal<-subset.sp$if.basal.spp
    if (i<20)
    {scen<-"Homogenious"} else if (i<35 && i>16) {scen<-"HeterogeneousMeso"} else if (i<40 && i>34) {scen<-"HeterogeneousOligo"} else {scen<-"HeterogeniousEutro"}
    
    vec<-data.frame(Scenario=scen,Nutrient.factor=i,Nutrient=Nutrients[p+1,], Web=w,Landscape=l,Patch=p,Species=k,Persist=persist, Growth1=mean.net.growth1, Growth2=mean.net.growth2, Plant=if.basal)
    
    persistance<-rbind(persistance,vec)
  
    
    }
   }
  }
  }
return(persistance)  
}
str(data.end)

persistance1<-data.end
write.table(persistance1,"data.out")
persistance1<-read.table("data.out")
persistence2<-read.table("../FullWeb/data.out")

for.plot<-ddply(persistance1,.(Scenario,Nutrient.factor,Nutrient,Web,Patch,Landscape, Plant), summarize, Persistance=sum(Persist))


plants<-subset(for.plot,for.plot$Plant==1)
animals<-subset(for.plot,for.plot$Plant==0)

pdf("FigPlantAnimal.pdf",width=12,height=8, onefile = T)
ggplot(plants, aes(x=log10(Nutrient),y=Persistance,color=Scenario))  + geom_point(size=0.2) + geom_smooth(size=0.5) + labs(x="Nutrient Supply log10[]",y=expression(bold(paste("Plant ",alpha,"-Diversity")))) +
  theme_classic() + theme(axis.text=element_text(size=22),text=element_text(size=22, face="bold"),
                          axis.title=element_text(size=22,face="bold")) 
ggplot(animals, aes(x=log10(Nutrient),y=Persistance,color=Scenario))   + geom_point(size=0.2) + geom_smooth(size=0.5) + labs(x="Nutrient Supply log10[]",y=expression(bold(paste("Animal ",alpha,"-Diversity")))) +
  theme_classic() + theme(axis.text=element_text(size=22),text=element_text(size=22, face="bold"),
                          axis.title=element_text(size=22,face="bold")) 

dev.off()
for.plot1<-ddply(persistance1,.(Scenario,Nutrient.factor,Nutrient,Web,Patch,Landscape), summarize, Persistance=sum(Persist))

pdf("Fig4.pdf",width=12,height=8)
ggplot(for.plot1, aes(x=log10(Nutrient),y=Persistance,color=Scenario))  + geom_count(position = position_jitter(w = 0.05, h = 0),shape=1, alpha=1, stroke = 0.8) + geom_smooth(size=1) + labs(x="Local Nutrient Supply log10[]",y=expression(bold(paste(alpha,"-Diversity")))) +
  theme_classic() + theme(axis.text=element_text(size=22),text=element_text(size=22, face="bold"),
                         axis.title=element_text(size=22,face="bold")) 

dev.off()

for.plot1.sub<-subset(for.plot1,for.plot1$Scenario!="Homogenious")
pdf("Fig4add.pdf", width=12, height=4, onefile = T)
ggplot(for.plot1, aes(y=log10(Nutrient),x=Scenario,color=Scenario)) + geom_violin( draw_quantiles = c(0.5), na.rm=T) + stat_summary(fun.y=mean, geom="point", shape=23, size=2) +coord_flip() + theme_classic() + theme(axis.text=element_text(size=22),text=element_text(size=22, face="bold"),
                                                                                                                                                             axis.title=element_text(size=22,face="bold")) 
ggplot(for.plot1, aes(y=log10(Nutrient),x=Scenario,color=Scenario)) + geom_violin( draw_quantiles = c(0.5), na.rm=T)  +coord_flip() + theme_classic()+ theme(axis.text=element_text(size=22),text=element_text(size=22, face="bold"),
                                                                                                                                                                 axis.title=element_text(size=22,face="bold")) 

dev.off()

for.plot2<-ddply(persistance1,.(Scenario,Nutrient.factor,Nutrient,Web,Landscape), summarize, Persistance=mean(Persist))
ggplot(for.plot2, aes(x=as.factor(Nutrient.factor),y=Persistance,color=as.factor(Scenario))) + geom_violin(trim=T, draw_quantiles = c(0.05, 0.5, 0.95), na.rm=T, scale="width")+ ggtitle("Gamma Diversity")

######
data.end<-foreach (i=c(21:23),.combine = "rbind") %dopar% {
  
  persistance<-data.frame(Scenario=numeric(0),Nutrient.factor=numeric(0),Nutrient=numeric(0),Web=numeric(0),Landscape=numeric(0),Patch=numeric(0), Species=numeric(0), Persist=numeric(0), Growth1=numeric(0),Growth2=numeric(0), Plant=numeric(0))
  
  local<-read.table("../Simulations/FullWeb/SummaryMass.out", sep=",", header=T)
  
  subset.nutrient<-subset(local,local$nut==i)
  Nutrients<- read.table(paste("../Simulations/FullWeb/Nutrients/Nutrient_", i, ".out",sep=""),header=F)
  
  for (w in 1)
  {
    subset.web<-subset(subset.nutrient,subset.nutrient$web==w)
    
    for (l in 1:5)
    {
      subset.landscape<-subset(subset.web,subset.web$landscape==l)
      
      for (p in 0:49)
      {
        subset.patch<- subset(subset.landscape, subset.landscape$patch==p)
        
        for (k in 0:39)
        {
          subset.sp<-subset(subset.patch,subset.patch$species==k)
          
          if (subset.sp$biomass.tend.10k==0 && subset.sp$biomass.tend.20k==0 && subset.sp$biomass.tend!=0) 
          {
            persist<-NA
            mean.net.growth1<-NA
            mean.net.growth2<-NA
            
          } else {
            persist<-(sum(subset.sp$biomass.tend)>0)*1
            mean.net.growth1<-subset.sp$Mean.net.growth1
            mean.net.growth2<-subset.sp$Mean.net.growth2
          }
          if.basal<-subset.sp$if.basal.spp
          if (i>16)
          {scen<-"Heterogenious"} else {scen<-"Homogenious"}
          
          vec<-data.frame(Scenario=scen,Nutrient.factor=i,Nutrient=Nutrients[p+1,], Web=w,Landscape=l,Patch=p,Species=k,Persist=persist, Growth1=mean.net.growth1, Growth2=mean.net.growth2, Plant=if.basal)
          
          persistance<-rbind(persistance,vec)
          
          
        }
      }
    }
  }
  return(persistance)  
}


persistance2<-data.end



for.plot<-ddply(persistance2,.(Scenario,Nutrient.factor,Nutrient,Web,Patch,Landscape, Plant), summarize, Persistance=mean(Persist))

plants<-subset(for.plot,for.plot$Plant==1)
animals<-subset(for.plot,for.plot$Plant==0)
pdf("PlotsHeterogenious.pdf", onefile=T)
ggplot(plants, aes(x=log10(Nutrient),y=Persistance,color=as.factor(Nutrient.factor)))  + geom_point(size=0.2) + geom_smooth(size=0.5) + ggtitle("Plants")
ggplot(animals, aes(x=log10(Nutrient),y=Persistance,color=as.factor(Nutrient.factor)))  + geom_point(size=0.2) + geom_smooth(size=0.5) + ggtitle("Animals")

for.plot1<-ddply(persistance2,.(Scenario,Nutrient.factor,Nutrient,Web,Patch,Landscape), summarize, Persistance=mean(Persist))

ggplot(for.plot1, aes(x=log10(Nutrient),y=Persistance,color=as.factor(Nutrient.factor)))  + geom_point(size=0.2) + geom_smooth(size=0.5) + ggtitle("All Species")


for.plot2<-ddply(persistance2,.(Scenario,Nutrient.factor,Nutrient,Web,Landscape), summarize, Persistance=mean(Persist))
ggplot(for.plot2, aes(x=as.factor(Nutrient.factor),y=Persistance,color=as.factor(Nutrient.factor))) + geom_violin(trim=T, draw_quantiles = c(0.05, 0.5, 0.95), na.rm=T, scale="width") + ggtitle("Gamma Diversity")
dev.off()


