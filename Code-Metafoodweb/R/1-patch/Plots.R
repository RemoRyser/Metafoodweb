require(gtools)
library(ggplot2)
library(reshape2)
library(plyr)
library(sfsmisc)
library(gridExtra)
library(grid)
library(EMD)



bm<-read.table("Nutrients/Nutrient.txt",header = T, sep=",")
keeps <- c("q","b3.min", "b3.max")
bm<-bm[ , (names(bm) %in% keeps)]
bm<-melt(bm, id.vars="q")
#bm$group[bm$variable %in% c("b1.min","b1.max")]="Plant"
#bm$group[bm$variable %in% c("b2.min","b2.max")]="Herbivore"
bm$group[bm$variable %in% c("b3.min","b3.max")]="Carnivore"
bm$group<- as.factor(bm$group)

pbm<-ggplot(bm,aes(q,value, colour=group))+ 
  geom_point(size=0.3) +
  facet_wrap(~group, nrow=3, scales = "free_y") +
  theme_bw() +theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none") +
  labs(x="Nutrient Supply",y="Biomass min-max") 


bme<-read.table("Emigration/Emigration.txt",header = T, sep=",")
keeps <- c("q","b3.min", "b3.max")
bme<-bme[ , (names(bme) %in% keeps)]
bme<-melt(bme, id.vars="q")
#bme$group[bme$variable %in% c("b1.min","b1.max")]="Plant"
#bme$group[bme$variable %in% c("b2.min","b2.max")]="Herbivore"
bme$group[bme$variable %in% c("b3.min","b3.max")]="Carnivore"
bme$group<- as.factor(bme$group)

pbme<-ggplot(bme,aes(q,value, colour=group))+ 
  geom_point(size=0.3) +
  facet_wrap(~group, nrow=3, scales = "free_y") +
  theme_bw() +theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none") +
  labs(x="Emigration rate",y="Biomass min-max") 

pdf("Figure2.pdf", width=6, height=3)
grid.arrange(pbm, pbme, nrow=2)
dev.off()
