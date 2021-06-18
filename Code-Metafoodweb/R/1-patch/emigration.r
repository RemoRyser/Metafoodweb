## this script produces a bifurcation diagramm for the 3-species food chain


rm(list=ls())

if (!"EMD" %in% installed.packages()) install.packages("EMD")
library(EMD)



# bifurcation data
bdd <- data.frame(
  q = numeric(),
  b1.min = numeric(),
  b1.max = numeric(),
  b2.min = numeric(),
  b2.max = numeric(),
  b3.min = numeric(),
  b3.max = numeric()
  )


nutrient<-seq(0,0.149,0.001)
for (j in 1:5)
{
for (i in 1:length(nutrient))
{
  out2 <-read.table(paste("../../Simulations/Emigration/timeseries_",i,"_",j,".out",sep=""),header=T)

  
  ### searching mins and maxs
  ## b1min
  ifelse(
    extrema(out2$B0P0)$nextreme<2,
    b1.min <- rep(min(out2$B0P0),25), # if no extrema have been found, take the minimum value 25 times (note: minima are in this case identical to the fixed point value!)
    b1.min <- out2$B0P0[extrema(out2$B0P0)$minindex[1:25]] # otherwise, take the first 25 minima that have been found
  )
  ## b1max
  ifelse(
    extrema(out2$B0P0)$nextreme<2,
    b1.max <- rep(min(out2$B0P0),25),
    b1.max <- out2$B0P0[extrema(out2$B0P0)$maxindex[1:25]]
  )  
  ## b2min
  ifelse(
    extrema(out2$B1P0)$nextreme<2,
    b2.min <- rep(min(out2$B1P0),25),
    b2.min <- out2$B1P0[extrema(out2$B1P0)$minindex[1:25]]
  )
  ## b2max
  ifelse(
    extrema(out2$B1P0)$nextreme<2,
    b2.max <- rep(min(out2$B1P0),25),
    b2.max <- out2$B1P0[extrema(out2$B1P0)$maxindex[1:25]]
  )
  ## b3min
  ifelse(
    extrema(out2$B2P0)$nextreme<2,
    b3.min <- rep(min(out2$B2P0),25),
    b3.min <- out2$B2P0[extrema(out2$B2P0)$minindex[1:25]]
  )
  ## b3max
  ifelse(
    extrema(out2$B2P0)$nextreme<2,
    b3.max <- rep(min(out2$B2P0),25),
    b3.max <- out2$B2P0[extrema(out2$B2P0)$maxindex[1:25]]
  )  
  bddx <- data.frame(
    q = rep(nutrient[i],25),
    b1.min = b1.min,
    b1.max = b1.max,
    b2.min = b2.min,
    b2.max = b2.max,
    b3.min = b3.min,
    b3.max = b3.max
    )
  bdd <- rbind(bdd,bddx)
  print(i)
 
}
}

write.table(bdd,"Emigration.txt", sep=",", row.names = F)

