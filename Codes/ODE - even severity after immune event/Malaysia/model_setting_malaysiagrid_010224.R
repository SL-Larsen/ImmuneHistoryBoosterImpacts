
source("ODE_GeneralModel_InfectionVaccination_gradients_v4_Malaysia_Evax.R")

# population distribution
sweep_fun <- function(t = t){
  
  library(tidyverse)
  library(plyr)
  ###Example code for running the sweep with vaccination and time constant NPIs##
  #cvax=c(0,1) #child vaccination on or off, assume that if it's 'on' it could be completed in ~180 days
  cvax=0
  boost.child=0
  boost.adult=c(0, 1)
  child.max=0.9
  #child.boost=c(0,1)
  #boost.adult=c(0, 1)
  #max.boost=0.9
  max.boost.child=0.9
  #max.boost=seq(from=0, to=1, by=0.1)
  max.boost_h=0.9
  max.boost_l=0.91
  #max.boost.child=c(0.1, 0.9)
  cp=seq(from=0, to=0.9, by=0.15)
  #cp=0.579
  #relax=seq(from=0, to=1, by=0.1)
  relax=c(0.1, 0.3)
  base.imm=seq(from=0, to=1, by=0.1)
  #base.imm=0.649
  sd=0.925
  wanevar=c(1/300)
  #boost.start=c(0, 60)
  boost.start=60
  boost.speed='fast'
  #boost.speed=c('fast', 'slow')
  dWVH<-10
  k_h<-5.5
  k_l<-4.97
  base.wane<-c(0.5, 0.9)
  
  check0<-as.data.frame(expand.grid(boost.child=boost.child, cvax=cvax, boost.adult=boost.adult, max.boost_h=max.boost_h, 
                                    max.boost_l=max.boost_l, max.boost.child=max.boost.child, 
                                    wanevar=wanevar, child.max=child.max, cp=cp, base.imm=base.imm, sd=sd, relax=relax, boost.start=boost.start,
                                    boost.speed=boost.speed, dWVH=dWVH, k_h=k_h, k_l=k_l))
  
  #Remove redundant rollout conditions for no booster scenario
  check.boost<-subset(check0, check0$boost.adult==1)
  check.noboost<-subset(check0, check0$boost.adult==0)
  check.noboost$imm_sd_cond<-paste(check.noboost$cp, check.noboost$base.imm, check.noboost$relax, sep='-')
  check.noboost1<-ddply(check.noboost, .(imm_sd_cond), function(x) head(x, 1))
  
  #Add VE variables
  #Only need one set for the no booster scenario (i.e., VE doesn't matter when we don't give a booster)
  check1<-check.boost
  #(1-0.526)/(1-0.571)=(1-vei_HIM)
  check1$vei2<-0.5256
  #(1-0.657)/(1-0.571)=0.799
  check1$veh2<-0.657
  check2<-check.boost
  check2$vei2<-0.8*0.925
  check2$veh2<-0.925
  
  check3<-check.noboost1
  check3$vei2<-0.5256
  check3$veh2<-0.657

  check.allb<-rbind(check1, check2)
  check.all<-rbind(check.allb, check3[,c(1:17, 19:20)])
  #subset(check, check$cvax==1)
  #sweep0<-rbind(boost, noboost.nochild, noboost.child)
  
  #sweep<-check
  sweep<-check.all
  
  num_sweep <- nrow(sweep)
  return(list(sp = sweep, num_sweep = num_sweep))
}

