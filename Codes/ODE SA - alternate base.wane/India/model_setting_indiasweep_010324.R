
source("ODE_GeneralModel_InfectionVaccination_gradients_v4_India_Evax.R")

# population distribution
sweep_fun <- function(t = t){
  
  library(tidyverse)
  
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
  max.boost_h=0.35
  max.boost_l=0.35
  #max.boost.child=c(0.1, 0.9)
  #cp=c(0.2946, 0.3818, 0.3159) #old estimates
  cp=0.436
  #relax=seq(from=0, to=1, by=0.1)
  relax=c(0, 0.1, 0.3)
  base.imm=0.526
  sd=0.288
  wanevar=c(1/300)
  #pboost=c(0, 0.218, 0.27)
  #pboost_h=c(0, 0.019, 0.0727)
  #pboost_l=c(0, 0.0115, 0.0279)
  boost.speed=c('fast', 'slow')
  boost.start=c(0, 60)
  dWVH<-10
  k_h<-2.83
  k_l<-2.54
  
  #imm<-data.frame(cp=cp, base.imm=base.imm, sd0=sd)
  #boost1<-as.data.frame('boost'=pboost)
  imm<-data.frame(cp=cp, base.imm=base.imm, sd0=sd)
  #boost1<-as.data.frame('boost'=pboost)
  check0<-as.data.frame(expand.grid(boost.child=boost.child, cvax=cvax, boost.adult=boost.adult, max.boost_h=max.boost_h, 
                                    max.boost_l=max.boost_l, max.boost.child=max.boost.child, 
                                    wanevar=wanevar, child.max=child.max, cp=cp, base.imm=base.imm, sd=sd, relax=relax, boost.start=boost.start,
                                    boost.speed=boost.speed, dWVH=dWVH, k_h=k_h, k_l=k_l))
  
  #Remove redundant rollout conditions for no booster scenario
  check.boost<-subset(check0, check0$boost.adult==1)
  check.noboost<-subset(check0, check0$boost.adult==0)
  check0.2<-rbind(check.boost, check.noboost[c(1:3),])
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
  
  check3<-check.noboost[c(1:3),]
  check3$vei2<-0.5256
  check3$veh2<-0.657
  
  check.allb<-rbind(check1, check2)
  check.all<-rbind(check.allb, check3)
  #subset(check, check$cvax==1)
  #sweep0<-rbind(boost, noboost.nochild, noboost.child)
  
  #sweep<-check
  sweep<-check.all
  
  num_sweep <- nrow(sweep)
  return(list(sp = sweep, num_sweep = num_sweep))
}

