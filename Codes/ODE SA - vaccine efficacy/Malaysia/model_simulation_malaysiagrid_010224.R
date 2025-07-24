# set initial condition and run the ode model
library(deSolve)
library(tidyverse)
require(plyr)

rm(list=ls())

tmax<-195
source("model_setting_malaysiagrid_010224.R")

num_sweep <- sweep_fun(t = (0:tmax))$num_sweep

model_sims <- function(i){
  
  library(tidyverse)
  #i<-1
  t = (1:tmax)
  #t = (1:365)
  
  source("model_setting_malaysiagrid_010224.R")
  spst <- sweep_fun(t = t)
  sweep <- spst$sp
  
  ##Initial state values##
  #i<-1
  #pboost<-sweep$pboost[i]
  boost.start<<-sweep$boost.start[i]
  boost.speed<-sweep$boost.speed[i]
  dWVH<-sweep$dWVH[i]
  cvax<<-sweep$cvax[i]  
  boost.adult<<-sweep$boost.adult[i]
  boost.child<<-sweep$boost.child[i]
  
  k_h<-sweep$k_h[i]
  k_l<-sweep$k_l[i]
  
  max.boost_h<<-sweep$max.boost_h[i]
  max.boost_l<<-sweep$max.boost_l[i]
  
  WVH_h<-19.92
  WVH_l<-21.28
  
  #Define vaccine rate function
  type_III_resp <- function (t, Vmax, WVH, k) {
    
    return(Vmax*(t^k)/((WVH^k) + (t^k)))
  }
  
  #### check type III values for starting coverage####
  ## Slow: 
  slow_high <- type_III_resp(60/7, max.boost_h, WVH_h, k_h)
  slow_low <- type_III_resp(60/7,max.boost_l, WVH_l, k_l)
  
  ## Fast: 
  fast_high <- type_III_resp(60/7,max.boost_h, WVH_h-dWVH, k_h)
  fast_low <- type_III_resp(60/7,max.boost_l, WVH_l-dWVH, k_l)
  
  
  if(boost.adult==0){
    pboost_h<-0
    pboost_l<-0
    #won't use, but adding so the loop will work
    WVH_h<-19.92 
    WVH_l<-21.28
  }
  if(boost.adult==1){
    if(boost.start==0){
      pboost_h<-0
      pboost_l<-0
      if(boost.speed=='slow'){
        WVH_h<-19.92
        WVH_l<-21.28 
      }
      if(boost.speed=='fast'){
        WVH_h<-19.92-dWVH
        WVH_l<-21.28-dWVH
      }}
    if(boost.start==60){
      if(boost.speed=='slow'){
        pboost_h<-slow_high
        pboost_l<-slow_low
        WVH_h<-19.92
        WVH_l<-21.28
      }
      if(boost.speed=='fast'){
        pboost_h<-fast_high
        pboost_l<-fast_low
        WVH_h<-19.92-dWVH
        WVH_l<-21.28-dWVH
      }
    }
  }

  #Define vaccine rate function
  type_III_resp <- function (t, Vmax, WVH, k) {
    
    return(Vmax*(t^k)/((WVH^k) + (t^k)))
  }
  
  lses<-0.5
  hses<-0.5
  #Basic problem: too many unvaccinated people are being put in the infected or recovered classes,
  #leading to 
  dist<-c(0.362*lses, 0.362*hses, 0.59*lses, 0.59*hses, 0.048*lses, 0.048*hses)
  
  start.Ns<-dist*(28552712)
  distI<-c(0.367*lses, 0.367*hses, 0.578*lses, 0.578*hses, 0.055*lses, 0.055*hses) #
  distH<-c(0.058*lses, 0.058*hses, 0.611*lses, 0.611*hses, 0.331*lses, 0.331*hses) #Age distribution for hospitalizations, 
  #obtained by applying hospitalization rates to distI and assuming equal probability by risk group within age strata
  
  base.imm<-sweep$base.imm[i]
  #base.imm<-0.649
  base.wane<-0.5
  
  #To get Is, I subtracted the number of deaths in the past week and didn't take any recoveries away
  report<-(0.178) #Reporting rate, lower bound on reporting with upper bound on CI
  
  #pvax<-c(0.88, 0.65, 0.97, 0.72, 0.97, 0.72) #https://www.pwc.com/my/en/publications/2021/pwc-annual-report-2021/client/supporting-the-malaysian-national-covid-19-immunisation-programm.html
  pvax<-c(0.97, 0.72, 0.97, 0.72, 0.97, 0.72)
  #child vax numbers: 0.89/0.98 for 12-17 year olds vs. adults, assume similar ratio achieved for other child age groups
  #(0.89/0.98)*0.97
  #(0.89/0.98)*0.72
  
  boostvec<-c(0, 0, pboost_h, pboost_l, pboost_h, pboost_l)
  start.Nvax<-start.Ns*pvax*(1-boostvec)
  start.Nboost<-start.Ns*pvax*boostvec
  start.Nunvax<-start.Ns-(start.Nvax+start.Nboost)
  
  #start.Is<-distI*(12600)*(1/report) #To get current Is, did sum(I.new, t=0:7) with rounding 
  start.Is<-distI*(11236)*(1/report) #To get current Is, did sum(I.new, t=0:7) with rounding 
  start.Iunvax<-start.Is/(1+((0.1*start.Nvax)/start.Nunvax))
  start.Ivax<-(start.Is-start.Iunvax)*(1-boostvec)
  start.Iboost<-(start.Is-start.Iunvax)*(boostvec)
  
  start.As<-(1/0.4)*(1-0.4)*start.Is #1/p gives total infections and (1-p)*infections gives asymptomatic infections 
  start.Aunvax<-start.As/(1+((0.1*start.Nvax)/start.Nunvax))
  start.Avax<-(start.As-start.Aunvax)*(1-boostvec)
  start.Aboost<-(start.As-start.Aunvax)*(boostvec)
  
  start.Es<-start.Is+start.As
  start.Eunvax<-start.Es/(1+((0.1*start.Nvax)/start.Nunvax))
  start.Evax<-(start.Es-start.Eunvax)*(1-boostvec)
  start.Eboost<-(start.Es-start.Eunvax)*(boostvec)
  
  BaseR<-base.imm*(1-base.wane)
  start.Runvax<-c(rep(BaseR, 6))*start.Nunvax #Upper bound based on all infections with no waning
  start.Rvax<-c(rep(BaseR, 6))*start.Nvax #Upper bound based on all infections with no waning
  start.Rboost<-c(rep(BaseR, 6))*start.Nboost #Upper bound based on all infections with no waning
  
  BaseSwane<-base.imm*base.wane
  start.Sunvax<-rep(NA, 6)
  start.Sunvax2<-rep(NA, 6)
  for(l in 1:6){
    if(start.Nunvax[l]-(start.Iunvax[l]+start.Eunvax[l]+start.Aunvax[l]+start.Runvax[l])-(BaseSwane)*start.Nunvax[l]<0){
      start.Sunvax[l]<-0
      start.Sunvax2[l]<-start.Nunvax[l]-(start.Iunvax[l]+start.Eunvax[l]+start.Aunvax[l]+start.Runvax[l])
    }
    if(start.Nunvax[l]-(start.Iunvax[l]+start.Eunvax[l]+start.Aunvax[l]+start.Runvax[l])-(BaseSwane)*start.Nunvax[l]>=0){
      start.Sunvax[l]<-start.Nunvax[l]-(start.Iunvax[l]+start.Eunvax[l]+start.Aunvax[l]+start.Runvax[l])-(BaseSwane)*start.Nunvax[l]
      start.Sunvax2[l]<-BaseSwane*start.Nunvax[l]
    }}
  
  start.Svax<-rep(NA, 6)
  start.Svax2<-rep(NA, 6)
  for(j in 1:6){
    if(start.Nvax[j]-(start.Ivax[j]+start.Evax[j]+start.Avax[j]+start.Rvax[j])-(BaseSwane)*start.Nvax[j]<0){
      start.Svax[j]<-0
      start.Svax2[j]<-start.Nvax[j]-(start.Ivax[j]+start.Evax[j]+start.Avax[j]+start.Rvax[j])
    }
    if(start.Nvax[j]-(start.Ivax[j]+start.Evax[j]+start.Avax[j]+start.Rvax[j])-(BaseSwane)*start.Nvax[j]>=0){
      start.Svax[j]<-start.Nvax[j]-(start.Ivax[j]+start.Evax[j]+start.Avax[j]+start.Rvax[j])-(BaseSwane)*start.Nvax[j]
      start.Svax2[j]<-BaseSwane*start.Nvax[j]
    }}
  
  start.Sboost<-rep(NA, 6)
  start.Sboost2<-rep(NA, 6)
  for(k in 1:6){
    if(start.Nboost[k]-(start.Iboost[k]+start.Eboost[k]+start.Aboost[k]+start.Rboost[k])-(BaseSwane)*start.Nboost[k]<0){
      start.Sboost[k]<-0
      start.Sboost2[k]<-start.Nboost[k]-(start.Iboost[k]+start.Eboost[k]+start.Aboost[k]+start.Rboost[k])
    }
    if(start.Nboost[k]-(start.Iboost[k]+start.Eboost[k]+start.Aboost[k]+start.Rboost[k])-(BaseSwane)*start.Nboost[k]>=0){
      start.Sboost[k]<-start.Nboost[k]-(start.Iboost[k]+start.Eboost[k]+start.Aboost[k]+start.Rboost[k])-(BaseSwane)*start.Nboost[k]
      start.Sboost2[k]<-BaseSwane*start.Nboost[k]
    }}
  
  #sum(start.Sunvax, start.Sunvax2, start.Svax, start.Svax2, start.Sboost, start.Sboost2, start.Runvax, start.Rvax, start.Rboost,
  #    start.Eunvax, start.Evax, start.Eboost, start.Aunvax, start.Avax, start.Aboost, start.Iunvax, start.Ivax, start.Iboost)
  
  start = c(Sch=start.Sunvax[1], Scl=start.Sunvax[2], Sah=start.Sunvax[3], Sal=start.Sunvax[4], Seh=start.Sunvax[5], Sel=start.Sunvax[6], #S
            Ech=start.Eunvax[1], Ecl=start.Eunvax[2], Eah=start.Eunvax[3], Eal=start.Eunvax[4], Eeh=start.Eunvax[5], Eel=start.Eunvax[6], #E
            Ach=start.Aunvax[1], Acl=start.Aunvax[2], Aah=start.Aunvax[3], Aal=start.Aunvax[4], Aeh=start.Aunvax[5], Ael=start.Aunvax[6], #A
            Ich=start.Iunvax[1], Icl=start.Iunvax[2], Iah=start.Iunvax[3], Ial=start.Iunvax[4], Ieh=start.Iunvax[5], Iel=start.Iunvax[6], #I
            Rch=start.Runvax[1], Rcl=start.Runvax[2], Rah=start.Runvax[3], Ral=start.Runvax[4], Reh=start.Runvax[5], Rel=start.Runvax[6], #R
            Dch=0, Dcl=0, Dah=0, Dal=0, Deh=0, Del=0, #D
            
            Sch2=start.Sunvax2[1], Scl2=start.Sunvax2[2], Sah2=start.Sunvax2[3], Sal2=start.Sunvax2[4], Seh2=start.Sunvax2[5], Sel2=start.Sunvax2[6], #S
            Ech2=0, Ecl2=0, Eah2=0, Eal2=0, Eeh2=0, Eel2=0, #E
            Ach2=0, Acl2=0, Aah2=0, Aal2=0, Aeh2=0, Ael2=0, #A
            Ich2=0, Icl2=0, Iah2=0, Ial2=0, Ieh2=0, Iel2=0, #I
            Rch2=0, Rcl2=0, Rah2=0, Ral2=0, Reh2=0, Rel2=0, #R
            Dch2=0, Dcl2=0, Dah2=0, Dal2=0, Deh2=0, Del2=0, 
            
            SVch=start.Svax[1], SVcl=start.Svax[2], SVah=start.Svax[3], SVal=start.Svax[4], SVeh=start.Svax[5], SVel=start.Svax[6], #SV
            EVch=start.Evax[1], EVcl=start.Evax[2], EVah=start.Evax[3], EVal=start.Evax[4], EVeh=start.Evax[5], EVel=start.Evax[6], #EV
            AVch=start.Avax[1], AVcl=start.Avax[2], AVah=start.Avax[3], AVal=start.Avax[4], AVeh=start.Avax[5], AVel=start.Avax[6], #AV
            IVch=start.Ivax[1], IVcl=start.Ivax[2], IVah=start.Ivax[3], IVal=start.Ivax[4], IVeh=start.Ivax[5], IVel=start.Ivax[6], #IV
            RVch=start.Rvax[1], RVcl=start.Rvax[2], RVah=start.Rvax[3], RVal=start.Rvax[4], RVeh=start.Rvax[5], RVel=start.Rvax[6], #RV
            DVch=0, DVcl=0, DVah=0, DVal=0, DVeh=0, DVel=0,
            
            SVch2=start.Svax2[1], SVcl2=start.Svax2[2], SVah2=start.Svax2[3], SVal2=start.Svax2[4], SVeh2=start.Svax2[5], SVel2=start.Svax2[6], #S
            EVch2=0, EVcl2=0, EVah2=0, EVal2=0, EVeh2=0, EVel2=0, #E
            AVch2=0, AVcl2=0, AVah2=0, AVal2=0, AVeh2=0, AVel2=0, #A
            IVch2=0, IVcl2=0, IVah2=0, IVal2=0, IVeh2=0, IVel2=0, #I
            RVch2=0, RVcl2=0, RVah2=0, RVal2=0, RVeh2=0, RVel2=0, #R
            DVch2=0, DVcl2=0, DVah2=0, DVal2=0, DVeh2=0, DVel2=0,
            
            SBch=start.Sboost[1], SBcl=start.Sboost[2], SBah=start.Sboost[3], SBal=start.Sboost[4], SBeh=start.Sboost[5], SBel=start.Sboost[6], #SV
            EBch=start.Eboost[1], EBcl=start.Eboost[2], EBah=start.Eboost[3], EBal=start.Eboost[4], EBeh=start.Eboost[5], EBel=start.Eboost[6], #EV
            ABch=start.Aboost[1], ABcl=start.Aboost[2], ABah=start.Aboost[3], ABal=start.Aboost[4], ABeh=start.Aboost[5], ABel=start.Aboost[6], #AV
            IBch=start.Iboost[1], IBcl=start.Iboost[2], IBah=start.Iboost[3], IBal=start.Iboost[4], IBeh=start.Iboost[5], IBel=start.Iboost[6], #IV
            RBch=start.Rboost[1], RBcl=start.Rboost[2], RBah=start.Rboost[3], RBal=start.Rboost[4], RBeh=start.Rboost[5], RBel=start.Rboost[6], #RV
            DBch=0, DBcl=0, DBah=0, DBal=0, DBeh=0, DBel=0,
            
            SBch2=start.Sboost2[1], SBcl2=start.Sboost2[2], SBah2=start.Sboost2[3], SBal2=start.Sboost2[4], SBeh2=start.Sboost2[5], SBel2=start.Sboost2[6], #SV
            EBch2=0, EBcl2=0, EBah2=0, EBal2=0, EBeh2=0, EBel2=0, #EV
            ABch2=0, ABcl2=0, ABah2=0, ABal2=0, ABeh2=0, ABel2=0, #AV
            IBch2=0, IBcl2=0, IBah2=0, IBal2=0, IBeh2=0, IBel2=0, #IV
            RBch2=0, RBcl2=0, RBah2=0, RBal2=0, RBeh2=0, RBel2=0, #RV
            DBch2=0, DBcl2=0, DBah2=0, DBal2=0, DBeh2=0, DBel2=0) #DV
  
  #Get vax rate
  
  #sum(start)
  #sum(start.Ns)

  #R0.var <- sweep$base.r0[i]
  
  #Get bl
  #bl.new<-(R0.var/130.931)
  #bl<-(10/130.931) #Largest eigenvalue from mathematica workbook
  # match ABM
  # ABM_scale <- 1.78
  # 
  # bl<-(1.79/26.14)*(18.266/26.14)*ABM_scale 
  # bh<-bl/2.63*ABM_scale
  
  bl <- 0.25*0.427*0.7*1.25
  bh <- 0.25*0.427*0.4*1.25
  
  #sd0<-0.75 #AK changed #Does this apply?
  #Severity
  m<-0.6  #Mortality multiplier, impacts phi (hospitalization rate)
  m2<-0.65
  
  #Escape
  wanevar<-sweep$wanevar[i]
  
  #Vaccine efficacy
  #vei<-1-sweep$vei[i]
  #veic<-1-sweep$veic[i]
  vei<-1-0
  veic<-1-0
  vei2<-1-sweep$vei2[i]
  ves<-1-0
  #veh<-1-sweep$veh[i] #now applied to rho
  veh<-1-0.7 #Changed this here
  veh2<-1-sweep$veh2[i] #now applied to rho
  

  cp<<-sweep$cp[i]
  max.boost_h<<-sweep$max.boost_h[i]
  max.boost_l<<-sweep$max.boost_l[i]
  max.boost.child<<-sweep$max.boost.child[i]
  child.max<<-sweep$child.max[i]
  
  #boost.adult<<-0
  SES.sev<-c(3.1, 2.4, 1) 
  lowCFR<-0.2*0.5 # cut CFR in half
  sd0<-sweep$sd[i]
  relax<-sweep$relax[i]
  
  #### Vaccination ####
  type_III_resp <- function (t, Vmax, WVH, k) {
    
    return(Vmax*(t^k)/((WVH^k) + (t^k)))
  }
  
  vacc_weekly <- function (Vmax, WVH, k){
    datOut <- tibble(Week = 1:60, Vac = type_III_resp(1:60, Vmax,WVH, k)) %>%
      mutate(dailyVac =  c(Vac[1],diff(Vac)),
             VacCumSum = cumsum(dailyVac)) %>% 
      arrange(Week) %>%
      mutate(Day = Week*7-7+1, dailyVac=dailyVac/7)
    
    vaccF <- splinefun(x = datOut$Day, y = datOut$dailyVac)
    return(vaccF)
  }
  
  #if(boost.adult==0){max.boost<-0}
  vacc_H<<-vacc_weekly(max.boost_h, WVH_h, k_h)
  vacc_L<<-vacc_weekly(max.boost_l, WVH_l, k_l)
  
  params <<- c('bI'=1,'bA'=1, 
               'nu_Ech'=0.4, 'sigma_Ech'=1/4, 'nu_Ecl'=0.4, 'sigma_Ecl'=1/4, 'nu_Eah'=0.4, 'sigma_Eah'=1/4, 'nu_Eal'=0.4, 'sigma_Eal'=1/4, 'nu_Eeh'=0.4, 'sigma_Eeh'=1/4, 'nu_Eel'=0.4, 'sigma_Eel'=1/4, #E
               'gA_Ach'=1/7, 'gA_Acl'=1/7, 'gA_Aah'=1/7, 'gA_Aal'=1/7, 'gA_Aeh'=1/7, 'gA_Ael'=1/7, #A
               'gI_Ich'=1/7, 'gI_Icl'=1/7, 'gI_Iah'=1/7, 'gI_Ial'=1/7, 'gI_Ieh'=1/7, 'gI_Iel'=1/7, 
               'rho_Ich'=0.0022*lowCFR,'rho_Icl'=0.0022*SES.sev[1]*lowCFR, 'rho_Iah'=0.0178*lowCFR,  'rho_Ial'=0.0178*SES.sev[2]*lowCFR,  'rho_Ieh'=0.173*lowCFR,  'rho_Iel'=0.173*SES.sev[3]*lowCFR, 
               'nuv_EVch'=ves*0.4, 'sigmav_EVch'=1/4, 'nuv_EVcl'=ves*0.4, 'sigmav_EVcl'=1/4, 'nuv_EVah'=ves*0.4, 'sigmav_EVah'=1/4, 'nuv_EVal'=ves*0.4, 'sigmav_EVal'=1/4, 'nuv_EVeh'=ves*0.4, 
               'sigmav_EVeh'=1/4, 'nuv_EVel'=ves*0.4, 'sigmav_EVel'=1/4, 
               'gAV_AVch'=1/7, 'gAV_AVcl'=1/7, 'gAV_AVah'=1/7, 'gAV_AVal'=1/7, 'gAV_AVeh'=1/7, 'gAV_AVel'=1/7, 
               'gIV_IVch'=1/7, 'gIV_IVcl'=1/7, 'gIV_IVah'=1/7, 'gIV_IVal'=1/7, 'gIV_IVeh'=1/7, 'gIV_IVel'=1/7,
               'rhov_IVch'=veh*0.0022*SES.sev[1]*lowCFR, 'rhov_IVcl'=veh*0.0022*SES.sev[1]*lowCFR,  'rhov_IVah'=veh*0.0178*SES.sev[2]*lowCFR, 'rhov_IVal'=veh*0.0178*SES.sev[2]*lowCFR, 'rhov_IVeh'=veh*0.173*SES.sev[3]*lowCFR, 'rhov_IVel'=veh*0.173*SES.sev[3]*lowCFR,
               'rhov_IVch2'=veh2*0.0022*SES.sev[1]*lowCFR, 'rhov_IVcl2'=veh2*0.0022*SES.sev[1]*lowCFR,  'rhov_IVah2'=veh2*0.0178*SES.sev[2]*lowCFR, 'rhov_IVal2'=veh2*0.0178*SES.sev[2]*lowCFR, 'rhov_IVeh2'=veh2*0.173*SES.sev[3]*lowCFR, 'rhov_IVel2'=veh2*0.173*SES.sev[3]*lowCFR,
               'omega'=wanevar, 'vei'=vei, 'veic'=veic, 'vei2'=vei2, 'betah'=bh, 'betal'=bl, 'sd'=sd0*(1+relax)) 
 
  
  model_out <- as.data.frame(ode(y = start, times = seq(from = 0, to = tmax, by = 1), 
                                 fun = COVID_vocBoost_hybrid_v4, parms = params, method='ode45'))
  
  return(list(mod_op = model_out, start = start, params = params))
}


# =============================================================================


# ========================
model_out_list <- list()
model_smr_list <- list()
cvec_list <- list()
sweep <- sweep_fun(t = 0:195)$sp # ()
#sweep <- sweep_fun(t = 0:500)$sp # ()
sweep$sweepnum<-c(1:nrow(sweep))


#each run is ~2 to 3 seconds
#i<-29000
#test<-model_sims(i)
#system.time(model_sims(i))
#out_test<-test$mod_op
#head(out_test)

for (i in 1:num_sweep) {
  #i<-1
  test<-model_sims(i)
  model_out_list[[i]] <- test$mod_op
  model_out_list[[i]]$sweepnum <-i
}

# format of output
model_out_dat <- bind_rows(model_out_list) %>% left_join(sweep)
summary(model_out_dat)
save.image('../../../../Data/ODE_SA_Efficacy/sweep_010224_malaysiagrid_Evax.RData')

end<-subset(model_out_dat, model_out_dat$time==tmax)

end$DeathsAll<-end$Dch+end$Dch2+end$DVch+end$DVch2+end$DBch+end$DBch2+
  end$Dah+end$Dah2+end$DVah+end$DVah2+end$DBah+end$DBah2+
  end$Deh+end$Deh2+end$DVeh+end$DVeh2+end$DBeh+end$DBeh2+
  end$Dcl+end$Dcl2+end$DVcl+end$DVcl2+end$DBcl+end$DBcl2+
  end$Dal+end$Dal2+end$DVal+end$DVal2+end$DBal+end$DBal2+
  end$Del+end$Del2+end$DVel+end$DVel2+end$DBel+end$DBel2

saveRDS(end, '../../../../Data/ODE_SA_Efficacy/endfile_sweep010224_malaysiagrid_Evax.RDS')

