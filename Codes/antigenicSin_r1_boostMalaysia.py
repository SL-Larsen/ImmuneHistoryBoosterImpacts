########################################
class Unbuffered(object):
   def __init__(self, stream):
       self.stream = stream
   def write(self, data):
       self.stream.write(data)
       self.stream.flush()
   def writelines(self, datas):
       self.stream.writelines(datas)
       self.stream.flush()
   def __getattr__(self, attr):
       return getattr(self.stream, attr)

import numpy as np
import random as rnd
import pandas as pnd
from scipy import stats
import math
import sys
import time
sys.stdout = Unbuffered(sys.stdout)
start_time = time.time()

#grid = pnd.read_csv("C:/Users/sophi/Box/WHO_project/Codes/Cluster_Jobs_2023/May/leap_edits16/Grid_Param.csv")
#indexID = 1 ## update for BioCluster

grid = pnd.read_csv('Grid_Param.csv')
indexID = int(sys.argv[1])
print(indexID)


########### set counts of immune histories ###########
## no boosters 
count_Delta = 0
count_Delta_vaccine = 0
count_two = 0
count_three = 0
count_four = 0
count_Omicron = 0
count_Omicron_star = 0
count_Omicron_star_vaccine = 0
count_Omicron_vaccine = 0
count_Vaccine = 0
count_WT = 0
count_WT_vaccine = 0 
## with boosters 
count_WT_vaccine_boost = 0
count_Omicron_vaccine_boost = 0
count_Delta_vaccine_boost = 0
count_Omicron_star_vaccine_boost = 0
count_Vaccine_boost = 0

## track number of people vaccinated
numVaccinatedH = 0
numVaccinatedL = 0
## track number of people boosted
numBoostedH = 0
numBoostedL = 0

############# Set Parameters  ####################
multPop = 1
N_child = grid['child'][indexID-1]*multPop ## 334165328 make pop smaller by 10^2
N_adult = grid['adult'][indexID-1]*multPop ## 898812536 make pop smaller by 10^2
N_old = grid['old'][indexID-1]*multPop ## 156659582 make pop smaller by 10^2
N_L = (N_child + N_adult + N_old)
N_H = (N_child + N_adult + N_old)

count_Naive = N_L + N_H

# proportion that are not children
nonchild = grid['nonchild'][indexID-1]

# contact base
b_LL = grid['c_LL'][indexID-1]  ## low contacts low
b_LH = grid['c_LH'][indexID-1] ## low contacts high
b_HL = grid['c_HL'][indexID-1] ## high contacts low
b_HH = grid['c_HH'][indexID-1] ## high contacts high  

## stringency
stringency = grid["string_wild"][indexID-1]
## multiply all contact by stringency
c_LL = b_LL*stringency  ## low contacts low
c_LH = b_LH*stringency ## low contacts high
c_HL = b_HL*stringency ## high contacts low
c_HH = b_HH*stringency ## high contacts high  

epsilon = 1/4 ## (4 day exposed period )
gamma = 1/7 ## asymptomatic infectious period
w = grid['w'][indexID-1] ## waning immunity
b = grid['b'][indexID-1] ## modify waning for vaccinated individuals
mu_child = 0.189 ## SAR for WT
mu_adult = 0.189 ## SAR for WT
mu_old = 0.189 ## SAR for WT
lowCFR = 0.2 ## reduce deaths in range with estimates
alphaL1_child =   0.4*0.0069*lowCFR ## death rate
alphaL1_adult =   0.4*0.0426*lowCFR ## death rate
alphaL1_old =   0.4*0.173*lowCFR ## death rate
alphaH1_child =  0.4*0.0022*lowCFR  ## death rate
alphaH1_adult =  0.4*0.0178*lowCFR  ## death rate
alphaH1_old =  0.4*0.173*lowCFR  ## death rate

si = 0.27 ## scale infectiousness of secondary infections

countInfected_H = 0
countInfected_L = 0

## determine whether to boost
boost = grid['boost'][indexID-1]
## booster shape
boost_shape = grid['boost_shape'][indexID-1]
## constant rates 
vax_constant_H = 0 #grid['vax_constant_H'][indexID-1]
vax_constant_L = 0 #grid['vax_constant_L'][indexID-1]
## vaccinate during omicron*?
vaccine_star = grid['vaccine_star'][indexID-1]
## maximum for boosters
Vm_H_boost = grid['MH_boost'][indexID-1]
Vm_L_boost = grid['ML_boost'][indexID-1]
## timing for boosters
Wh_H_boost = grid['WH_boost'][indexID-1]
Wh_L_boost = grid['WL_boost'][indexID-1]
## shape for boosters
kH_boost = grid["kH_boost"][indexID-1]
kL_boost = grid["kL_boost"][indexID-1]


######## Serotype parameters during wild-type covid (these will be applied to severe disease, infection will be modified) ##########
k_Delta = 0.61
k_Delta_vaccine = 0.11

k_two = 0.16
k_three = 0.1
k_four = 0.050

k_Omicron = 0.85
k_Omicron_star = 0.95
k_Omicron_star_vaccine = 0.15
k_Omicron_vaccine = 0.05
k_Vaccine = 0.58
k_WT = 0.41
k_WT_vaccine = 0.32 


## WANED serotype parameters - how much to dampen their protection when waned 
dampening1 = 0.4 ## percentage of protection left after waning - 1 immune event
dampening2 = 0.7 ## percentage of protection left after waning - 2 immune events
dampening3 = 0.85 ## percentage of protection left after waning - 3 or more immune events

## infection penalty - reduction in protection against infection compared to severe disease
inf_pen = grid['inf_vs_sev'][indexID-1]

## vaccination parameters
Vm_H = grid['MH'][indexID-1]
Wh_H = grid['WH'][indexID-1]
h_H = grid['hH'][indexID-1]
Vm_L = grid['ML'][indexID-1]
Wh_L = grid['WL'][indexID-1]
h_L = grid['hL'][indexID-1]

############# Function transition rates ############################

## vaccination function (functional) (time in days)
def vaxTrend(Vm, Wh, h, t, timeplus, unVax):
  vaxperc_top = Vm*((t/7)**h)
  vaxperc_bottom = Wh**h + ((t/7)**h)
  vaxperc = vaxperc_top/vaxperc_bottom
  lead_vaxperc = (Vm*(((t+timeplus)/7)**h))/(Wh**h + (((t+timeplus)/7)**h))
  vaxtotal = lead_vaxperc - vaxperc
  ## if at least one person should be vaccinated, do a rounding to get an integer
  if vaxtotal*unVax >= 0:
    roundvax = math.ceil(vaxtotal*unVax)
  else: 
    roundvax = 0
  return(int(roundvax))

## vaccination function (constant)
def vaxConstant(vax_constant, Vm, numVax, pop_t): ## pop_t is the population we are vaccinating out of, numVax is the total that have received this vaccine
  ## if at least one person should be vaccinated, do a rounding to get an integer
  if vax_constant*pop_t >= 0 and numVax/pop_t < Vm: ## if there are unvaccinated left, and the max hasn't been reached
    roundvax = math.ceil(vax_constant*pop_t)
  else: 
    roundvax = 0
  return(int(roundvax))

######### Do Tau function #########

## events: "contLL","contLH","contHL","contHH","infectiousL","infectiousH","death_or_recL","death_or_recH","waneFL","waneFH","waneSL","waneSH"

def do_Tau_event(event,countInfected_H, countInfected_L, variant):
  if event == 'contLL_sus':
    h = rnd.choice(Sus_L)
    g = rnd.choice(Pop_L)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist:
        h.contact_event_L()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_L()
  elif event == 'contLH_sus':
    h = rnd.choice(Sus_L)
    g = rnd.choice(Pop_H)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist:
        h.contact_event_L()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_L()
  elif event == 'contHH_sus':
    h = rnd.choice(Sus_H)
    g = rnd.choice(Pop_H)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist: ## len(g.immunityhist) == 0:
        h.contact_event_H()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_H()
  elif event == 'contHL_sus':
    h = rnd.choice(Sus_H)
    g = rnd.choice(Pop_L)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist: ## len(g.immunityhist) == 0:
        h.contact_event_H()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_H()
  elif event == 'contLL_nat':
    h = rnd.choice(Rec_nat_L)
    g = rnd.choice(Pop_L)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist:
        h.contact_event_L()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_L()
  elif event == 'contLH_nat':
    h = rnd.choice(Rec_nat_L)
    g = rnd.choice(Pop_H)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist:
        h.contact_event_L()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_L()
  elif event == 'contHH_nat':
    h = rnd.choice(Rec_nat_H)
    g = rnd.choice(Pop_H)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist: ## len(g.immunityhist) == 0:
        h.contact_event_H()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_H()
  elif event == 'contHL_nat':
    h = rnd.choice(Rec_nat_H)
    g = rnd.choice(Pop_L)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist: ## len(g.immunityhist) == 0:
        h.contact_event_H()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_H()
  elif event == 'contLL_vax':
    h = rnd.choice(Rec_vax_L)
    g = rnd.choice(Pop_L)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist:
        h.contact_event_L()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_L()
  elif event == 'contLH_vax':
    h = rnd.choice(Rec_vax_L)
    g = rnd.choice(Pop_H)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist:
        h.contact_event_L()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_L()
  elif event == 'contHH_vax':
    h = rnd.choice(Rec_vax_H)
    g = rnd.choice(Pop_H)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist: ## len(g.immunityhist) == 0:
        h.contact_event_H()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_H()
  elif event == 'contHL_vax':
    h = rnd.choice(Rec_vax_H)
    g = rnd.choice(Pop_L)
    if g.infectious == True:
      if len(g.immunityhist) == 0: ## 'Booster' not in g.immunityhist and 'Vaccine' not in g.immunityhist: ## len(g.immunityhist) == 0:
        h.contact_event_H()
      else: 
        ## scale down infectiousness for secondary infections
        if np.random.binomial(1,si)==1:
          h.contact_event_H()
  elif event == 'infectiousL':
    h = rnd.choice(Exp_L)
    h.can_infect_others_L()
    countInfected_L +=1
  elif event == 'infectiousH':
    h = rnd.choice(Exp_H)
    h.can_infect_others_H()
    countInfected_H += 1
  elif event == 'death_or_recL':
    h = rnd.choice(Inf_L)
    if h.age == 'child':
      h.death_or_recover_L(alphaL1_child, variant)
    elif h.age == 'adult':
      h.death_or_recover_L(alphaL1_adult, variant)
    elif h.age == 'old':
      h.death_or_recover_L(alphaL1_old,variant)
  elif event == 'death_or_recH':
    h = rnd.choice(Inf_H)
    if h.age == 'child':
      h.death_or_recover_H(alphaH1_child, variant)
    elif h.age == 'adult':
      h.death_or_recover_H(alphaH1_adult,variant)
    elif h.age == 'old':
      h.death_or_recover_H(alphaH1_old,variant)
  elif event == 'waneFL':
    h = rnd.choice(Rec_nat_L)
    h.has_waned_nat_L()
  elif event == 'waneFH':
    h = rnd.choice(Rec_nat_H)
    h.has_waned_nat_H()
  elif event == 'waneSL':
    h = rnd.choice(Rec_vax_L)
    h.has_waned_vax_L()
  else: 
    h = rnd.choice(Rec_vax_H)
    h.has_waned_vax_H()
  return(countInfected_H, countInfected_L)


################## Tau leaping algorithm ######################

# It returns time_to_next_event,event

n_c = 20 # critical number for tau-leaping 
error = 0.08 # error allowed in tau-leaping

crit_tot = 0
total_rate = 0

############# BEGIN "choose event tau" function ############

def choose_event_tau(c_LL, c_LH, c_HL, c_HH, epsilon, gamma, 
S_L, S_H, E_L, E_H, I_L, I_H, b, w, R_nat_L, R_nat_H, R_vax_L, R_vax_H, countInfected_H, countInfected_L, variant): 

## reset all lists, variables

  no_critlist = []
  no_critlist_string = []
  mu_list = []
  sigma_sq_list = []
  critlist = []
  critlist_string = []
  crit_tot = 0
  total_rate = 0
  do_nothing = []

  ## we will use this function below to calculate rates of each of our reactions given a parameter and a state variable
  def calc_rate(rate, variable, total_rate, crit_tot, name):
    full_rate = rate*variable
    total_rate += full_rate
    if full_rate == 0:
      do_nothing.append(full_rate) ## if the rate is 0, we are not going to do anything with it; it will never be one of the minimums
    elif full_rate > 0 and variable > n_c:
      no_critlist.append(full_rate)
      no_critlist_string.append(name)
      mu_i = 1/(rate*full_rate)
      sigma_i = 1/(full_rate*(rate**2))
      mu_list.append(mu_i)
      sigma_sq_list.append(sigma_i)
    else: 
      critlist.append(full_rate)
      critlist_string.append(name)
      crit_tot = crit_tot + full_rate
    return full_rate, total_rate, crit_tot

  ## The rule is: if the rate is zero (so, no chance of firing) or the number of agents
  ## is high enough that firing would not exhaust the remaining agents, the event is noncritical 
  ## and is allowed to fire

  ## we go sequentially by each of the below, adding to the total rate and critital rxns at each step

  ## IMPORTANT: the contact in a deterministic SEIR model is mu*c*S*I/N. Here, a contact function is 
  ## implemented as c*S for how much contact a susceptible has with the whole population,
  #  and later, we choose a random infected individual from the population (I/N)
  # and run the get_exposed function to accomplish c*I/N. While the result of either method
  # will be equivalent, it is mathematically MUCH easier to implement this way in terms of the Tau leap, because
  # of the presence of partial derivatives. 

  ####### contact: fully susceptibles #########
  contLL_sus, total_rate, crit_tot = calc_rate(c_LL,S_L, total_rate, crit_tot, "contLL_sus")
  contLH_sus, total_rate, crit_tot = calc_rate(c_LH,S_L, total_rate, crit_tot, "contLH_sus")
  contHL_sus, total_rate, crit_tot = calc_rate(c_HL,S_H, total_rate, crit_tot, "contHL_sus")
  contHH_sus, total_rate, crit_tot = calc_rate(c_HH,S_H, total_rate, crit_tot, "contHH_sus")

####### contact: recovered from natural infection #########
  contLL_nat, total_rate, crit_tot = calc_rate(c_LL,R_nat_L, total_rate, crit_tot, "contLL_nat")
  contLH_nat, total_rate, crit_tot = calc_rate(c_LH,R_nat_L, total_rate, crit_tot, "contLH_nat")
  contHL_nat,total_rate, crit_tot = calc_rate(c_HL,R_nat_H, total_rate, crit_tot, "contHL_nat")
  contHH_nat,total_rate, crit_tot = calc_rate(c_HH,R_nat_H, total_rate, crit_tot, "contHH_nat")

####### contact: recovered from vaccine #########
  contLL_vax, total_rate, crit_tot = calc_rate(c_LL,R_vax_L, total_rate, crit_tot, "contLL_vax")
  contLH_vax, total_rate, crit_tot = calc_rate(c_LH,R_vax_L, total_rate, crit_tot, "contLH_vax")
  contHL_vax,total_rate, crit_tot = calc_rate(c_HL,R_vax_H, total_rate, crit_tot, "contHL_vax")
  contHH_vax,total_rate, crit_tot = calc_rate(c_HH,R_vax_H, total_rate, crit_tot, "contHH_vax")

######### infectious, exposed, recovered ##########
  infectiousL, total_rate, crit_tot = calc_rate(epsilon,E_L, total_rate, crit_tot, "infectiousL")
  infectiousH, total_rate, crit_tot = calc_rate(epsilon,E_H, total_rate, crit_tot, "infectiousH")

  death_or_recL, total_rate, crit_tot = calc_rate(gamma,I_L, total_rate, crit_tot, "death_or_recL")
  death_or_recH, total_rate, crit_tot = calc_rate(gamma,I_H, total_rate, crit_tot, "death_or_recH")

  waneFL, total_rate, crit_tot = calc_rate(w,R_nat_L, total_rate, crit_tot, "waneFL") # natural wane
  waneFH, total_rate, crit_tot = calc_rate(w,R_nat_H, total_rate, crit_tot, "waneFH") # natural  wane

  waneSL, total_rate, crit_tot = calc_rate(b*w,R_vax_L, total_rate, crit_tot, "waneSL") # vaccinated wane
  waneSH, total_rate, crit_tot = calc_rate(b*w,R_vax_H, total_rate, crit_tot, "waneSH") # vaccinated wane
    
## If the noncritical list has no members,
#  then nothing is allowed to fire by tau-leaping (essentially, proceed by Gillespie)
# Otherwise, calculate the min that should fire

  ## tau options
  ## if both the lists are nonempty
  if len(no_critlist) > 0 and len(critlist) > 0:
    ## get the minimum of all the possible values
    error_tot = error*total_rate
    mins1 = min([i*error_tot for i in mu_list])
    mins2 = min([i*(error_tot**2) for i in sigma_sq_list])
    ## get tau prime
    tau_prime = min(mins1, mins2)
    ## get tau 2prime
    tau_2prime = np.random.exponential(1/crit_tot)
  elif len(no_critlist)>0 and len(critlist) == 0:
    ## get the minimum of all the possible values
    error_tot = error*total_rate
    mins1 = min([i*error_tot for i in mu_list])
    mins2 = min([i*(error_tot**2) for i in sigma_sq_list])
    ## get tau prime
    tau_prime = min(mins1, mins2)
    ## since length of crit list is 0, tau2prime is infinity
    tau_2prime = tau_prime + 1 ## just set it to tau_prime + 1, because tau2prime is really infinity since crit_tot = 0
  elif len(no_critlist) == 0 and len(critlist) > 0:
    tau_2prime = np.random.exponential(1/crit_tot)
    tau_prime = tau_2prime + 1 ## not going to have a tau_prime so it will defautl to tau_2prime

  ## now do the actual tau selection
  if tau_prime < tau_2prime:
    tau = tau_prime
    ## we are not going to touch the criticals, i.e. k = 0  
    ## these are the noncritical k values
    k_noncrit = [np.random.poisson(i*tau,1) for i in no_critlist]   
    k_crit = [0 for i in critlist]
  else: 
    tau = tau_2prime
    ## number of times to fire noncritical
    k_noncrit = [np.random.poisson(i*tau,1) for i in no_critlist]   
    ## don't fire most criticals,
    k_crit = [0 for i in critlist]
    ## except, we will pick a single one of them to fire
    crit_event = rnd.choices(critlist, weights = critlist, k=1)[0]
    ## go along the list of critical events, and the one that equals the critical event we chose, it gets reassigned kcrit = 1
    for i in range(len(k_crit)):
      if critlist[i] == crit_event:
        k_crit[i] = 1

  ## fire the non_critical reactions
  for i in range(len(no_critlist_string)):
    q=0
    if k_noncrit[i]>0:
      while q < k_noncrit[i]:
        countInfected_H, countInfected_L = do_Tau_event(no_critlist_string[i],countInfected_H, countInfected_L, variant)
        q = q+1
        
  ## fire the critical reactions
  for i in range(len(critlist_string)):
    if k_crit[i]>0:
      q = 0
      while q < k_crit[i]:
        countInfected_H, countInfected_L = do_Tau_event(critlist_string[i],countInfected_H, countInfected_L, variant)
        q = q+1
  return(tau, countInfected_H, countInfected_L)

######### END choose event tau function #################

########## function to update the dataframe #########
def update_rows(): 
  rows.append([indexID, t, len(host_pop), len(Inf_L), len(Inf_H), len(Dead_L), len(Dead_H), countInfected_L, countInfected_H,
  (countInfected_L+countInfected_H)/(len(Pop_L)+len(Pop_H)), len(Exp_L), len(Exp_H), len(Sus_L), len(Sus_H), len(Rec_nat_L), len(Rec_nat_H), len(Rec_vax_L), len(Rec_vax_H), numVaccinatedH, numVaccinatedL, 
  numBoostedH, numBoostedL,
  count_Naive, count_Delta, count_Delta_vaccine, count_two, count_three, count_four, 
count_Omicron, count_Omicron_star, count_Omicron_star_vaccine, count_Omicron_vaccine, count_Vaccine,
count_WT, count_WT_vaccine, count_WT_vaccine_boost, count_Omicron_vaccine_boost,
count_Delta_vaccine_boost, count_Omicron_star_vaccine_boost, count_Vaccine_boost
  ])

############## Class Host - Low SES ###########################

class HostL(object): ## Host class
  def __init__(self, age):
    self.age = age
    self.immunityhist = []
    self.exposed= False
    self.infectious = False
    self.recovered = False

################ Class functions Low SES #####################

  def param_assign_L(self): ## grid_add
    ## naive
    if len(self.immunityhist) == 0:
      k = 1
    ## vaccine only
    elif self.immunityhist == ["Vaccine"]:
      k = k_Vaccine
    elif self.immunityhist == ["Vaccine", "Booster"]:
      k = k_Vaccine_boost 
    ## wild type
    elif self.immunityhist == ["WT"]:
      k = k_WT 
    elif (self.immunityhist == ["WT", "Vaccine"] or self.immunityhist == ["Vaccine", "WT"]):
      k = k_WT_vaccine
    elif (self.immunityhist == ["WT", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "WT", "Booster"]):
      k = k_WT_vaccine_boost
    ## delta
    elif self.immunityhist == ["Delta"]:
      k = k_Delta
    elif (self.immunityhist == ["Delta", "Vaccine"] or self.immunityhist == ["Vaccine", "Delta"]):
      k = k_Delta_vaccine
    elif (self.immunityhist == ["Delta", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Delta", "Booster"]):
      k = k_Delta_vaccine_boost
    ## omicron
    elif self.immunityhist == ["Omicron"]:
      k = k_Omicron
    elif (self.immunityhist == ["Omicron", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron"]):
      k = k_Omicron_vaccine
    elif (self.immunityhist == ["Omicron", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron"]):
      k = k_Omicron_vaccine_boost
    ## omicron*
    elif self.immunityhist == ["Omicron*"]:
      k = k_Omicron_star
    elif (self.immunityhist == ["Omicron*", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron*"]):
      k = k_Omicron_star_vaccine
    elif (self.immunityhist == ["Omicron*", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron*", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron*"]):
      k = k_Omicron_star_vaccine_boost
    ## if not included above and length < 4
    elif len(self.immunityhist) == 2:
      k = k_two
    elif len(self.immunityhist) == 3:
      k = k_three
    ## four and above histories
    else: 
      k = k_four
    return k

## infection
  def can_get_infected_L(self, mu): ## boolean answer. 
    answer = False
    ## if they already have a pathogen, can't get infected
    if self.exposed == True or self.infectious == True: 
      return answer
    else:
      k = self.param_assign_L()
      ## since this is for infections, do the reduction in protection (1-k) 
      k = 1-(inf_pen * (1-k))
      ## if they are recently recovered
      if self.recovered == True:
        if np.random.binomial(1, k*mu) == 1:
          answer = True
      ## if they are not recently recovered
      else:
        ## if no prior infection, they are naive so k is just 1 
        if len(self.immunityhist) == 0:
          if np.random.binomial(1, k*mu) == 1:
            answer = True
        ## if one prior infection, use the dampening1 parameter on their protection (1-k)
        elif len(self.immunityhist) == 1:
          k = 1-(dampening1*(1-k))
          if np.random.binomial(1, k*mu) == 1:
            answer = True
        ## if two prior infections, use the dampening2 parameter on their protection (1-k)
        elif len(self.immunityhist) == 2:
          k = 1-(dampening2*(1-k))
          if np.random.binomial(1, k*mu) == 1:
            answer = True
        ## if three or more infections, use dampening3
        else:
          k = 1-(dampening3*(1-k))
          if np.random.binomial(1, k*mu) == 1:
            answer = True
      return answer

## get exposed
  def get_exposed_L(self):
    self.exposed = True
    if self in Rec_nat_L:
      Rec_nat_L.remove(self)
    elif self in Rec_vax_L:
      Rec_vax_L.remove(self)
    else: #Sus_L:
      Sus_L.remove(self)
    Exp_L.append(self)

## contact event
  def contact_event_L(self):
    if self.age == "child":
      if self.can_get_infected_L(mu_child): 
        self.get_exposed_L()
    elif self.age == "adult":
      if self.can_get_infected_L(mu_adult): 
        self.get_exposed_L()
    else: 
      if self.can_get_infected_L(mu_old): 
        self.get_exposed_L()

## exposed becomes infectious event
  def can_infect_others_L(self):
    self.infectious = True
    self.exposed = False
    Exp_L.remove(self)
    Inf_L.append(self)
    ## if they are someone without a vaccine, get them out of the vax-eligible class while they are infected
    if 'Vaccine' not in self.immunityhist:
      unVaxLow.remove(self)
    ## if they are someone without a booster, get them out of the booster-eligible class while they are infected
    if self in unBoostLow:
      unBoostLow.remove(self)

## function that kills the host
  def death_L(self):
    self.subtract_histories_L()
    host_pop.remove(self)
    Inf_L.remove(self)
    Pop_L.remove(self)
    Dead_L.append(self)

  def add_histories_L(self):
    ## no boosters 
    global count_Naive
    global count_Delta
    global count_Delta_vaccine
    global count_two
    global count_three
    global count_four
    global count_Omicron
    global count_Omicron_star
    global count_Omicron_star_vaccine
    global count_Omicron_vaccine
    global count_Vaccine
    global count_WT
    global count_WT_vaccine 
    ## with boosters 
    global count_WT_vaccine_boost
    global count_Omicron_vaccine_boost
    global count_Delta_vaccine_boost
    global count_Omicron_star_vaccine_boost
    global count_Vaccine_boost
    if len(self.immunityhist) == 0:
      count_Naive += 1
    ## vaccine only
    elif self.immunityhist == ["Vaccine"]:
      count_Vaccine += 1
    elif self.immunityhist == ["Vaccine", "Booster"]:
      count_Vaccine_boost += 1
    ## wild type
    elif self.immunityhist == ["WT"]:
      count_WT += 1
    elif (self.immunityhist == ["WT", "Vaccine"] or self.immunityhist == ["Vaccine", "WT"]):
      count_WT_vaccine += 1
    elif (self.immunityhist == ["WT", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "WT", "Booster"]):
      count_WT_vaccine_boost += 1
    ## delta
    elif self.immunityhist == ["Delta"]:
      count_Delta += 1
    elif (self.immunityhist == ["Delta", "Vaccine"] or self.immunityhist == ["Vaccine", "Delta"]):
      count_Delta_vaccine += 1
    elif (self.immunityhist == ["Delta", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Delta", "Booster"]):
      count_Delta_vaccine_boost  += 1
    ## omicron
    elif self.immunityhist == ["Omicron"]:
      count_Omicron  += 1
    elif (self.immunityhist == ["Omicron", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron"]):
      count_Omicron_vaccine  += 1
    elif (self.immunityhist == ["Omicron", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron", "Booster"] or self.immunityhist == ["Vaccine","Booster","Omicron"]):
      count_Omicron_vaccine_boost  += 1
    ## omicron*
    elif self.immunityhist == ["Omicron*"]:
      count_Omicron_star += 1
    elif (self.immunityhist == ["Omicron*", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron*"]):
      count_Omicron_star_vaccine += 1
    elif (self.immunityhist == ["Omicron*", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron*", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron*"]):
      count_Omicron_star_vaccine_boost +=1
    ## if not included above and length < 4
    elif len(self.immunityhist) == 2:
      count_two  += 1
    elif len(self.immunityhist) == 3:
      count_three  += 1
    ## four and above histories
    else: 
      count_four  += 1

  def subtract_histories_L(self):
    ## no boosters 
    global count_Naive
    global count_Delta
    global count_Delta_vaccine
    global count_two
    global count_three
    global count_four
    global count_Omicron
    global count_Omicron_star
    global count_Omicron_star_vaccine
    global count_Omicron_vaccine
    global count_Vaccine
    global count_WT
    global count_WT_vaccine 
    ## with boosters 
    global count_WT_vaccine_boost
    global count_Omicron_vaccine_boost
    global count_Delta_vaccine_boost
    global count_Omicron_star_vaccine_boost
    global count_Vaccine_boost
    if len(self.immunityhist) == 0:
      count_Naive -= 1
    ## vaccine only
    elif self.immunityhist == ["Vaccine"]:
      count_Vaccine -= 1
    elif self.immunityhist == ["Vaccine", "Booster"]:
      count_Vaccine_boost -= 1
    ## wild type
    elif self.immunityhist == ["WT"]:
      count_WT -= 1
    elif (self.immunityhist == ["WT", "Vaccine"] or self.immunityhist == ["Vaccine", "WT"]):
      count_WT_vaccine -= 1
    elif (self.immunityhist == ["WT", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "WT", "Booster"]):
      count_WT_vaccine_boost -= 1
    ## delta
    elif self.immunityhist == ["Delta"]:
      count_Delta -= 1
    elif (self.immunityhist == ["Delta", "Vaccine"] or self.immunityhist == ["Vaccine", "Delta"]):
      count_Delta_vaccine -= 1
    elif (self.immunityhist == ["Delta", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Delta", "Booster"]):
      count_Delta_vaccine_boost  -= 1
    ## omicron
    elif self.immunityhist == ["Omicron"]:
      count_Omicron  -= 1
    elif (self.immunityhist == ["Omicron", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron"]):
      count_Omicron_vaccine  -= 1
    elif (self.immunityhist == ["Omicron", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron"]):
      count_Omicron_vaccine_boost  -= 1
    ## omicron*
    elif self.immunityhist == ["Omicron*"]:
      count_Omicron_star -= 1
    elif (self.immunityhist == ["Omicron*", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron*"]):
      count_Omicron_star_vaccine -= 1
    elif (self.immunityhist == ["Omicron*", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron*", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron*"]):
      count_Omicron_star_vaccine_boost -=1
    ## if not included above and length < 4
    elif len(self.immunityhist) == 2:
      count_two  -= 1
    elif len(self.immunityhist) == 3:
      count_three  -= 1
    ## four and above histories
    else: 
      count_four  -= 1

## function that recovers the host
  def recover_L(self, variant):
    self.recovered = True
    Rec_nat_L.append(self)
    self.infectious = False
    Inf_L.remove(self)
    self.subtract_histories_L() ## take them off their old history count
    self.immunityhist.append(variant) ## immune history of having a variant
    self.add_histories_L() ## add them to new history count
    self.infecting_pathogen = None
    ## if they are someone without a vaccine, add them back to eligibility
    if 'Vaccine' not in self.immunityhist:
      unVaxLow.append(self)
    ## if they are someone without a booster, add them back to eligibility
    if 'Booster' not in self.immunityhist and 'Vaccine' in self.immunityhist and self.age != "child":
      unBoostLow.append(self)

  ## death or recovery event
  def death_or_recover_L(self, alpha, variant):
      k = self.param_assign_L()
      ## since this is for deaths, use the base parameter 
      ## if they are recently recovered
      if self.recovered == True:
        if np.random.binomial(1, k*alpha) == 1: self.death_L()
        else: self.recover_L(variant)
      ## if they are not recently recovered
      else:
        ## if no prior infection, they are naive so k is just 1 
        if len(self.immunityhist) == 0:
          if np.random.binomial(1, k*alpha) == 1: self.death_L()
          else: self.recover_L(variant)
        ## if one prior infection, use the dampening1 parameter on their protection (1-k)
        elif len(self.immunityhist) == 1:
          k = 1-(dampening1*(1-k))
          if np.random.binomial(1, k*alpha) == 1: self.death_L()
          else: self.recover_L(variant)
        ## if two prior infections, use the dampening2 parameter on their protection (1-k)
        elif len(self.immunityhist) == 2:
          k = 1-(dampening2*(1-k))
          if np.random.binomial(1, k*alpha) == 1: self.death_L()
          else: self.recover_L(variant)
        ## if three or more infections, use dampening3
        else:
          k = 1-(dampening3*(1-k))
          if np.random.binomial(1, k*alpha) == 1: self.death_L()
          else: self.recover_L(variant)

  ## waning immunity event
  def has_waned_vax_L(self):
    Rec_vax_L.remove(self)
    Sus_L.append(self)
    self.recovered = False
  
  def has_waned_nat_L(self):
    Rec_nat_L.remove(self)
    Sus_L.append(self)
    self.recovered = False
    
  ## vaccination event 
  def vaccinate_L(self):
    self.subtract_histories_L()
    self.immunityhist.append('Vaccine')
    self.add_histories_L()
    self.recovered = True
    ## they can't already be in Rec Vax since this is their first vaccine event
    Rec_vax_L.append(self)
    if self in Exp_L:
      Exp_L.remove(self)
      self.exposed = False
    elif self in Sus_L:
      Sus_L.remove(self)
    elif self in Rec_nat_L:
      Rec_nat_L.remove(self)
    unVaxLow.remove(self)
    if self.age == "adult" or self.age == "old":
    	unBoostLow.append(self)

  def boost_L(self):
    self.subtract_histories_L()
    self.immunityhist.append('Booster')
    self.add_histories_L()
    self.recovered = True
    if self in Exp_L:
      Exp_L.remove(self)
      self.exposed = False
    elif self in Sus_L:
      Sus_L.remove(self)
    elif self in Rec_nat_L:
      Rec_nat_L.remove(self)
    elif self in Rec_vax_L:
      Rec_vax_L.remove(self)
    Rec_vax_L.append(self)
    unBoostLow.remove(self)

############## Class Host - High SES ###########################
class HostH(object): ## host class
  def __init__(self, age):
    self.age = age
    self.immunityhist = []
    self.exposed = False
    self.infectious = False
    self.recovered = False

############## Class functions - High SES ###########################
  def param_assign_H(self): ## grid_add
    ## naive
    if len(self.immunityhist) == 0:
      k = 1
    ## vaccine only
    elif self.immunityhist == ["Vaccine"]:
      k = k_Vaccine
    elif self.immunityhist == ["Vaccine", "Booster"]:
      k = k_Vaccine_boost 
    ## wild type
    elif self.immunityhist == ["WT"]:
      k = k_WT 
    elif (self.immunityhist == ["WT", "Vaccine"] or self.immunityhist == ["Vaccine", "WT"]):
      k = k_WT_vaccine
    elif (self.immunityhist == ["WT", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "WT", "Booster"]):
      k = k_WT_vaccine_boost
    ## delta
    elif self.immunityhist == ["Delta"]:
      k = k_Delta
    elif (self.immunityhist == ["Delta", "Vaccine"] or self.immunityhist == ["Vaccine", "Delta"]):
      k = k_Delta_vaccine
    elif (self.immunityhist == ["Delta", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Delta", "Booster"]):
      k = k_Delta_vaccine_boost
    ## omicron
    elif self.immunityhist == ["Omicron"]:
      k = k_Omicron
    elif (self.immunityhist == ["Omicron", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron"]):
      k = k_Omicron_vaccine
    elif (self.immunityhist == ["Omicron", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron"]):
      k = k_Omicron_vaccine_boost
    ## omicron*
    elif self.immunityhist == ["Omicron*"]:
      k = k_Omicron_star
    elif (self.immunityhist == ["Omicron*", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron*"]):
      k = k_Omicron_star_vaccine
    elif (self.immunityhist == ["Omicron*", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron*", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron*"]):
      k = k_Omicron_star_vaccine_boost
    ## if not included above and length < 4
    elif len(self.immunityhist) == 2:
      k = k_two
    elif len(self.immunityhist) == 3:
      k = k_three
    ## four and above histories
    else: 
      k = k_four
    return k

## infection
  def can_get_infected_H(self, mu): ## boolean answer. 
    answer = False
    ## if they already have a pathogen, can't get infected
    if self.exposed == True or self.infectious == True: 
      return answer
    else:
      k = self.param_assign_H()
      ## since this is for infections, do the reduction in protection (1-k) 
      k = 1-(inf_pen * (1-k))
      ## if they are recently recovered
      if self.recovered == True:
        if np.random.binomial(1, k*mu) == 1:
          answer = True
      ## if they are not recently recovered
      else:
        ## if no prior infection, they are naive so k is just 1 
        if len(self.immunityhist) == 0:
          if np.random.binomial(1, k*mu) == 1:
            answer = True
        ## if one prior infection, use the dampening1 parameter on their protection (1-k)
        elif len(self.immunityhist) == 1:
          k = 1-(dampening1*(1-k))
          if np.random.binomial(1, k*mu) == 1:
            answer = True
        ## if two prior infections, use the dampening2 parameter on their protection (1-k)
        elif len(self.immunityhist) == 2:
          k = 1-(dampening2*(1-k))
          if np.random.binomial(1, k*mu) == 1:
            answer = True
        ## if three or more infections, use dampening3
        else:
          k = 1-(dampening3*(1-k))
          if np.random.binomial(1, k*mu) == 1:
            answer = True
      return answer

## get exposed
  def get_exposed_H(self):
    self.exposed = True
    if self in Rec_nat_H:
      Rec_nat_H.remove(self)
    elif self in Rec_vax_H:
      Rec_vax_H.remove(self)
    else: #Sus_H:
      Sus_H.remove(self)
    Exp_H.append(self)

## contact event
  def contact_event_H(self):
    if self.age == "child":
      if self.can_get_infected_H(mu_child): 
        self.get_exposed_H()
    elif self.age == "adult":
      if self.can_get_infected_H(mu_adult): 
        self.get_exposed_H()
    else: 
      if self.can_get_infected_H(mu_old): 
        self.get_exposed_H()

## exposed becomes infectious event
  def can_infect_others_H(self):
    self.infectious = True
    self.exposed = False
    Exp_H.remove(self)
    Inf_H.append(self)
    ## if they are someone without a vaccine, get them out of the booster-eligible class while they are infected
    if 'Vaccine' not in self.immunityhist:
      unVaxHigh.remove(self)
    ## if they are someone without a booster, get them out of the booster-eligible class while they are infected
    if self in unBoostHigh:
      unBoostHigh.remove(self)

## function that kills the host
  def death_H(self):
    self.subtract_histories_H()
    host_pop.remove(self)
    Inf_H.remove(self)
    Pop_H.remove(self)
    Dead_H.append(self)
  
## counting histories
  def add_histories_H(self):
    ## no boosters 
    global count_Naive
    global count_Delta
    global count_Delta_vaccine
    global count_two
    global count_three
    global count_four
    global count_Omicron
    global count_Omicron_star
    global count_Omicron_star_vaccine
    global count_Omicron_vaccine
    global count_Vaccine
    global count_WT
    global count_WT_vaccine 
    ## with boosters 
    global count_WT_vaccine_boost
    global count_Omicron_vaccine_boost
    global count_Delta_vaccine_boost
    global count_Omicron_star_vaccine_boost
    global count_Vaccine_boost
    if len(self.immunityhist) == 0:
      count_Naive += 1
    ## vaccine only
    elif self.immunityhist == ["Vaccine"]:
      count_Vaccine += 1
    elif self.immunityhist == ["Vaccine", "Booster"]:
      count_Vaccine_boost += 1
    ## wild type
    elif self.immunityhist == ["WT"]:
      count_WT += 1
    elif (self.immunityhist == ["WT", "Vaccine"] or self.immunityhist == ["Vaccine", "WT"]):
      count_WT_vaccine += 1
    elif (self.immunityhist == ["WT", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "WT", "Booster"]):
      count_WT_vaccine_boost += 1
    ## delta
    elif self.immunityhist == ["Delta"]:
      count_Delta += 1
    elif (self.immunityhist == ["Delta", "Vaccine"] or self.immunityhist == ["Vaccine", "Delta"]):
      count_Delta_vaccine += 1
    elif (self.immunityhist == ["Delta", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Delta", "Booster"]):
      count_Delta_vaccine_boost  += 1
    ## omicron
    elif self.immunityhist == ["Omicron"]:
      count_Omicron  += 1
    elif (self.immunityhist == ["Omicron", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron"]):
      count_Omicron_vaccine  += 1
    elif (self.immunityhist == ["Omicron", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron"]):
      count_Omicron_vaccine_boost  += 1
    ## omicron*
    elif self.immunityhist == ["Omicron*"]:
      count_Omicron_star += 1
    elif (self.immunityhist == ["Omicron*", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron*"]):
      count_Omicron_star_vaccine += 1
    elif (self.immunityhist == ["Omicron*", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron*", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron*"]):
      count_Omicron_star_vaccine_boost +=1
    ## if not included above and length < 4
    elif len(self.immunityhist) == 2:
      count_two  += 1
    elif len(self.immunityhist) == 3:
      count_three  += 1
    ## four and above histories
    else: 
      count_four  += 1

  def subtract_histories_H(self):
    ## no boosters 
    global count_Naive
    global count_Delta
    global count_Delta_vaccine
    global count_two
    global count_three
    global count_four
    global count_Omicron
    global count_Omicron_star
    global count_Omicron_star_vaccine
    global count_Omicron_vaccine
    global count_Vaccine
    global count_WT
    global count_WT_vaccine 
    ## with boosters 
    global count_WT_vaccine_boost
    global count_Omicron_vaccine_boost
    global count_Delta_vaccine_boost
    global count_Omicron_star_vaccine_boost
    global count_Vaccine_boost
    if len(self.immunityhist) == 0:
      count_Naive -= 1
    ## vaccine only
    elif self.immunityhist == ["Vaccine"]:
      count_Vaccine -= 1
    elif self.immunityhist == ["Vaccine", "Booster"]:
      count_Vaccine_boost -= 1
    ## wild type
    elif self.immunityhist == ["WT"]:
      count_WT -= 1
    elif (self.immunityhist == ["WT", "Vaccine"] or self.immunityhist == ["Vaccine", "WT"]):
      count_WT_vaccine -= 1
    elif (self.immunityhist == ["WT", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "WT", "Booster"]):
      count_WT_vaccine_boost -= 1
    ## delta
    elif self.immunityhist == ["Delta"]:
      count_Delta -= 1
    elif (self.immunityhist == ["Delta", "Vaccine"] or self.immunityhist == ["Vaccine", "Delta"]):
      count_Delta_vaccine -= 1
    elif (self.immunityhist == ["Delta", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Delta", "Booster"]):
      count_Delta_vaccine_boost  -= 1
    ## omicron
    elif self.immunityhist == ["Omicron"]:
      count_Omicron  -= 1
    elif (self.immunityhist == ["Omicron", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron"]):
      count_Omicron_vaccine  -= 1
    elif (self.immunityhist == ["Omicron", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron"]):
      count_Omicron_vaccine_boost  -= 1
    ## omicron*
    elif self.immunityhist == ["Omicron*"]:
      count_Omicron_star -= 1
    elif (self.immunityhist == ["Omicron*", "Vaccine"] or self.immunityhist == ["Vaccine", "Omicron*"]):
      count_Omicron_star_vaccine -= 1
    elif (self.immunityhist == ["Omicron*", "Vaccine", "Booster"] or self.immunityhist == ["Vaccine", "Omicron*", "Booster"] or self.immunityhist == ["Vaccine", "Booster", "Omicron*"]):
      count_Omicron_star_vaccine_boost -=1
    ## if not included above and length < 4
    elif len(self.immunityhist) == 2:
      count_two  -= 1
    elif len(self.immunityhist) == 3:
      count_three  -= 1
    ## four and above histories
    else: 
      count_four  -= 1

## function that recovers the host
  def recover_H(self, variant):
    self.recovered = True
    Rec_nat_H.append(self)
    self.infectious = False
    Inf_H.remove(self)
    self.subtract_histories_H() ## take them off their old history count
    self.immunityhist.append(variant) ## immune history of having a variant
    self.add_histories_H() ## add them to new history count
    self.infecting_pathogen = None
    ## if they are someone without a vaccine, add them back to eligibility
    if 'Vaccine' not in self.immunityhist:
      unVaxHigh.append(self)
    ## if they are someone without a booster, and not a child, add them back to eligibility
    if 'Booster' not in self.immunityhist and 'Vaccine' in self.immunityhist and self.age != "child":
      unBoostHigh.append(self)

  ## death or recovery event
  def death_or_recover_H(self, alpha, variant):
      k = self.param_assign_H()
      ## since this is for deaths, use the base parameter 
      ## if they are recently recovered
      if self.recovered == True:
        if np.random.binomial(1, k*alpha) == 1: self.death_H()
        else: self.recover_H(variant)
      ## if they are not recently recovered
      else:
        ## if no prior infection, they are naive so k is just 1 
        if len(self.immunityhist) == 0:
          if np.random.binomial(1, k*alpha) == 1: self.death_H()
          else: self.recover_H(variant)
        ## if one prior infection, use the dampening1 parameter on their protection (1-k)
        elif len(self.immunityhist) == 1:
          k = 1-(dampening1*(1-k))
          if np.random.binomial(1, k*alpha) == 1: self.death_H()
          else: self.recover_H(variant)
        ## if two prior infections, use the dampening2 parameter on their protection (1-k)
        elif len(self.immunityhist) == 2:
          k = 1-(dampening2*(1-k))
          if np.random.binomial(1, k*alpha) == 1: self.death_H()
          else: self.recover_H(variant)
        ## if three or more infections, use dampening3
        else:
          k = 1-(dampening3*(1-k))
          if np.random.binomial(1, k*alpha) == 1: self.death_H()
          else: self.recover_H(variant)

  ## waning immunity event
  def has_waned_vax_H(self):
    Rec_vax_H.remove(self)
    Sus_H.append(self)
    self.recovered = False
  
  def has_waned_nat_H(self):
    Rec_nat_H.remove(self)
    Sus_H.append(self)
    self.recovered = False
    
  ## vaccination event 
  def vaccinate_H(self):
    self.subtract_histories_H()
    self.immunityhist.append('Vaccine')
    self.add_histories_H()
    self.recovered = True
    ## they can't already be in Rec Vax since this is their first vaccine event
    Rec_vax_H.append(self)
    if self in Exp_H:
      Exp_H.remove(self)
      self.exposed = False
    elif self in Sus_H:
      Sus_H.remove(self)
    elif self in Rec_nat_H:
      Rec_nat_H.remove(self)
    unVaxHigh.remove(self)
    if self.age == "adult" or self.age == "old":
    	unBoostHigh.append(self)

  def boost_H(self):
    self.subtract_histories_H()
    self.immunityhist.append('Booster')
    self.add_histories_H()
    self.recovered = True
    if self in Exp_H:
      Exp_H.remove(self)
      self.exposed = False
    elif self in Sus_H:
      Sus_H.remove(self)
    elif self in Rec_nat_H:
      Rec_nat_H.remove(self)
    elif self in Rec_vax_H:
      Rec_vax_H.remove(self)
    Rec_vax_H.append(self)
    unBoostHigh.remove(self)

########## Initialize model ##########
########## Start rows ########## 
rows = []
########## Set host population ########## 

host_Lchild = [HostL("child") for i in range(N_child)]
host_Ladult = [HostL("adult") for i in range(N_adult)]
host_Lold = [HostL("old") for i in range(N_old)]
Pop_L = host_Lchild.copy() + host_Ladult.copy() + host_Lold.copy()

host_Hchild = [HostH("child") for i in range(N_child)]
host_Hadult = [HostH("adult") for i in range(N_adult)]
host_Hold = [HostH("old") for i in range(N_old)]
Pop_H = host_Hchild.copy() + host_Hadult.copy() + host_Hold.copy()

host_pop = Pop_H.copy() + Pop_L.copy()

t = 0.0
variant = 'WT'

############ Split into Classes ################
Sus_L = Pop_L.copy()
Sus_H = Pop_H.copy()
Exp_L = []
Exp_H = []
Inf_L = []
Inf_H = []
Rec_nat_L = []
Rec_nat_H = []
Rec_vax_L = []
Rec_vax_H = []
Dead_L = []
Dead_H = []
unVaxLow = Pop_L.copy()
unVaxHigh = Pop_H.copy()
unBoostLow = []
unBoostHigh = []

# ## add exposed:
for j in range(25):
  h = rnd.choice(Sus_L + Rec_nat_L + Rec_vax_L)
  h.get_exposed_L()
  g = rnd.choice(Sus_H + Rec_nat_H + Rec_vax_H)
  g.get_exposed_H()

################ STAGE 1 - WT #################
## run for 420 days of WT with vaccination starting at t = 320
print('WT')
print('--- %s seconds ---' % (time.time() - start_time))
countT = 1

while t < 420.0: 
  time_added, countInfected_H, countInfected_L = choose_event_tau(c_LL, c_LH, c_HL, c_HH, epsilon, gamma,
len(Sus_L), len(Sus_H), len(Exp_L), len(Exp_H), len(Inf_L), len(Inf_H), b, w, len(Rec_nat_L), len(Rec_nat_H), len(Rec_vax_L), len(Rec_vax_H),
countInfected_H, countInfected_L, variant)
  if(t > countT):
  ##  vaccinate
    if t>320:
      to_vax_H = vaxTrend(Vm_H, Wh_H, h_H, countT-320, 1, len(Pop_H))
      to_vax_L = vaxTrend(Vm_L, Wh_L, h_L, countT-320, 1, len(Pop_L))
      for j in range(to_vax_H):
        if len(unVaxHigh)>0:
          h = rnd.choice(unVaxHigh)
          h.vaccinate_H()
          numVaccinatedH = numVaccinatedH + 1
      for j in range(to_vax_L):
        if len(unVaxLow)>0:
          h = rnd.choice(unVaxLow)
          h.vaccinate_L()
          numVaccinatedL = numVaccinatedL + 1
    update_rows()
    countT = countT + 1
  t+=time_added

################ STAGE 2 - DELTA AND VACCINATION #################
## this is the current variant
variant = 'Delta'
mu_child = 0.297 ## SAR for Delta
mu_adult = 0.297 ## SAR for Delta
mu_old = 0.297 ## SAR for Delta

######## Serotype parameters during Delta covid (these will be applied to severe disease, infection will be modified) ##########
k_Delta = 0.33
k_Delta_vaccine = 0.30

k_two = 0.33
k_three = 0.19
k_four = 0.050

k_Omicron = 0.95
k_Omicron_star = 1
k_Omicron_star_vaccine = 0.32
k_Omicron_vaccine = 0.22
k_Vaccine = 0.75
k_WT = 0.62
k_WT_vaccine = 0.48

## stringency
stringency = grid["string_delta"][indexID-1]
## multiply all contact by stringency
c_LL = b_LL*stringency  ## low contacts low
c_LH = b_LH*stringency ## low contacts high
c_HL = b_HL*stringency ## high contacts low
c_HH = b_HH*stringency ## high contacts high  

# ## add pathogens:
for j in range(25):
  h = rnd.choice(Sus_L + Rec_nat_L + Rec_vax_L)
  h.get_exposed_L()
  g = rnd.choice(Sus_H + Rec_nat_H + Rec_vax_H)
  g.get_exposed_H()

## run for another 7 months with vaccination 
print('Delta')
print('--- %s seconds ---' % (time.time() - start_time))
while t < 630.0:
  time_added, countInfected_H, countInfected_L = choose_event_tau(c_LL, c_LH, c_HL, c_HH, epsilon, gamma,
len(Sus_L), len(Sus_H), len(Exp_L), len(Exp_H), len(Inf_L), len(Inf_H), b, w, len(Rec_nat_L), len(Rec_nat_H), len(Rec_vax_L), len(Rec_vax_H),
countInfected_H, countInfected_L, variant)
  if(t > countT):
    to_vax_H = vaxTrend(Vm_H, Wh_H, h_H, countT-320, 1, len(Pop_H))
    to_vax_L = vaxTrend(Vm_L, Wh_L, h_L, countT-320, 1, len(Pop_L))
    for j in range(to_vax_H):
      if len(unVaxHigh)>0:
        h = rnd.choice(unVaxHigh)
        h.vaccinate_H()
        numVaccinatedH = numVaccinatedH + 1
    for j in range(to_vax_L):
      if len(unVaxLow)>0:
        h = rnd.choice(unVaxLow)
        h.vaccinate_L()
        numVaccinatedL = numVaccinatedL + 1
    update_rows()
    countT = countT + 1
  t+=time_added

################ STAGE 3 - OMICRON VARIANT #################
## this is the current variant
variant = 'Omicron'

## Update IFR

alphaL1_child =   alphaL1_child*0.5 ## death rate
alphaL1_adult =   alphaL1_adult*0.5 ## death rate
alphaL1_old =   alphaL1_old*0.5 ## death rate
alphaH1_child =   alphaH1_child*0.5 ## death rate
alphaH1_adult =   alphaH1_adult*0.5 ## death rate
alphaH1_old =   alphaH1_old*0.5 ## death rate

## update secondary attack rate
mu_child = 0.427 ## SAR for Omicron
mu_adult = 0.427 ## SAR for Omicron
mu_old = 0.427 ## SAR for Omicron

######## Serotype parameters during Omicron covid (these will be applied to severe disease, infection will be modified) ##########
k_Delta = 0.87
k_Delta_vaccine = 0.62

k_two = 0.60
k_three = 0.32
k_four = 0.050

k_Omicron = 0.4
k_Omicron_star = 0.5
k_Omicron_star_vaccine = 0.56
k_Omicron_vaccine = 0.46
k_Vaccine = 0.83
k_WT = 0.92
k_WT_vaccine = 0.75

## with boosters 
k_WT_vaccine_boost = grid['WT_vaccine_Om'][indexID-1]
k_Omicron_vaccine_boost = grid['Omicron_vaccine_Om'][indexID-1]
k_Delta_vaccine_boost = grid['Delta_vaccine_Om'][indexID-1]
k_Omicron_star_vaccine_boost = grid['Omicron_star_vaccine_Om'][indexID-1]
k_Vaccine_boost = grid['Vaccine_Om'][indexID-1]

## stringency
stringency = grid["string_omicron"][indexID-1]
## multiply all contact by stringency
c_LL = b_LL*stringency  ## low contacts low
c_LH = b_LH*stringency ## low contacts high
c_HL = b_HL*stringency ## high contacts low
c_HH = b_HH*stringency ## high contacts high  

# ## add pathogens:
for j in range(25):
  h = rnd.choice(Sus_L + Rec_nat_L + Rec_vax_L)
  h.get_exposed_L()
  g = rnd.choice(Sus_H + Rec_nat_H + Rec_vax_H)
  g.get_exposed_H()

## run for 7 months with BA1
print('Omicron')
print('--- %s seconds ---' % (time.time() - start_time))
while t < 900.0:
  time_added, countInfected_H, countInfected_L = choose_event_tau(c_LL, c_LH, c_HL, c_HH, epsilon, gamma,
len(Sus_L), len(Sus_H), len(Exp_L), len(Exp_H), len(Inf_L), len(Inf_H), b, w, len(Rec_nat_L), len(Rec_nat_H), len(Rec_vax_L), len(Rec_vax_H),
countInfected_H, countInfected_L, variant)
  if(t > countT):
## vaccines
    to_vax_H = vaxTrend(Vm_H, Wh_H, h_H, countT-320, 1, len(Pop_H))
    to_vax_L = vaxTrend(Vm_L, Wh_L, h_L, countT-320, 1, len(Pop_L))
    for j in range(to_vax_H):
      if len(unVaxHigh)>0:
        h = rnd.choice(unVaxHigh)
        h.vaccinate_H()
        numVaccinatedH = numVaccinatedH + 1
    for j in range(to_vax_L):
      if len(unVaxLow)>0:
        h = rnd.choice(unVaxLow)
        h.vaccinate_L()
        numVaccinatedL = numVaccinatedL + 1
## boosters
    if boost == "yes" and t > 840.0:
      if boost_shape == "functional":
        to_boost_H = vaxTrend(Vm_H_boost, Wh_H_boost, kH_boost, countT-840, 1, nonchild*numVaccinatedH)
        to_boost_L = vaxTrend(Vm_L_boost, Wh_L_boost, kL_boost, countT-840, 1, nonchild*numVaccinatedL)
      else:
        to_boost_H = vaxConstant(vax_constant_H, Vm_H_boost, numBoostedH, nonchild*numVaccinatedH)
        to_boost_L = vaxConstant(vax_constant_L, Vm_L_boost, numBoostedL, nonchild*numVaccinatedL)
        # print(len(unBoostHigh), len(unBoostLow))
      # print(len(unBoostHigh), len(unBoostLow))
      for j in range(to_boost_H):
        ## pick a random choice from the unvaccinated host population
        if len(unBoostHigh)>0:
          h = rnd.choice(unBoostHigh)
          h.boost_H()
          numBoostedH = numBoostedH + 1
    ## boost low SES  (starting 840 days after time 0)
      for j in range(to_boost_L):
        ## pick a random choice from the unvaccinated host population
        if len(unBoostLow)>0:
          h = rnd.choice(unBoostLow)
          h.boost_L()
          numBoostedL = numBoostedL + 1
    update_rows()
    countT = countT + 1
  t+=time_added ## add time

################ STAGE 4 - Omicron* #################
## this is the current variant - a new Omicron lineage
variant = 'Omicron*'
## update secondary attack rate
mu_child = grid['SAR'][indexID-1]
mu_adult = grid['SAR'][indexID-1]
mu_old = grid['SAR'][indexID-1]
######## Serotype parameters during Omicron_star covid (these will be applied to severe disease, infection will be modified) ##########
k_Delta = 0.97
k_Delta_vaccine = 0.72

k_two = 0.65
k_three = 0.35
k_four = 0.050

k_Omicron = 0.5
k_Omicron_star = 0.4
k_Omicron_star_vaccine = 0.46
k_Omicron_vaccine = 0.56
k_Vaccine = 0.93
k_WT = 1
k_WT_vaccine = 0.85

## boosters
k_WT_vaccine_boost = grid['WT_vaccine'][indexID-1]
k_Omicron_vaccine_boost = grid['Omicron_vaccine'][indexID-1]
k_Delta_vaccine_boost = grid['Delta_vaccine'][indexID-1]
k_Omicron_star_vaccine_boost = grid['Omicron_star_vaccine'][indexID-1]
k_Vaccine_boost = grid['Vaccine'][indexID-1]

## stringency
stringency = grid["string_om2"][indexID-1]
## multiply all contact by stringency
c_LL = b_LL*stringency  ## low contacts low
c_LH = b_LH*stringency ## low contacts high
c_HL = b_HL*stringency ## high contacts low
c_HH = b_HH*stringency ## high contacts high  

# ## add exposed:
for j in range(25):
  h = rnd.choice(Sus_L + Rec_nat_L + Rec_vax_L)
  h.get_exposed_L()
  g = rnd.choice(Sus_H + Rec_nat_H + Rec_vax_H)
  g.get_exposed_H()

## NEW VARIANT SCENARIO - 195 days
print('Omicron*')
print('--- %s seconds ---' % (time.time() - start_time))
update_rows()
while t < 1095.0:
## call Gillespie algorithm
  time_added, countInfected_H, countInfected_L = choose_event_tau(c_LL, c_LH, c_HL, c_HH, epsilon, gamma,
len(Sus_L), len(Sus_H), len(Exp_L), len(Exp_H), len(Inf_L), len(Inf_H), b, w, len(Rec_nat_L), len(Rec_nat_H), len(Rec_vax_L), len(Rec_vax_H),
countInfected_H, countInfected_L, variant)
## boost high SES (starting 840 days after time 0)
  ## for j from 1 to the number of people we should vaccinate at this step
  if(t > countT):
    # vaccines - ONLY IN THE SUPPLEMENT SCENARIO
    if vaccine_star == "yes":
      to_vax_H = vaxTrend(Vm_H, Wh_H, h_H, countT-320, 1, len(Pop_H))
      to_vax_L = vaxTrend(Vm_L, Wh_L, h_L, countT-320, 1, len(Pop_L))
      for j in range(to_vax_H):
        if len(unVaxHigh)>0:
          h = rnd.choice(unVaxHigh)
          h.vaccinate_H()
          numVaccinatedH = numVaccinatedH + 1
      for j in range(to_vax_L):
        if len(unVaxLow)>0:
          h = rnd.choice(unVaxLow)
          h.vaccinate_L()
          numVaccinatedL = numVaccinatedL + 1
## boosters
    if boost == "yes":
      if boost_shape == "functional":
        to_boost_H = vaxTrend(Vm_H_boost, Wh_H_boost, kH_boost, countT-840, 1, nonchild*numVaccinatedH)
        to_boost_L = vaxTrend(Vm_L_boost, Wh_L_boost, kL_boost, countT-840, 1, nonchild*numVaccinatedL)
      else:
        to_boost_H = vaxConstant(vax_constant_H, Vm_H_boost, numBoostedH, nonchild*numVaccinatedH)
        to_boost_L = vaxConstant(vax_constant_L, Vm_L_boost, numBoostedL, nonchild*numVaccinatedL)
        # print(len(unBoostHigh), len(unBoostLow))
      for j in range(to_boost_H):
        ## pick a random choice from the unvaccinated host population
        if len(unBoostHigh)>0:
          h = rnd.choice(unBoostHigh)
          h.boost_H()
          numBoostedH = numBoostedH + 1
    ## boost low SES  (starting 840 days after time 0)
      for j in range(to_boost_L):
        ## pick a random choice from the unvaccinated host population
        if len(unBoostLow)>0:
          h = rnd.choice(unBoostLow)
          h.boost_L()
          numBoostedL = numBoostedL + 1
    update_rows()
    countT = countT + 1
  t+=time_added

testing = count_Naive + count_Delta + count_Delta_vaccine + count_two + count_three + count_four + count_Omicron + count_Omicron_star + count_Omicron_star_vaccine + count_Omicron_vaccine + count_Vaccine + count_WT + count_WT_vaccine + count_WT_vaccine_boost + count_Omicron_vaccine_boost + count_Delta_vaccine_boost + count_Omicron_star_vaccine_boost + count_Vaccine_boost
classes = len(Sus_L) + len(Sus_H) + len(Exp_H) + len(Exp_L) + len(Inf_H) + len(Inf_L) + len(Rec_nat_H) + len(Rec_nat_L) + len(Rec_vax_H) + len(Rec_vax_L)
print("------------ population ---------- ")
print(len(Pop_H) + len(Pop_L))
print("----- sum of counts --------- ")
print(testing)
print("--------- sum of classes ------------ ")
print(classes)


df = pnd.DataFrame(rows, columns=['indexID','time','N','Inf_L','Inf_H','Dead_L','Dead_H','cum_cases_L','cum_cases_H','cum_sum_per_N',
                      'Exposed_L', "Exposed_H", "Sus_L", "Sus_H", "Rec_nat_L", "Rec_nat_H", "Rec_vax_L", "Rec_vax_H",
                      "vaccinatedH", "vaccinatedL", "boostedH", "boostedL", 
                      "count_Naive", "count_Delta", "count_Delta_vaccine", "count_two", "count_three", "count_four", 
                      "count_Omicron", "count_Omicron_star", "count_Omicron_star_vaccine", 
                      "count_Omicron_vaccine", "count_Vaccine",
                      "count_WT", "count_WT_vaccine", "count_WT_vaccine_boost", "count_Omicron_vaccine_boost",
                      "count_Delta_vaccine_boost", "count_Omicron_star_vaccine_boost", "count_Vaccine_boost"])
df.to_csv('output_files_AS/output_' + str(multPop) + 'x_ID_' + str(indexID) + '.csv')
