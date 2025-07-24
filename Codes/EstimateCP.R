library(tidyverse)


cpest<-function(wt, delta, omicron, omicronstar, length2, length3, length4){
  totinfect=sum(wt, delta, omicron, omicronstar, length2, length3, length4)
  cp=(wt/totinfect)*0+((delta/totinfect)*0.03)+((omicron/totinfect)*0.5)+((omicronstar/totinfect)*0.6)+
    ((length2/totinfect)*0.35)+((length3/totinfect)*0.65)+((length4/totinfect)*0.95)
  return(cp)
}

base.imm <- function(naive, vaccine){
  totinfect=100-sum(naive, vaccine)
  return(totinfect)
}


#### MAIN ####

data <- read_csv("../../Data/ABM Main Runs/dataHistories899.csv") %>% filter(name == "Mean")

malaysia_est <- cpest(wt=0.07, delta=0.31, omicron=2.65, omicronstar=0, length2=4.23, length3=33.14, length4=22.60)
ecuador_est <- cpest(wt = 0.11, delta = 0.21, omicron = 4.36, omicronstar = 0, length2 = 3.56, length3 = 24.04, length4 = 8.37)
india_est <- cpest(wt = 1.20, delta = 5.01, omicron = 1.56, omicronstar = 0, length2 = 2.90, length3 = 9.31, length4 = 1.38)

basemalaysia_est <- base.imm(naive = 0.44, vaccine = 8.85)
baseecuador_est <- base.imm(naive = 1.20, vaccine = 20.62)
baseindia_est <- base.imm(naive = 6.98, vaccine = 40.47)

#### T CELL SENSITIVITY ANALYSIS ####

data <- read_csv("../../Data/ABM Supplement/TCell_Severity/concatenated data/dataHistories899.csv") %>% filter(name == "Mean")

malaysia_est_tcell <- cpest(wt = 0.07, delta = 0.30, omicron = 2.65, omicronstar = 0, length2 = 4.22, length3 = 33.15, length4 = 22.61)
ecuador_est_tcell <- cpest(wt = 0.11, delta = 0.21, omicron = 4.35, omicronstar = 0, length2 = 3.56, length3 = 24.05, length4 = 8.40)
india_est_tcell <- cpest(wt = 1.23, delta = 5.10, omicron = 1.47, omicronstar = 0, length2 = 2.87, length3 = 9.23, length4 = 1.36)

basemalaysia_est_tcell <- base.imm(naive = 0.44, vaccine = 8.86)
baseecuador_est_tcell <- base.imm(naive = 1.20, vaccine = 20.62)
baseindia_est_tcell <- base.imm(naive = 6.99, vaccine = 40.54)




