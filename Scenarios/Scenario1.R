##Scenario 1: default model

setwd("C:/Users/Tim/Desktop/Epidemiology Thesis/Epidemiology-thesis")
source("Epidemiology_Functions.R")
library(lme4)

# To reproduce the results used in the thesis
set.seed(8736)

### Input parameters ##########################################################
## population stats ##
nsire <- 100
ndams <- nsire * 100
nherds <- 100
nSiresPerHerd <- 5
nOffspringPerHerd <- 20 # number of offspring of one sire in herd

## Trait stats ##
vAsus <- 0.25 # variation in suseptibility
vEsus <- 0.25 # Environmental variation in suseptibility
vAinf <- 0.25 # variation in infectivity
vEinf <- 0.25 # Environmental variation in infectivity

## Infection stats ##
alpha <- 0.02
contactrate <- 0.06
timepoints <- 0:20 * 14

### Main script ###############################################################

## Simulate pedigree and add traits and herds ## 
pedigree <- get_pedigree(nsire, ndams)
pedigree <- add_trait_to_pedigree("susceptibility", vAsus, vEsus, pedigree, 
                                  SireBVFile = "BVsus.csv")
pedigree <- add_trait_to_pedigree("infectivity", vAinf, vEinf, pedigree, 
                                  SireBVFile = "BVinf.csv")
pedigree <- set_herd(pedigree, nherds, nOffspringPerHerd, nSiresPerHerd)

## Simulate infection for every herd ##
InfectedPedigree <- data.frame()
events <- data.frame()
for(herd in levels(pedigree$herd)){
  herddata <- pedigree[pedigree$herd == herd,]
  repeat{
    output <- simulate_infection(herddata, alpha, contactrate, max(timepoints),
                                 infectivity = herddata$infectivity,
                                 susceptibility = herddata$susceptibility,
                                 initialstate = sample(c(rep("S", 99), "I")))
    if(nrow(output[[2]]) > 2){
      break
    }
  }
  InfectedPedigree <- rbind(InfectedPedigree, output[[1]])
  events <- rbind(events, output[[2]])
}
rm(pedigree, output, herddata)
## Generate timeseries barplot ##
InfectedPedigree <- Generate_time_series_data(timepoints, events, InfectedPedigree)

Plot_time_series(InfectedPedigree, timepoints)

Plot_infected_fraction(events, 1:nherds, max(timepoints))
rm(events)
Write_infectivity_file_for_SIRE(InfectedPedigree, "Scenario1.txt", timepoints)

GLMM_Data <- generate_GLMM_Data(InfectedPedigree, timepoints)

for (row in 1:nrow(GLMM_Data)){
  rows = which(GLMM_Data$Herd == GLMM_Data$Herd[row] & GLMM_Data$t == GLMM_Data$t[row])
  count = 1
  for (sire in GLMM_Data$Sire[rows]){
    GLMM_Data[rows, paste0("sire", count)] = paste(sire)
    GLMM_Data[rows, paste0("I", count)] = GLMM_Data[rows[rows %in% which(GLMM_Data$Sire == sire)], "Isire"]
    count = count + 1
  }
}

#small check
sum(sapply(1:nrow(GLMM_Data), function(row){
  sum(GLMM_Data[row, c("I1","I2","I3","I4","I5")]) != GLMM_Data$I[row]
}))

#(I1 | sire1) + (I2 | sire2) + (I3 | sire3) + (I4 | sire4) + (I5 | sire5)

model = glmer(data = GLMM_Data,
              cbind(C, S-C) ~  + (1 | Sire) + (1 | Herd), 
              offset = log(GLMM_Data$I/nherds * GLMM_Data$DeltaT),
              family = binomial(link="cloglog"))
summary(model)
ranef(model)

