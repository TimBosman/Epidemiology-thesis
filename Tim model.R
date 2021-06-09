setwd("C:/Users/Tim/Desktop/Epidemiology Thesis/Epidemiology-thesis")
source("Epidemiology_Functions.R")

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
contactrate <- 0.03
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
Write_infectivity_file_for_SIRE(InfectedPedigree, "SIRE.txt", timepoints)

GLMM_Data <- generate_GLMM_Data(InfectedPedigree, timepoints)

saveRDS(GLMM_Data, file = "dif.RDS")

GLMM_Data = readRDS(file = "dif.RDS")
GLMM_Data$Herd <- as.factor(GLMM_Data$Herd)
GLMM_Data$Sire <- as.factor(GLMM_Data$Sire)
GLMM_Data$S <- as.integer(GLMM_Data$S)
GLMM_Data$C <- as.integer(GLMM_Data$C)
GLMM_Data$I <- as.integer(GLMM_Data$I)
GLMM_Data$DeltaT <- as.numeric(GLMM_Data$DeltaT)

library(lme4)
model = glmer(data = GLMM_Data,
              cbind(C, S-C) ~ (1 | Sire) + (1 | Herd), 
              offset = log(GLMM_Data$I/nherds * GLMM_Data$DeltaT),
              family = binomial(link="cloglog"))
summary(model)
