source("Epidemiology_Functions.R")

### Input parameters ##########################################################
## population stats ##
nsire <- 102
ndams <- nsire * 102
nherds <- 102
nSiresPerHerd <- 6
nOffspringPerHerd <- 17 # number of offspring of one sire in herd

## Trait stats ##
vAsus <- 0.5 # variation in suseptibility
vEsus <- 0.5 # Environmental variation in suseptibility
vAinf <- 0.5 # variation in infectivity
vEinf <- 0.5 # Environmental variation in infectivity

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
  output <- simulate_infection(herddata, alpha, contactrate, max(timepoints),
                               infectivity = herddata$infectivity,
                               susceptibility = herddata$susceptibility)
  InfectedPedigree <- rbind(InfectedPedigree, output[[1]])
  events <- rbind(events, output[[2]])
}

## Generate timeseries barplot ##
pedigree <- Generate_time_series_data(timepoints, events, InfectedPedigree)

Plot_time_series(pedigree, timepoints)

Plot_infected_fraction(events, 1:102)

Write_infectivity_file_for_SIRE(pedigree, "SIRE.txt", timepoints)
