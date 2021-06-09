source("Epidemiology_Functions.R")

### Input parameters ##########################################################
## population stats ##
nsire <- 100
ndams <- nsire * 100
nherds <- 100
nSiresPerHerd <- 5
nOffspringPerHerd <- 20 # number of offspring of one sire in herd

## Trait stats ##
vArec <- 0.2 # variation in recoverability
vErec <- 0.2 # Environmental variation in recoverability

## Infection stats ##
alpha <- 0.1
contactrate <- 0.15
timepoints <- c(0, 280, 1400)

### Main script ###############################################################

## Simulate pedigree and add traits and herds ## 
pedigree <- get_pedigree(nsire, ndams)
pedigree <- add_trait_to_pedigree("recoverability", vArec, vErec, pedigree)
pedigree <- set_herd(pedigree, nherds, nOffspringPerHerd, nSiresPerHerd)


## Simulate infection for every herd ##
InfectedPedigree <- data.frame()
events <- data.frame()
for(herd in levels(pedigree$herd)){
  herddata <- pedigree[pedigree$herd == herd,]
  output <- simulate_infection(herddata, alpha, contactrate, max(timepoints),
                               recoverability = herddata$recoverability, 
                               model =  "SIS")
  InfectedPedigree <- rbind(InfectedPedigree, output[[1]])
  events <- rbind(events, output[[2]])
}

## Generate timeseries barplot ##
Plot_infected_fraction(events, 1:nherds, 1400)

