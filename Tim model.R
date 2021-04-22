#### Imports #################################################################
#No imports

### Functions ################################################################

get_pedigree <- function(nsires, ndams = NA, OffspringPerMate = 1, 
                         ndamsPerSire = NA){
  #If number of dams per sire is set overrule the ndams
  if(is.na(ndamsPerSire)){
    ndamsPerSire <- ndams / nsires
  }
  #Calculate the total number of offspring
  noff <- nsires*ndamsPerSire * OffspringPerMate
  #Get numbers for individuals
  sires = rep(1 : nsires, each=ndamsPerSire)
  dams = (nsires + 1) : ((nsires * ndamsPerSire) + nsires)
  offspring = (dams[length(dams)] + 1) : (dams[length(dams)] + noff)
  ped <- as.data.frame(cbind(offspring, sires, dams))
  ped$sires <- as.factor(ped$sires)
  ped$dams <- as.factor(ped$dams)
  ped$offspring <- as.factor(ped$offspring)
  return(ped)
}

add_trait_to_pedigree <- function(name, vA, vE, pedigree){
  nsires <- length(levels(pedigree$sires))
  mates <- nrow(pedigree) / nsires
  Adams <- rnorm(nrow(pedigree),0, sqrt(vA))
  Asire <- rnorm(nsires,0, sqrt(vA))
  vMS <- 0.5*vA
  MS <- rnorm(nrow(pedigree),0,sqrt(vMS))
  Arec <- 0.5* c(sapply(Asire, function(sire){rep(sire, mates)})) + 0.5 * Adams + MS
  Erec <- rnorm(nrow(pedigree),0,sqrt(vE))
  trait <- exp(Arec+Erec)
  pedigree <- cbind(pedigree, trait)
  colnames(pedigree)[colnames(pedigree)=="trait"] <- name
  return(pedigree)
}

set_herd <- function(pedigree, nherds, off_sire_g, nsires_g){
  #Based on script of Dries Hulst
  #Download on: https://doi.org/10.25386/genetics.13090076 (S4)
  samples_sg <- matrix(nrow = nherds, ncol = nsires_g)
  samples_sg[,1] <- sample(1:nherds, nherds, replace = FALSE)
  for (s in 2:nsires_g){
    repeat{
      samples_sg[,s]<- sample(1:nherds, nherds, replace = FALSE)
      samples_sgt<-samples_sg[,1:s]
      if (sum(apply(samples_sgt,MARGIN = 1, function(x) length(unique(x))))==(nrow(samples_sg)*s)) break
    }
  }
  offspring <- matrix(nrow = 0, ncol = ncol(pedigree)+1)
  colnames(offspring) <- c(colnames(pedigree), "group")
  for(s in 1:max(levels(pedigree$sires))){
    offspring_sel <- pedigree[pedigree$sires==s,]
    offspring_sel$herd <- rep(samples_sg[s,1:nsires_g],off_sire_g)
    offspring <- rbind(offspring, offspring_sel)
  }
  offspring$herd <- as.factor(offspring$herd)
  return(offspring)
}

simulateinfection <- function(herddata, alpha, beta, R0, totaltime){
  # return(herddata)
}

### Input parameters #########################################################
## population stats ##
nsire = 102
ndams = nsire * 102
nherds = 102
nSiresPerHerd = 6
nOffspringPerHerd = 17 # number of offspring of one sire in herd

## Trait stats ##
vAsus = 0.5 # variation in suseptibility
vEsus = 0.5 # Environmental variation in suseptibility
vAinf = 0.5 # variation in infectivity
vEinf = 0.5 # Environmental variation in infectivity

## Infection stats ##
R0 = 1.5
alpha = 0.02
beta = 0.03
totaltime = 1

### Main script ############################################################## 

## Simulate pedigree and add traits and herds ## 
pedigree <- get_pedigree(nsire, ndams)
pedigree <- add_trait_to_pedigree("suseptibility", vAsus, vEsus, pedigree)
pedigree <- add_trait_to_pedigree("infectivity", vAinf, vEinf, pedigree)
pedigree <- set_herd(pedigree, nherds, nOffspringPerHerd, nSiresPerHerd)

## Simulate the infection for all the herds
for(herd in levels(pedigree$herd)){
  herddata <- pedigree[pedigree$herd == herd,]
}

Prevalence<- 1-1/R0
herddata$initialstate <- c("S", "I")[rbinom(nrow(herddata),1,Prevalence) + 1] 
herddata$currentstate <- herddata$initialstate
time = 0
happenings = c(0, "-", "-")
while(time < totaltime & herddata$currentstate != rep("R", nrow(herddata))){
  #Calculate chances of a happening
  Rinf = (herddata$currentstate == "S") * herddata$suseptibility * 
    mean((herddata$currentstate == "I") * herddata$infectivity)
  Rrec = (herddata$currentstate == "I") * alpha
  totalR = sum(Rinf) + sum(Rrec)
  #Calculate next time point
  time =+ rexp(1, Rtotal)
}
