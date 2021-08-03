##Scenario 2: Less offspring

setwd("C:/Users/Tim/Desktop/Epidemiology Thesis/Epidemiology-thesis")
source("Epidemiology_Functions.R")

# To reproduce the results used in the thesis
set.seed(4965) 

### Input parameters ##########################################################
## population stats ##
nsire <- 100
ndams <- nsire * 50
nherds <- 100
nSiresPerHerd <- 5
nOffspringPerHerd <- 10 # number of offspring of one sire in herd

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
## Simulate pedigree and add traits and herds
pedigree <- get_pedigree(nsire, ndams)
pedigree <- add_trait_to_pedigree("susceptibility", vAsus, vEsus, pedigree, 
                                  SireBVFile = "BVsusScen2.csv")
pedigree <- add_trait_to_pedigree("infectivity", vAinf, vEinf, pedigree, 
                                  SireBVFile = "BVinfScen2.csv")
pedigree <- set_herd(pedigree, nherds, nOffspringPerHerd, nSiresPerHerd)
## Simulate infection for every herd
InfectedPedigree <- data.frame()
events <- data.frame()
for(herd in levels(pedigree$herd)){
  herddata <- pedigree[pedigree$herd == herd,]
  repeat{
    output <- simulate_infection(herddata, alpha, contactrate, max(timepoints),
                                 infectivity = herddata$infectivity,
                                 susceptibility = herddata$susceptibility,
                                 initialstate = sample(c(rep("S", nrow(herddata) - 1), "I")))
    if(nrow(output[[2]]) > 2){
      break
    }
  }
  InfectedPedigree <- rbind(InfectedPedigree, output[[1]])
  events <- rbind(events, output[[2]])
}
rm(pedigree, output, herddata)
## Generate timeseries data
InfectedPedigree <- Generate_time_series_data(timepoints, events, InfectedPedigree)
##Plot infection information
Plot_time_series(InfectedPedigree, timepoints)
Plot_infected_fraction(events, 1:nherds, max(timepoints))
rm(events)
##Write file for SIRE
Write_infectivity_file_for_SIRE(InfectedPedigree, "Scenario2.txt", timepoints)
##Prepare data for GLMM
GLMM_Data <- generate_GLMM_Data(InfectedPedigree, timepoints)
GLMM_Data <- addInfectedSiresInHerd(GLMM_Data)
##Make Infectivity effect matrix
infectivityeffects = matrix(0, ncol = length(levels(GLMM_Data$Sire)), nrow = nrow(GLMM_Data))
colnames(infectivityeffects) <- as.character(1:length(levels(GLMM_Data$Sire)))
for(i in 1:nrow(infectivityeffects)){
  for(sire in 1:nSiresPerHerd){
    infectivityeffects[i, GLMM_Data[i, paste0("sire", sire)]] <- as.numeric(GLMM_Data[i, paste0("I", sire)] / GLMM_Data[i,"I"])
  }
}
infectivityeffects <- infectivityeffects[,sort(colnames(infectivityeffects))]
#Add dummy variable to replace with the infectivityeffects
GLMM_Data$InfSire <- GLMM_Data$Sire
#GLMM model
m2 <- glmer_multimemb(formula = cbind(C, S-C) ~  (1 | InfSire) + (1 | Sire) +  Herd,
                      data = GLMM_Data,
                      memb_mat = list(InfSire = infectivityeffects),
                      herdsize = nSiresPerHerd * nOffspringPerHerd)
summary(m2)
#plot Susceptiblity information
data = GetEstimations("BVsusScen2.csv", "SireOutput/Scenario2Sire.txt", m2, "b_g.", "Sire")
SIREplot = ggplot(data) + geom_point(aes(x = BV.for.susceptibility, y = SIRE), col = "#F8766D") +
  ggtitle("Susceptibility estimation by SIRE") +
  geom_abline(slope = 1) +
  geom_text(x = -0.75, y = 0.5,  label = paste("SIRE r2 =", round(cor(data$BV.for.susceptibility, data$SIRE)^2, 3)))
GLMMplot = ggplot(data) + geom_point(aes(x = GLMM, y = BV.for.susceptibility), col = "#F8766D") +
  ggtitle("Susceptibility estimation by GLMM") +
  geom_text(x = -0.45, y = 1,  label = paste("GLMM r2 =", round(cor(data$BV.for.susceptibility, data$GLMM)^2, 3)))
SIREplot + GLMMplot  + plot_layout(nrow = 2)
#plot infectivity information
data = GetEstimations("BVinfScen2.csv", "SireOutput/Scenario2Sire.txt", m2, "b_f.", "InfSire")
SIREplot = ggplot(data) + geom_point(aes(x = BV.for.infectivity, y = SIRE), col = "#00BA38") +
  ggtitle("Infectivity estimation by SIRE") +
  geom_abline(slope = 1) +
  geom_text(x = -0.9, y = 1.5,  label = paste("SIRE r2 =", round(cor(data$BV.for.infectivity, data$SIRE)^2, 3)))
GLMMplot = ggplot(data) + geom_point(aes(x = GLMM, y = BV.for.infectivity), col = "#00BA38") +
  ggtitle("Infectivity estimation by GLMM") +
  geom_text(x = -0.015, y = 1,  label = paste("GLMM r2 =", round(cor(data$BV.for.infectivity, data$GLMM)^2, 3)))
SIREplot + GLMMplot  + plot_layout(nrow = 2)
