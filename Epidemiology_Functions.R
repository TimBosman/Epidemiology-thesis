#### Imports ##################################################################
library(tidyr)
library(plyr)
library(ggplot2)
library(patchwork)

### Functions #################################################################

get_pedigree <- function(nsires, ndams = NA, OffspringPerMate = 1, 
                         ndamsPerSire = NA){
  #If number of dams per sire is set overrule the ndams
  if (is.na(ndamsPerSire)) {
    ndamsPerSire <- ndams / nsires
  }
  #Calculate the total number of offspring
  noff <- nsires*ndamsPerSire * OffspringPerMate
  #Get numbers for individuals
  sires <- rep(1:nsires, each = ndamsPerSire)
  dams <- (nsires + 1):((nsires * ndamsPerSire) + nsires)
  offspring <- (dams[length(dams)] + 1):(dams[length(dams)] + noff)
  ped <- as.data.frame(cbind(offspring, sires, dams))
  ped$sires <- as.factor(ped$sires)
  ped$dams <- as.factor(ped$dams)
  ped$offspring <- as.factor(ped$offspring)
  return(ped)
}

add_trait_to_pedigree <- function(name, vA, vE, pedigree, SireBVFile = NA, 
                                  DamBVFile = NA){
  # nsires <- length(levels(pedigree$sires))
  # mates <- nrow(pedigree) / nsires
  Adams <- rnorm(length(levels(pedigree$dams)), 0, sqrt(vA))
  Asire <- rnorm(length(levels(pedigree$sires)), 0, sqrt(vA))
  vMS <- 0.5 * vA
  MS <- rnorm(nrow(pedigree), 0, sqrt(vMS))
  Arec <- 0.5 * Asire[pedigree$sires] + 0.5 * Adams[pedigree$dams] + MS
  Erec <- rnorm(nrow(pedigree), 0, sqrt(vE))
  trait <- exp(Arec + Erec)
  pedigree <- cbind(pedigree, trait, Arec)
  names(pedigree)[names(pedigree) %in%  c("trait", "Arec")] <- 
    c(name, paste("breeding value of", name))
  
  if (!is.na(SireBVFile)) {
    df <- unique.data.frame(cbind(pedigree$sires, Asire[pedigree$sires]))
    colnames(df) <- c("Sire ID", paste("BV for", name))
    write.csv(df, SireBVFile, row.names = FALSE)
  }
  if (!is.na(DamBVFile)){
    df <- unique.data.frame(cbind(pedigree$dams, Adams[pedigree$dams]))
    colnames(df) <- c("Dam ID", paste("BV for", name))
    write.csv(df, DamBVFile, row.names = FALSE)
  }
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
      if (sum(apply(samples_sgt, MARGIN = 1, function(x) length(unique(x))))
          == (nrow(samples_sg)*s)) break
    }
  }
  offspring <- matrix(nrow = 0, ncol = ncol(pedigree)+1)
  colnames(offspring) <- c(colnames(pedigree), "group")
  for(s in 1:length(levels(pedigree$sires))){
    offspring_sel <- pedigree[pedigree$sires==s,]
    offspring_sel$herd <- rep(samples_sg[s,1:nsires_g], off_sire_g)
    offspring <- rbind(offspring, offspring_sel)
  }
  offspring$herd <- as.factor(offspring$herd)
  return(offspring)
}

simulate_infection <- function(herddata, alpha, contactrate, totaltime, 
                               susceptibility = rep(1, nrow(herddata)), 
                               infectivity = rep(1, nrow(herddata)), 
                               recoverability = rep(1, nrow(herddata)), 
                               model = "SIR", 
                               initialstate = c("S", "I")[rbinom(nrow(herddata), 1, 1-1/(contactrate / alpha)) + 1]){
  herddata$initialstate <- initialstate
  herddata$currentstate <- herddata$initialstate
  time <- 0
  FractionInfected <- sum(herddata$initialstate == "I") / nrow(herddata)
  events <- data.frame(0, NA, NA, NA, NA, FractionInfected)
  colnames(events) <- c("Time", "Event", "Cow ID", "herd", "Status after event", "Fraction infected")
  while(time < totaltime & sum(herddata$currentstate == "I") >= 1){
    #Calculate beta
    beta <- contactrate * 
      mean(infectivity[which(herddata$currentstate == "I")])
    #Calculate chances of a happening
    Rinf <- beta * sum(herddata$currentstate == "I") * 
      susceptibility * as.numeric(herddata$currentstate == "S") / 
      length(herddata$currentstate)
    Rrec <- recoverability * as.numeric(herddata$currentstate == "I") * alpha
    totalR <- sum(Rinf) + sum(Rrec)
    #take cumsums
    cuminf <- cumsum(Rinf / totalR)
    cumrec <- cumsum(Rrec / totalR) + cuminf[length(cuminf)] 
    randomNumber <- runif(1, 0, 1)
    if (randomNumber <= cuminf[length(cuminf)]) {
      event <- "Infection"
      FractionInfected <- FractionInfected + 1 / nrow(herddata)
      IndI <- min(which(cuminf >= randomNumber))
      IndID <- as.character(herddata$offspring[IndI])
      newstatus = "I"
      herddata$currentstate[IndI] <- newstatus
    } else {
      event <- "Recovery"
      FractionInfected <- FractionInfected - 1 / nrow(herddata)
      IndI <- min(which(cumrec >= randomNumber))
      IndID <- as.character(herddata$offspring[IndI])
      if(model == "SIR"){
        newstatus = "R"
        herddata$currentstate[IndI] <- newstatus
      } else if (model == "SIS") {
        newstatus <- "S"
        herddata$currentstate[IndI] <- newstatus
      }
    }
    #Calculate time point
    time <- time + rexp(1, totalR)
    events <- rbind(events, c(time, event, IndID, herd, newstatus, FractionInfected))
  }
  events$Time <- as.numeric(events$Time)
  events <- events[!is.na(events$Event), ]
  return(list(herddata, events))
}

Generate_time_series_data <- function(timepoints, events, pedigree, model = "SIR"){
  events$Time <- as.numeric(events$Time)
  for(timepoint in timepoints) {
    temp_events <- events[events$Time <= timepoint, ]
    temp_events <- ddply(temp_events, .(`Cow ID`), subset, subset = Time == 
                           max(Time), select = c(`Cow ID`, Time, Event))
    recovered <- temp_events$`Cow ID`[temp_events$Event == "Recovery"]
    infected <- temp_events$`Cow ID`[temp_events$Event == "Infection"]
    pedigree[, as.character(timepoint)] <- pedigree$initialstate
    pedigree[as.character(pedigree$offspring) %in% infected, 
             as.character(timepoint)] <- "I"
    if (model == "SIR"){
      pedigree[as.character(pedigree$offspring) %in% recovered, 
               as.character(timepoint)] <- "R"
    } else {
      pedigree[as.character(pedigree$offspring) %in% recovered, 
               as.character(timepoint)] <- "S"
    }
  }
  return(pedigree)
}

Plot_time_series <- function(TimeSeries, timepoints){
  temp <- pivot_longer(TimeSeries, as.character(timepoints), "Time")
  temp$Time <- as.factor(as.numeric(temp$Time))
  temp$value[temp$value == "R"] <- "Recovered"
  temp$value[temp$value == "I"] <- "Infected"
  temp$value[temp$value == "S"] <-  "Susceptible"
  ggplot(temp) + geom_bar(aes(x= Time, fill = value), stat = "count") 
}

Plot_infected_fraction <- function(events, herds, timelimit){
  events <- events[events$herd %in% herds, ]
  points <- ggplot(events, aes(x= Time, y = as.numeric(`Fraction infected`))) +
    geom_bin2d() + ylim(0, 1) + scale_fill_gradient(low="Grey", high="Black") +
    xlim(0, timelimit) + ylab("Fraction infected animals") + ggtitle("A")
  line <-  ggplot(events, aes(x= Time, y = as.numeric(`Fraction infected`))) +
    geom_smooth(se = TRUE, method = "gam")  + ylim(0, 1) + xlim(0, timelimit) +
    ylab("Fraction infected animals") + ggtitle("B")

  points + line + plot_layout(nrow = 2)
}

Write_infectivity_file_for_SIRE <- function(pedigree, filename, timepoints){
  pedigree$type <- "Contact"
  pedigree$type[pedigree$`0` == "I"] <- "Seeder"
  pedigree <- pivot_longer(pedigree, as.character(timepoints), "time")
  write.table(pedigree, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}