#### Imports ##################################################################
library(tidyr)
library(plyr)
library(ggplot2)
library(patchwork)
library(lme4)

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
  #If the BV for the sires are required later save them in a file
  if (!is.na(SireBVFile)) {
    df <- unique.data.frame(cbind(pedigree$sires, Asire[pedigree$sires]))
    colnames(df) <- c("Sire ID", paste("BV for", name))
    write.csv(df, SireBVFile, row.names = FALSE)
  }
  #Same for dams
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
  #WARINING only works for nherds == nsires
  if(nherds != length(levels(pedigree$sires))) stop("This algorithm is not prepared for situations where nherds != nsires")
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
  for(s in as.numeric(levels(pedigree$sires))){
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
    events <- rbind(events, c(time, event, IndID, herddata$herd[1], newstatus, FractionInfected))
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
  temp$value = factor(temp$value, levels = c("Susceptible", "Infected", "Recovered"))
  ggplot(temp) + geom_bar(aes(x= Time, fill = value), stat = "count") +
  ylab("Number of infected individuals")
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
  pedigree <- pedigree[, c("offspring", "sires", "herd", "time", "type", "value")]
  write.table(pedigree, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}

generate_GLMM_Data <- function(InfectedPedigree, timepoints){
  dif = data.frame(matrix(NA, ncol = 9))
  names(dif)<- c("Herd", "Sire", "S", "C", "DeltaT", "I", "t", "Isire" , "R")
  row = 1
  for (herd in levels(InfectedPedigree$herd)){
    herdrows = which(InfectedPedigree$herd == herd)
    for(timepoint in as.character(timepoints[2:length(timepoints)])){
      previoustime = as.character(timepoints[which(timepoints == timepoint)-1])
      I = sum(InfectedPedigree[herdrows, previoustime] == "I")
      if(I < 1){break}
      DeltaT = as.numeric(timepoint) - as.numeric(previoustime)
      for(sire in unique(InfectedPedigree$sires[herdrows])){
        t1_value <- InfectedPedigree[herdrows[herdrows %in% which(InfectedPedigree$sires == sire)], previoustime]
        t2_value <- InfectedPedigree[herdrows[herdrows %in% which(InfectedPedigree$sires == sire)], timepoint]
        Isire = sum(t1_value == "I")
        S = sum(t1_value == "S")
        C = S - sum(t2_value == "S")
        R = sum(t2_value == "R") - sum(t1_value == "R")
        dif[row, ] = c(herd, sire, S, C, DeltaT, I, timepoint, Isire, R)
        row = row + 1
      }
    }
  }
  dif$Herd <- as.factor(dif$Herd)
  dif$Sire <- as.factor(dif$Sire)
  dif$S <- as.integer(dif$S)
  dif$C <- as.integer(dif$C)
  dif$I <- as.integer(dif$I)
  dif$DeltaT <- as.numeric(dif$DeltaT)
  dif$Isire <- as.integer(dif$Isire)
  dif$R <- as.integer(dif$R)
  
  return(dif)
}

addInfectedSiresInHerd <- function(GLMM_Data){
  for (row in 1:nrow(GLMM_Data)){
    rows = which(GLMM_Data$Herd == GLMM_Data$Herd[row] & GLMM_Data$t == GLMM_Data$t[row])
    count = 1
    for (sire in GLMM_Data$Sire[rows]){
      GLMM_Data[rows, paste0("sire", count)] = paste(sire)
      GLMM_Data[rows, paste0("I", count)] = GLMM_Data[rows[rows %in% which(GLMM_Data$Sire == sire)], "Isire"]
      count = count + 1
    }
  }
  return(GLMM_Data)
}

# quickRun <- function(nsire, ndams, nherds,nSiresPerHerd, nOffspringPerHerd, vAsus, vEsus, vAinf,vEinf, alpha, contactrate, timepoints){
#   pedigree <- get_pedigree(nsire, ndams)
#   pedigree <- add_trait_to_pedigree("susceptibility", vAsus, vEsus, pedigree, 
#                                     SireBVFile = "BVsus.csv")
#   pedigree <- add_trait_to_pedigree("infectivity", vAinf, vEinf, pedigree, 
#                                     SireBVFile = "BVinf.csv")
#   pedigree <- set_herd(pedigree, nherds, nOffspringPerHerd, nSiresPerHerd)
#   
#   ## Simulate infection for every herd ##
#   InfectedPedigree <- data.frame()
#   events <- data.frame()
#   for(herd in levels(pedigree$herd)){
#     herddata <- pedigree[pedigree$herd == herd,]
#     repeat{
#       output <- simulate_infection(herddata, alpha, contactrate, max(timepoints),
#                                    infectivity = herddata$infectivity,
#                                    susceptibility = herddata$susceptibility,
#                                    initialstate = sample(c(rep("S", 99), "I")))
#       if(nrow(output[[2]]) > 2){
#         break
#       }
#     }
#     InfectedPedigree <- rbind(InfectedPedigree, output[[1]])
#     events <- rbind(events, output[[2]])
#   }
#   rm(pedigree, output, herddata)
#   ## Generate timeseries barplot ##
#   InfectedPedigree <- Generate_time_series_data(timepoints, events, InfectedPedigree)
#   
#   # Plot_time_series(InfectedPedigree, timepoints)
#   # 
#   # Plot_infected_fraction(events, 1:nherds, max(timepoints))
#   rm(events)
#   # Write_infectivity_file_for_SIRE(InfectedPedigree, "Scenario1.txt", timepoints)
#   
#   GLMM_Data <- generate_GLMM_Data(InfectedPedigree, timepoints)
#   
#   GLMM_Data <- addInfectedSiresInHerd(GLMM_Data)
#   return(GLMM_Data)
# }

GetEstimations <- function(bvFile, Sirefile, GLMM_model, sireCode, glmmterm){
  #Get real BVs
  data <- read.csv(bvFile)
  
  #Get SIRE estimates
  sireEstimates = read.table(Sirefile, header = TRUE)
  SIREestimates = sapply(c(2:100), function(i){
    mean(sireEstimates[,paste0(sireCode, i)])
  })
  SIREestimates <- c(0, SIREestimates)
  data$SIRE <- SIREestimates
  
  #Get GLMM estimates
  data2 <-data.frame(ranef(GLMM_model)[glmmterm])
  data$GLMM <- data2[match(data$Sire.ID, rownames(data2)),1]
  return(data)
}

glmer_multimemb <- function(formula, data, memb_mat=list(), herdsize, ...) {
  data$herdsize <- herdsize
  mnms <- names(memb_mat)
  fb <- findbars(formula)
  gvars <- vapply(fb, function(x) deparse(x[[3]]), character(1))
  Ztlist <- list()
  for (i in seq_along(fb)) {
    fbnm <- deparse(fb[[i]])
    ## find corresponding random-effects term
    w <- which(mnms==gvars[i])
    if (length(w)>0) {
      M <- Matrix::Matrix(memb_mat[[w]])
      ## extract LHS (effect)
      form <- as.formula(substitute(~z,list(z=fb[[i]][[2]])))
      ## construct model matrix & compute Khatri-Rao product
      X <- model.matrix(form,data=data)
      Zt <- Matrix::KhatriRao(t(M),t(X),make.dimnames=TRUE)
      ## FIXME: mess with names?
      Ztlist[[fbnm]] <- Zt
    } ## if  (length(w)>0)
  } ## for i in seq(fb)
  lmod <- lFormula(formula,data=data)
  ## substitute new Ztlist elements
  # m = names(Ztlist)
  for (m in names(Ztlist)) {
    lmod$reTrms$Ztlist[[m]] <- Ztlist[[m]]
  }
  lmod$reTrms$Zt <- do.call(rbind,lmod$reTrms$Ztlist)
  ## finish fitting
  m1 <- glmer(lmod,
              family = binomial(link="cloglog"),
              data = data,
              offset = log(data$I / herdsize * data$DeltaT))
  return(m1)
}
