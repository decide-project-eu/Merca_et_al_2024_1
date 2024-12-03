### Functions for the Salmon study

#libraries
library(dplyr)
library(lubridate)
library(splines)
library(gdata)
library('MuMIn')
library(rlist)
library(stringr)


################################################################################################################
### CREATE LEARNING AND TEST SETS ###

# Function to create Learning and Test sets for DLMs ----
# - df is the original data set
# - relevant.names are the names of the variables to be co-modeled
# - N is the time step where the datasets are divided into learning and test sets
learning.test.sets <- function(df, relevant.names, N){
  
  # - see which nseqs != 0 are cut and move them to the closest side
  set_nseqs <- data.frame("nseq"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "start"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "end"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "set"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])))
  loop=0
  for(i in unique(df$nseq)[unique(df$nseq)!=0]){
    
    dates <- unique(df[order(df[,"date"]), "date"])
    df_nseq <- subset(df, df$nseq == i)
    min_nseq <- df_nseq[1,"date"]
    max_nseq <- df_nseq[dim(df_nseq)[1],"date"]
    
    range <- which(dates==min_nseq):which(dates==max_nseq)
    a <- intersect(range,N+1)
    
    loop=loop+1
    set_nseqs[loop, "nseq"] <- i
    set_nseqs[loop, "start"] <- which(dates==min_nseq)
    set_nseqs[loop, "end"] <- which(dates==max_nseq)
    
    if(isTRUE(length(a)>0)){ #what to do with the nseqs that are cut
      dist_min <- N - which(dates==min_nseq)
      dist_max <- which(dates==max_nseq) - N
      
      if(dist_min < dist_max){
        set_nseqs[loop, "set"] <- "Test"
      }
      if(dist_max < dist_min){
        set_nseqs[loop, "set"] <- "Learning"
      }
      if(dist_max == dist_min){ #if it's in the middle goes to Learning.set
        set_nseqs[loop, "set"] <- "Learning"
      }
      
    }else{#put the others in learning or test set
      
      if(isTRUE(which(dates==max_nseq) <= N)){
        set_nseqs[loop, "set"] <- "Learning"
      }else{
        set_nseqs[loop, "set"] <- "Test"
      }
    }
  }
  sum(set_nseqs$set=="Learning")/dim(set_nseqs)[1]*100   
  sum(set_nseqs$set=="Test")/dim(set_nseqs)[1]*100       
  
  learning.nseqs <- set_nseqs$nseq[set_nseqs$set== "Learning"]
  test.nseqs <- set_nseqs$nseq[set_nseqs$set== "Test"]
  
  # - for nseqs=0, use the first 3/4 of time (dates) for Learning.set
  df_0 <- subset(df, df$nseq==0)
  learning.dates.0 <- unique(df_0[order(df_0[,"date"]), "date"])[1:N]
  Learning.set.0 <- subset(df_0, df_0$date %in% learning.dates.0)
  Test.set.0 <- subset(df_0, df_0$date %in% learning.dates.0==FALSE)
  
  # - create final Learning and Test sets
  Learning.set.cut <- subset(df, df$nseq  %in% learning.nseqs)
  Learning.set <- rbind(Learning.set.0, Learning.set.cut) 
  Test.set.cut <- subset(df, df$nseq  %in% test.nseqs) 
  Test.set <- rbind(Test.set.0, Test.set.cut)            
  # - order a by site and date
  Learning.set <- Learning.set %>% 
    arrange(site, date)
  Test.set <- Test.set %>% 
    arrange(site, date)
  
  
  #remove sites in Learning.set with log0.mortality.rel.20/log0.mortality.rel.10 all equal to NA, 
  #farms with only 1 nseq (!=0) in Learning.set and 
  #all nseqs in the Learning.set with less than 6 observations != NA
  sites.exclude.all.na <- c()
  sites.exclude.1.nseq <- c()
  nseq.to.exclue.6.obs <- c()
  var <- relevant.names[length(relevant.names)]
  
  for(farm in unique(Learning.set$site)){ 
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    n.nseqs <- unique(D.farm$nseq[D.farm$nseq!=0])
    
    if(sum(is.na(D.farm[,var]))==dim(D.farm)[1]){
      a <- unique(D.farm$site)
      sites.exclude.all.na <- c(a, sites.exclude.all.na)
      
    }else{
      
      for(i.nseq in n.nseqs){ 
        nseq.farm <- subset(D.farm, D.farm$nseq == i.nseq)
        n.not.na.nseq <- dim(nseq.farm)[1] - sum(is.na(nseq.farm[,var]))
        
        if(n.not.na.nseq<6){ #if has less than 6 observations != NA
          nseq.to.exclue.6.obs <- c(i.nseq, nseq.to.exclue.6.obs)
        }
      }
      f.nseqs <- setdiff(n.nseqs, nseq.to.exclue.6.obs)
      if(length(f.nseqs)<2){ #if has only 1 nseq - remove farm
        b <- unique(D.farm$site)
        sites.exclude.1.nseq <- c(b, sites.exclude.1.nseq)
      }
    }
  }
  # - exclude nseqs with less than 6 observations != NA on the Learning.set
  Learning.set <- subset(Learning.set, Learning.set$nseq %in% nseq.to.exclue.6.obs == FALSE) 
  
  # - exclude farms in the Learning.set with log0.mortality.rel.20/log0.mortality.rel.10 always equal to NA and 
  # - farms with only 1 nseq in Learning.set  
  all.sites.exclude <- c(sites.exclude.all.na, sites.exclude.1.nseq) 
  Learning.set <- subset(Learning.set, Learning.set$site %in% all.sites.exclude == FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% all.sites.exclude == FALSE)             
  
  
  #remove the farms that all environmental data is NA - none
  table.is.na <- c()
  for(farm in unique(Learning.set$site)){
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    a <- c()
    relevant.names <- c("d1.temp", "d9.temp", "log.max.daily.range.temp",
                        "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal",
                        "log1.d9.phyc", "log.d9.chl",
                        "log.d1.do", "log.d9.do", "log.max.daily.range.do",
                        "log.d9.prep", "log.d9.dino", "log.d9.diato",
                        "log.d9.nano", "log.d9.pico",
                        "d1.ph", "d9.ph", "max.daily.range.ph",
                        "log.d9.no3", "log.max.daily.range.no3",
                        "log0.mortality.rel.20")
    for (name in relevant.names[1:length(relevant.names)-1]){
      
      total <- dim(D.farm)[1]
      assign(paste0("is.na.", name, sep=""), sum(is.na(D.farm[,name])))
      a <- c(a, get(paste0("is.na.", name, sep="")))
      vector <- c(farm, total, a)
    }
    table.is.na <- rbind(table.is.na, vector)
    table.is.na <- as.data.frame(table.is.na)
    colnames(table.is.na)[1:2] <- c("site", "total")
    colnames(table.is.na)[3:length(colnames(table.is.na))] <- relevant.names[1:length(relevant.names)-1]
    rownames(table.is.na) <- NULL
  }
  exclude.sites.i <- which(table.is.na$total==table.is.na$d1.temp) #if 1 environmental variable is all missing, all variables are
  exclude.sites <- table.is.na$site[exclude.sites.i]
  Learning.set <- subset(Learning.set, Learning.set$site %in% exclude.sites== FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% exclude.sites == FALSE)           
  
  #have the same sites in Learning and Test sets
  sites.learn.not.test <- setdiff(Learning.set$site,Test.set$site) #sites that are in Learning.set but not in Test.set
  sites.test.not.learn <- setdiff(Test.set$site,Learning.set$site) #sites in Test.set that are not in Learning.set
  Learning.set <- subset(Learning.set, Learning.set$site %in% sites.learn.not.test == FALSE)
  Test.set <- subset(Test.set, Test.set$site %in% sites.test.not.learn == FALSE)             
  
  return(list(Learning.set=Learning.set, Test.set=Test.set))
}


# Function to create Learning and Test sets for the EM algorithm ----
# - df is the original data set
# - relevant.names are the names of the variables to be co-modeled
# - N is the time step where the datasets are divided into learning and test sets
## Difference is that doesn't remove farms on the Learning set with only 1 nseq (!=0)
learning.test.sets.EM <- function(df, relevant.names, N){
  
  # - see which nseqs != 0 are cut and move them to the closest side
  set_nseqs <- data.frame("nseq"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "start"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "end"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])),
                          "set"=rep(NA, length(unique(df$nseq)[unique(df$nseq)!=0])))
  loop=0
  for(i in unique(df$nseq)[unique(df$nseq)!=0]){
    
    dates <- unique(df[order(df[,"date"]), "date"])
    df_nseq <- subset(df, df$nseq == i)
    min_nseq <- df_nseq[1,"date"]
    max_nseq <- df_nseq[dim(df_nseq)[1],"date"]
    
    range <- which(dates==min_nseq):which(dates==max_nseq)
    a <- intersect(range,N+1)
    
    loop=loop+1
    set_nseqs[loop, "nseq"] <- i
    set_nseqs[loop, "start"] <- which(dates==min_nseq)
    set_nseqs[loop, "end"] <- which(dates==max_nseq)
    
    if(isTRUE(length(a)>0)){ 
      dist_min <- N - which(dates==min_nseq)
      dist_max <- which(dates==max_nseq) - N
      
      if(dist_min < dist_max){
        set_nseqs[loop, "set"] <- "Test"
      }
      if(dist_max < dist_min){
        set_nseqs[loop, "set"] <- "Learning"
      }
      if(dist_max == dist_min){ #if it's in the middle goes to Learning.set
        set_nseqs[loop, "set"] <- "Learning"
      }
      
    }else{#put the others in learning or test set
      
      if(isTRUE(which(dates==max_nseq) <= N)){
        set_nseqs[loop, "set"] <- "Learning"
      }else{
        set_nseqs[loop, "set"] <- "Test"
      }
    }
  }
  sum(set_nseqs$set=="Learning")/dim(set_nseqs)[1]*100   
  sum(set_nseqs$set=="Test")/dim(set_nseqs)[1]*100       
  
  learning.nseqs <- set_nseqs$nseq[set_nseqs$set== "Learning"]
  test.nseqs <- set_nseqs$nseq[set_nseqs$set== "Test"]
  
  # - for nseqs=0, use the first 3/4 of time (dates) for Learning.set
  df_0 <- subset(df, df$nseq==0)
  learning.dates.0 <- unique(df_0[order(df_0[,"date"]), "date"])[1:N]
  Learning.set.0 <- subset(df_0, df_0$date %in% learning.dates.0)
  Test.set.0 <- subset(df_0, df_0$date %in% learning.dates.0==FALSE)
  
  # - create final Learning and Test sets
  Learning.set.cut <- subset(df, df$nseq  %in% learning.nseqs)
  Learning.set <- rbind(Learning.set.0, Learning.set.cut) 
  Test.set.cut <- subset(df, df$nseq  %in% test.nseqs) 
  Test.set <- rbind(Test.set.0, Test.set.cut)             
  # - order a by site and date
  Learning.set <- Learning.set %>% 
    arrange(site, date)
  Test.set <- Test.set %>% 
    arrange(site, date)
  
  
  #remove sites in Learning.set with mortality.rel.20/log0.mortality.rel.10 all equal to NA, 
  #farms with only 1 nseq (!=0) in Learning.set and 
  #all nseqs in the Learning.set with less than 6 observations != NA
  sites.exclude.all.na <- c()
  sites.exclude.1.nseq <- c()
  nseq.to.exclue.6.obs <- c()
  var <- relevant.names[length(relevant.names)]
  
  for(farm in unique(Learning.set$site)){ 
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    n.nseqs <- unique(D.farm$nseq[D.farm$nseq!=0])
    
    if(sum(is.na(D.farm[,var]))==dim(D.farm)[1]){
      a <- unique(D.farm$site)
      sites.exclude.all.na <- c(a, sites.exclude.all.na)
      
    }else{
      
      for(i.nseq in n.nseqs){
        nseq.farm <- subset(D.farm, D.farm$nseq == i.nseq)
        n.not.na.nseq <- dim(nseq.farm)[1] - sum(is.na(nseq.farm[,var]))
        
        if(n.not.na.nseq<6){ #if has less than 6 observations != NA
          nseq.to.exclue.6.obs <- c(i.nseq, nseq.to.exclue.6.obs)
        }
      }
      f.nseqs <- setdiff(n.nseqs, nseq.to.exclue.6.obs)
      if(length(f.nseqs)<2){ #if has only 1 nseq - don't remove farm
        b <- unique(D.farm$site)
        sites.exclude.1.nseq <- c(b, sites.exclude.1.nseq)
      }
    }
  }
  # - exclude nseqs with less than 6 observations != NA on the Learning.set
  Learning.set <- subset(Learning.set, Learning.set$nseq %in% nseq.to.exclue.6.obs == FALSE) 
  
  # - exclude farms in the Learning.set with mortality.rel.20/log0.mortality.rel.10 always equal to NA
  Learning.set <- subset(Learning.set, Learning.set$site %in% sites.exclude.all.na == FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% sites.exclude.all.na == FALSE)             
  
  
  #remove the farms that all environmental data is NA - none
  table.is.na <- c()
  for(farm in unique(Learning.set$site)){
    D.farm <- subset(Learning.set, Learning.set$site == farm)
    a <- c()
    relevant.names <- c("d1.temp", "d9.temp", "log.max.daily.range.temp",
                        "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal",
                        "log1.d9.phyc", "log.d9.chl",
                        "log.d1.do", "log.d9.do", "log.max.daily.range.do",
                        "log.d9.prep", "log.d9.dino", "log.d9.diato",
                        "log.d9.nano", "log.d9.pico",
                        "d1.ph", "d9.ph", "max.daily.range.ph",
                        "log.d9.no3", "log.max.daily.range.no3",
                        "log0.mortality.rel.20")
    for (name in relevant.names[1:length(relevant.names)-1]){
      
      total <- dim(D.farm)[1]
      assign(paste0("is.na.", name, sep=""), sum(is.na(D.farm[,name])))
      a <- c(a, get(paste0("is.na.", name, sep="")))
      vector <- c(farm, total, a)
    }
    table.is.na <- rbind(table.is.na, vector)
    table.is.na <- as.data.frame(table.is.na)
    colnames(table.is.na)[1:2] <- c("site", "total")
    colnames(table.is.na)[3:length(colnames(table.is.na))] <- relevant.names[1:length(relevant.names)-1]
    rownames(table.is.na) <- NULL
  }
  exclude.sites.i <- which(table.is.na$total==table.is.na$d1.temp) #if 1 environmental variable is all missing, all variables are
  exclude.sites <- table.is.na$site[exclude.sites.i]
  Learning.set <- subset(Learning.set, Learning.set$site %in% exclude.sites== FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% exclude.sites == FALSE)            
  
  #have the same sites in Learning and Test sets
  sites.learn.not.test <- setdiff(Learning.set$site,Test.set$site) #sites that are in Learning.set but not in Test.set
  sites.test.not.learn <- setdiff(Test.set$site,Learning.set$site) #sites in Test.set that are not in Learning.set
  Learning.set <- subset(Learning.set, Learning.set$site %in% sites.learn.not.test == FALSE) 
  Test.set <- subset(Test.set, Test.set$site %in% sites.test.not.learn == FALSE)             
  
  return(list(Learning.set=Learning.set, Test.set=Test.set))
}



# Function to get the final Learning and Test sets ----
# - df is the original data set
# - relevant.names are the names of the variables to be co-modeled
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not
# - N is the time step where the datasets are divided into learning and test sets
## Use only the sites that will be there after dividing the Learning set twice
## (so that the sites in the train.set used for the EM algorithm will be the same as in the Learning.set)
get.learning.test.sets <- function(df, relevant.names, hierarchical, N){
  
  df_0 <- df[-c(which(df$nseq==0)),] 
  
  sets <- learning.test.sets(df_0, relevant.names, N)
  Learning.set_0 <- sets[["Learning.set"]]
  Test.set_0 <- sets[["Test.set"]]
  
  N.EM <- round(3*(length(unique(Learning.set_0[order(Learning.set_0[,"date"]), "date"]))/4))
  
  sets.EM <- learning.test.sets.EM(df=Learning.set_0, relevant.names, N=N.EM)
  train.set.EM_0 <- sets.EM[["Learning.set"]]
  test.set.EM_0 <- sets.EM[["Test.set"]]
  
  sites.to.remove <- setdiff(Learning.set_0$site, train.set.EM_0$site)
  sites.to.keep <- unique(train.set.EM_0$site)
  
  if(hierarchical==FALSE){
    Learning.set.final <- subset(Learning.set_0, Learning.set_0$site %in% sites.to.remove == FALSE)
    Test.set.final <- subset(Test.set_0, Test.set_0$site %in% sites.to.remove == FALSE)
  }
  if(hierarchical==TRUE){
    sets <- learning.test.sets(df, relevant.names, N)
    Learning.set <- sets[["Learning.set"]]
    Test.set <- sets[["Test.set"]]
    
    Learning.set.final <- subset(Learning.set, Learning.set$site %in% sites.to.keep == TRUE)
    Test.set.final <- subset(Test.set, Test.set$site %in% sites.to.keep == TRUE)
  }
  return(list(Learning.set=Learning.set.final, Test.set=Test.set.final))
}


################################################################################################################
### WHEN WORKING WITH HARMONIC WAVES ###

# Function for defining the LM for a data set as the sum of some trend and a number of harmonic waves ----
# - D is the data set, containing the time series we want to model
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on
# - N.w is the number of harmonic waves you want to include in your model
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
# - relevant.name is the name of the variables to be modeled
# - time.var is the variable that describes the observation time
# - stratify.by is the column name which gives the names of the unit of observation, e.g. the site name, the herd name, the pen number, the pen_batch ID, etc
# - parameters_per_stratify.by, if is TRUE the DLM parameters will be created per stratify.by per relevant.name, if is FALSE they are going o be created only per relevant.name
# - plot.it (TRUE of FALSE) determines whether plots should be made to show the fit of the model to the data of some unique identified units
# - remove.zeros whether observations in the time series with the value 0 should be removed. This is not usually recommended, unless 0s are stand-ins for NAs
# - round.by is the number of decimals you want to round the estimates
# - silent determines if the function should print the r^2 values of the fitted model to the console
# - useInteractions (TRUE of FALSE) whether you want to use time as an interaction 
make.LM.wHarmonics <- function(D, time.in.year, N.w, trend, relevant.name, time.var, stratify.by, parameters_per_stratify.by, plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE, useInteractions=FALSE){
  
  w <- (2*pi)/time.in.year 
  
  # Start by plotting the average level per period and over longer time scales
  if(length(which(is.na(D[,relevant.name]))) > 0){
    D <- D[-which(is.na(D[,relevant.name])),]
  }
  
  par(mfrow=c(1,2))
  
  # - per period
  period <- (2*pi)/w
  D$time.short <- D[,time.var]%%period
  agg.short <- aggregate(D[,relevant.name], by=list(D$time.short), FUN=mean)
  plot(agg.short, type='l', xlab='Time in one period', ylab=relevant.name)
  
  # - over longer time
  agg.long <- aggregate(D[,relevant.name], by=list(D[,time.var]), FUN=mean)
  plot(agg.long, type='l', xlab=time.var, ylab=relevant.name)
  
  par(mfrow=c(1,1))
  
  # Required packages
  library('MuMIn')
  lmer <- lme4::lmer 
  
  # Remove zero-values, if relevant
  if(remove.zeros == TRUE){
    remove.i <- which(D[,relevant.name] == 0)
    D <- D[-remove.i,]
  }
  
  # We make an empty vector, into which we will add the relevant text for the linear function
  lm.vector <- c()
  
  # We also make an empty vector, into which we will add the names for the parameter vector
  names.vector <- c()
  
  # We add the target variable (the variable we model with this collection of harmonics)
  lm.vector <- c(lm.vector, relevant.name, '~')
  names.vector <- c(names.vector, 'Level')
  
  # Make the linear trend component, if relevant
  if(trend == TRUE){
    lm.vector <- c(lm.vector, time.var)
    names.vector <- c(names.vector, 'Trend +')
  }
  
  # Now we add the relevant number of harmonics, one at a time
  if(N.w > 0){
    
    if(trend == TRUE){
      lm.vector <- c(lm.vector, '+')
    }
    
    for(i in 1:N.w){
      if(i == 1){
        a <- paste('cos(', i, '*w*', time.var, ') + sin(', i, '*w*', time.var, ')')
        lm.vector <- c(lm.vector, a)
        names.vector <- c(names.vector, paste('wave.', i, '.cos', sep=''), paste('wave.', i, '.sin', sep=''))
      }else{
        a <- paste('+ cos(', i, '*w*', time.var, ') + sin(', i, '*w*', time.var, ')')
        lm.vector <- c(lm.vector, a)
        names.vector <- c(names.vector, paste('wave.', i, '.cos', sep=''), paste('wave.', i, '.sin', sep=''))
      }
    }
  }
  
  # We can also add the interactions between the waves and the time
  # - this should NOT be used for the DLM, but may be relevant to assess the fit of the model!
  if(useInteractions == TRUE){ #use time as an interaction
    if(N.w > 0){
      for(i in 1:N.w){
        if(i < N.w){
          a <- paste('+ cos(', i, '*w*', time.var, '):', time.var, '+' , sep='')
          b <- paste('sin(', i, '*w*', time.var, '):', time.var, sep='')
        }else{
          a <- paste('cos(', i, '*w*', time.var, '):', time.var, '+' , sep='')
          b <- paste('sin(', i, '*w*', time.var, '):', time.var , sep='')
        }
        
        lm.vector <- c(lm.vector, a, b)
      }
    }
  }
  
  # We can now make the final name vector
  names.vector <- paste(relevant.name, '_', names.vector, sep='')
  
  # Make sure we can handle the case with no trend AND no waves
  if(length(lm.vector) == 2){
    lm.vector <- c(lm.vector, 1)
  }
  
  # Lastly, we add the random effect variable
  a <- paste( ' + (1|', stratify.by, ')', sep='')
  lm.vector_LME <- c(lm.vector, a)
  
  # Now we make it into a linear function - without random effect
  b1 <- as.formula(paste(lm.vector, collapse=' '))
  lm.1 <- lm(b1, data = D)
  S1 <- summary(lm.1)
  r.squared <- S1$r.squared 
  adj.r.squared <- S1$adj.r.squared
  
  # Now we make it into a linear function - with random effect
  if(isTRUE(parameters_per_stratify.by)){
    b2 <- as.formula(paste(lm.vector_LME, collapse=' '))
    lm.2 <- lmer(b2, data = D)
    S2 <- summary(lm.2)
    r.sqrd.out <- r.squaredGLMM(lm.2) 
  }else{
    r.sqrd.out <- cbind()
  }
  
  # - now we directly have the initial parameter vector
  if(isTRUE(parameters_per_stratify.by)){
    mu0 <- matrix(S2$coefficients[,'Estimate'])
  }else{
    mu0 <- matrix(S1$coefficients[,'Estimate'])
  }
  mu0 <- matrix(mu0[1:length(names.vector)])
  rownames(mu0) <- names.vector
  
  # - and we directly have the initial prior variance matrix
  if(isTRUE(parameters_per_stratify.by)){
    C0 <- as.matrix(vcov(lm.2))
  }else{
    C0 <- as.matrix(vcov(lm.1))
  }
  C0 <- C0[1:length(names.vector),1:length(names.vector)]
  
  # Present the fit of this model
  r.sqrd.out <- cbind(r.sqrd.out, 'r.squared'=r.squared, 'adj.r.squared'=adj.r.squared)
  r.sqrd.out <- round(r.sqrd.out, round.by)
  if(!silent){
    print(r.sqrd.out)
  }
  
  if(plot.it == TRUE){
    N.plots <- min(4,length(unique(D[,stratify.by])))
    par(mfrow=c( max( (N.plots/2), 1) , 2))
    set.seed(42)
    random.IDs <- sample(x = unique(D[,stratify.by]), size = min(4,length(unique(D[,stratify.by]))), replace = FALSE)
    for(ID in random.IDs){ #random.IDs=1
      ID.set <- subset(D, D[,stratify.by] == ID)
      pred <- predict(object = lm.1, newdata=ID.set)
      plot(ID.set[,relevant.name]~ID.set[,time.var], xlab=time.var, ylab = relevant.name, main=ID)
      lines(pred~ID.set[,time.var], col='red')
      
      Residuals <- ID.set[,relevant.name] - pred
      hist(Residuals)
      abline(v=0, col='red', lwd=2)
    }
  }
  
  # Return the mu0 and C0 for use in the DLM
  return(list('mu0'=mu0,
              'C0'=C0,
              'adj.r.squared'=round(adj.r.squared, round.by),
              'LM'=lm.1))
  
}



# A function to see the best number of waves for each variable ----
# - D is the data set, containing the time series we want to model
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on. 
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
# - relevant.name is the name of the variables to be modeled
# - time.var is the variable that describes the observation time
# - stratify.by is the column name which gives the names of the unit of observation, e.g. the site name, the herd name, the pen number, the pen_batch ID, etc
# - parameters_per_stratify.by, if is TRUE the DLM parameters will be created per stratify.by per relevant.name, if is FALSE they are going o be created only per relevant.name
# - plot.it (TRUE of FALSE) determines whether plots should be made to show the fit of the model to the data of some unique identified units
# - remove.zeros whether observations in the time series with the value 0 should be removed. This is not usually recommended, unless 0s are stand-ins for NAs
# - round.by is the number of decimals you want to round the estimates
get.best.number.waves <- function(D, time.in.year=12, trend = FALSE, relevant.name, time.var, stratify.by, parameters_per_stratify.by, plot.it=TRUE, remove.zeros=FALSE, round.by=2){
  
  w <- (2*pi)/time.in.year
  
  # Find the optimal number of waves for the model
  start.time <- Sys.time()
  N.w = -1 
  adj.r.squared.best <- -1
  no.better <- FALSE
  while(no.better == FALSE){
    N.w <- N.w + 1
    out <- make.LM.wHarmonics(D, time.in.year, N.w, trend, relevant.name, time.var, stratify.by, parameters_per_stratify.by, plot.it=TRUE, remove.zeros=FALSE, round.by=2, silent=TRUE, useInteractions=FALSE)
    adj.r.squared <- out$adj.r.squared
    if(adj.r.squared > adj.r.squared.best){
      adj.r.squared.best <- adj.r.squared
      N.w.best <- N.w
    }else{
      no.better <- TRUE
    }
  }
  return(list('Best N.w'=N.w.best,
              'adj.r.squared'=adj.r.squared))
}


################################################################################################################
### MODELLING ###

# Function to estimate the Gt (system matrix) ----
# - D is the data set, containing the time series we want to model
# - i is the time step you are modelling in the DLM, you can also give it equal to NA if you want to use dates instead
# - Date is the date you are currently modelling in the DLM
# - Date_1 is the previously date to what you are currently modelling in the DLM
# - time.var is the variable that describes the observation time
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on. 
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not
# - level2 is the second level (here correspond to the region)
# - level3 is the third level (here correspond to the farm/site)
# - best.n.waves is a vector with the best number of waves for each variable
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
# - relevant.names are the names of the variables to be co-modeled
# - Spline.list is a list with the splines calculated for each farm
getGt <- function(D, i, Date=NA, Date_1=NA, time.var, time.in.year, hierarchical, level2="local.authority", 
                  level3="site", best.n.waves, trend, relevant.names, Spline.list){
  
  if(hierarchical==FALSE){
    
    D.full <- D
    
    # Prepare for making Gt matrix
    loop=0
    end <- c()
    n.runs <- 1
    n.zeros <- length(which(best.n.waves==0))
    n.na <- length(which(is.na(best.n.waves)))
    n.nyquist <- length(which((best.n.waves==time.in.year/2)))
    n.rows.cols <- sum(best.n.waves*2 + 1, na.rm=T) - n.nyquist + n.zeros + n.na*2
    
    time.var.all <- time.var
    
    # Make the diagonal of the Gt matrix
    Gt <- as.data.frame(diag(1, ncol = n.rows.cols, nrow = n.rows.cols))
    
    for(N.w in best.n.waves){
      
      loop=loop+1
      
      if(isTRUE(is.na(N.w))){
        
        #get the name we are referring to
        name <- relevant.names[loop]
        
        #get the time.var for the corresponding name
        time.var=time.var.all[loop]
        
        #where to add the trend on the Gt matrix
        if(loop==1){
          start = n.runs
        }else{
          start=end + 1
        }
        end <- start + 1
        
        #get the relevant spline function
        spline.name <- paste(name, '_Spline', sep='') #1 spline per relevant.name
        Spline <- Spline.list[[which(names(Spline.list)==spline.name)]]
        
        time <- D[i,time.var] #The observation time of the current observations
        time_1 <- D[i-1,time.var] #The observation time of the previous observations
        if(length(time_1) == 0){
          time_1 <- 0
        }
        
        #get trend for Gt matrix
        pred <- predict(object = Spline, x = time)$y
        pred_1 <- predict(object = Spline, x = time_1)$y
        Trend <- pred - pred_1
        
        Gt[start,end] <- Trend
        
      }else{
        
        if(N.w!=0){
          
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #define the base frequency
          w <- (2*pi)/time.in.year
          
          #get Gt matrix using HW for name
          if(trend==TRUE){
            G <- diag(1, 2+2*N.w)
            G[1,2] <- 1
            wave.i <- 3
          }else{
            G <- diag(1, 1+2*N.w)
            wave.i <- 2
          }
          count <- 0
          wave.elements <- c(cos(w), sin(w), -sin(w), cos(w))
          for(n in 1:N.w){
            G[wave.i,wave.i] <- cos(n*w)
            G[wave.i,wave.i+1] <- sin(n*w)
            G[wave.i+1,wave.i] <- -sin(n*w)
            G[wave.i+1,wave.i+1] <- cos(n*w)
            wave.i <- wave.i + 2
          }
          
          #add the Gt matrix for name to the final Gt matrix
          if(N.w!=time.in.year/2){
            if(loop==1){
              start = n.runs
            }else{
              start=end + 1
            }
            end <- start + 2*N.w
          }
          if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
            
            G <- G[1:(dim(G)[1]-1), 1:(dim(G)[1]-1)]
            G[dim(G)[1],dim(G)[1]] <- -1
            
            if(loop==1){
              start = n.runs
            }else{
              start=end + 1
            }
            end <- start + 2*N.w-1
          }
          
          Gt[start:end,start:end] <- G
        }
        
        
        if(N.w==0){
          
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #get the time.var for the corresponding name
          time.var=time.var.all[loop]
          
          #where to add the trend on the Gt matrix
          if(loop==1){
            start = n.runs
          }else{
            start=end + 1
          }
          end <- start + 1
          
          #get the relevant spline function
          spline.name <- paste(name, '_Spline', sep='') 
          Spline <- Spline.list[[which(names(Spline.list)==spline.name)]]
          
          time <- D[i,time.var] #The observation time of the current observations
          time_1 <- D[i-1,time.var] #The observation time of the previous observations
          if(length(time_1) == 0){
            time_1 <- 0
          }
          
          #get trend for Gt matrix
          pred <- predict(object = Spline, x = time)$y
          pred_1 <- predict(object = Spline, x = time_1)$y
          Trend <- pred - pred_1
          
          Gt[start,end] <- Trend
        }
      }
    }
  }
  
  if(hierarchical==TRUE){
    
    D.full <- D
    
    # Prepare for making Gt matrix
    loop=0
    end <- c()
    n.runs <- 1
    n.zeros <- length(which(best.n.waves==0))
    n.na <- length(which(is.na(best.n.waves)))
    n.nyquist <- length(which((best.n.waves==time.in.year/2)))
    n.rows.cols <- sum(best.n.waves*2 + 1, na.rm=T) - n.nyquist + n.zeros + n.na*2
    #since is hierarchical:
    n.rows.cols <- n.rows.cols + (n.rows.cols*length(unique(D.full[,level2]))) + (n.rows.cols*length(unique(D.full[,level3])))
    
    time.var.all <- time.var
    
    # Make the diagonal of the Gt matrix
    Gt <- as.data.frame(diag(1, ncol = n.rows.cols, nrow = n.rows.cols))
    
    
    #### FOR COUNTRY
    
    for(N.w in best.n.waves){
      
      loop=loop+1
      
      if(isTRUE(is.na(N.w))){
        
        #get the name we are referring to
        name <- relevant.names[loop] 
        
        #get the time.var for the corresponding name
        time.var=time.var.all[loop]
        
        #where to add the trend on the Gt matrix
        if(loop==1){
          start = n.runs
        }else{
          start=end + 1
        }
        end <- start + 1
        
        Gt[start,end] <- 1 #because it doesn't make sense to use splines for predicting the mortality of the country
        
      }else{
        
        if(N.w!=0){
          
          #get the name we are referring to
          name <- relevant.names[loop] 
          
          #define the base frequency
          w <- (2*pi)/time.in.year  
          
          #get Gt matrix using HW for name
          if(trend==TRUE){
            G <- diag(1, 2+2*N.w)
            G[1,2] <- 1
            wave.i <- 3
          }else{
            G <- diag(1, 1+2*N.w)
            wave.i <- 2
          }
          count <- 0
          wave.elements <- c(cos(w), sin(w), -sin(w), cos(w))
          for(n in 1:N.w){
            G[wave.i,wave.i] <- cos(n*w)
            G[wave.i,wave.i+1] <- sin(n*w)
            G[wave.i+1,wave.i] <- -sin(n*w)
            G[wave.i+1,wave.i+1] <- cos(n*w)
            wave.i <- wave.i + 2
          }
          
          #add the Gt matrix for name to the final Gt matrix
          if(N.w!=time.in.year/2){
            if(loop==1){
              start = n.runs
            }else{
              start=end + 1
            }
            end <- start + 2*N.w
          }
          if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
            
            G <- G[1:(dim(G)[1]-1), 1:(dim(G)[1]-1)]
            G[dim(G)[1],dim(G)[1]] <- -1
            
            if(loop==1){
              start = n.runs
            }else{
              start=end + 1
            }
            end <- start + 2*N.w-1
          }
          
          Gt[start:end,start:end] <- G
        }
        
        if(N.w==0){   
          
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #where to add 1 on the Gt matrix
          if(loop==1){
            start = n.runs
          }else{
            start=end + 1
          }
          end <- start + 1
          
          Gt[start,end] <- 1}
      }
    }
    
    #### FOR REGION
    
    for(region in sort(unique(D.full[,level2]))){ 
      
      D.region <- subset(D.full, D.full[,level2] == region)
      loop=0
      n.runs <- end + 1
      
      for(N.w in best.n.waves){
        
        loop=loop+1
        
        if(isTRUE(is.na(N.w))){
          
          #get the right data
          D <- D.region
          
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #get the time.var for the corresponding name
          time.var=time.var.all[loop]
          
          #where to add the trend on the Gt matrix
          if(loop==1){
            start = n.runs
          }else{
            start=end + 1
          }
          end <- start + 1
          
          Gt[start,end] <- 1 #because it doesn't make sense to use splines for predicting the mortality of the region
          
        }else{
          
          if(N.w!=0){
            
            #get the name we are referring to
            name <- relevant.names[loop] 
            
            #define the base frequency
            w <- (2*pi)/time.in.year  
            
            #get Gt matrix using HW for name
            if(trend==TRUE){
              G <- diag(1, 2+2*N.w)
              G[1,2] <- 1
              wave.i <- 3
            }else{
              G <- diag(1, 1+2*N.w)
              wave.i <- 2
            }
            count <- 0
            wave.elements <- c(cos(w), sin(w), -sin(w), cos(w))
            for(n in 1:N.w){
              G[wave.i,wave.i] <- cos(n*w)
              G[wave.i,wave.i+1] <- sin(n*w)
              G[wave.i+1,wave.i] <- -sin(n*w)
              G[wave.i+1,wave.i+1] <- cos(n*w)
              wave.i <- wave.i + 2
            }
            
            #add the Gt matrix for name to the final Gt matrix
            if(N.w!=time.in.year/2){
              if(loop==1){
                start = n.runs
              }else{
                start=end + 1
              }
              end <- start + 2*N.w
            }
            if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
              
              G <- G[1:(dim(G)[1]-1), 1:(dim(G)[1]-1)]
              G[dim(G)[1],dim(G)[1]] <- -1
              
              if(loop==1){
                start = n.runs
              }else{
                start=end + 1
              }
              end <- start + 2*N.w-1
            }
            
            Gt[start:end,start:end] <- G
          }
          
          if(N.w==0){
            
            #get the name we are referring to
            name <- relevant.names[loop]
            
            #where to add 1 on the Gt matrix
            if(loop==1){
              start = n.runs
            }else{
              start=end + 1
            }
            end <- start + 1
            
            Gt[start,end] <- 1
          }
        }
      }
    }
    
    #### FOR SITE
    
    for(farm in sort(unique(D.full[,level3]))){
      
      D.farm <- subset(D.full, D.full[,level3] == farm)
      loop=0
      n.runs <- end + 1
      
      for(N.w in best.n.waves){
        
        loop=loop+1
        
        if(isTRUE(is.na(N.w))){
          
          #get the right data
          D <- D.farm
          
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #get the time.var for the corresponding name
          time.var=time.var.all[loop]
          
          #where to add the trend on the Gt matrix
          if(loop==1){
            start = n.runs
          }else{
            start=end + 1
          }
          end <- start + 1
          
          #get the relevant spline function
          spline.name <- paste(farm, name, 'Spline', sep='_') 
          Spline <- Spline.list[[which(names(Spline.list)==spline.name)]]
          
          i <- which(unique(D$date) == Date)
          i_1 <- which(unique(D$date) == Date_1)
          
          if(isTRUE(length(i)!=0)){
            
            time <- D[i,time.var] #The observation time of the current observations
            time_1 <- D[i_1,time.var] #The observation time of the previous observations
            if(length(time_1) == 0){
              time_1 <- 0
            }
            
            #get trend for Gt matrix
            if(!is.na(time)){
              pred <- predict(object = Spline, x = time)$y
              if(isTRUE(time != time_1)){
                pred_1 <- predict(object = Spline, x = time_1)$y
                Trend <- pred - pred_1
              }else{
                Trend <- pred 
              }
            }else{
              Trend <- 0
            }
            
            Gt[start,end] <- Trend
          }else{
            Gt[start,end] <- 0
          }
          
        }else{
          
          if(N.w!=0){
            
            #get the name we are referring to
            name <- relevant.names[loop] 
            
            #define the base frequency
            w <- (2*pi)/time.in.year
            
            #get Gt matrix using HW for name
            if(trend==TRUE){
              G <- diag(1, 2+2*N.w)
              G[1,2] <- 1
              wave.i <- 3
            }else{
              G <- diag(1, 1+2*N.w)
              wave.i <- 2
            }
            count <- 0
            wave.elements <- c(cos(w), sin(w), -sin(w), cos(w))
            for(n in 1:N.w){
              G[wave.i,wave.i] <- cos(n*w)
              G[wave.i,wave.i+1] <- sin(n*w)
              G[wave.i+1,wave.i] <- -sin(n*w)
              G[wave.i+1,wave.i+1] <- cos(n*w)
              wave.i <- wave.i + 2
            }
            
            #add the Gt matrix for name to the final Gt matrix
            if(N.w!=time.in.year/2){
              if(loop==1){
                start = n.runs
              }else{
                start=end + 1
              }
              end <- start + 2*N.w
            }
            if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
              
              G <- G[1:(dim(G)[1]-1), 1:(dim(G)[1]-1)] #last trend is only -1
              G[dim(G)[1],dim(G)[1]] <- -1
              
              if(loop==1){
                start = n.runs
              }else{
                start=end + 1
              }
              end <- start + 2*N.w-1
            }
            
            Gt[start:end,start:end] <- G
          }
          
          if(N.w==0){
            
            #get the right data
            D <- D.farm
            
            #get the name we are referring to
            name <- relevant.names[loop]
            
            #get the time.var for the corresponding name
            time.var=time.var.all[loop]
            
            #where to add the trend on the Gt matrix
            if(loop==1){
              start = n.runs
            }else{
              start=end + 1
            }
            end <- start + 1
            
            #get the relevant spline function
            spline.name <- paste(farm, name, 'Spline', sep='_')
            Spline <- Spline.list[[which(names(Spline.list)==spline.name)]]
            
            i <- which(unique(D$date) == Date)
            i_1 <- which(unique(D$date) == Date_1)
            
            if(isTRUE(length(i)!=0)){
              
              time <- D[i,time.var] #The observation time of the current observations
              time_1 <- D[i_1,time.var] #The observation time of the previous observations
              if(length(time_1) == 0){
                time_1 <- 0
              }
              
              #get trend for Gt matrix
              if(!is.na(time)){
                pred <- predict(object = Spline, x = time)$y
                if(isTRUE(time != time_1)){
                  pred_1 <- predict(object = Spline, x = time_1)$y
                  Trend <- pred - pred_1
                }else{
                  Trend <- pred 
                }
              }else{
                Trend <- 0
              }
              
              Gt[start,end] <- Trend
            }else{
              Gt[start,end] <- 0
            }
          }
        }
      }
    }
  }
  return(as.matrix(Gt))
}



# Function to estimate the Ft (design matrix) ----
# - D is the data set, containing the time series we want to model
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on. 
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not
# - level2 is the second level (here correspond to the region)
# - level3 is the third level (here correspond to the farm/site)
# - best.n.waves is a vector with the best number of waves for each variable
# - relevant.names are the names of the variables to be co-modeled
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
getFt <- function(D, time.in.year, hierarchical, level2, level3, best.n.waves, relevant.names, trend){
  
  if(hierarchical==FALSE){
    
    # Prepare for making Ft matrix
    loop=0
    end <- c()
    n.runs <- 1
    n.zeros <- length(which(best.n.waves==0))
    n.na <- length(which(is.na(best.n.waves)))
    n.nyquist <- length(which((best.n.waves==time.in.year/2)))
    n.rows.cols <- sum(best.n.waves*2 + 1, na.rm=T) - n.nyquist + n.zeros + n.na*2
    
    # Make the diagonal of the Ft matrix
    Ft <- as.data.frame(matrix(0, ncol = length(relevant.names), nrow = n.rows.cols))    
    
    for(N.w in best.n.waves){
      
      loop=loop+1
      
      if(isTRUE(is.na(N.w))){
        
        #get the name we are referring to
        name <- relevant.names[loop]
        
        #get Fti matrix for name
        Fti <- c(1,0)
        
        #add Fti to Ft matrix
        if(loop==1){
          start=n.runs
        }else{
          start=end + 1
        }
        end <- start + 1
        
        Ft[start:end,loop] <- Fti
        colnames(Ft)[loop] <- name
        
      }else{
        
        if(N.w!=0){
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #get Fti matrix for name
          if(trend == TRUE){
            Fti <- c(1,0)
          }else{
            Fti <- c(1)
          }
          for(n in 1:N.w){
            Fti <- c(Fti, c(1,0))
          }
          
          #add Fti to Ft matrix
          if(N.w!=time.in.year/2){
            
            if(loop==1){
              start = n.runs
            }else{
              start = end + 1
            }
            end <- start + 2*N.w
          }
          if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
            
            Fti <- Fti[1:(length(Fti)-1)]
            
            if(loop==1){
              start = n.runs
            }else{
              start = end + 1
            }
            end <- start + 2*N.w-1
          }
          
          Ft[start:end,loop] <- Fti
          colnames(Ft)[loop] <- name
        }
        
        if(N.w==0){
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #get Fti matrix for name
          Fti <- c(1,0)
          
          #add Fti to Ft matrix
          if(loop==1){
            start=n.runs
          }else{
            start=end + 1
          }
          end <- start + 1
          
          Ft[start:end,loop] <- Fti
          colnames(Ft)[loop] <- name
        }
      }
    }
  }
  
  
  if(hierarchical==TRUE){
    
    D.full <- D
    
    # Prepare for making Ft matrix
    loop=0
    end <- c()
    n.runs <- 1
    n.runs.cols = 0
    n.zeros <- length(which(best.n.waves==0))
    n.na <- length(which(is.na(best.n.waves)))
    n.nyquist <- length(which((best.n.waves==time.in.year/2)))
    n.rows <- sum(best.n.waves*2 + 1, na.rm=T) - n.nyquist + n.zeros + n.na*2
    #since is hierarchical:
    n.rows <- n.rows + (n.rows*length(unique(D.full[,level2]))) + (n.rows*length(unique(D.full[,level3])))
    n.cols <- length(relevant.names)+(length(relevant.names)*length(unique(D.full[,level2])))+(length(relevant.names)*length(unique(D.full[,level3])))
    
    # Make the diagonal of the Ft matrix
    Ft <- as.data.frame(matrix(0, ncol = n.cols, nrow = n.rows))
    
    #### FOR COUNTRY
    
    for(N.w in best.n.waves){
      
      loop=loop+1
      n.runs.cols <- n.runs.cols+1 
      
      if(isTRUE(is.na(N.w))){
        
        #get the name we are referring to
        name <- relevant.names[loop]
        
        #get Fti matrix for name
        Fti <- c(1,0)
        
        #add Fti to Ft matrix
        if(loop==1){
          start=n.runs
        }else{
          start=end + 1
        }
        end <- start + 1
        
        Ft[start:end,n.runs.cols] <- Fti
        colnames(Ft)[n.runs.cols] <-  paste("Country", name, sep='.')
        
      }else{
        
        if(N.w!=0){
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #get Fti matrix for name
          if(trend == TRUE){
            Fti <- c(1,0)
          }else{
            Fti <- c(1)
          }
          for(n in 1:N.w){
            Fti <- c(Fti, c(1,0))
          }
          
          #add Fti to Ft matrix
          if(N.w!=time.in.year/2){
            
            if(loop==1){
              start = n.runs
            }else{
              start = end + 1
            }
            end <- start + 2*N.w
          }
          if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
            
            Fti <- Fti[1:(length(Fti)-1)] 
            
            if(loop==1){
              start = n.runs
            }else{
              start = end + 1
            }
            end <- start + 2*N.w-1
          }
          
          Ft[start:end,n.runs.cols] <- Fti
          colnames(Ft)[n.runs.cols] <- paste("Country", name, sep='.')
        }
        
        if(N.w==0){
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #get Fti matrix for name
          Fti <- c(1,0)
          
          #add Fti to Ft matrix
          if(loop==1){
            start=n.runs
          }else{
            start=end + 1
          }
          end <- start + 1
          
          Ft[start:end,n.runs.cols] <- Fti
          colnames(Ft)[n.runs.cols] <-  paste("Country", name, sep='.')
        }
      }
    }
    
    #### FOR REGION
    
    for(region in sort(unique(D.full[,level2]))){
      
      loop=0
      n.runs <- end + 1
      
      for(N.w in best.n.waves){
        
        loop=loop+1
        n.runs.cols <- n.runs.cols+1 
        
        if(isTRUE(is.na(N.w))){
          
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #get Fti matrix for name
          Fti <- c(1,0)
          
          #add Fti to Ft matrix
          if(loop==1){
            start=n.runs
          }else{
            start=end + 1
          }
          end <- start + 1
          
          Ft[start:end,n.runs.cols] <- Fti
          colnames(Ft)[n.runs.cols] <- paste(region, name, sep='.')
          
        }else{
          
          if(N.w!=0){
            #get the name we are referring to
            name <- relevant.names[loop]
            
            #get Fti matrix for name
            if(trend == TRUE){
              Fti <- c(1,0)
            }else{
              Fti <- c(1)
            }
            for(n in 1:N.w){
              Fti <- c(Fti, c(1,0))
            }
            
            #add Fti to Ft matrix
            if(N.w!=time.in.year/2){
              
              if(loop==1){
                start = n.runs
              }else{
                start = end + 1
              }
              end <- start + 2*N.w
            }
            if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
              
              Fti <- Fti[1:(length(Fti)-1)]
              
              if(loop==1){
                start = n.runs
              }else{
                start = end + 1
              }
              end <- start + 2*N.w-1
            }
            
            Ft[start:end,n.runs.cols] <- Fti
            colnames(Ft)[n.runs.cols] <- paste(region, name, sep='.')
          }
          
          if(N.w==0){
            #get the name we are referring to
            name <- relevant.names[loop]
            
            #get Fti matrix for name
            Fti <- c(1,0)
            
            #add Fti to Ft matrix
            if(loop==1){
              start=n.runs
            }else{
              start=end + 1
            }
            end <- start + 1
            
            Ft[start:end,n.runs.cols] <- Fti
            colnames(Ft)[n.runs.cols] <- paste(region, name, sep='.')
          }
        }
      }
    }
    
    #### FOR SITE
    
    for(farm in sort(unique(D.full[,level3]))){
      
      loop=0
      n.runs <- end + 1
      
      for(N.w in best.n.waves){
        
        loop=loop+1
        n.runs.cols <- n.runs.cols+1 
        
        if(isTRUE(is.na(N.w))){
          
          #get the name we are referring to
          name <- relevant.names[loop]
          
          #get Fti matrix for name
          Fti <- c(1,0)
          
          #add Fti to Ft matrix
          if(loop==1){
            start=n.runs
          }else{
            start=end + 1
          }
          end <- start + 1
          
          Ft[start:end,n.runs.cols] <- Fti
          colnames(Ft)[n.runs.cols] <- paste(farm, name, sep='.')
          
        }else{
          
          if(N.w!=0){
            #get the name we are referring to
            name <- relevant.names[loop]
            
            #get Fti matrix for name
            if(trend == TRUE){
              Fti <- c(1,0)
            }else{
              Fti <- c(1)
            }
            for(n in 1:N.w){
              Fti <- c(Fti, c(1,0))
            }
            
            #add Fti to Ft matrix
            if(N.w!=time.in.year/2){
              
              if(loop==1){
                start = n.runs
              }else{
                start = end + 1
              }
              end <- start + 2*N.w
            }
            if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
              
              Fti <- Fti[1:(length(Fti)-1)]
              
              if(loop==1){
                start = n.runs
              }else{
                start = end + 1
              }
              end <- start + 2*N.w-1
            }
            
            Ft[start:end,n.runs.cols] <- Fti
            colnames(Ft)[n.runs.cols] <- paste(farm, name, sep='.')
          }
          
          if(N.w==0){
            #get the name we are referring to
            name <- relevant.names[loop]
            
            #get Fti matrix for name
            Fti <- c(1,0)
            
            #add Fti to Ft matrix
            if(loop==1){
              start=n.runs
            }else{
              start=end + 1
            }
            end <- start + 1
            
            Ft[start:end,n.runs.cols] <- Fti
            colnames(Ft)[n.runs.cols] <- paste(farm, name, sep='.')
          }
        }
      }
    }
  }
  return(as.matrix(Ft))
}



# Function to update mu0 ----
## when a dataset does not start in January we have to update the variables that use harmonic waves in mu0 for the correspondent month
# - mu is the initial parameter vector
# - Gt is the  system matrix
# - start.time is the month when the data starts
update.mu0 <- function(mu, Gt, start.time){
  for(i in 1:(start.time-1)){
    mu <- Gt %*% mu
  }
  return(mu)
}



# Function to estimate the VSum for EM algorithm (non-hierarchical) ----
# - D is the data set, containing the time series we want to model
# - n is the time step
## Get a matrix with 1 in the cells where the observation 
## contributes to the observation variance-covariance matrix.
## Other cells are 0 (When it's NA). Used by the EM-algorithm
## for non-hierarchical model
getVSumElement = function(D, n) {
  row = D[n, ]
  rem = c()
  V = matrix(1, nrow=length(relevant.names), ncol=length(relevant.names))
  for(cont.var in relevant.names){
    if (is.na(row[cont.var])) { 
      cont.var.index <- which(relevant.names == cont.var)
      V[cont.var.index, ] = 0
      V[, cont.var.index] = 0
    }
  }
  return(V)
}



# Function for moving average ----
# - x is a vector with the observations
# - n is the moving window
# - FUN is the function to use (usually=mean)
moving.function <- function(x, n, FUN){
  start <- floor(n/2)+1
  N <- length(x) - floor(n/2)
  out <- c()
  if(N > start){
    for(i in start:N){
      obs <- x[(i-floor(n/2)):(i+floor(n/2))]
      res <- FUN(na.omit(obs))
      out <- c(out,res)
    }
    NAs <- rep(NA, floor(n/2))
    out <- c(NAs, out, NAs)
  }else{
    out <- rep(NA, length(x))
  }
  return(out)
}



# Function for the DLM for non-hierarchical models ----
# - D is the data set, containing the time series we want to model
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - V is the observation variance
# - W is the system variance
# - time.var is the variable that describes the observation time
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on. 
# - best.n.waves is a vector with the best number of waves for each variable
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
# - relevant.names are the names of the variables to be co-modeled
# - Spline.list is a list with the splines calculated for each farm
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not, here is hierarchical=FALSE
# - level2 is the second level (correspond to the region, here is NULL because is for non-hierarchical models)
# - level3 is the third level (correspond to the farm/site, here is NULL because is for non-hierarchical models)
FUNFILTERG <- function(D, mu0, C0, V, W, time.var, time.in.year, best.n.waves, trend, relevant.names, 
                       Spline.list, hierarchical, level2=NULL, level3=NULL){
  
  n <- nrow(D)
  
  Yt.list <- list()
  at.list <- list()		 
  Rt.list <- list()
  ft.list <- list()
  Qt.list <- list()
  At.list <- list()
  et.list <- list()
  ut.list <- list()
  mt.list <- list()
  Ct.list <- list()
  Ft.list <- list()
  VSE.list <- list()
  Gt.list <- list()
  fullFt.list <- list()
  
  # Update mu0 for the variables that have harmonic waves if the production cycle does not start in January
  if(sum(!is.na(best.n.waves) & best.n.waves!=0)>0){
    if(month(D$date[1])!=1){
      m <- month(D$date[1])
      Gt <- getGt(D, i=1, Date=NA, Date_1=NA, time.var, time.in.year, hierarchical, level2, 
                  level3, best.n.waves, trend, relevant.names, Spline.list)
      mu0.all <- update.mu0(mu = mu0, Gt = Gt, start.time = m)
      rownames(mu0.all) <- rownames(mu0)
      
      # - which rownames in mu0 should be replaced (harmonic waves) 
      names <- relevant.names[which(!is.na(best.n.waves) & best.n.waves!=0)]
      r.names.all <- c()
      for (name in names){
        a <- grepl(name, rownames(mu0.all))
        r.names <- rownames(mu0.all)[a]
        r.names.all <- c(r.names.all, r.names)
      }
      
      # - replace mu0 of the harmonic waves variables
      mu0[r.names.all,] <- mu0.all[r.names.all,]
    }
  }
  
  mt <- mu0				
  Ct <- C0
  
  # Make sure Ct is symmetrical
  Ct <- (Ct + t(Ct))/2
  
  for(i in (1:n)){ 
    
    # Get V-sum-element, to be used in the EM-algorithm
    VSE <- getVSumElement(D, i)
    
    # Get the observation vector
    Yt <- t(as.matrix(D[i,relevant.names]))
    colnames(Yt) <- NULL
    
    # Get the observational variance (Vt)  
    Vt <- V
    
    # Get the system variance (Wt)
    Wt <- W
    
    # Get Gt - independent of the current Yt
    Gt <- getGt(D, i, Date=NA, Date_1=NA, time.var, time.in.year, hierarchical, level2, 
                level3, best.n.waves, trend, relevant.names, Spline.list)
    
    # Get Ft given the current Yt 
    Ft <- getFt(D, time.in.year, hierarchical, level2, level3, best.n.waves, relevant.names, trend)
    
    # Remove missing values from Yt, Ft and Vt
    missing = which(is.na(Yt))
    fullFt <- Ft
    fullVt <- Vt
    
    # If length of missing is > 0 there is at least one missing
    if (length(missing) > 0) {
      # remove from Yt
      Yt = as.matrix(Yt[-missing])
      rownames(Yt) <- relevant.names[-missing]
      
      # Remove from Ft
      Ft = as.matrix(Ft[, -missing])
      colnames(Ft) <- relevant.names[-missing]
      
      # Remove from Vt
      Vt = Vt[-missing, -missing]
    }
    
    # if the only observation is NA (univariate case) or all observations are NA (multivariate case)
    if (length(Yt) == 0) { 
      at <- Gt %*%  mt		                  
      rownames(at) <- rownames(mt)
      
      Rt <- Gt %*%  Ct %*% t(Gt) + Wt        
      mt=at
      Ct=Rt
      ft <- t(fullFt) %*% at
      ft2=ft
      Qt <- t(fullFt) %*% Rt %*% fullFt + fullVt
      Qt2=Qt
      At=NA
      et=NA
      ut=NA
      
    }else{
      
      # Run the Kalman filter - only if we observe at least one variable!
      at <- Gt %*%  mt		                           
      rownames(at) <- rownames(mt)
      
      Rt <- Gt %*%  Ct %*% t(Gt) + Wt               
      
      ft <- t(Ft) %*% at		      	                
      ft2 <- t(fullFt) %*% at                       
      
      Qt <- t(Ft) %*% Rt %*%  Ft + Vt                
      Qt2 <- t(fullFt) %*% Rt %*% fullFt + fullVt   
      
      At <- Rt %*% Ft %*% solve(Qt)                 
      et <- Yt  - ft	                             
      ut <- et / sqrt(diag(Qt))                    
      
      # - update the parameter vector and variance matrix
      mt <- at + At %*% et                         
      Ct <- Rt - At  %*% Qt %*% t(At)	             
    }
    
    # Make sure Ct is symmetrical
    Ct <- (Ct + t(Ct))/2
    
    
    # Save the values in lists
    Yt.list[[i]] <- Yt
    at.list[[i]] <- at
    Rt.list[[i]] <- Rt
    ft.list[[i]] <- ft2
    Qt.list[[i]] <- Qt2
    At.list[[i]] <- At
    et.list[[i]] <- et
    ut.list[[i]] <- ut
    mt.list[[i]] <- mt
    Ct.list[[i]] <- Ct
    Ft.list[[i]] <- t(Ft)
    VSE.list[[i]] <- VSE
    Gt.list[[i]] <- Gt
    fullFt.list[[i]] <- t(fullFt)
    
  }
  
  return(list(
    Yt=Yt.list,
    at=at.list,
    Rt=Rt.list,
    ft=ft.list,
    Qt=Qt.list,
    At=At.list,
    et=et.list,
    ut=ut.list,
    mt=mt.list,
    Ct=Ct.list,
    F=Ft.list,
    vse=VSE.list,
    Gt.list=Gt.list,
    fullFt.list=fullFt.list
  ))
}



# Function for the EM algorithm for non-hierarchical models ----
# - Des is the data set for the EM algorithm
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - V0 is the observation variance
# - W0 is the system variance
# - time.var is the variable that describes the observation time
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on. 
# - best.n.waves is a vector with the best number of waves for each variable
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
# - relevant.names are the names of the variables to be co-modeled
# - Spline.list is a list with the splines calculated for each farm
# - steps is the number of iterations
# - silent is to put equal to TRUE of FALSE whether you want to print the step you are in or not
# - DLM.version is the name of the function that runs the DLM
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not, here is hierarchical=FALSE
# - level2 is the second level (correspond to the region, here is NULL because is for non-hierarchical models)
# - level3 is the third level (correspond to the farm/site, here is NULL because is for non-hierarchical models)
runEM = function(Des, mu0, C0, V0, W0, time.var, time.in.year, best.n.waves, trend, 
                 relevant.names, Spline.list, steps = 1, silent = TRUE, 
                 DLM.version = FUNFILTERG, hierarchical=hierarchical, level2=NULL, level3=NULL) {
  Vs = list()
  Ws = list()
  Vn = V0
  Wn = W0
  Cn <- diag(0, length(mu0))
  mu0n <- rep(0, length(mu0))
  
  #Choose wich version of DLM to use
  runDLM <- DLM.version
  
  # Iterate over steps
  for (s in 1:steps) {
    print(paste("Step", s))
    
    # Set sums and counts to zero
    # Count structure for observation variance
    sumV = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
    # Sum for observation variance
    sumObs = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
    # Count for system variance
    sumW = 0
    # Sum for observation variance
    sumSys = matrix(0, nrow=length(mu0), ncol=length(mu0))
    # Iterate over Des
    for (b in 1:length(Des)) {
      
      # Expectation step - filter and smooth
      if(silent == FALSE){
        print(b)
      }
      
      res <- runDLM(D = Des[[b]], mu0 = mu0, C0 = C0, V = Vn, W = Wn, time.var, time.in.year, 
                    best.n.waves, trend, relevant.names, Spline.list, hierarchical=FALSE, level2, level3)
      smot = runSmoother(res)
      
      # Get the smoothed C0 from this one
      Cn <- Cn + smot$Cts[,,1]
      mu0n <-  mu0n + smot$mts[,,1]
      
      # Set contributions to sums and counts to zero for this D
      bSumV = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
      sumCountV = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
      bSumW = matrix(0, nrow=length(mu0), ncol=length(mu0))
      sumCountW = 0
      
      # Sum contributions
      n = length(res$mt)
      
      # Iterate over time within D 
      for (t in 1:n) {
        
        # Get the Gt 
        Gt <- res$Gt.list[[t]]
        
        # Only if observations at all
        if (length(res$Yt[[t]]) > 0) {
          
          # Observation variance
          
          # Find the contribution to the sum even though it does not have the correct dimension
          Vcont = res$F[[t]]%*%smot$Cts[,,t]%*%t(res$F[[t]]) + 
            (res$Y[[t]] - res$F[[t]]%*%smot$mts[,,t])%*%t((res$Y[[t]] - res$F[[t]]%*%smot$mts[,,t]))
          # Get the pointer matrix
          vse = res$vse[[t]]
          # Create a full matrix
          Vfull = matrix(0, nrow = length(relevant.names), ncol = length(relevant.names))
          # Find out which cells to enter
          ind = c()
          
          for (i in 1:dim(Vcont)[1]) { 
            if (vse[i,i] > 0) {
              ind = c(ind, i)
            }
          }
          
          # Enter the contribution into the right cells of the 3 x 3 matrix
          Vfull[ind,ind] = Vcont[ind,ind]
          # Add the resulting 3 x 3 matrix
          bSumV = bSumV + Vfull
          # Adjust the counts matrix
          sumCountV = sumCountV + vse          
          
          
          # System variance
          if (t > 1) {
            # Find the contribution - the dimension is always correct
            bSumW = bSumW + smot$Lt[,,t] + 
              (smot$mts[,,t] - Gt%*%smot$mts[,,t-1])%*%t(smot$mts[,,t] - Gt%*%smot$mts[,,t-1])            
            sumCountW = sumCountW + 1
          }
        }
      }
      
      # Check for negative variances
      ignore = FALSE
      for (j in 1:length(relevant.names)) {
        if (bSumV[j, j] < 0) {
          # Adjust to 0
          bSumV[j, j] = 0
          bSumV[j, ] = 0
          bSumV[, j] = 0
          # Print a comment
          if (! silent) print(paste("Negative contribution to observation variance", j, "for D", b))
        }
      }
      
      # Check for negative variances
      ignore = FALSE
      for (j in 1:length(mu0)) {
        if (bSumW[j, j] < 0) {
          # Adjust to 0
          bSumW[j, j] = 0
          bSumW[, j] = 0
          bSumW[j, ] = 0
          # Print a message
          if (! silent) print(paste("Negative contribution to system variance", j, "for D", b))
        }
      }
      
      # This will never happen (used for debugging)
      if (ignore) {
        print(paste("Contribution to observation variance ignored from D", b))
        for (i in 1:length(relevant.names)) {
          for (j in 1:length(relevant.names)) {
            bSumV[i, j] = 0
          }
        }
        sumCountV = matrix(0, nrow=length(relevant.names), ncol=length(relevant.names))
      }
      # Add the contribution from the D to the total
      sumObs = sumObs + bSumV
      sumV = sumV + sumCountV
      
      
      # This will never happen (used for debugging)
      if (ignore) {
        print(paste("Contribution to system variance ignored from D", b))
        for (i in 1:length(mu0)) {
          for (j in 1:length(mu0)) {
            bSumW[i, j] = 0
          }
        }
        sumCountW = 0
      }
      # Add the contribution from the D to the total
      sumSys = sumSys + bSumW
      sumW = sumW + sumCountW
      
    }
    # Normalize by counts
    Vn = sumObs/sumV
    Wn = sumSys/sumW
    Cn <- Cn/n
    mu0n <-  mu0n/n
    # Make sure they are symmetric
    Vn = (Vn + t(Vn))/2
    Wn = (Wn + t(Wn))/2
    Cn <- (Cn + t(Cn))/2
    Vn <<- Vn
    Wn <<- Wn
    # Mke sure they have the right names
    colnames(Vn) <- relevant.names
    rownames(Vn) <- relevant.names
    
    colnames(Wn) <- rownames(mu0)
    rownames(Wn) <- rownames(mu0)
    
    colnames(Cn) <- rownames(mu0)
    rownames(Cn) <- rownames(mu0)
    
    mu0n <- as.matrix(mu0n)
    rownames(mu0n) <- rownames(mu0)
    
    # Save them in lists
    Vs[[s]] = Vn
    Ws[[s]] = Wn
  }
  return(list(V=Vs, W=Ws, Cn=Cn, mu0n=mu0n, smot=smot))
}



# Function of the Kalman Smoother ----
# - res is the result returned from the filter (DLM)
runSmoother <- function(res) {
  
  n = length(res$mt)
  p = length(res$mt[[1]])
  mts <- array(NA,dim=c(p,1,n));
  Cts <- array(NA,dim=c(p,p,n));
  
  # Put last value equal to filtered
  mts[,,n] <- res$mt[[n]]
  Cts[,,n] <- res$Ct[[n]]
  
  # These are useful
  Bt <- array(NA,dim=c(p,p,n))
  Lt <- array(NA,dim=c(p,p,n));  
  
  # Iterate backwards over days
  for(i in ((n-1):1))   {
    
    # Get Gt
    Gt <- res$Gt.list[[i+1]]
    
    res$R[[i+1]] <- as.matrix(res$R[[i+1]])
    
    Bt[,,i] <- as.matrix( res$Ct[[i]] %*% t(Gt) %*% solve(res$R[[i+1]]) )
    mts[,,i] <- res$mt[[i]] + Bt[,,i] %*% (mts[,,i+1] - res$a[[i+1]])
    Cts[,,i] <- as.matrix( res$C[[i]] + Bt[,,i] %*% (Cts[,,i+1] - res$R[[i+1]]) %*% t(Bt[,,i]) )
    
    mts[,,i] <- as.matrix(mts[,,i])
    
  }
  # give names
  rownames(mts) <- rownames(res$mt[[i]])
  colnames(mts) <- "mts"
  rownames(Cts) <- rownames(res$Ct[[i]])
  colnames(Cts) <- colnames(res$Ct[[i]])
  
  # Now when we are at it: Find L and store it for the EM algorithm
  for(i in ((n):2))  {
    Lt[,,i] <- Cts[,,i] + Gt%*%Cts[,,i-1]%*%t(Gt) - Cts[,,i]%*%t(Bt[,,i-1]) - Bt[,,i-1]%*%Cts[,,i]
  }
  rownames(Lt) <- rownames(res$Ct[[i]])
  colnames(Lt) <- colnames(res$Ct[[i]])
  
  return(list(mts=mts,
              Cts=Cts,
              Lt=Lt,
              D=res$D));
}



# Function to extract the relevant information from the DLM and the Smoother ----
# - res is the result returned from the filter (DLM)
# - smot is the result returned from the smoother 
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not
# - D is the data set, containing the time series we want to model
extract.res <- function(res, smot, hierarchical, D){
  
  if(hierarchical==FALSE){
    
    df <- c()
    
    if(length(relevant.names)>1){
      
      # Get the raw observations
      # -- if there are missing values in some variable, change to NA
      for(i in 1:length(res$Yt)){ 
        
        if(length(res$Yt[[i]]) < length(relevant.names)){
          matrix <- matrix(nrow = length(relevant.names), ncol = 1)
          rownames(matrix) <- relevant.names
          colnames(matrix) <- colnames(res$Yt[[i]])
          
          for(a in rownames(res$Yt[[i]])){ 
            matrix[a,1] <- res$Yt[[i]][a,1]
          }
          res$Yt[[i]] <- matrix
        }
      }
      df_Yt <- data.frame(t(sapply(res$Yt,c)))
      colnames(df_Yt) <- relevant.names
      df <- df_Yt
      
      # Get the filtered mean (multiply mt * Ft to get the right values for mt of the harmonic waves variables)
      df_mt <- data.frame()
      
      for(i in 1:length(res$mt)){ 
        
        val <- res$fullFt.list[[i]] %*% res$mt[[i]]
        df_val_mt <- data.frame(t(sapply(val,c)))
        df_mt <- rbind(df_mt, df_val_mt)
      }
      colnames(df_mt) <- paste("mt_", rownames(val), sep="")
      df <- cbind(df, df_mt)
      
      # Get the forecasts
      for(i in 1:length(res$ft)){ 
        
        if(length(res$ft[[i]]) < length(relevant.names)){
          matrix <- matrix(nrow = length(relevant.names), ncol = 1)
          rownames(matrix) <- relevant.names
          colnames(matrix) <- colnames(res$ft[[i]])
          
          for(a in rownames(res$ft[[i]])){ 
            matrix[a,1] <- res$ft[[i]][a,1]
          }
          res$ft[[i]] <- matrix
        }
      }
      df_ft <- data.frame(t(sapply(res$ft,c)))
      colnames(df_ft) <- paste('ft_', relevant.names,sep='')
      df <- cbind(df, df_ft)
      
      # Get the raw forecasts errors
      for(i in 1:length(res$et)){ 
        
        if(length(res$et[[i]]) < length(relevant.names)){
          matrix <- matrix(nrow = length(relevant.names), ncol = 1)
          rownames(matrix) <- relevant.names
          colnames(matrix) <- colnames(res$et[[i]])
          
          for(a in rownames(res$et[[i]])){ 
            matrix[a,1] <- res$et[[i]][a,1]
          }
          res$et[[i]] <- matrix
        }
      }
      df_et <- data.frame(t(sapply(res$et,c)))
      colnames(df_et) <- paste('et_', rownames(res$et[[1]]),sep='')
      df <- cbind(df, df_et)
      
      # Get the standardized forecasts errors
      for(i in 1:length(res$ut)){ 
        
        if(length(res$ut[[i]]) < length(relevant.names)){
          matrix <- matrix(nrow = length(relevant.names), ncol = 1)
          rownames(matrix) <- relevant.names
          colnames(matrix) <- colnames(res$ut[[i]])
          
          for(a in rownames(res$ut[[i]])){ 
            matrix[a,1] <- res$ut[[i]][a,1]
          }
          res$ut[[i]] <- matrix
        }
      }
      df_ut <- data.frame(t(sapply(res$ut,c)))
      colnames(df_ut) <- paste('ut_', rownames(res$ut[[1]]),sep='')
      df <- cbind(df, df_ut)
      
      # Get the Ct of the DLM (multiply t(Ft) * Ct * Ft to get the right values for Ct of the harmonic waves variables))
      for(i in 1:length(res$Ct)){ 
        
        mat <- res$fullFt.list[[i]] %*% res$Ct[[i]] %*% t(res$fullFt.list[[i]])
        
        for(name in relevant.names){ 
          df[i, paste("Ct_", name, sep="")] <- mat[name, name] 
        }
      }
      
      # Get the Qt of the DLM
      for(i in 1:length(res$Qt)){ 
        
        for(name in colnames(res$Qt[[i]])){
          df[i, paste("Qt_", name, sep="")] <- res$Qt[[i]][name, name]
        }
      }
      
      if (!is.null(smot)) {
        
        # Get the smoothed mean (multiply mts * Ft to get the right values for mt of the harmonic waves variables)
        df_mts <- data.frame()
        
        # Get the mts from the smoother
        for(i in 1:dim(smot$mts)[3]){ 
          
          val2 <- res$fullFt.list[[i]] %*% smot$mts[,,i]
          df_val_mts <- data.frame(t(sapply(val2,c)))
          df_mts <- rbind(df_mts, df_val_mts)
        }
        colnames(df_mts) <- paste("mts_", rownames(val2), sep="")
        df <- cbind(df, df_mts)
        
        # Get the smoothed variance (multiply t(Ft) * Cts * Ft to get the right values for Cts of the harmonic waves variables))
        for(i in 1:dim(smot$Cts)[3]){ 
          
          mat2 <- res$fullFt.list[[i]] %*% smot$Cts[,,i] %*% t(res$fullFt.list[[i]])
          
          for(name in relevant.names){ 
            df[i, paste("Cts_", name, sep="")] <- mat2[name, name]
          }
        }
      }
    }
    
    if(length(relevant.names)==1){
      
      # Get the raw observations
      res$Yt[lengths(res$Yt) == 0] <- NA  
      df_Yt <- data.frame(sapply(res$Yt,c))
      colnames(df_Yt) <- relevant.names
      df <- df_Yt
      
      # Get the filtered means 
      df_mt <- data.frame(t(sapply(res$mt,c)))
      colnames(df_mt) <- c(paste('mt_level_', relevant.names, sep=''),
                           paste('mt_trend_', relevant.names, sep=''))
      df <- cbind(df, df_mt)
      
      # Get the forecasts
      df_ft <- data.frame(sapply(res$ft,c))
      colnames(df_ft) <- paste('ft_', relevant.names, sep='')
      df <- cbind(df, df_ft)
      
      # Get the raw forecasts errors
      df_et <- data.frame(sapply(res$et,c))
      colnames(df_et) <- paste('et_', relevant.names, sep='')
      df <- cbind(df, df_et)
      
      # Get the standardized forecasts errors
      df_ut <- data.frame(sapply(res$ut,c))
      colnames(df_ut) <- paste('ut_', relevant.names, sep='')
      df <- cbind(df, df_ut)
      
      # Get the Cts for the level of the DLM
      for(i in 1:length(res$Ct)){ 
        
        # - get variables names just for the level on Ct
        all.names <- rownames(res$Ct[[i]])
        library(gdata)
        name <- all.names[! startsWith(all.names, "d.")]  
        
        df[i, paste("Ct_", name, sep="")] <- res$Ct[[i]][name, name]
      }
      
      # Get the Qts of the DLM
      for(i in 1:length(res$Qt)){ 
        df[i, paste("Qt_", relevant.names, sep="")] <- res$Qt[[i]]
      }
      
      if (!is.null(smot)) {
        
        # Get the mts from the smoother
        for(i in 1:dim(smot$mts)[3]){ 
          
          # - get variables names just for the level on mts
          all.names <- rownames(smot$mts)
          library(gdata)
          name <- all.names[! startsWith(all.names, "d.")] 
          
          # - get mts and add it to df
          df[i, paste("mts_", name, sep="")] <- smot$mts[,,i][name]
        }
        
        # Get the Cts from the smoother
        for(i in 1:dim(smot$Cts)[3]){ 
          
          # - get variables names just for the level on Cts
          all.names <- rownames(smot$Cts)
          library(gdata)
          name <- all.names[! startsWith(all.names, "d.")] 
          
          # - get Cts and add it to df
          df[i, paste("Cts_", name, sep="")] <- smot$Cts[,,i][name, name]
        }
      }
    }
  }
  
  if(hierarchical==TRUE){
    
    if(length(relevant.names)>1){
      
      # Get the raw observations
      # -- if there are missing values in some variable, change to NA
      for(i in 1:length(res$Yt)){ 
        
        if(length(res$Yt[[i]]) < length(relevant.names)){
          matrix <- matrix(nrow = length(relevant.names), ncol = 1)
          rownames(matrix) <- relevant.names
          colnames(matrix) <- colnames(res$Yt[[i]])
          
          for(a in rownames(res$Yt[[i]])){ 
            matrix[a,1] <- res$Yt[[i]][a,1]
          }
          res$Yt[[i]] <- matrix
        }
      }
      df_Yt <- data.frame(t(sapply(res$Yt,c)))
      colnames(df_Yt) <- relevant.names
      df <- df_Yt
      
      # Get the filtered mean (multiply mt * Ft to get the right values for mt of the harmonic waves variables)
      df_mt <- data.frame()
      
      for(i in 1:length(res$mt)){ 
        
        val <- res$fullFt.list[[i]] %*% res$mt[[i]]
        df_val_mt <- data.frame(t(sapply(val,c)))
        df_mt <- rbind(df_mt, df_val_mt)
      }
      colnames(df_mt) <- paste("mt_", rownames(val), sep="")
      df <- cbind(df, df_mt)
      
      # Get the forecasts
      for(i in 1:length(res$ft)){ 
        
        if(length(res$ft[[i]]) < length(relevant.names)){
          matrix <- matrix(nrow = length(relevant.names), ncol = 1)
          rownames(matrix) <- relevant.names
          colnames(matrix) <- colnames(res$ft[[i]])
          
          for(a in rownames(res$ft[[i]])){ 
            matrix[a,1] <- res$ft[[i]][a,1]
          }
          res$ft[[i]] <- matrix
        }
      }
      df_ft <- data.frame(t(sapply(res$ft,c)))
      colnames(df_ft) <- paste('ft_', relevant.names,sep='')
      df <- cbind(df, df_ft)
      
      # Get the raw forecasts errors
      for(i in 1:length(res$et)){ 
        
        if(length(res$et[[i]]) < length(relevant.names)){
          matrix <- matrix(nrow = length(relevant.names), ncol = 1)
          rownames(matrix) <- relevant.names
          colnames(matrix) <- colnames(res$et[[i]])
          
          for(a in rownames(res$et[[i]])){ 
            matrix[a,1] <- res$et[[i]][a,1]
          }
          res$et[[i]] <- matrix
        }
      }
      df_et <- data.frame(t(sapply(res$et,c)))
      colnames(df_et) <- paste('et_', rownames(res$et[[1]]),sep='')
      df <- cbind(df, df_et)
      
      # Get the standardized forecasts errors
      for(i in 1:length(res$ut)){ 
        
        if(length(res$ut[[i]]) < length(relevant.names)){
          matrix <- matrix(nrow = length(relevant.names), ncol = 1)
          rownames(matrix) <- relevant.names
          colnames(matrix) <- colnames(res$ut[[i]])
          
          for(a in rownames(res$ut[[i]])){ 
            matrix[a,1] <- res$ut[[i]][a,1]
          }
          res$ut[[i]] <- matrix
        }
      }
      df_ut <- data.frame(t(sapply(res$ut,c)))
      colnames(df_ut) <- paste('ut_', rownames(res$ut[[1]]),sep='')
      df <- cbind(df, df_ut)
      
      # Get the Ct of the DLM (multiply t(Ft) * Ct * Ft to get the right values for Ct of the harmonic waves variables))
      for(i in 1:length(res$Ct)){ 
        
        mat <- res$fullFt.list[[i]] %*% res$Ct[[i]] %*% t(res$fullFt.list[[i]])
        
        for(name in relevant.names){ 
          df[i, paste("Ct_", name, sep="")] <- mat[name, name]
        }
      }
      
      # Get the Qt of the DLM
      for(i in 1:length(res$Qt)){ 
        
        for(name in colnames(res$Qt[[i]])){
          df[i, paste("Qt_", name, sep="")] <- res$Qt[[i]][name, name]
        }
      }
      
      if (!is.null(smot)) {
        
        # Get the smoothed mean (multiply mts * Ft to get the right values for mt of the harmonic waves variables)
        df_mts <- data.frame()
        
        # Get the mts from the smoother
        for(i in 1:dim(smot$mts)[3]){ 
          
          val2 <- res$fullFt.list[[i]] %*% smot$mts[,,i]
          df_val_mts <- data.frame(t(sapply(val2,c)))
          df_mts <- rbind(df_mts, df_val_mts)
        }
        colnames(df_mts) <- paste("mts_", rownames(val2), sep="")
        df <- cbind(df, df_mts)
        
        # Get the smoothed variance (multiply t(Ft) * Cts * Ft to get the right values for Cts of the harmonic waves variables)
        for(i in 1:dim(smot$Cts)[3]){ 
          
          mat2 <- res$fullFt.list[[i]] %*% smot$Cts[,,i] %*% t(res$fullFt.list[[i]])
          
          for(name in relevant.names){ 
            df[i, paste("Cts_", name, sep="")] <- mat2[name, name] 
          }
        }
      }
    }
    
    if(length(relevant.names)==1){
      
      # Get the raw observations
      df <- data.frame("date"=sort(unique(D$date)))
      
      # - columns for the country
      df[,paste("Country", relevant.names, sep=".")] <- rep(NA, length(unique(D$date)))
      
      # - get columns for the regions
      for (region in sort(unique(D$local.authority))){ 
        df[,paste(region, relevant.names, sep=".")] <- rep(NA, length(unique(D$date)))
      }
      
      # - get columns for the farms
      for (farm in sort(unique(D$site))){
        df[,paste(farm, relevant.names, sep=".")] <- rep(NA, length(unique(D$date)))
      }
      
      for(i in 1:length(res$Yt)){ 
        for(a in rownames(res$Yt[[i]])){
          df[i,a] <- res$Yt[[i]][a,1]
        }
      }
      
      # Get the filtered mean (multiply mt * Ft to get the right values for mt (without the trend))
      df_mt <- data.frame()
      
      for(i in 1:length(res$mt)){
        
        val <- res$fullFt.list[[i]] %*% res$mt[[i]]
        df_val_mt <- data.frame(t(sapply(val,c)))
        df_mt <- rbind(df_mt, df_val_mt)
      }
      colnames(df_mt) <- paste("mt_", rownames(val), sep="")
      df <- cbind(df, df_mt)
      
      # Get the forecasts
      # - get columns for the county
      df[,paste("ft_", "Country.", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      
      # - get columns for the regions
      for (region in sort(unique(D$local.authority))){
        df[,paste("ft_", region, ".", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      }
      
      # - get columns for the farms
      for (farm in sort(unique(D$site))){ 
        df[,paste("ft_", farm, ".", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      }
      
      for(i in 1:length(res$ft)){ 
        for(a in rownames(res$ft[[i]])){ 
          df[i, paste("ft_", a, sep="")] <- res$ft[[i]][a,1]
        }
      }
      
      # Get the raw forecasts errors
      # - get columns for the county
      df[,paste("et_", "Country.", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      
      # - get columns for the regions
      for (region in sort(unique(D$local.authority))){ 
        df[,paste("et_", region, ".", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      }
      
      # - get columns for the farms
      for (farm in sort(unique(D$site))){ 
        df[,paste("et_", farm, ".", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      }
      
      for(i in 1:length(res$et)){ 
        for(a in rownames(res$et[[i]])){
          df[i, paste("et_", a, sep="")] <- res$et[[i]][a,1]
        }
      }
      
      # Get the standardized forecasts errors
      # - get columns for the county
      df[,paste("ut_", "Country.", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      
      # - get columns for the regions
      for (region in sort(unique(D$local.authority))){ 
        df[,paste("ut_", region, ".", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      }
      
      # - get columns for the farms
      for (farm in sort(unique(D$site))){ 
        df[,paste("ut_", farm, ".", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      }
      
      for(i in 1:length(res$ut)){ 
        for(a in rownames(res$ut[[i]])){ 
          df[i, paste("ut_", a, sep="")] <- res$ut[[i]][a,1]
        }
      }
      
      # Get the Ct of the DLM (multiply t(Ft) * Ct * Ft to get the right values for Ct 
      # (only variance of the relevant.names))
      for(i in 1:length(res$Ct)){ 
        
        mat <- res$fullFt.list[[i]] %*% res$Ct[[i]] %*% t(res$fullFt.list[[i]])
        names <- row.names(mat)
        
        for(n in names){
          df[i, paste("Ct_", n, sep="")] <- mat[n, n] #diagonal values
        }
      }
      
      # Get the Qt of the DLM
      # - get columns for the county
      df[,paste("Qt_", "Country.", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      
      # - get columns for the regions
      for (region in sort(unique(D$local.authority))){ 
        df[,paste("Qt_", region, ".", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      }
      
      # - get columns for the farms
      for (farm in sort(unique(D$site))){ 
        df[,paste("Qt_", farm, ".", relevant.names, sep="")] <- rep(NA, length(unique(D$date)))
      }
      
      for(i in 1:length(res$Qt)){ 
        for(a in rownames(res$Qt[[i]])){ 
          df[i, paste("Qt_", a, sep="")] <- res$Qt[[i]][a,a]
        }
      }
      
      if (!is.null(smot)) {
        
        df_mts <- data.frame()
        
        # Get the mts from the smoother
        # multiply mts * Ft to get the right values for mts (without trend)
        for(i in 1:dim(smot$mts)[3]){ 
          
          val2 <- res$fullFt.list[[i]] %*% smot$mts[,,i]
          df_val_mts <- data.frame(t(sapply(val2,c)))
          df_mts <- rbind(df_mts, df_val_mts)
        }
        colnames(df_mts) <- paste("mts_", rownames(val2), sep="")
        df <- cbind(df, df_mts)
        
        # Get the Cts from the smoother
        # multiply t(Ft) * Ct * Ft to get the right values for Ct (only variance of the relevant.names)
        for(i in 1:dim(smot$Cts)[3]){ 
          
          mat2 <- res$fullFt.list[[i]] %*% smot$Cts[,,i] %*% t(res$fullFt.list[[i]])
          names <- row.names(mat2)
          
          for(n in names){
            df[i, paste("Cts_", n, sep="")] <- mat2[n, n]
          }
        }
      }
    }
  }
  return(df)
}



# Function to get all initial parameters needed for the DLM (with splines or harmonic waves) ----
# - D is the data set, containing the time series we want to model
# - identifyer is the variable that identifies the production cycle/batch
# - time.var is the variable that describes the observation time
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on. 
# - best.n.waves is a vector with the best number of waves for each variable
# - relevant.names are the names of the variables to be co-modeled
# - metadata.names are the variables' names we need to give context to our time series (specification/information about the data)
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - C.all is a list with all Ct matrices
# - for.V are the parameters needed to create V
# - level2 is the second level (here correspond to the region)
# - level3 is the third level (here correspond to the farm/site)
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not
# - level is the first level (here correspond to the country)
# - Spline.list is a list with the splines calculated for each farm
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
# - w is the is the base frequency
# - loop is the iteration number (no need to define it)
# - n.runs is the number of runs (no need to define it)
# - end is the end of each iteration (no need to define it)
# - expected.start.time is the observation time when we expect the model to start
# - plot.it (TRUE of FALSE) determines whether plots should be made to show the fit of the model to the data of some unique identified units
# - remove.zeros whether observations in the time series with the value 0 should be removed. This is not usually recommended, unless 0s are stand-ins for NAs
# - round.by is the number of decimals you want to round the estimates
# - silent determines if the function should print the r^2 values of the fitted model to the console
get.DLM.parameters <- function(D, identifyer, time.var, time.in.year, best.n.waves, relevant.names, metadata.names, 
                               mu0=mu0, C0=C0, C.all=C.all, for.V=for.V, level2="local.authority", level3="site",
                               hierarchical=hierarchical, level=level, Spline.list=Spline.list,
                               trend, w=w, loop=loop, n.runs=n.runs, end=end, expected.start.time, 
                               plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE){
  
  #get original data for later
  D.full <- D
  time.var.all <- time.var
  
  for(N.w in best.n.waves){
    
    loop=loop+1
    
    if(isTRUE(is.na(N.w))){
      
      #get original data
      D <- D.full
      D.A <- D
      
      #get the name we are referring to
      name <- relevant.names[loop]
      
      #get the time.var for the corresponding name
      time.var=time.var.all[loop]
      
      print('###############################################################')
      print(paste('DEFINING PARAMETERS RELATED TO relevant.name', '=', name))
      start.V <- Sys.time()
      
      # Remove outliers, based on overall mean and moving SD
      par(mfrow=c(1,2))
      ylim <- range(na.omit(D[,name]))
      Mean <- mean(na.omit(D[,name]))
      SD <- sd(na.omit(D[,name]))
      upper <- Mean + 3*SD #99.7% of data occurs within 3 standard deviations of the mean within a normal distribution
      lower <- Mean - 3*SD
      remove.i <- which(D[,name] > upper | D[,name] < lower)
      D[remove.i,name] <- NA
      par(mfrow=c(1,1))
      
      
      # Elements for the initial parameter vector (mu0)
      D.B <- subset(D, !is.na(D[,name]))
      x <- D.B[,time.var]
      y <- D.B[,name]
      Spline = smooth.spline(x = x, y = y)
      
      plot(y~x, xlab=time.var, ylab=name, main = paste("Spline", '=', name))
      lines(Spline, col='red', lwd=3)
      
      # - Save the spline function - we will need it for the Gt matrix
      if(isFALSE(hierarchical)){
        spline.name <- paste(name, '_Spline', sep='')
      }
      if(isTRUE(hierarchical)){
        spline.name <- paste(level, name, 'Spline', sep='_')
      }
      
      Spline.list[[length(Spline.list)+1]] <- Spline
      names(Spline.list)[length(Spline.list)] <- spline.name
      
      # - Finalize and save mu
      pred0 <- predict(object = Spline, x = (expected.start.time-1))$y
      mu <- c(pred0, 1)
      
      #give names to columns and rows
      if(isFALSE(hierarchical)){
        mu.names <- c(name, paste('d.', name, sep=''))      
      }
      if(isTRUE(hierarchical)){
        mu.names <- c(paste(level, name, sep="."),
                      paste('d', level, sep=".", name))      
      }
      
      mu <- matrix(mu)
      row.names(mu) <- mu.names
      
      # - Initial parameter vector for all variables (mu0)
      mu0 <- rbind(mu0, mu)
      colnames(mu0) <- "mu0"
      
      
      # Elements for the prior variance (C0)
      # - Diff for the prior variance on the initial rate of change (difference between the present observation and the previous observation)
      Diff <- diff(D[,name])
      
      # - Remove outliers, as these likely represent transitions from 1 production cycle to the other
      Q <- quantile(x = na.omit(Diff), probs = c(0.025, 0.975))
      Diff[which(Diff <= Q[1] | Diff >= Q[2])] <- NA
      D.A[,paste('d.', name, sep='')] <- c(NA,Diff)
      
      # - Re-order the columns
      D.A <- D.A[,c(metadata.names, name, paste('d', sep=".", name))]
      #only consider the start
      D.A <- subset(D.A, D.A[,time.var] %in% (expected.start.time):(expected.start.time+6))
      
      # - Finalize and save C0
      D.A <- na.omit(D.A)
      C <- cov(D.A[,c(name, paste('d', sep=".", name))])
      
      rm(D.A)
      
      if(loop==1){
        start = n.runs
      }else{
        start=end + 1
      }
      end <- start + 1
      
      # - Prior variance for all variables (C0)
      C0[start:end,start:end] <- C
      rownames(C0)[start:end] <- mu.names
      colnames(C0)[start:end] <- mu.names   
      C0 <- as.matrix(C0)
      
      # - Save all C matrices to see where there is NA values and replace them afterwards
      C.all[[paste(level, name, sep=".")]] <- C
      
      
      # Elements for the observational matrix (V)
      resid.ID.all <- c()
      
      for(ID in unique(D[,identifyer])){
        
        ID.set <- subset(D, D[,identifyer] == ID)
        
        # use a two-sided moving average to estimate V
        y <- ID.set[,name]
        x <- ID.set[,time.var]
        MA <- moving.function(x = y, n = 5, FUN = mean)
        
        # we find the residuals between the observed and filtered/estimated values
        resid.ID <- ID.set[,name] - MA
        
        # we combine the residuals for all individual IDs to estimate an overall V
        resid.ID.all <- c(resid.ID.all, resid.ID)
      }
      
      #get the variance per name
      var.resid <- var(na.omit(resid.ID.all))
      
      if(isFALSE(hierarchical)){
        a <- paste(name, sep='')
      }
      if(isTRUE(hierarchical)){
        a <- paste(level, name, sep='.')
      }
      
      for.V <- cbind(for.V, var.resid)
      colnames(for.V)[ncol(for.V)] <- a
      
      
      print(paste('--Done! It took', difftime(Sys.time(), start.V, units = 'min'), 'minutes'))
      
      gc()
      
    }else{
      
      if(N.w!=0){
        
        #get the name we are referring to
        name <- relevant.names[loop] 
        
        #get the time.var for the corresponding name
        time.var=time.var.all[loop]
        
        print('###############################################################')
        print(paste('DEFINING PARAMETERS RELATED TO relevant.name', '=', name))
        start.V <- Sys.time()
        
        # Start by plotting the average level per period and over longer time scales
        if(length(which(is.na(D[,name]))) > 0){
          D <- D[-which(is.na(D[,name])),]
        }
        
        par(mfrow=c(1,2))
        
        # - per period
        period <- (2*pi)/w
        D$time.short <- D[,time.var]%%period
        agg.short <- aggregate(D[,name], by=list(D$time.short), FUN=mean)
        plot(agg.short, type='l', xlab='Time in one period', ylab=name)
        
        # - over longer time
        agg.long <- aggregate(D[,name], by=list(D[,time.var]), FUN=mean)
        plot(agg.long, type='l', xlab=time.var, ylab=name)
        
        par(mfrow=c(1,1))
        
        # Required packages
        library('MuMIn')
        lmer <- lme4::lmer #(linear regression model)
        
        # Remove zero-values, if relevant
        if(remove.zeros == TRUE){
          remove.i <- which(D[,name] == 0)
          D <- D[-remove.i,]
        }
        
        # We make an empty vector, into which we will add the relevant text for the linear function
        lm.vector <- c()
        
        # We also make an empty vector, into which we will add the names for the parameter vector
        names.vector <- c()
        
        # We add the target variable (the variable we model with this collection of harmonics)
        lm.vector <- c(lm.vector, name, '~')
        names.vector <- c(names.vector, 'Level')
        
        # Make the linear trend component, if relevant
        if(trend == TRUE){
          lm.vector <- c(lm.vector, time.var)
          names.vector <- c(names.vector, 'Trend +')
        }
        
        # Now we add the relevant number of harmonics, one at a time
        if(N.w > 0){
          
          if(trend == TRUE){
            lm.vector <- c(lm.vector, '+')
          }
          
          for(i in 1:N.w){
            if(i == 1){
              a <- paste('cos(', i, '*w*', time.var, ') + sin(', i, '*w*', time.var, ')')
              lm.vector <- c(lm.vector, a)
              names.vector <- c(names.vector, paste('wave.', i, '.cos', sep=''), paste('wave.', i, '.sin', sep=''))
            }else{
              a <- paste('+ cos(', i, '*w*', time.var, ') + sin(', i, '*w*', time.var, ')')
              lm.vector <- c(lm.vector, a)
              names.vector <- c(names.vector, paste('wave.', i, '.cos', sep=''), paste('wave.', i, '.sin', sep=''))
            }
          }
        }
        
        # We can now make the final name vector
        names.vector <- paste(name, '_', names.vector, sep='')
        
        # Make sure we can handle the case with no trend AND no waves
        if(length(lm.vector) == 2){
          lm.vector <- c(lm.vector, 1)
        }
        
        # Now we make it into a linear function - without random effect
        b1 <- as.formula(paste(lm.vector, collapse=' '))
        lm.1 <- lm(b1, data = D )
        S1 <- summary(lm.1)
        r.squared <- S1$r.squared 
        adj.r.squared <- S1$adj.r.squared
        r.sqrd.out <- cbind()
        
        
        # - now we directly have the initial parameter vector
        mu <- matrix(S1$coefficients[,'Estimate'])
        mu <- matrix(mu[1:length(names.vector)])
        
        #give names to columns and rows
        if(isFALSE(hierarchical)){
          table.names <- names.vector
        }
        if(isTRUE(hierarchical)){
          table.names <- c(paste(level, names.vector, sep="."))
        }
        
        if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
          
          mu <- as.matrix(mu[1:(dim(mu)[1]-1)])
          
          #give names to columns and rows
          if(isFALSE(hierarchical)){
            table.names <- names.vector[1:(length(names.vector)-1)]
          }
          if(isTRUE(hierarchical)){
            table.names <- c(paste(level, names.vector[1:(length(names.vector)-1)], sep="."))
          }
        }
        
        rownames(mu) <- table.names
        colnames(mu) <- "mu0"
        
        #initial parameter vector for all variables (mu0)
        mu0 <- rbind(mu0, mu)
        
        
        # - and we directly have the initial prior variance matrix
        C <- as.matrix(vcov(lm.1))
        C <- C[1:length(names.vector),1:length(names.vector)]
        
        if(N.w!=time.in.year/2){
          if(loop==1){
            start = n.runs
          }else{
            start=end + 1
          }
          end <- start + 2*N.w
        }
        if(N.w==time.in.year/2){ #if it's the last harmonic wave possible to create (nyquist harmonic)
          
          C <- C[1:length(table.names),1:length(table.names)]
          
          if(loop==1){
            start = n.runs
          }else{
            start=end + 1
          }
          end <- start + 2*N.w-1
        }
        
        #prior variance for all variables (C0)
        C0[start:end,start:end] <- C 
        rownames(C0)[start:end] <- table.names  
        colnames(C0)[start:end] <- table.names
        C0 <- as.matrix(C0)
        
        # - Save all C matrices to see where there is NA values and replace them afterwards
        C.all[[paste(level, name, sep=".")]] <- C
        
        
        # - and also the observational variance
        e.all <- lm.1$residuals
        
        #get the variance per name
        var.resid <- var(na.omit(e.all))
        
        if(isFALSE(hierarchical)){
          a <- paste(name, sep='')
        }
        if(isTRUE(hierarchical)){
          a <- paste(level, name, sep='.')
        }
        
        for.V <- cbind(for.V, var.resid)
        colnames(for.V)[ncol(for.V)] <- a
        
        
        # Present the fit of this model
        r.sqrd.out <- cbind(r.sqrd.out, 'r.squared'=r.squared, 'adj.r.squared'=adj.r.squared)
        r.sqrd.out <- round(r.sqrd.out, round.by)
        if(!silent){
          print(r.sqrd.out)
        }
        
        if(plot.it == TRUE){
          if(!is.null(level3)){
            N.plots <- min(4,length(unique(D[,level3])))
            par(mfrow=c( max( (N.plots/2), 1) , 2))
            set.seed(42)
            random.IDs <- sample(x = unique(D[,level3]), size = min(4,length(unique(D[,level3]))), replace = FALSE)
            for(ID in random.IDs){ 
              ID.set <- subset(D, D[,level3] == ID)
              pred <- predict(object = lm.1, newdata=ID.set)
              plot(ID.set[,name]~ID.set[,time.var], xlab=time.var, ylab = name, main=ID)
              lines(pred~ID.set[,time.var], col='red')
              Residuals <- ID.set[,name] - pred
              hist(Residuals)
              abline(v=0, col='red', lwd=2)
            }
          }
          if(is.null(level3)){
            next
          }
        }
        
        print(paste('--Done! It took', difftime(Sys.time(), start.V, units = 'min'), 'minutes'))
      }
      
      
      if(N.w==0){
        
        #get original data
        D <- D.full
        D.A <- D
        
        #get the name we are referring to
        name <- relevant.names[loop]
        
        #get the time.var for the corresponding name
        time.var=time.var.all[loop]
        
        print('###############################################################')
        print(paste('DEFINING PARAMETERS RELATED TO relevant.name', '=', name))
        start.V <- Sys.time()
        
        # Remove outliers, based on overall mean and moving SD
        par(mfrow=c(1,2))
        ylim <- range(na.omit(D[,name]))
        Mean <- mean(na.omit(D[,name]))
        SD <- sd(na.omit(D[,name]))
        upper <- Mean + 3*SD #99.7% of data occurs within 3 standard deviations of the mean within a normal distribution
        lower <- Mean - 3*SD
        remove.i <- which(D[,name] > upper | D[,name] < lower)
        D[remove.i,name] <- NA
        par(mfrow=c(1,1))
        
        
        # Elements for the initial parameter vector (mu0)
        D.B <- subset(D, !is.na(D[,name]))
        x <- D.B[,time.var]
        y <- D.B[,name]
        Spline = smooth.spline(x = x, y = y)
        
        plot(y~x, xlab=time.var, ylab=name, main = paste("Spline", '=', name))
        lines(Spline, col='red', lwd=3)
        
        # - Save the spline function - we will need it for the Gt matrix
        if(isFALSE(hierarchical)){
          spline.name <- paste(name, '_Spline', sep='')
        }
        if(isTRUE(hierarchical)){
          spline.name <- paste(level, name, 'Spline', sep='_')
        }
        
        Spline.list[[length(Spline.list)+1]] <- Spline
        names(Spline.list)[length(Spline.list)] <- spline.name
        
        # - Finalize and save mu
        pred0 <- predict(object = Spline, x = (expected.start.time-1))$y
        mu <- c(pred0, 1)
        
        #give names to columns and rows
        if(isFALSE(hierarchical)){
          mu.names <- c(name, paste('d.', name, sep=''))      
        }
        if(isTRUE(hierarchical)){
          mu.names <- c(paste(level, name, sep="."),
                        paste('d', level, sep=".", name))}
        
        mu <- matrix(mu)
        row.names(mu) <- mu.names
        
        # - Initial parameter vector for all variables (mu0)
        mu0 <- rbind(mu0, mu)
        colnames(mu0) <- "mu0"
        
        
        # Elements for the prior variance (C0)
        # - Diff for the prior variance on the initial rate of change (difference between the present observation and the previous observation)
        Diff <- diff(D[,name])
        
        # - Remove outliers
        Q <- quantile(x = na.omit(Diff), probs = c(0.025, 0.975))
        Diff[which(Diff <= Q[1] | Diff >= Q[2])] <- NA
        D.A[,paste('d.', name, sep='')] <- c(NA,Diff)
        
        # - Re-order the columns
        D.A <- D.A[,c(metadata.names, name, paste('d', sep=".", name))]
        #only consider the start
        D.A <- subset(D.A, D.A[,time.var] %in% (expected.start.time):(expected.start.time+6))
        
        # - Finalize and save C0
        D.A <- na.omit(D.A)
        C <- cov(D.A[,c(name, paste('d', sep=".", name))])
        
        rm(D.A)
        
        if(loop==1){
          start = n.runs
        }else{
          start=end + 1
        }
        end <- start + 1
        
        # - Prior variance for all variables (C0)
        C0[start:end,start:end] <- C
        rownames(C0)[start:end] <- mu.names
        colnames(C0)[start:end] <- mu.names   
        C0 <- as.matrix(C0)
        
        # - Save all C matrices to see where there is NA values and replace them afterwards
        C.all[[paste(level, name, sep=".")]] <- C
        
        
        # Elements for the observational matrix (V)
        resid.ID.all <- c()
        
        for(ID in unique(D[,identifyer])){
          
          ID.set <- subset(D, D[,identifyer] == ID)
          
          # use a two-sided moving average to estimate V
          y <- ID.set[,name]
          x <- ID.set[,time.var]
          MA <- moving.function(x = y, n = 5, FUN = mean)
          
          # we find the residuals between the observed and filtered/estimated values
          resid.ID <- ID.set[,name] - MA
          
          # we combine the residuals for all individual IDs to estimate an overall V
          resid.ID.all <- c(resid.ID.all, resid.ID)
        }
        
        #get the variance per name
        var.resid <- var(na.omit(resid.ID.all))
        
        if(isFALSE(hierarchical)){
          a <- paste(name, sep='')
        }
        if(isTRUE(hierarchical)){
          a <- paste(level, name, sep='.')
        }
        
        for.V <- cbind(for.V, var.resid)
        colnames(for.V)[ncol(for.V)] <- a
        
        print(paste('--Done! It took', difftime(Sys.time(), start.V, units = 'min'), 'minutes'))
        
        gc()
      }
    }
  }
  
  return(list('mu0'=mu0,
              'C0'=C0,
              'C.all'=C.all,
              'for.V'=for.V,
              'Spline.list'=Spline.list,
              'end'=end))
}



# Function that does everything with nested get.DLM.parameters function ----
# - D is the data set, containing the time series we want to model
# - identifyer is the variable that identifies the production cycle/batch
# - time.var is the variable that describes the observation time
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on. 
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not
# - level2 is the second level (here correspond to the region)
# - level3 is the third level (here correspond to the farm/site)
# - best.n.waves is a vector with the best number of waves for each variable
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
# - expected.start.time is the observation time when we expect the model to start
# - no.better.limit defines how many iterations the EM algorithm must run after not having improved the first time
# - relevant.names are the names of the variables to be co-modeled
# - metadata.names are the variables' names we need to give context to our time series (specification/information about the data)
# - plot.it (TRUE of FALSE) determines whether plots should be made to show the fit of the model to the data of some unique identified units
# - remove.zeros whether observations in the time series with the value 0 should be removed. This is not usually recommended, unless 0s are stand-ins for NAs
# - round.by is the number of decimals you want to round the estimates
# - silent determines if the function should print the r^2 values of the fitted model to the console
define.DLM.parameters <- function(D, identifyer, time.var, time.in.year, hierarchical, 
                                  level2="local.authority", level3="site", 
                                  best.n.waves, trend, expected.start.time, no.better.limit, 
                                  relevant.names, metadata.names, plot.it=FALSE, remove.zeros=FALSE, 
                                  round.by=3, silent=TRUE){
  
  # Standardize the data (force the data to be standard normal distributed)
  SDs <- c()
  Means <- c()
  for(name in relevant.names){
    SDs <- c(SDs, sd(na.omit(D[,name])))
    Means <- c(Means, mean(na.omit(D[,name])))
  }
  
  for(name in relevant.names){
    Mean.name <- Means[which(relevant.names == name)]
    SD.name <- SDs[which(relevant.names == name)]
    D[,name] <- (D[,name] - Mean.name)/SD.name
  }
  
  Standarized.factors <- cbind(as.data.frame(Means), as.data.frame(SDs))
  rownames(Standarized.factors) <- relevant.names
  
  # Define the base frequency
  w <- (2*pi)/time.in.year  
  
  # Save original dataset and prepare empty list for results
  D.full <- D
  mu0.list <- list()
  C0.list <- list()
  V.list <- list()
  W.list <- list()
  Spline.list <- list()
  results.list <- list()
  
  if(hierarchical==FALSE){
    
    # Prepare for making mu0, C0, V and Splines (if needed)
    loop=0
    end <- c()
    n.runs <- 1
    mu0 <- c()
    n.zeros <- length(which(best.n.waves==0))
    n.na <- length(which(is.na(best.n.waves)))
    n.nyquist <- length(which((best.n.waves==time.in.year/2)))
    n.rows.cols <- sum(best.n.waves*2 + 1, na.rm=T) - n.nyquist + n.zeros + n.na*2
    C0 <- as.data.frame(matrix(0, ncol = n.rows.cols, nrow = n.rows.cols))
    e.all <- c()
    for.V <- cbind()
    time.var.all <- time.var
    C.all <- list()
    
    parameters_non_H <- get.DLM.parameters(D=D.full, identifyer, time.var=time.var.all, time.in.year=12,
                                           best.n.waves, relevant.names, metadata.names,
                                           mu0=mu0, C0=C0,  C.all=C.all, for.V=for.V, level2="local.authority", level3="site",
                                           hierarchical=FALSE, level=NULL, Spline.list=Spline.list,
                                           trend, w=w, loop=loop, n.runs=n.runs, end=end, expected.start.time, 
                                           plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE)
    
    # Get the parameters from parameters_non_H
    mu0 <- parameters_non_H$mu0
    C0 <- parameters_non_H$C0
    for.V <- parameters_non_H$for.V
    Spline.list <- parameters_non_H$Spline.list
    
    # Estimate the observational variance matrix for all variables (V), and save it in a list
    V <- diag(as.vector(for.V), nrow=length(for.V))
    colnames(V) <- colnames(for.V)
    rownames(V) <- colnames(for.V)
    
    # Arbitrarily define the initial W as C0/10, and save it in a list
    W <- C0/10
    
    
    ### Run the EM algorithm
    print('-------------------------------------------')
    divide.start <- Sys.time()
    D <- D.full
    
    # Randomly select a training set and a test set to determine when to stop the EM algorithm
    print('Dividing the data into a training and test set for EM ...')
    
    #get Learning and Test sets for EM algorithm
    N.EM <- round(3*(length(unique(D[order(D[,"date"]), "date"]))/4))
    sets.EM <- learning.test.sets.EM(D, relevant.names, N=N.EM)
    train.set.EM <- sets.EM[["Learning.set"]]
    test.set.EM <- sets.EM[["Test.set"]]
    
    # - we need all the training ID's to be elements of a list (that is how the EM algorithm takes them)
    start.makeList <- Sys.time()
    Des <- list()
    for(ID in sort(unique(train.set.EM[,identifyer]))){ 
      ID.set <- subset(train.set.EM, train.set.EM[,identifyer] == ID)
      Des[[length(Des)+1]] <- ID.set
    }
    
    train.set <- Des
    
    print(paste('--Done! It took', difftime(Sys.time(), divide.start, units = 'min'), 'minutes'))
    
    # Run the EM algorithm until it no longer improves the performance of the DLM
    print('Running the EM algorithm ...')
    
    time.var <- time.var.all
    no.better <- 0
    V = V
    W = W
    C <- C0
    mu <- mu0
    RMSE.best <- 9999999999
    total.steps <- 0
    V.best <- V
    W.best <- W
    C0.best <- C
    mu0.best <- mu
    while(no.better <= no.better.limit){
      
      # Run it on the training set
      # - it runs until converges on the best initial parameters
      a <- runEM(Des = train.set, mu0 = mu, C0 = C, V0 = V, W0 = W, time.var, time.in.year, 
                 best.n.waves, trend, relevant.names, Spline.list, steps = 1, silent = TRUE, 
                 DLM.version = FUNFILTERG, hierarchical=FALSE, level2=level2, level3=level3)
      V <- a$V[[1]]
      W <- a$W[[1]]
      
      # Try running the DLM with these variance parameters on all IDs in the test set
      results.list.all <- cbind()
      
      for(ID in sort(unique(test.set.EM[,identifyer]))){ 
        test.set.EM.i <- subset(test.set.EM, test.set.EM[,identifyer] == ID)
        res <- FUNFILTERG(D = test.set.EM.i, mu, C0 = C, V, W, time.var, time.in.year, 
                          best.n.waves, trend, relevant.names, Spline.list, 
                          hierarchical=FALSE, level2=level2, level3=level3)
        
        # get all results from res
        results.list <- extract.res(res, smot=NULL, hierarchical=FALSE, D = test.set.EM.i)
        results.list.all <- rbind(results.list.all, results.list)
      }
      
      et.all <- results.list.all[, grep("^[et_]", names(results.list.all), value=TRUE)]
      et.all <- as.data.frame(et.all)
      row.names(et.all) <- NULL
      
      # Calculate the RMSE and check if it is still improving
      RMSE <- round(sqrt(mean(unlist(na.omit(et.all))^2)), 4)
      
      if(RMSE >= RMSE.best){
        no.better <- no.better + 1
      }else{
        no.better <- 0
        RMSE.best <- RMSE
        V.best <- V
        W.best <- W
        C0.best <- C
        mu0.best <- mu
      }
      
      total.steps <- total.steps+1
      print(paste('Total steps:', total.steps, '| Current RMSE:', round(RMSE,4), '| Best RMSE:', round(RMSE.best,4)))
    }
    print('EM algorithm finished running')
    
    # Save the estimated variance components
    mu0.list <- mu0.best
    C0.list <- C0.best
    V.list <- V.best
    W.list <- W.best
    
    # Check the distribution of the standardized forecast errors on the test set
    # - it should be standard normally distributed!
    results.list.all <- cbind()
    
    for(CP in sort(unique(test.set.EM[,identifyer]))){
      test.set.EM.i <- subset(test.set.EM, test.set.EM[,identifyer] == CP)
      res <- FUNFILTERG(D = test.set.EM.i, mu0=mu0.best, C0=C0.best, V=V.best, W=W.best, time.var, time.in.year, 
                        best.n.waves, trend, relevant.names, Spline.list, 
                        hierarchical=FALSE, level2=level2, level3=level3)
      
      # get all results from res
      results.list <- extract.res(res, smot=NULL, hierarchical=FALSE, D = test.set.EM.i)
      results.list.all <- rbind(results.list.all, results.list)
    }
    
    ut.all <- results.list.all[, grep("^[ut_]", names(results.list.all), value=TRUE)]
    ut.all <- as.data.frame(ut.all)
    row.names(ut.all) <- NULL
    colnames(ut.all) <- relevant.names
    
    # - plot it and get the percentage outside of the 95 % CI
    par(mfrow=c(1,1))
    for(name in relevant.names){
      A <- na.omit(ut.all[,name])
      percent.outside.CI <- round(length(which(A < -1.96 | A > 1.96))/length(A),3) * 100
      hist(ut.all[,name],100, main=paste(name, '|', percent.outside.CI, '%'), xlim=c(-4,4), xlab="Standarized forecast errors")
      abline(v=0, col='red')
      abline(v=-1.96, col='blue')
      abline(v=1.96, col='blue')
    }
  }
  
  
  if(hierarchical==TRUE){
    
    # Prepare for making mu0, C0, V and Splines (if needed)
    loop=0
    end <- c()
    n.runs <- 1
    mu0 <- c()
    n.zeros <- length(which(best.n.waves==0))
    n.na <- length(which(is.na(best.n.waves)))
    n.nyquist <- length(which((best.n.waves==time.in.year/2)))
    n.rows.cols <- sum(best.n.waves*2 + 1, na.rm=T) - n.nyquist + n.zeros + n.na*2
    #since is hierarchical:
    n.rows.cols <- n.rows.cols + (n.rows.cols*length(unique(D.full[,level2]))) + (n.rows.cols*length(unique(D.full[,level3])))
    C0 <- as.data.frame(matrix(0, ncol = n.rows.cols, nrow = n.rows.cols))
    e.all <- c()
    for.V <- cbind()
    time.var.all <- time.var
    C.all <- list()
    
    
    print(paste('#### FOR COUNTRY'))
    
    parameters_H_Country <- get.DLM.parameters(D=D.full, identifyer, time.var=time.var.all, time.in.year=12,
                                               best.n.waves, relevant.names, metadata.names,
                                               mu0=mu0, C0=C0, C.all=C.all, for.V=for.V, level2="local.authority", level3="site",
                                               hierarchical=TRUE, level="Country", Spline.list=Spline.list,
                                               trend, w=w, loop=loop, n.runs=n.runs, end=end, expected.start.time, 
                                               plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE)
    
    # Get the parameters from parameters_H_Country
    mu0 <- parameters_H_Country$mu0
    C0 <- parameters_H_Country$C0
    for.V <- parameters_H_Country$for.V
    Spline.list <- parameters_H_Country$Spline.list
    end <- parameters_H_Country$end
    
    
    print(paste('#### FOR REGION'))
    
    for(region in sort(unique(D.full[,level2]))){
      
      D.region <- subset(D.full, D.full[,level2] == region)
      loop=0
      n.runs <- end + 1
      
      parameters_H_Region <- get.DLM.parameters(D=D.region, identifyer, time.var=time.var.all, time.in.year=12,
                                                best.n.waves, relevant.names, metadata.names,
                                                mu0=mu0, C0=C0,  C.all=C.all, for.V=for.V, level2="local.authority", level3="site",
                                                hierarchical=TRUE, level=region, Spline.list=Spline.list,
                                                trend, w=w, loop=loop, n.runs=n.runs, end=end, expected.start.time, 
                                                plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE)
      
      # Get the parameters from parameters_H_Region
      mu0 <- parameters_H_Region$mu0
      C0 <- parameters_H_Region$C0
      for.V <- parameters_H_Region$for.V
      Spline.list <- parameters_H_Region$Spline.list
      end <- parameters_H_Region$end
    }
    
    
    print(paste('#### FOR SITE'))
    
    C.all <- list()
    
    for(farm in sort(unique(D.full[,level3]))){
      
      D.farm <- subset(D.full, D.full[,level3] == farm)
      loop=0
      n.runs <- end + 1
      parameters_H_Farm <- get.DLM.parameters(D=D.farm, identifyer, time.var=time.var.all, time.in.year=12,
                                              best.n.waves, relevant.names, metadata.names,
                                              mu0=mu0, C0=C0,  C.all=C.all, for.V=for.V, level2="local.authority", level3="site",
                                              hierarchical=TRUE, level=farm, Spline.list=Spline.list,
                                              trend, w=w, loop=loop, n.runs=n.runs, end=end, expected.start.time, 
                                              plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE)
      
      # Get the parameters from parameters_H_Farm
      mu0 <- parameters_H_Farm$mu0
      C0 <- parameters_H_Farm$C0
      for.V <- parameters_H_Farm$for.V
      Spline.list <- parameters_H_Farm$Spline.list
      end <- parameters_H_Farm$end
      C.all <- parameters_H_Farm$C.all
    }
    
    # Replace the C0 of relevant names in the farms that are NA by the C0 of its region
    library(rlist)
    library(stringr)
    
    #just works for the variables without harmonic waves
    if(length(names(list.filter(C.all, anyNA(unlist(.))))) > 0){
      
      #farms and relevant names with NA
      names.na <- names(list.filter(C.all, anyNA(unlist(.))))
      
      for(i in names.na){
        
        #get region each farm belongs to
        farm.name <- sub("\\..*", "", i)
        D.farm <- subset(D.full, D.full$site == farm.name)
        region <- unique(D.farm$local.authority)
        
        #get relevant name with NA
        rel.name.na <- str_extract_all(i, relevant.names)
        rel.name.na <- as.character(rel.name.na[lengths(rel.name.na) != 0L])
        
        #get C0 for that region and relevant name
        region.rel.name <- paste(region, rel.name.na, sep=".")
        col.names.region <- c(region.rel.name,
                              paste('d', region.rel.name, sep=".")) 
        
        C0.region <- C0[col.names.region, col.names.region]
        
        # Make sure Ct is symmetrical
        C0.region <- (C0.region + t(C0.region))/2
        
        #put C0 of the region in the farms with NA, for the specific relevant name
        col.names <- c(i,
                       paste('d', i, sep="."))
        
        n.cols <- which(colnames(C0) %in% col.names)
        C0[n.cols,n.cols] <- C0.region
      }
    }
    
    # Estimate the observational variance matrix for all variables (V), and save it in a list
    V <- diag(as.vector(for.V), nrow=length(for.V))
    colnames(V) <- colnames(for.V)
    rownames(V) <- colnames(for.V)
    
    # Arbitrarily define the initial W as C0/10, and save it in a list
    W <- C0/10
    
    
    print(paste('ALL PARAMETERS CALCULATED, MOVE TO EM ALGORITHM'))
    
    ### Run the EM algorithm
    print('-------------------------------------------')
    divide.start <- Sys.time()
    D <- D.full
    
    # Randomly select a training set and a test set to determine when to stop the EM algorithm
    print('Dividing the data into a training and test set for EM ...')
    
    #get Learning and Test sets for EM algorithm
    N.EM <- round(3*(length(unique(D[order(D[,"date"]), "date"]))/4))
    sets.EM <- learning.test.sets.EM(D, relevant.names, N=N.EM)
    train.set.EM <- sets.EM[["Learning.set"]]
    test.set.EM <- sets.EM[["Test.set"]]
    
    #order the data frames per date
    train.set.EM <- train.set.EM[order(as.POSIXct(train.set.EM$date)),]
    test.set.EM <- test.set.EM[order(as.POSIXct(test.set.EM$date)),]
    
    print(paste('--Done! It took', difftime(Sys.time(), divide.start, units = 'min'), 'minutes'))
    
    # Run the EM algorithm until it no longer improves the performance of the DLM
    print('Running the EM algorithm ...')
    divide.start <- Sys.time()
    
    time.var <- time.var.all
    no.better <- 0
    V = V
    W = W
    C <- C0
    mu <- mu0
    RMSE.best <- 9999999999
    total.steps <- 0
    V.best <- V
    W.best <- W
    C0.best <- C
    mu0.best <- mu
    while(no.better <= no.better.limit){
      
      # Run it on the training set
      # - it runs until converges on the best initial parameters
      a <- runEM_H(Des = train.set.EM, mu0 = mu, C0 = C, V0 = V, W0 = W, time.var, time.in.year, 
                   best.n.waves, trend, relevant.names, Spline.list, steps = 1, silent = TRUE, 
                   DLM.version = FUNFILTERG_H, hierarchical=TRUE, level2=level2, level3=level3)
      
      Vn <- a$V[[1]]
      Wn <- a$W[[1]]
      
      # Try running the DLM with these variance parameters on all IDs in the test set
      res <- FUNFILTERG_H(D = test.set.EM, mu, C0 = C, V=try.V, W=try.W, time.var, time.in.year, 
                          best.n.waves, trend, relevant.names, Spline.list, 
                          hierarchical=TRUE, level2=level2, level3=level3)
      
      results.list.all <- extract.res(res, smot=NULL, hierarchical=TRUE, D = test.set.EM)
      
      et.all <- results.list.all[, grep("^[et_]", names(results.list.all), value=TRUE)]
      et.all <- as.data.frame(et.all)
      row.names(et.all) <- NULL
      
      
      # Calculate the RMSE and check if it is still improving
      RMSE <- round(sqrt(mean(unlist(na.omit(et.all))^2)), 4)
      
      if(RMSE >= RMSE.best){
        no.better <- no.better + 1
      }else{
        no.better <- 0
        RMSE.best <- RMSE
        V.best <- V
        W.best <- W
        C0.best <- C
        mu0.best <- mu
      }
      
      total.steps <- total.steps+1
      print(paste('Total steps:', total.steps, '| Current RMSE:', round(RMSE,4), '| Best RMSE:', round(RMSE.best,4)))
    }
    print(paste('--Done with the EM algorithm! It took', difftime(Sys.time(), divide.start, units = 'hours'), 'hours'))
    
    # Save the estimated variance components
    mu0.list <- mu0.best
    C0.list <- C0.best
    V.list <- V.best
    W.list <- W.best
    
    # Check the distribution of the standardized forecast errors on the test set
    # - it should be standard normally distributed!
    res <- FUNFILTERG_H(D = test.set.EM.i, mu0=mu0.best, C0=C0.best, V=V.best, W=W.best, time.var, time.in.year, 
                        best.n.waves, trend, relevant.names, Spline.list, 
                        hierarchical=TRUE,level2=level2, level3=level3)
    
    # get all results from res
    results.list.all <- extract.res(res, smot=NULL, hierarchical=TRUE, D = test.set.EM)
    
    ut.all <- results.list.all[, grep("^[ut_]", names(results.list.all), value=TRUE)]
    ut.all <- as.data.frame(ut.all)
    row.names(ut.all) <- NULL
    colnames(ut.all) <- relevant.names
    
    # - plot it and get the percentage outside of the 95 % CI
    par(mfrow=c(1,1))
    for(name in relevant.names){
      A <- na.omit(ut.all[,name])
      percent.outside.CI <- round(length(which(A < -1.96 | A > 1.96))/length(A),3) * 100
      hist(ut.all[,name],100, main=paste(name, '|', percent.outside.CI, '%'), xlim=c(-4,4), xlab="Standarized forecast errors")
      abline(v=0, col='red')
      abline(v=-1.96, col='blue')
      abline(v=1.96, col='blue')
    }
  }
  print('Finished')
  return(list(
    mu0.list=mu0.list,
    C0.list=C0.list,
    V.list=V.list,
    W.list=W.list,
    Spline.list=Spline.list,
    Standarized.factors=Standarized.factors,
    results.list.all = results.list.all))
}


################################################################################################################
### WHEN USING HIERARCHICAL MODELS ###

# Function to estimate the Yt (Observation vector) for hierarchical models ----
# - Date is the date you are currently modelling in the DLM
# - D is the data set, containing the time series we want to model
# - relevant.names are the names of the variables to be co-modeled
get.Yt.hierarchical <- function(Date, D, relevant.names){
  
  sites <- sort(unique(D$site))
  local.authorities <- sort(unique(D$local.authority))
  Date.set <- subset(D, D$date == Date)
  
  Yt <- c()
  Yt.names <- c()
  
  for (name in relevant.names) {
    mean.country <- mean(Date.set[,name], na.rm = T)
    Yt <- c(Yt, mean.country)
    name.country <- paste("Country", name, sep=".")
    Yt.names <- c(Yt.names, name.country)
  }
  
  for (region in local.authorities){
    
    if(region %in% Date.set$local.authority){
      D.region <- subset(Date.set, Date.set$local.authority == region)
      
      for (name in relevant.names) {
        mean.region <- mean(D.region[,name], na.rm = T)
        Yt <- c(Yt, mean.region)
        name.region <- paste(region, name, sep=".")
        Yt.names <- c(Yt.names, name.region)
      }
    }else{
      for (name in relevant.names) {
        mean.region <- NA
        Yt <- c(Yt, mean.region)
        name.region <- paste(region, name, sep=".")
        Yt.names <- c(Yt.names, name.region)
      }
    }
  }
  
  for (farm in sites){
    
    if(farm %in% Date.set$site){
      D.farm <- subset(Date.set, Date.set$site == farm)
      
      for (name in relevant.names) {
        mean.farm <- D.farm[,name]
        Yt <- c(Yt, mean.farm)
        name.farm <- paste(farm, name, sep=".")
        Yt.names <- c(Yt.names, name.farm)
      }
    }else{
      for (name in relevant.names) {
        mean.farm <- NA
        Yt <- c(Yt, mean.farm)
        name.farm <- paste(farm, name, sep=".")
        Yt.names <- c(Yt.names, name.farm)
      }
    }
  }
  
  Yt <- as.matrix(Yt)
  rownames(Yt) <- Yt.names
  
  return(Yt)
}



# Function to estimate the initalize W (system variance) in hierarchical models ----
# - D is the data set, containing the time series we want to model
# - W is the system variance
# - Date is the date you are currently modelling in the DLM
# - name is/are the names of the variables to be co-modeled
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not
# - level3 is the third level (here correspond to the farm/site)
## Get the system variance (Wt), replace by 0 for the farm-dates when
## there is no information about the farm.
getWt <- function(D, W, Date=NA, name=relevant.names[length(relevant.names)], hierarchical, level3){
  
  initial.W <- W
  
  if(hierarchical==FALSE){
    Wt <- initial.W
  }
  
  if(hierarchical==TRUE){
    Wt <- initial.W
    
    for(farm in sort(unique(D[,level3]))){
      
      D.farm <- subset(D, D[,level3] == farm)
      
      #get row number
      i <- which(D.farm$date == Date)
      
      if(Date %in% D.farm$date){
        
        Wt <- Wt
        
      }else{
        
        #replace W for that farm on that time by a matrix of 0's when there is no information for that farm-date
        rows.cols.names <- c(paste(farm, name, sep="."), paste("d", farm, name, sep="."))
        Wt[rows.cols.names, rows.cols.names] <- matrix(nrow=2, ncol=2, 0)
      }
    }
  }
  return(as.matrix(Wt))
}



# Function to initialized mt for hierarchical models ----
## when a new production cycle starts the mt for that farm should be initialized (=mu0) 
# - mt is the filtered mean
# - mu0 is the initial parameter vector
# - D is the data set, containing the time series we want to model
# - Date is the date you are currently modelling in the DLM
# - var are the names of the variables to be co-modeled
# - level3 is the third level (here correspond to the farm/site)
get.mt.hierarchical <- function(mt, mu0, D, Date, var=relevant.names[length(relevant.names)], level3){
  
  #get right time.var
  i=which(relevant.names==var)
  time=time.var[i]
  
  #dataset with the observation only for that Date
  Date.set <- subset(D, D$date == Date)
  
  #dataset and farms that started a new production cycle on that Date
  Date.first.set <- subset(Date.set, Date.set[, time]==0)
  farms.first <- Date.first.set[,level3]
  
  mt.h <- mt
  
  if(length(farms.first)!=0){
    
    for(farm in farms.first){ #farm=farms.first[1]
      
      #get mu0 for those farms and var
      farm.name.mu0.l <- paste(farm, var, sep=".")
      farm.name.mu0.t <- paste("d", farm, var, sep=".")
      
      #replace mt for those farms and var by mu0 (initialize)
      mt.h[c(farm.name.mu0.l, farm.name.mu0.t),] <- as.matrix(mu0[c(farm.name.mu0.l, farm.name.mu0.t),])
    }
  }
  return(mt.h)
}



# Function to initialized Ct for hierarchical models ----
## when a new production cycle starts the Ct for that farm should be initialized (=C0) and
## covariances should have the same correlation but different magnitudes
# - Ct is the filtered variance
# - C0 is the prior variance
# - D is the data set, containing the time series we want to model
# - Date is the date you are currently modelling in the DLM
# - var are the names of the variables to be co-modeled
# - level3 is the third level (here correspond to the farm/site)
get.Ct.hierarchical <- function(Ct, C0, D, Date, var=relevant.names[length(relevant.names)], level3){
  
  #get right time.var
  i=which(relevant.names==var)
  time=time.var[i]
  
  #dataset with the observation only for that Date
  Date.set <- subset(D, D$date == Date)
  
  #dataset and farms that started a new production cycle on that Date (or started having information on that day)
  Date.first.set <- subset(Date.set, Date.set[, time]==0)
  farms.first <- Date.first.set[,level3]
  
  Ct.h <- Ct
  
  if(length(farms.first)!=0){
    
    for(farm in farms.first){
      
      #get C0 for those farms and var
      farm.name.C0.l <- paste(farm, var, sep=".")
      farm.name.C0.t <- paste("d", farm, var, sep=".")
      
      #replace Ct for those farms and var by C0 (initialize)
      # - replace 4 "diagonal" values
      Ct.h[c(farm.name.C0.l, farm.name.C0.t),c(farm.name.C0.l, farm.name.C0.t)] <- 
        as.matrix(C0[c(farm.name.C0.l, farm.name.C0.t),c(farm.name.C0.l, farm.name.C0.t)])
      
      #replace covariances
      # - get variance of C0 and Ct
      C0.a <- as.matrix(C0[farm.name.C0.l, farm.name.C0.l]) #variance of level
      C0.b <- as.matrix(C0[farm.name.C0.t, farm.name.C0.t]) #variance of trend
      Ct.a <- as.matrix(Ct[farm.name.C0.l, farm.name.C0.l]) #variance of level
      Ct.b <- as.matrix(Ct[farm.name.C0.t, farm.name.C0.t]) #variance of trend
      
      # - start replacing the rows and the columns that have the level (a)
      # - first rows and columns upper left side (1)
      cov.a1 <- t(as.matrix(Ct.h[farm.name.C0.l, 1:((which(colnames(Ct.h)==farm.name.C0.l))-1)]))
      rownames(cov.a1) <- farm.name.C0.l
      
      for(i in 1:length(cov.a1)){
        Ct.h[farm.name.C0.l,i] <- (cov.a1[1,i] %*% sqrt(C0.a)) / sqrt(Ct.a)
        Ct.h[i,farm.name.C0.l] <- (cov.a1[1,i] %*% sqrt(C0.a)) / sqrt(Ct.a)
      }
      
      # - second rows and columns down right side (2)
      if(farm.name.C0.l != rownames(C0)[(dim(C0)[1])-1]){ #only if is not the last farm (nothing after it)
        cov.a2 <- t(as.matrix(Ct.h[farm.name.C0.l, ((which(colnames(Ct.h)==farm.name.C0.t))+1):dim(Ct.h)[1]]))
        rownames(cov.a2) <- farm.name.C0.l
        
        for(i in 1:length(cov.a2)){
          row.name <- colnames(cov.a2)[i] 
          Ct.h[farm.name.C0.l,row.name] <- (cov.a2[1,i] %*% sqrt(C0.a)) / sqrt(Ct.a)
          Ct.h[row.name,farm.name.C0.l] <- (cov.a2[1,i] %*% sqrt(C0.a)) / sqrt(Ct.a)
        }
      }
      
      # - now replace the rows and the columns that have the trend (b)
      # - first rows and columns upper left side (1)
      cov.b1 <- t(as.matrix(Ct.h[farm.name.C0.t, 1:((which(colnames(Ct.h)==farm.name.C0.l))-1)]))
      rownames(cov.b1) <- farm.name.C0.t
      
      for(i in 1:length(cov.b1)){
        Ct.h[farm.name.C0.t,i] <- (cov.b1[1,i] %*% sqrt(C0.b)) / sqrt(Ct.b)
        Ct.h[i,farm.name.C0.t] <- (cov.b1[1,i] %*% sqrt(C0.b)) / sqrt(Ct.b)
      }
      
      # - second rows and columns down right side (2)
      if(farm.name.C0.l != rownames(C0)[(dim(C0)[1])-1]){ #only if is not the last farm (nothing after it)
        cov.b2 <- t(as.matrix(Ct.h[farm.name.C0.t, ((which(colnames(Ct.h)==farm.name.C0.t))+1):dim(Ct.h)[1]]))
        rownames(cov.b2) <- farm.name.C0.t
        
        for(i in 1:length(cov.b2)){
          row.name <- colnames(cov.b2)[i] 
          Ct.h[farm.name.C0.t,row.name] <- (cov.b2[1,i] %*% sqrt(C0.b)) / sqrt(Ct.b)
          Ct.h[row.name,farm.name.C0.t] <- (cov.b2[1,i] %*% sqrt(C0.b)) / sqrt(Ct.b)
        }
      }
    }
  }
  return(Ct.h)
}


# Function to estimate the VSum for EM algorithm (hierarchical) ----
# - Yt is the observation vector
## Get a matrix with 1 in the cells where the observation 
## contributes to the observation variance-covariance matrix.
## Other cells are 0 (When it's NA). Used by the EM-algorithm
## for hierarchical model
getVSumElement_H = function(Yt) {
  VSE = matrix(1, nrow=length(Yt), ncol=length(Yt))
  for(i in 1:length(Yt)){
    if (is.na(Yt[i])) { 
      VSE[i, ] = 0
      VSE[, i] = 0
    }
  }
  return(VSE)
}



# Function for the DLM for hierarchical models ----
# - D is the data set, containing the time series we want to model
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - V is the observation variance
# - W is the system variance
# - time.var is the variable that describes the observation time
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on. 
# - best.n.waves is a vector with the best number of waves for each variable
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
# - relevant.names are the names of the variables to be co-modeled
# - Spline.list is a list with the splines calculated for each farm
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not, here is hierarchical=TRUE
# - level2 is the second level (correspond to the region, here is level2="local.authority" because is a hierarchical model)
# - level3 is the third level (correspond to the farm/site, here is level3="site" because is a hierarchical model)
FUNFILTERG_H <- function(D, mu0, C0, V, W, time.var, time.in.year, best.n.waves, trend, relevant.names, 
                         Spline.list, hierarchical, level2="local.authority", level3="site"){
  
  n <- 1:length(unique(D$date))
  first.date <- as.Date(D$date)[1]
  
  Yt.list <- list()
  at.list <- list()		 
  Rt.list <- list()
  ft.list <- list()
  Qt.list <- list()
  At.list <- list()
  et.list <- list()
  ut.list <- list()
  mt.list <- list()
  Ct.list <- list()
  Ft.list <- list()
  VSE.list <- list()
  Gt.list <- list()
  fullFt.list <- list()
  
  mt <- mu0				
  Ct <- C0
  
  # Make sure Ct is symmetrical
  Ct <- (Ct + t(Ct))/2
  
  for(i in n){
    
    Date <- sort(unique(D$date))[i]
    Date_1 <- sort(unique(D$date))[i-1]
    if(length(Date_1)==0){
      Date_1 <- 0
    }
    
    # Get the observation vector
    Yt <- get.Yt.hierarchical(Date, D, relevant.names)
    
    # Get V-sum-element, to be used in the EM-algorithm
    VSE <- getVSumElement_H(Yt)
    
    # Get the observational variance (Vt)
    Vt <- V
    
    # Get the system variance (Wt), replace by 0 for the dates where there is no information about the farm
    Wt <- getWt(D, W, Date, name=relevant.names[length(relevant.names)], hierarchical, level3)
    
    # Get Gt - independent of the current Yt
    Gt <- getGt(D, i=NA, Date, Date_1, time.var, time.in.year, hierarchical, level2, 
                level3, best.n.waves, trend, relevant.names, Spline.list)
    
    # Get Ft given the current Yt 
    Ft <- getFt(D, time.in.year, hierarchical, level2, level3, best.n.waves, relevant.names, trend)
    
    # Remove missing values from Yt, Ft and Vt
    missing = which(is.na(Yt))
    fullFt <- Ft
    fullVt <- Vt
    
    # If length of missing is > 0 there is at least one missing
    if (length(missing) > 0) {
      
      # remove from Yt
      a <- rownames(Yt)[-missing]
      Yt = as.matrix(Yt[-missing])
      rownames(Yt) <- a
      
      # Remove from Ft
      Ft = as.matrix(Ft[, -missing])
      
      # Remove from Vt
      Vt = Vt[-missing, -missing]
    }
    
    #if the only observation is NA (univariate case) or all observations are NA (multivariate case)
    if (length(Yt) == 0) {
      
      #Initialize mt and Ct when a new production cycle starts or when starts to have information for one farm
      if(Date > first.date){
        mt <- get.mt.hierarchical(mt, mu0, D, Date, var=relevant.names[length(relevant.names)], level3)
        Ct <- get.Ct.hierarchical(Ct, C0, D, Date, var=relevant.names[length(relevant.names)], level3)
      }
      
      at <- Gt %*%  mt		                   
      rownames(at) <- rownames(mt)
      Rt <- Gt %*%  Ct %*% t(Gt) + Wt       
      mt=at
      Ct=Rt
      ft <- t(fullFt) %*% at
      Qt <- t(fullFt) %*% Rt %*% fullFt + fullVt
      At=NA
      et=NA
      ut=NA
      
    }else{
      
      #Initialize mt and Ct when a new production cycle starts or when starts to have information for one farm
      if(Date > first.date){
        mt <- get.mt.hierarchical(mt, mu0, D, Date, var=relevant.names[length(relevant.names)], level3)
        Ct <- get.Ct.hierarchical(Ct, C0, D, Date, var=relevant.names[length(relevant.names)], level3)
      }
      
      # Run the Kalman filter 
      at <- Gt %*%  mt		                 
      rownames(at) <- rownames(mt)
      
      Rt <- Gt %*%  Ct %*% t(Gt) + Wt       
      
      ft <- t(Ft) %*% at		      	        
      Qt <- t(Ft) %*% Rt %*%  Ft + Vt   
      
      At <- Rt %*% Ft %*% solve(Qt)       
      et <- Yt  - ft	                      
      ut <- et / sqrt(diag(Qt))          
      
      # - update the parameter vector and variance matrix
      mt <- at + At %*% et                  
      Ct <- Rt - At  %*% Qt %*% t(At)	      
    }
    
    # Make sure Ct is symmetrical
    Ct <- (Ct + t(Ct))/2
    
    # Save the values in lists
    Yt.list[[i]] <- Yt
    at.list[[i]] <- at
    Rt.list[[i]] <- Rt
    ft.list[[i]] <- ft
    Qt.list[[i]] <- Qt
    At.list[[i]] <- At
    et.list[[i]] <- et
    ut.list[[i]] <- ut
    mt.list[[i]] <- mt
    Ct.list[[i]] <- Ct
    Ft.list[[i]] <- t(Ft)
    VSE.list[[i]] <- VSE
    Gt.list[[i]] <- Gt
    fullFt.list[[i]] <- t(fullFt)
  }
  
  return(list(
    Yt=Yt.list,
    at=at.list,
    Rt=Rt.list,
    ft=ft.list,
    Qt=Qt.list,
    At=At.list,
    et=et.list,
    ut=ut.list,
    mt=mt.list,
    Ct=Ct.list,
    F=Ft.list,
    vse=VSE.list,
    Gt.list=Gt.list,
    fullFt.list=fullFt.list
    
  ))
}



# Function for the EM algorithm for hierarchical models ----
# - Des is the data set for the EM algorithm
# - mu0 is the initial parameter vector
# - C0 is the prior variance
# - V0 is the observation variance
# - W0 is the system variance
# - time.var is the variable that describes the observation time
# - time.in.year is the base frequency. For a model related to 24 hour periods, time.in.year=24, for a model related to one year (365 days) periods, time.in.year=365, and so on. 
# - best.n.waves is a vector with the best number of waves for each variable
# - trend is whether or not (TRUE or FALSE) the model should include a linear trend for the mean level over time
# - relevant.names are the names of the variables to be co-modeled
# - Spline.list is a list with the splines calculated for each farm
# - steps is the number of iterations
# - silent is to put equal to TRUE of FALSE whether you want to print the step you are in or not
# - DLM.version is the name of the function that runs the DLM
# - hierarchical is to put equal to TRUE of FALSE whether it is a multi-level model or not, here is hierarchical=TRUE
# - level2 is the second level (correspond to the region, here is level2="local.authority" because is a hierarchical model)
# - level3 is the third level (correspond to the farm/site, here is level3="site" because is a hierarchical model)
runEM_H = function(Des, mu0, C0, V0, W0, time.var, time.in.year, best.n.waves, trend, 
                   relevant.names, Spline.list, steps = 1, silent = TRUE, 
                   DLM.version = FUNFILTERG_H, hierarchical=hierarchical, level2="local.authority", level3="site") {
  Vs = list()
  Ws = list()
  Vn = V0
  Wn = W0
  Cn <- diag(0, length(mu0))
  mu0n <- rep(0, length(mu0))
  
  #Choose wich version of DLM to use
  runDLM <- DLM.version
  
  # Iterate over steps
  for (s in 1:steps) {
    print(paste("Step", s))
    
    # Set sums and counts to zero
    # Count structure for observation variance
    n.rows <- length(relevant.names)+length(relevant.names)*length(unique(Des[,level2]))+
      length(relevant.names)*length(unique(Des[,level3]))
    sumV = matrix(0, nrow=n.rows, ncol=n.rows)
    # Sum for observation variance
    sumObs = matrix(0, nrow=n.rows, ncol=n.rows)
    # Count for system variance
    sumW = 0
    # Sum for system variance
    sumSys = matrix(0, nrow=length(mu0), ncol=length(mu0))
    
    # Apply Em to Des
    
    # Expectation step - filter and smooth
    if(silent == FALSE){
      print(b)
    }
    res <- runDLM(D = Des, mu0 = mu0, C0 = C0, V = Vn, W = Wn, time.var, time.in.year, 
                  best.n.waves, trend, relevant.names, Spline.list, hierarchical, level2, level3)
    smot = runSmoother(res)
    
    # Get the smoothened C0 from this one
    Cn <- Cn + smot$Cts[,,1]
    mu0n <-  mu0n + smot$mts[,,1]
    
    # Set contributions to sums and counts to zero for this D
    bSumV = matrix(0, nrow=n.rows, ncol=n.rows)
    sumCountV = matrix(0, nrow=n.rows, ncol=n.rows)
    bSumW = matrix(0, nrow=length(mu0), ncol=length(mu0))
    sumCountW = 0
    
    # Sum contributions
    n = length(res$mt)
    
    # Iterate over time within D 
    for (t in 1:n) {
      
      # Get the Gt 
      Gt <- res$Gt.list[[t]]
      
      # Only if observations at all
      if (length(res$Yt[[t]]) > 0) {
        
        # Observation variance
        
        # Find the contribution to the sum even though it does not have the correct dimension
        Vcont = res$F[[t]]%*%smot$Cts[,,t]%*%t(res$F[[t]]) + 
          (res$Y[[t]] - res$F[[t]]%*%smot$mts[,,t])%*%t((res$Y[[t]] - res$F[[t]]%*%smot$mts[,,t]))
        # Get the pointer matrix
        vse = res$vse[[t]]
        # Create a full matrix
        Vfull = matrix(0, nrow=n.rows, ncol=n.rows)
        rownames(Vfull) <- rownames(V)
        colnames(Vfull) <- colnames(V)
        
        # Find out which cells to enter
        ind = c()
        for (i in 1:dim(Vfull)[1]) { 
          if (vse[i,i] > 0) {
            ind = c(ind, i)
          }
        }
        
        # Enter the contribution into the right cells
        a <- rownames(Vfull)[ind]
        Vfull[a,a] = Vcont[a,a]
        
        # Add the resulting
        bSumV = bSumV + Vfull
        # Adjust the counts matrix
        sumCountV = sumCountV + vse          
        
        
        # System variance
        if (t > 1) {
          # Find the contribution - the dimension is always correct
          bSumW = bSumW + smot$Lt[,,t] + 
            (smot$mts[,,t] - Gt%*%smot$mts[,,t-1])%*%t(smot$mts[,,t] - Gt%*%smot$mts[,,t-1])            
          sumCountW = sumCountW + 1
        }
      }
    }
    
    # Check for negative variances
    ignore = FALSE
    for (j in 1:n.rows) {
      if (bSumV[j, j] < 0) {
        # Adjust to 0
        bSumV[j, j] = 0
        bSumV[j, ] = 0
        bSumV[, j] = 0
        # Print a comment
        if (! silent) print(paste("Negative contribution to observation variance", j, "for D", b))
      }
    }
    
    # Check for negative variances
    ignore = FALSE
    for (j in 1:length(mu0)) {
      if (bSumW[j, j] < 0) {
        # Adjust to 0
        bSumW[j, j] = 0
        bSumW[, j] = 0
        bSumW[j, ] = 0
        # Print a message
        if (! silent) print(paste("Negative contribution to system variance", j, "for D", b))
      }
    }
    
    # This will never happen (used for debugging)
    if (ignore) {
      print(paste("Contribution to observation variance ignored from D", b))
      for (i in 1:n.rows) {
        for (j in 1:n.rows) {
          bSumV[i, j] = 0
        }
      }
      sumCountV = matrix(0, nrow=n.rows, ncol=n.rows)
    }
    # Add the contribution from the D to the total
    sumObs = sumObs + bSumV
    sumV = sumV + sumCountV
    
    # This will never happen (used for debugging)
    if (ignore) {
      print(paste("Contribution to system variance ignored from D", b))
      for (i in 1:length(mu0)) {
        for (j in 1:length(mu0)) {
          bSumW[i, j] = 0
        }
      }
      sumCountW = 0
    }
    # Add the contribution from the D to the total
    sumSys = sumSys + bSumW
    sumW = sumW + sumCountW
    
    # Normalize by counts
    Vn = sumObs/sumV
    Vnn = Vn
    # - if there is NaN on Vn: some farms don't exist at the same time so covariance is unknown
    if(sum(is.na(Vnn))!=0){
      na.all <- which(is.nan(Vnn), arr.ind=TRUE)
      names <- unique(sub("^[^.]+.","", rownames(na.all))) 
      
      for(name in names){
        names.na <- rownames(na.all)[endsWith(rownames(na.all), name)] 
        if(length(names)>1){ 
          na <- na.all[-which(rownames(na.all) == setdiff(rownames(na.all), names.na)),]
        }else{
          na <- na.all
        }
        
        #remove repeated/opposite rows and cols from the matrix
        naa <- na
        for(i in 1:dim(na)[1]){
          try({cel <- naa[i,]}, silent=TRUE) 
          if(length(which(naa[,1]==cel[2] & naa[,2]==cel[1]))!=0){
            naa <- naa[-which(naa[,1]==cel[2] & naa[,2]==cel[1]),]
          } else {
            naa <- naa
          }
        }
        
        for(i in 1:dim(naa)[1]){
          
          cel <- naa[i,]
          farm.r.i <- rownames(Vnn)[cel[1]]
          farm.c.i <- rownames(Vnn)[cel[2]]
          farm.r <- sub("\\..*", "", farm.r.i)
          farm.c <- sub("\\..*", "", farm.c.i)
          region.r <- unique(subset(Des,Des[,level3]==farm.r)[,level2])
          region.c <- unique(subset(Des,Des[,level3]==farm.c)[,level2])
          
          if(isTRUE(region.r==region.c)){
            
            #farms on the same region as farm.c and farm.r that the covariance 
            #with farm.c and farm.r is not NaN and are not farm.c/farm.r
            farms.same.region <- unique(subset(Des,Des[,level2]==region.c)[,level3])
            farms.same.region.nan.farm.c.i <- names(which(is.nan(Vnn[paste(farms.same.region, name, sep="."), farm.c.i])))
            farms.same.region.nan.farm.r.i <- names(which(is.nan(Vnn[paste(farms.same.region, name, sep="."), farm.r.i])))
            farms.same.region.nan.farm.c <- sub("\\..*", "", farms.same.region.nan.farm.c.i)
            farms.same.region.nan.farm.r <- sub("\\..*", "", farms.same.region.nan.farm.r.i)
            farms.same.region.not.nan.farm.c <- farms.same.region[-which(farms.same.region %in% farms.same.region.nan.farm.c)]
            farms.same.region.not.nan.farm.c <- farms.same.region.not.nan.farm.c[-which(farms.same.region.not.nan.farm.c == farm.c)] #remove farm.c
            farms.same.region.not.nan.farm.r <- farms.same.region[-which(farms.same.region %in% farms.same.region.nan.farm.r)]
            farms.same.region.not.nan.farm.r <- farms.same.region.not.nan.farm.r[-which(farms.same.region.not.nan.farm.r == farm.r)] #remove farm.r
            
            save.corr.farms.same.region.farm.c <- list() 
            loop=0
            #loop over all farms that are from the same region as farm.c and farm.r and the covariance with farm.c is not NaN and are not farm.c
            #to get all the correlations of the covariances between farm.c and all other farms from the same region != NaN
            for(farm.not.nan.farm.c in farms.same.region.not.nan.farm.c){ 
              loop=loop+1
              farm.not.nan.farm.c.i <- paste(farm.not.nan.farm.c, name, sep=".")
              
              # - calculate the correlation of the covariance of those farms (farm.not.nan.farm.c and farm.c)
              corr.farm.same.region.farm.c <- Vnn[farm.not.nan.farm.c.i, farm.c.i]/(sqrt(Vnn[farm.not.nan.farm.c.i,farm.not.nan.farm.c.i])*sqrt(Vnn[farm.c.i,farm.c.i]))
              save.corr.farms.same.region.farm.c[[loop]] <- corr.farm.same.region.farm.c
            }
            
            save.corr.farms.same.region.farm.r <- list()
            loop=0
            #loop over all farms that are from the same region as farm.r and farm.c and the covariance with farm.r is not NaN and are not farm.r
            #to get all the correlations of the covariances between farm.r and all other farms from the same region != NaN
            for(farm.not.nan.farm.r in farms.same.region.not.nan.farm.r){ 
              loop=loop+1
              farm.not.nan.farm.r.i <- paste(farm.not.nan.farm.r, name, sep=".")
              
              # - calculate the correlation of the covariance of those farms (farm.not.nan.farm.r and farm.r)
              corr.farm.same.region.farm.r <- Vnn[farm.not.nan.farm.r.i, farm.r.i]/(sqrt(Vnn[farm.not.nan.farm.r.i,farm.not.nan.farm.r.i])*sqrt(Vnn[farm.r.i,farm.r.i]))
              save.corr.farms.same.region.farm.r[[loop]] <- corr.farm.same.region.farm.r
            }
            
            #get mean of all correlations of all covariances 
            #(all farms from the same region with farm.c and all farms from the same region with farm.r)
            corr.final <- mean(c(unlist(save.corr.farms.same.region.farm.c), unlist(save.corr.farms.same.region.farm.r)))
            
            #calculate the covariance with the right correlation between farm.r and farm.c
            cov <- corr.final*sqrt(Vnn[farm.r.i,farm.r.i])*sqrt(Vnn[farm.c.i,farm.c.i])
            
            # - replace NaN farm values in Vn for cov
            Vn[cel[1], cel[2]] <- cov
            Vn[cel[2], cel[1]] <- cov
            
          } else { #if the farms belong to different regions
            
            #farms on the same region as farm.c/farm.r that are not NaN (do it by columns) and are not farm.c/farm.r
            farms.same.region.c <- unique(subset(Des,Des[,level2]==region.c)[,level3])
            farms.same.region.r <- unique(subset(Des,Des[,level2]==region.r)[,level3])
            farms.region.c.nan.farm.r.i <- names(which(is.nan(Vnn[paste(farms.same.region.c, name, sep="."), farm.r.i])))
            farms.region.r.nan.farm.c.i <- names(which(is.nan(Vnn[paste(farms.same.region.r, name, sep="."), farm.c.i])))
            farms.region.c.nan.farm.r <- sub("\\..*", "", farms.region.c.nan.farm.r.i)
            farms.region.r.nan.farm.c <- sub("\\..*", "", farms.region.r.nan.farm.c.i)
            
            if(length(farms.region.c.nan.farm.r)!=0){
              farms.region.c.not.nan.farm.r <- farms.same.region.c[-which(farms.same.region.c %in% farms.region.c.nan.farm.r)]
            }
            if(length(farms.region.r.nan.farm.c)!=0){
              farms.region.r.not.nan.farm.c <- farms.same.region.r[-which(farms.same.region.r %in% farms.region.r.nan.farm.c)]
            }
            if(length(farms.region.c.nan.farm.r)==0){ 
              farms.region.c.not.nan.farm.r <- farms.same.region.c
            }
            if(length(farms.region.r.nan.farm.c)==0){
              farms.region.r.not.nan.farm.c <- farms.same.region.r
            }
            
            save.corr.farms.region.c.farm.r <- list() 
            loop=0
            
            #loop over all farms that are from region.c and the covariance with farm.r is not NaN
            #to get all the correlations of the covariances between farms for region.c and farm.r
            for(farm.region.c in farms.region.c.not.nan.farm.r){ #farm.region.c=farms.region.c.not.nan.farm.r[1]
              loop=loop+1
              farm.region.c.i <- paste(farm.region.c, name, sep=".")
              
              # - calculate the correlation of the covariance of those farms (farms.region.c and farm.r)
              corr.farm.region.c.farm.r <- Vnn[farm.region.c.i, farm.r.i]/(sqrt(Vnn[farm.region.c.i,farm.region.c.i])*sqrt(Vnn[farm.r.i,farm.r.i]))
              save.corr.farms.region.c.farm.r[[loop]] <- corr.farm.region.c.farm.r
            }
            
            save.corr.farms.region.r.farm.c <- list() 
            loop=0
            #loop over all farms that are from region.r and the covariance with farm.c is not NaN
            #to get all the correlations of the covariances between farms for region.r and farm.c
            for(farm.region.r in farms.region.r.not.nan.farm.c){ #farm.region.r=farms.region.r.not.nan.farm.c[14]
              loop=loop+1
              farm.region.r.i <- paste(farm.region.r, name, sep=".")
              
              # - calculate the correlation of the covariance of those farms (farms.region.r and farm.c)
              corr.farm.region.r.farm.c <- Vnn[farm.region.r.i, farm.c.i]/(sqrt(Vnn[farm.region.r.i,farm.region.r.i])*sqrt(Vnn[farm.c.i,farm.c.i]))
              save.corr.farms.region.r.farm.c[[loop]] <- corr.farm.region.r.farm.c
            }
            
            #get mean of all correlations of all covariances 
            #(all farms from region.c with farm.r and all farms from region.r with farm.c)
            corr.final <- mean(c(unlist(save.corr.farms.region.c.farm.r), unlist(save.corr.farms.region.r.farm.c)))
            
            #calculate the covariance with the right correlation between farm.r and farm.c
            cov <- corr.final*sqrt(Vnn[farm.r.i,farm.r.i])*sqrt(Vnn[farm.c.i,farm.c.i])
            
            # - replace NaN farm values in Vn for cov
            Vn[cel[1], cel[2]] <- cov
            Vn[cel[2], cel[1]] <- cov
          }
        }
      }
    }
    
    # Normalize by counts
    Wn = sumSys/sumW
    Cn <- Cn/n
    mu0n <-  mu0n/n
    # Make sure they are symmetric
    Vn = (Vn + t(Vn))/2
    Wn = (Wn + t(Wn))/2
    Cn <- (Cn + t(Cn))/2
    Vn <<- Vn
    Wn <<- Wn
    # Make sure they have the right names
    colnames(Vn) <- colnames(V)
    rownames(Vn) <- rownames(V)
    
    colnames(Wn) <- rownames(mu0)
    rownames(Wn) <- rownames(mu0)
    
    colnames(Cn) <- rownames(mu0)
    rownames(Cn) <- rownames(mu0)
    
    mu0n <- as.matrix(mu0n)
    rownames(mu0n) <- rownames(mu0)
    
    # Save them in lists
    Vs[[s]] = Vn
    Ws[[s]] = Wn
  }
  return(list(V=Vs, W=Ws, Cn=Cn, mu0n=mu0n, smot=smot))
}
