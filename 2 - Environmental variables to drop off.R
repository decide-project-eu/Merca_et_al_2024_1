#### Several non-hierarchical DLMs to see which environmental variables actually improve the model 

#libraries
library(MASS)

#get data
setwd("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/Datasets")
load("df_salmon_20.RData")
df <- df_salmon_20


### Define the constant relevant columns in the data set ----

#  -- identifyer is the variables that identifies production cycle or the site
identifyer =  "nseq"

# -- relevant.names are the names of the columns to be co-modeled
relevant.names <- c("d1.temp", "d9.temp", "log.max.daily.range.temp",
                    "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal",
                    "log1.d9.phyc", "log.d9.chl",
                    "log.d1.do", "log.d9.do", "log.max.daily.range.do",
                    "log.d9.prep", "log.d9.dino", "log.d9.diato",
                    "log.d9.nano", "log.d9.pico",
                    "d1.ph", "d9.ph", "max.daily.range.ph",
                    "log.d9.no3", "log.max.daily.range.no3",
                    "log0.mortality.rel.20")

#  -- expected.start.time is the observation time when we expect the model to start
expected.start.time=0

# -- no.better.limit defines how many iterations the EM algorithm must run after not having improved the first time
no.better.limit=0

# -- metadata.names are the variables' names we need to give context to our time series (specification/information about the data)
metadata.names <- c("date", "site", "local.authority", "sealochs", "lat", "lon", "nseq", "production.cycle", "production.year",
                    "months.since.start", "n.months")


### Create Learning and Test set ----

source("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/Functions for Salmon study.R")

#get Learning and Test sets
#use the first ~3/4 (that will be 70%) of time (dates) for Learning.set
N <- round(3*(length(unique(df[order(df[,"date"]), "date"]))/4))
sets <- get.learning.test.sets(df, relevant.names, hierarchical=FALSE, N=N)
Learning.set <- sets[["Learning.set"]]
Test.set <- sets[["Test.set"]]


## Create several DMLs, mortality with 1 environmental group ----

# Make sure you avoid scientific notation in your output - it will make things easier!
options(scipen=999)

# Get the functions needed
function.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio"
source(paste(function.dir, '/Functions for Salmon study.R', sep=''))

# Get variables
group.env.variables <- c("temp", "sal", "phy_types", "do", "prep", "ph", "no3")

initial.relevant.names <- relevant.names
RMSE.all.l <- cbind()
RMSE.all.t <- cbind()

for(var in group.env.variables){ 
  
  progress <- round(which(group.env.variables == var)/length(group.env.variables)*100,2)
  print(paste('Applying DLM to', var, ':', progress, '%'))
  
  #get relevant names
  relevant.names <- c(initial.relevant.names[grepl(var, initial.relevant.names)], "log0.mortality.rel.20")
  
  if(length(relevant.names)==1){
    relevant.names <- c("log1.d9.phyc", "log.d9.chl", "log.d9.dino", "log.d9.diato", 
                        "log.d9.nano", "log.d9.pico", "log0.mortality.rel.20")
  }
  if(var=="ph"){
    relevant.names <- c("d1.ph", "d9.ph", "max.daily.range.ph", "log0.mortality.rel.20")
  }
  
  # -- best.n.waves are the best number of waves for each variable (0 means that is not going to be modeled using HW)
  df_best_number_waves <- data.frame()
  loop=0
  for(i in relevant.names){
    loop=loop+1
    out <- get.best.number.waves(D=df, time.in.year=12, trend = FALSE, relevant.name=i, 
                                 time.var="month", stratify.by="site", plot.it=TRUE, 
                                 remove.zeros=FALSE, round.by=2, parameters_per_stratify.by=FALSE)
    out <- as.data.frame(out)
    df_best_number_waves <- rbind(df_best_number_waves, out)
    rownames(df_best_number_waves)[loop] <- i
  }
  
  if("log.max.daily.range.sal" %in% relevant.names){
    #as salinity does not have any seasonality, it has 0 harmonic waves
    df_best_number_waves["log.max.daily.range.sal", "Best.N.w"] <- 0 
  }
  
  best.n.waves <- as.vector(df_best_number_waves$Best.N.w)
  best.n.waves[length(best.n.waves)] <- NA #make sure mortality has NA waves - do it with splines (not 0 to be different from salinity)
  
  # -- time.var is the variable that describes the observation time
  time.var <- c(rep("month", length(relevant.names)))
  time.var[which(is.na(best.n.waves))] <- "months.since.start"  #months.since.start for mortality
  time.var[which(best.n.waves==0)] <- "n.months"                #n.months for salinity
  
  
  ### Learn all parameters for the DLM - parameters per name in relevant.names
  Start.time <- Sys.time()
  
  out <- define.DLM.parameters(D=Learning.set, identifyer, time.var, time.in.year=12, hierarchical=FALSE,
                               level2=NULL, level3=NULL, best.n.waves, trend=FALSE, expected.start.time=0, no.better.limit=0,
                               relevant.names, metadata.names, plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE)
  print(Sys.time()-Start.time)
  
  # - Save the relevant parameters based on the learning set
  model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Learning.set"
  file.name <- paste("mortality", "all", var, "non-H.RDS", sep="_")
  saveRDS(object = out, file = paste(model.dir, '/', file.name, sep=''))
  
  ### Get RMSE (applied to the Learning set)
  
  et.all.l <-  out$results.list.all[, paste("et_","log0.mortality.rel.20", sep="")]
  et.all.l <- na.omit(et.all.l)
  RMSE.l <- sqrt(mean(et.all.l^2)) 
  RMSE.all.l <- rbind(RMSE.all.l, RMSE.l)
  
  
  ### Apply to the Test set
  runDLM <- FUNFILTERG
  D <- Test.set
  
  # The same standardization factors apply to the test set
  Standarized.factors <- out$Standarized.factors
  for(name in relevant.names){
    Mean.name <- Standarized.factors[name, "Means"]
    SD.name <- Standarized.factors[name, "SDs"]
    D[,name] <- (D[,name] - Mean.name)/SD.name
  }
  
  # Get the relevant DLM parameters
  mu0 <- out$mu0.list
  C0 <- out$C0.list
  V <- out$V.list
  W <- out$W.list
  Spline.list <- out$Spline.list
  
  # Apply the model for each nseq (production cycle)
  start.time <- Sys.time()
  results.list.all <- cbind()
  
  for (ID in unique(D[,identifyer])){
    
    progress <- round(which(unique(D[,identifyer]) == ID)/length(unique(D[,identifyer]))*100,2)
    
    # Get the data for the current nseq
    ID.set <- subset(D, D[,identifyer] == ID)
    
    # run the DLM
    res <- runDLM(D = ID.set, mu0 = mu0, C0 = C0, V = V, W = W, time.var, time.in.year=12, 
                  best.n.waves, trend=FALSE, relevant.names = relevant.names, Spline.list,
                  hierarchical=FALSE, level2=NULL, level3=NULL)
    
    # get all results from res
    results.list <- extract.res(res, smot=NULL, hierarchical=FALSE, D = ID.set)
    results.list.all <- rbind(results.list.all, results.list) 
  }
  
  # add meta data to the results.list.all
  results.test.set <- cbind(D[,metadata.names], results.list.all)
  
  
  ### Get RMSE (applied to the Test set!)
  
  et.all.t <- results.test.set[, paste("et_","log0.mortality.rel.20", sep="")]
  et.all.t <- na.omit(et.all.t)
  RMSE.t <- sqrt(mean(et.all.t^2)) 
  RMSE.all.t <- rbind(RMSE.all.t, RMSE.t)
  
  # - Save the results and RMSE applied to the Test set
  model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set"
  file.name <- paste("mortality", "all", var, "non-H.RDS", sep="_")
  saveRDS(object = results.test.set, file = paste(model.dir, '/', file.name, sep=''))
}

# RMSE applied to the Learning set
row.names(RMSE.all.l) <- group.env.variables
colnames(RMSE.all.l) <- "RMSE"
RMSE.all.l <- rbind(RMSE.all.l, "mortality"=0.8568) #univariate non-hierarchical DLM (only with mortality)
RMSE.all.l <- as.data.frame(RMSE.all.l)
for(i in 1:c(dim(RMSE.all.l)[1]-1)){
  if(RMSE.all.l[i,"RMSE"] > RMSE.all.l[dim(RMSE.all.l)[1],"RMSE"]){
    RMSE.all.l[i, "Improved"] <- "FALSE"
  }else{
    RMSE.all.l[i, "Improved"] <- "TRUE"
  }
}
View(RMSE.all.l)

# RMSE applied to the Test set
row.names(RMSE.all.t) <- group.env.variables
colnames(RMSE.all.t) <- "RMSE"
RMSE.all.t <- rbind(RMSE.all.t, "mortality"=0.8602789) #univariate non-hierarchical DLM (only with mortality)
RMSE.all.t <- as.data.frame(RMSE.all.t)
for(i in 1:c(dim(RMSE.all.t)[1]-1)){
  if(RMSE.all.t[i,"RMSE"] > RMSE.all.t[dim(RMSE.all.t)[1],"RMSE"]){
    RMSE.all.t[i, "Improved"] <- "FALSE"
  }else{
    RMSE.all.t[i, "Improved"] <- "TRUE"
  }
}
View(RMSE.all.t)


### Save it 
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Learning.set"
saveRDS(object = RMSE.all.l, file = paste(model.dir, '/', "RMSE.all.l.RDS", sep=''))
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set"
saveRDS(object = RMSE.all.t, file = paste(model.dir, '/', "RMSE.all.t.RDS", sep=''))
#salinity was the best



## Create several DMLs, mortality+sal with 1 other environmental group ----

# Get variables
group.env.variables <- c("temp", "do", "phy_types", "prep", "ph", "no3") #without sal
RMSE.all.2.l <- cbind()
RMSE.all.2.t <- cbind()

for(var in group.env.variables){ 
  
  progress <- round(which(group.env.variables == var)/length(group.env.variables)*100,2)
  print(paste('Applying DLM to', var, ':', progress, '%'))
  
  #get relevant names
  relevant.names <- c(initial.relevant.names[grepl(var, initial.relevant.names)],
                      "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal",
                      "log0.mortality.rel.20")
  
  if(var=="phy_types"){
    relevant.names <- c("log1.d9.phyc", "log.d9.chl", "log.d9.dino", "log.d9.diato", 
                        "log.d9.nano", "log.d9.pico", 
                        "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal", "log0.mortality.rel.20")
  }
  
  if(var=="ph"){
    relevant.names <- c("d1.ph", "d9.ph", "max.daily.range.ph",
                        "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal", "log0.mortality.rel.20")
  }
  
  # -- best.n.waves are the best number of waves for each variable (0 means that is not going to be modeled using HW)
  df_best_number_waves <- data.frame()
  loop=0
  for(i in relevant.names){
    loop=loop+1
    out <- get.best.number.waves(D=df, time.in.year=12, trend = FALSE, relevant.name=i, 
                                 time.var="month", stratify.by="site", plot.it=TRUE, 
                                 remove.zeros=FALSE, round.by=2, parameters_per_stratify.by=FALSE)
    out <- as.data.frame(out)
    df_best_number_waves <- rbind(df_best_number_waves, out)
    rownames(df_best_number_waves)[loop] <- i
  }
  
  if("log.max.daily.range.sal" %in% relevant.names){
    #as salinity does not have any seasonality, it has 0 harmonic waves
    df_best_number_waves["log.max.daily.range.sal", "Best.N.w"] <- 0 
  }
  
  best.n.waves <- as.vector(df_best_number_waves$Best.N.w)
  best.n.waves[length(best.n.waves)] <- NA #make sure mortality has NA waves - do it with splines (not 0 to be different from salinity)
  
  # -- time.var is the variable that describes the observation time
  time.var <- c(rep("month", length(relevant.names)))
  time.var[which(is.na(best.n.waves))] <- "months.since.start"  #months.since.start for mortality
  time.var[which(best.n.waves==0)] <- "n.months"                #n.months for salinity
  
  
  ### Learn all parameters for the DLM - parameters per name in relevant.names
  Start.time <- Sys.time()
  
  out <- define.DLM.parameters(D=Learning.set, identifyer, time.var, time.in.year=12, hierarchical=FALSE,
                               level2=NULL, level3=NULL, best.n.waves, trend=FALSE, expected.start.time=0, no.better.limit=0,
                               relevant.names, metadata.names, plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE)
  print(Sys.time()-Start.time)
  
  # - Save the relevant parameters based on the learning set
  model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Learning.set"
  file.name <- paste("mortality+sal", "all", var, "non-H.RDS", sep="_")
  saveRDS(object = out, file = paste(model.dir, '/', file.name, sep=''))
  
  ### Get RMSE (applied to the Learning set)
  
  et.all.2.l <-  out$results.list.all[, paste("et_","log0.mortality.rel.20", sep="")]
  et.all.2.l <- na.omit(et.all.2.l)
  RMSE.2.l <- sqrt(mean(et.all.2.l^2)) 
  RMSE.all.2.l <- rbind(RMSE.all.2.l, RMSE.2.l)
  
  
  ### Apply to the Test set
  runDLM <- FUNFILTERG
  D <- Test.set
  
  # The same standardization factors apply to the test set
  Standarized.factors <- out$Standarized.factors
  for(name in relevant.names){
    Mean.name <- Standarized.factors[name, "Means"]
    SD.name <- Standarized.factors[name, "SDs"]
    D[,name] <- (D[,name] - Mean.name)/SD.name
  }
  
  # Get the relevant DLM parameters
  mu0 <- out$mu0.list
  C0 <- out$C0.list
  V <- out$V.list
  W <- out$W.list
  Spline.list <- out$Spline.list
  
  # Apply the model for each nseq (production cycle)
  start.time <- Sys.time()
  results.list.all <- cbind()
  
  for (ID in unique(D[,identifyer])){ 
    
    progress <- round(which(unique(D[,identifyer]) == ID)/length(unique(D[,identifyer]))*100,2)
    
    # Get the data for the current nseq
    ID.set <- subset(D, D[,identifyer] == ID)
    
    # run the DLM
    res <- runDLM(D = ID.set, mu0 = mu0, C0 = C0, V = V, W = W, time.var, time.in.year=12, 
                  best.n.waves, trend=FALSE, relevant.names = relevant.names, Spline.list,
                  hierarchical=FALSE, level2=NULL, level3=NULL)
    
    # get all results from res
    results.list <- extract.res(res, smot=NULL, hierarchical=FALSE, D = ID.set)
    results.list.all <- rbind(results.list.all, results.list) 
  }
  
  # add meta data to the results.list.all
  results.test.set <- cbind(D[,metadata.names], results.list.all)
  
  
  ### Get RMSE (applied to the Test set!)
  
  et.all.2.t <- results.test.set[, paste("et_","log0.mortality.rel.20", sep="")]
  et.all.2.t <- na.omit(et.all.2.t)
  RMSE.2.t <- sqrt(mean(et.all.2.t^2)) 
  RMSE.all.2.t <- rbind(RMSE.all.2.t, RMSE.2.t)
  
  # - Save the results and RMSE applied to the Test set
  model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set"
  file.name <- paste("mortality+sal", "all", var, "non-H.RDS", sep="_")
  saveRDS(object = results.test.set, file = paste(model.dir, '/', file.name, sep=''))
}

# RMSE applied to the Learning set
row.names(RMSE.all.2.l) <- group.env.variables
colnames(RMSE.all.2.l) <- "RMSE"
RMSE.all.2.l <- rbind(RMSE.all.2.l, "mortality+sal"=0.8572404)
RMSE.all.2.l <- as.data.frame(RMSE.all.2.l)
for(i in 1:c(dim(RMSE.all.2.l)[1]-1)){
  if(RMSE.all.2.l[i,"RMSE"] > RMSE.all.2.l[dim(RMSE.all.2.l)[1],"RMSE"]){
    RMSE.all.2.l[i, "Improved"] <- "FALSE"
  }else{
    RMSE.all.2.l[i, "Improved"] <- "TRUE"
  }
}
View(RMSE.all.2.l)

# RMSE applied to the Test set
row.names(RMSE.all.2.t) <- group.env.variables
colnames(RMSE.all.2.t) <- "RMSE"
RMSE.all.2.t <- rbind(RMSE.all.2.t, "mortality+sal"=0.8586170)
RMSE.all.2.t <- as.data.frame(RMSE.all.2.t)
for(i in 1:c(dim(RMSE.all.2.t)[1]-1)){
  if(RMSE.all.2.t[i,"RMSE"] > RMSE.all.2.t[dim(RMSE.all.2.t)[1],"RMSE"]){
    RMSE.all.2.t[i, "Improved"] <- "FALSE"
  }else{
    RMSE.all.2.t[i, "Improved"] <- "TRUE"
  }
}
View(RMSE.all.2.t) #did not improve


### Save it 
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Learning.set"
saveRDS(object = RMSE.all.2.l, file = paste(model.dir, '/', "RMSE.all.2.l.RDS", sep=''))
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set"
saveRDS(object = RMSE.all.2.t, file = paste(model.dir, '/', "RMSE.all.2.t.RDS", sep=''))



### Create several DMLs, what environmental variables from sal + mortality to keep ----

model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set"
RMSE.all.2.t <- readRDS(file = paste(model.dir, '/RMSE.all.2.t.RDS', sep=''))

#keep mortality + the group of variables with smaller RMSE
RMSE.all.2.t #mortality + sal
env.var.to.keep <- c("log.d1.sal", "log.d9.sal", "log.max.daily.range.sal")

#get all combinations of variables (mortality+sal)
v1 <- env.var.to.keep
vars.comb <- do.call("c", lapply(seq_along(v1), function(i) combn(v1, i, FUN = list)))
RMSE.sal.mort.all.l <- cbind()
RMSE.sal.mort.all.t <- cbind()

#run all the possible DLMs with mortality+sal variables
for(i in 1:length(vars.comb)){ 
  
  progress <- round(i/length(vars.comb)*100,2)
  print(paste("Progress:" ,progress, '%', sep=""))
  
  # -- get relevant names for a specific combination
  relevant.names <- c(vars.comb[[i]], "log0.mortality.rel.20")
  
  # get the correspondent number for the variables used
  ns <- which(env.var.to.keep %in% vars.comb[[i]])
  ns <- c(as.integer(ns), 4)
  
  # -- best.n.waves are the best number of waves for each variable (0 means that is not going to be modeled using HW)
  df_best_number_waves <- data.frame()
  loop=0
  for(i in relevant.names){
    loop=loop+1
    out <- get.best.number.waves(D=df, time.in.year=12, trend = FALSE, relevant.name=i, 
                                 time.var="month", stratify.by="site", plot.it=TRUE, 
                                 remove.zeros=FALSE, round.by=2, parameters_per_stratify.by=FALSE)
    out <- as.data.frame(out)
    df_best_number_waves <- rbind(df_best_number_waves, out)
    rownames(df_best_number_waves)[loop] <- i
  }
  if("log.max.daily.range.sal" %in% relevant.names){
    #as salinity does not have any seasonality, it has 0 harmonic waves
    df_best_number_waves["log.max.daily.range.sal", "Best.N.w"] <- 0 
  }
  
  best.n.waves <- as.vector(df_best_number_waves$Best.N.w)
  best.n.waves[length(best.n.waves)] <- NA #make sure mortality has NA waves - do it with splines (not 0 to be different from salinity)
  
  # -- time.var is the variable that describes the observation time
  time.var <- c(rep("month", length(relevant.names)))
  time.var[which(is.na(best.n.waves))] <- "months.since.start"  #months.since.start for mortality
  time.var[which(best.n.waves==0)] <- "n.months"                #n.months for salinity
  
  
  ### Learn all parameters for the DLM - parameters per name in relevant.names
  Start.time <- Sys.time()
  
  out <- define.DLM.parameters(D=Learning.set, identifyer, time.var, time.in.year=12, hierarchical=FALSE,
                               level2=NULL, level3=NULL, best.n.waves, trend=FALSE, expected.start.time=0, no.better.limit=0,
                               relevant.names, metadata.names, plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE)
  print(Sys.time()-Start.time)
  
  # - Save the relevant parameters based on the learning set
  model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Learning.set"
  file.name <- paste("mortality+sal", paste(ns, collapse="+"), "non-H.RDS", sep="_")
  saveRDS(object = out, file = paste(model.dir, '/', file.name, sep=''))
  
  
  ### Get RMSE (applied to the Learning set)
  
  et.sal.mort.l <-  out$results.list.all[, paste("et_","log0.mortality.rel.20", sep="")]
  et.sal.mort.l <- na.omit(et.sal.mort.l)
  RMSE.sal.mort.l <- sqrt(mean(et.sal.mort.l^2)) 
  RMSE.sal.mort.all.l <- rbind(RMSE.sal.mort.all.l, RMSE.sal.mort.l)
  
  ### Apply to the Test set
  runDLM <- FUNFILTERG
  D <- Test.set
  
  # The same standardization factors apply to the test set
  Standarized.factors <- out$Standarized.factors
  for(name in relevant.names){
    Mean.name <- Standarized.factors[name, "Means"]
    SD.name <- Standarized.factors[name, "SDs"]
    D[,name] <- (D[,name] - Mean.name)/SD.name
  }
  
  # Get the relevant DLM parameters
  mu0 <- out$mu0.list
  C0 <- out$C0.list
  V <- out$V.list
  W <- out$W.list
  Spline.list <- out$Spline.list
  
  # Apply the model for each nseq (production cycle)
  start.time <- Sys.time()
  results.list.all <- cbind()
  
  for (ID in unique(D[,identifyer])){
    
    progress <- round(which(unique(D[,identifyer]) == ID)/length(unique(D[,identifyer]))*100,2)
    
    # Get the data for the current nseq
    ID.set <- subset(D, D[,identifyer] == ID)
    
    # run the DLM
    res <- runDLM(D = ID.set, mu0 = mu0, C0 = C0, V = V, W = W, time.var, time.in.year=12, 
                  best.n.waves, trend=FALSE, relevant.names = relevant.names, Spline.list,
                  hierarchical=FALSE, level2=NULL, level3=NULL)
    
    # get all results from res
    results.list <- extract.res(res, smot=NULL, hierarchical=FALSE, D = ID.set)
    results.list.all <- rbind(results.list.all, results.list) 
  }
  
  # add meta data to the results.list.all
  results.test.set <- cbind(D[,metadata.names], results.list.all)
  
  
  ### Get RMSE (applied to the Test set!)
  
  et.sal.mort.t <- results.test.set[, paste("et_","log0.mortality.rel.20", sep="")]
  et.sal.mort.t <- na.omit(et.sal.mort.t)
  RMSE.sal.mort.t <- sqrt(mean(et.sal.mort.t^2)) 
  RMSE.sal.mort.all.t <- rbind(RMSE.sal.mort.all.t, RMSE.sal.mort.t)
  
  # - Save the results and RMSE applied to the Test set
  model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set"
  file.name <- paste("mortality+sal", paste(ns, collapse="+"), "non-H.RDS", sep="_")
  saveRDS(object = results.test.set, file = paste(model.dir, '/', file.name, sep=''))
}

# RMSE applied to the Learning set
vars.comb.all <- vars.comb
for(i in 1:length(vars.comb)){
  vars.comb.all[[i]] <- c(vars.comb[[i]], "log0.mortality.rel.20")
}
row.names(RMSE.sal.mort.all.l) <- vars.comb.all
colnames(RMSE.sal.mort.all.l) <- "RMSE"
RMSE.sal.mort.all.l <- as.data.frame(RMSE.sal.mort.all.l)
RMSE.sal.mort.all.l[,"What's the best"] <- rep(NA, dim(RMSE.sal.mort.all.l)[1])
RMSE.sal.mort.all.l[which(RMSE.sal.mort.all.l[,"RMSE"]==min(RMSE.sal.mort.all.l[,"RMSE"])), "What's the best"] <- "Best"

View(RMSE.sal.mort.all.l)

# RMSE applied to the Test set
row.names(RMSE.sal.mort.all.t) <- vars.comb.all
colnames(RMSE.sal.mort.all.t) <- "RMSE"
RMSE.sal.mort.all.t <- as.data.frame(RMSE.sal.mort.all.t)
RMSE.sal.mort.all.t[,"What's the best"] <- rep(NA, dim(RMSE.sal.mort.all.t)[1])
RMSE.sal.mort.all.t[which(RMSE.sal.mort.all.t[,"RMSE"]==min(RMSE.sal.mort.all.t[,"RMSE"])), "What's the best"] <- "Best"

View(RMSE.sal.mort.all.t)

### BEST COMBINATION ###
#"log.d9.sal", "log.max.daily.range.sal" and "log0.mortality.rel.20" (vars 2 , 3 and 4)

### Save it 
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Learning.set"
saveRDS(object = RMSE.sal.mort.all.l, file = paste(model.dir, '/', "RMSE.sal.mort.all.l.RDS", sep=''))
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set"
saveRDS(object = RMSE.sal.mort.all.t, file = paste(model.dir, '/', "RMSE.sal.mort.all.t.RDS", sep=''))



## Create several DMLs, log.d9.sal + log.max.daily.range.sal + log0.mortality.rel.20 with 1 var at a time ----

#remove mortality + log.d9.sal + log.max.daily.range.sal
relevant.names <- c("d1.temp", "d9.temp", "log.max.daily.range.temp",
                    "log.d1.sal", "log.d9.sal", "log.max.daily.range.sal",
                    "log1.d9.phyc", "log.d9.chl",
                    "log.d1.do", "log.d9.do", "log.max.daily.range.do",
                    "log.d9.prep", "log.d9.dino", "log.d9.diato",
                    "log.d9.nano", "log.d9.pico",
                    "d1.ph", "d9.ph", "max.daily.range.ph",
                    "log.d9.no3", "log.max.daily.range.no3",
                    "log0.mortality.rel.20")
vars.to.test <- relevant.names[-c(5,6,22)]
RMSE.mort.9sal.max.d.r.sal.all.l <- cbind()
RMSE.mort.9sal.max.d.r.sal.all.t <- cbind()

#run all the possible DLMs with log.d9.sal + log.max.daily.range.sal + log0.mortality.rel.20 + 1 var at a time
for(var in vars.to.test){
  
  progress <- round(which(var==vars.to.test)/length(vars.to.test)*100,2)
  print(paste("Progress: " ,progress, '%', sep=""))
  
  # -- get relevant names for a specific combination
  relevant.names <- c(var, "log.d9.sal", "log.max.daily.range.sal", "log0.mortality.rel.20")
  
  # get the correspondent number for the var used
  ns <- which(var==vars.to.test)
  
  # -- best.n.waves are the best number of waves for each variable (0 means that is not going to be modeled using HW)
  df_best_number_waves <- data.frame()
  loop=0
  for(i in relevant.names){
    loop=loop+1
    out <- get.best.number.waves(D=df, time.in.year=12, trend = FALSE, relevant.name=i, 
                                 time.var="month", stratify.by="site", plot.it=TRUE, 
                                 remove.zeros=FALSE, round.by=2, parameters_per_stratify.by=FALSE)
    out <- as.data.frame(out)
    df_best_number_waves <- rbind(df_best_number_waves, out)
    rownames(df_best_number_waves)[loop] <- i
  }
  
  if("log.max.daily.range.sal" %in% relevant.names){
    #as salinity does not have any seasonality, it has 0 harmonic waves
    df_best_number_waves["log.max.daily.range.sal", "Best.N.w"] <- 0 
  }
  
  best.n.waves <- as.vector(df_best_number_waves$Best.N.w)
  best.n.waves[length(best.n.waves)] <- NA #make sure mortality has NA waves - do it with splines (not 0 to be different from salinity)
  
  # -- time.var is the variable that describes the observation time
  time.var <- c(rep("month", length(relevant.names)))
  time.var[which(is.na(best.n.waves))] <- "months.since.start"  #months.since.start for mortality
  time.var[which(best.n.waves==0)] <- "n.months"                #n.months for salinity
  
  
  ### Learn all parameters for the DLM - parameters per name in relevant.names
  Start.time <- Sys.time()
  
  out <- define.DLM.parameters(D=Learning.set, identifyer, time.var, time.in.year=12, hierarchical=FALSE,
                               level2=NULL, level3=NULL, best.n.waves, trend=FALSE, expected.start.time=0, no.better.limit=0,
                               relevant.names, metadata.names, plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE)
  print(Sys.time()-Start.time)
  
  # - Save the relevant parameters based on the learning set
  model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Learning.set"
  file.name <- paste("mortality+9sal+max.d.r.sal+", paste(ns, collapse="+"), "_non-H.RDS", sep="")
  saveRDS(object = out, file = paste(model.dir, '/', file.name, sep=''))
  
  
  ### Get RMSE (applied to the Learning set)
  
  et.mort.9sal.max.d.r.sal.l <-  out$results.list.all[, paste("et_","log0.mortality.rel.20", sep="")]
  et.mort.9sal.max.d.r.sal.l <- na.omit(et.mort.9sal.max.d.r.sal.l)
  RMSE.mort.9sal.max.d.r.sal.l <- sqrt(mean(et.mort.9sal.max.d.r.sal.l^2)) 
  RMSE.mort.9sal.max.d.r.sal.all.l <- rbind(RMSE.mort.9sal.max.d.r.sal.all.l, RMSE.mort.9sal.max.d.r.sal.l)
  
  
  ### Apply to the Test set
  runDLM <- FUNFILTERG
  D <- Test.set
  
  # The same standardization factors apply to the test set
  Standarized.factors <- out$Standarized.factors
  for(name in relevant.names){
    Mean.name <- Standarized.factors[name, "Means"]
    SD.name <- Standarized.factors[name, "SDs"]
    D[,name] <- (D[,name] - Mean.name)/SD.name
  }
  
  # Get the relevant DLM parameters
  mu0 <- out$mu0.list
  C0 <- out$C0.list
  V <- out$V.list
  W <- out$W.list
  Spline.list <- out$Spline.list
  
  # Apply the model for each nseq (production cycle)
  start.time <- Sys.time()
  results.list.all <- cbind()
  
  for (ID in unique(D[,identifyer])){
    
    progress <- round(which(unique(D[,identifyer]) == ID)/length(unique(D[,identifyer]))*100,2)
    
    # Get the data for the current nseq
    ID.set <- subset(D, D[,identifyer] == ID)
    
    # run the DLM
    res <- runDLM(D = ID.set, mu0 = mu0, C0 = C0, V = V, W = W, time.var, time.in.year=12, 
                  best.n.waves, trend=FALSE, relevant.names = relevant.names, Spline.list,
                  hierarchical=FALSE, level2=NULL, level3=NULL)
    
    # get all results from res
    results.list <- extract.res(res, smot=NULL, hierarchical=FALSE, D = ID.set)
    results.list.all <- rbind(results.list.all, results.list) 
  }
  
  # add meta data to the results.list.all
  results.test.set <- cbind(D[,metadata.names], results.list.all)
  
  
  ### Get RMSE (applied to the Test set!)
  
  et.mort.9sal.max.d.r.sal.t <- results.test.set[, paste("et_","log0.mortality.rel.20", sep="")]
  et.mort.9sal.max.d.r.sal.t <- na.omit(et.mort.9sal.max.d.r.sal.t)
  RMSE.mort.9sal.max.d.r.sal.t <- sqrt(mean(et.mort.9sal.max.d.r.sal.t^2)) 
  RMSE.mort.9sal.max.d.r.sal.all.t <- rbind(RMSE.mort.9sal.max.d.r.sal.all.t, RMSE.mort.9sal.max.d.r.sal.t)
  
  # - Save the results and RMSE applied to the Test set
  model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set"
  file.name <- paste("mortality+9sal+max.d.r.sal+", paste(ns, collapse="+"), "_non-H.RDS", sep="")
  saveRDS(object = results.test.set, file = paste(model.dir, '/', file.name, sep=''))
}

# RMSE applied to the Learning set
row.names(RMSE.mort.9sal.max.d.r.sal.all.l) <- vars.to.test
colnames(RMSE.mort.9sal.max.d.r.sal.all.l) <- "RMSE"
RMSE.mort.9sal.max.d.r.sal.all.l <- rbind(RMSE.mort.9sal.max.d.r.sal.all.l, "mortality+9sal+max.d.r.sal"=0.8573477)
RMSE.mort.9sal.max.d.r.sal.all.l <- as.data.frame(RMSE.mort.9sal.max.d.r.sal.all.l)
for(i in 1:c(dim(RMSE.mort.9sal.max.d.r.sal.all.l)[1]-1)){
  if(RMSE.mort.9sal.max.d.r.sal.all.l[i,"RMSE"] > RMSE.mort.9sal.max.d.r.sal.all.l[dim(RMSE.mort.9sal.max.d.r.sal.all.l)[1],"RMSE"]){
    RMSE.mort.9sal.max.d.r.sal.all.l[i, "Improved"] <- "FALSE"
  }else{
    RMSE.mort.9sal.max.d.r.sal.all.l[i, "Improved"] <- "TRUE"
  }
}
View(RMSE.mort.9sal.max.d.r.sal.all.l)

# RMSE applied to the Test set
row.names(RMSE.mort.9sal.max.d.r.sal.all.t) <- vars.to.test
colnames(RMSE.mort.9sal.max.d.r.sal.all.t) <- "RMSE"
RMSE.mort.9sal.max.d.r.sal.all.t <- rbind(RMSE.mort.9sal.max.d.r.sal.all.t, "mortality+9sal+max.d.r.sal"=0.8585984)
RMSE.mort.9sal.max.d.r.sal.all.t <- as.data.frame(RMSE.mort.9sal.max.d.r.sal.all.t)
for(i in 1:c(dim(RMSE.mort.9sal.max.d.r.sal.all.t)[1]-1)){
  if(RMSE.mort.9sal.max.d.r.sal.all.t[i,"RMSE"] > RMSE.mort.9sal.max.d.r.sal.all.t[dim(RMSE.mort.9sal.max.d.r.sal.all.t)[1],"RMSE"]){
    RMSE.mort.9sal.max.d.r.sal.all.t[i, "Improved"] <- "FALSE"
  }else{
    RMSE.mort.9sal.max.d.r.sal.all.t[i, "Improved"] <- "TRUE"
  }
}
View(RMSE.mort.9sal.max.d.r.sal.all.t) #did not improved

## Best combination ##
#"log.d9.sal" + "log.max.daily.range.sal" + "log0.mortality.rel.20"

### Save it 
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Learning.set"
saveRDS(object = RMSE.mort.9sal.max.d.r.sal.all.l, file = paste(model.dir, '/', "RMSE.mort.9sal.max.d.r.sal.all.l.RDS", sep=''))
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set"
saveRDS(object = RMSE.mort.9sal.max.d.r.sal.all.t, file = paste(model.dir, '/', "RMSE.mort.9sal.max.d.r.sal.all.t.RDS", sep=''))



## Are the univariate DLM and the best multivariate DLM significantly different from each other? ----

# BEST MULTIVARIATE
# - get the model results
`mortality+sal_2+3+4_non-H` <- readRDS("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Learning.set/mortality+sal_2+3+4_non-H.RDS")
results.test.set <- readRDS("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Several non-hierarchical DLM/Final/Applied to the Test.set/mortality+sal_2+3+4_non-H.RDS")

# - RMSE
et.all <- results.test.set[, paste("et_","log0.mortality.rel.20", sep="")]
et.all <- na.omit(et.all)
RMSE <- sqrt(mean(et.all^2)) #0.8585984


# UNIVARIATE
# - get the model results
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Univariate non-hierarchical DLM"
results.test.set.uni <- readRDS(file = paste(model.dir, '/results_Test.set_uni_DLM_non-H.RDS', sep=''))

# - RMSE
et.all.uni <- results.test.set.uni[, grep("^[et_]", names(results.test.set.uni), value=TRUE)]
et.all.uni <- na.omit(et.all.uni)
RMSE <- sqrt(mean(et.all.uni^2)) #0.8602789


# T-TEST
t.test(et.all^2, et.all.uni^2, paired =TRUE)
#they are significantly different