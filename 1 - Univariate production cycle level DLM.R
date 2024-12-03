#### Univariate non-hierarchical DLM 

#libraries
library(MASS)

#get data
setwd("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/Datasets")
load("df_salmon_20.RData")
df <- df_salmon_20


### Define the relevant columns in the data set ----

#  -- identifyer is the variables that identifies production cycle or the site
identifyer =  "nseq"

# -- time.var is the variable that describes the observation time
time.var <- "months.since.start"

#  -- expected.start.time is the observation time when we expect the model to start
expected.start.time=0

# -- no.better.limit defines how many iterations the EM algorithm must run after not having improved the first time
no.better.limit=0

# -- relevant.names are the names of the columns to be co-modeled
relevant.names="log0.mortality.rel.20"

# -- best.n.waves are the best number of waves for each variable (NA or 0 means that is not going to be modeled using HW)
best.n.waves <- as.numeric(NA)

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


### Train Univariate non-hierarchical DLM ----
##Parameters per name in relevant.names

# Make sure you avoid scientific notation in your output - it will make things easier!
options(scipen=999)

# Get the functions needed
function.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio"
source(paste(function.dir, '/Functions for Salmon study.R', sep=''))

# Learn all parameters for the DLM - parameters per name in relevant.names
# define.DLM.parameters is a function that defines mu0, C0, V, W, and spline functions for the relevant variables, all from the given data
Start.time <- Sys.time()

out <- define.DLM.parameters(D=Learning.set, identifyer, time.var, time.in.year=12, hierarchical=FALSE,
                             level2=NULL, level3=NULL, best.n.waves, trend=FALSE, expected.start.time=0, no.better.limit=0,
                             relevant.names, metadata.names, plot.it=FALSE, remove.zeros=FALSE, round.by=3, silent=TRUE)

print(Sys.time()-Start.time)

# Save the relevant parameters
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Univariate non-hierarchical DLM"
saveRDS(object = out, file = paste(model.dir, '/Uni_non-H_mu0andC0_not_improved_by_EM.RDS', sep=''))


### How good the Univariate non-hierarchical DLM is? ----
#### Get the DLM parameters
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Univariate non-hierarchical DLM"
DLM_parameters <- readRDS(file = paste(model.dir, '/Uni_non-H_mu0andC0_not_improved_by_EM.RDS', sep=''))

# - Get RMSE (train set)
# - Plot a histogram of the standardized forecast errors
# - Q-Q plot 
#R^2 closer to 1 the better
#(intercept should be=0 and slope=1) - summary(ml)

RMSE.all <- data.frame()
loop=0

par(mfrow=c(1,2))
for(name in relevant.names){ 
  
  table <- DLM_parameters[["results.list.all"]]
  
  #RMSE
  et <- table[,paste("et_", name, sep="")]
  RMSE <- round(sqrt(mean(na.omit(et)^2)),4)
  
  RMSE.all <- rbind(RMSE.all, RMSE)
  rownames(RMSE.all)[loop] <- name
  colnames(RMSE.all)[1] <- "RMSE"
  
  #HIST
  ut <- table[,paste("ut_", name, sep="")]
  
  percent.outside.CI <- round(length(which(ut < -1.96 | ut > 1.96))/length(ut)*100,2)
  hist(ut, 50, xlab="Standardized forecast errors", main=paste(name, '|', percent.outside.CI, '%', '|', RMSE), xlim=c(-4,4), freq=FALSE)
  abline(v=0, col='blue', lwd=2)
  abline(v=-1.96, col='red', lty=2)
  abline(v=1.96, col='red', lty=2)
  curve(dnorm(x, mean=0, sd=1), 
        col="darkblue", lwd=2, add=TRUE, yaxt="n")
  
  #Q-Q plot
  qq <- qqnorm(ut, main = paste("Q-Q Plot Original for ", name, sep=""))
  qqline(ut)
  abline(h=0, col="blue", lty = "dashed")
  abline(v=0, col="blue", lty = "dashed")
  ml <- lm(qq$x~qq$y, data = as.data.frame(qq))
  S <- summary(ml)
  Coef.original <- paste(round(S$coefficients[,1],2), collapse=' | ')
  text(c(-1.5,6,9), max(qq$y, na.rm=T), paste("R^2 =", round(summary(ml)$r.squared,4)))
  text(x = -1.5, y =  (max(qq$y, na.rm=T)-0.2*max(qq$y, na.rm=T)), labels = Coef.original)
}
par(mfrow=c(1,1))



### Apply Univariate non-hierarchical DLM ----
#### Get the DLM parameters
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Univariate non-hierarchical DLM"
DLM_parameters <- readRDS(file = paste(model.dir, '/Uni_non-H_mu0andC0_not_improved_by_EM.RDS', sep=''))

runDLM <- FUNFILTERG
D <- Test.set

# The same standardization factors apply to the test set
Standarized.factors <- DLM_parameters$Standarized.factors
for(name in relevant.names){
  Mean.name <- Standarized.factors[name, "Means"]
  SD.name <- Standarized.factors[name, "SDs"]
  D[,name] <- (D[,name] - Mean.name)/SD.name
}

# Get the relevant DLM parameters
mu0 <- DLM_parameters$mu0.list
C0 <- DLM_parameters$C0.list
V <- DLM_parameters$V.list
W <- DLM_parameters$W.list
Spline.list <- DLM_parameters$Spline.list

# Apply the model for each nseq (production cycle)
start.time <- Sys.time()
results.list.all <- cbind()

for (ID in unique(D[,identifyer])){ 
  
  progress <- round(which(unique(D[,identifyer]) == ID)/length(unique(D[,identifyer]))*100,2)
  print(paste('Applying DLM to new data:', progress, '%'))
  
  # Get the data for the current nseq
  ID.set <- subset(D, D[,identifyer] == ID)
  
  # run the DLM
  res <- runDLM(D = ID.set, mu0 = mu0, C0 = C0, V = V, W = W, time.var, time.in.year=12, 
                best.n.waves, trend=FALSE, relevant.names = relevant.names, Spline.list,
                hierarchical=FALSE, level2=NULL, level3=NULL)
  
  # run the Smoother
  smot = runSmoother(res)
  
  # get all results from res
  results.list <- extract.res(res, smot, hierarchical=FALSE, D = ID.set)
  results.list.all <- rbind(results.list.all, results.list) 
}

# add meta data to the results.list.all
results.test.set <- cbind(D[,metadata.names], results.list.all)

#### save it
results.path <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Univariate non-hierarchical DLM"
saveRDS(object = results.test.set, file = paste(results.path, '/results_Test.set_uni_DLM_non-H.RDS', sep=''))


# - Get RMSE (applied to the Test set!)
#### Get the results
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Univariate non-hierarchical DLM"
results.test.set <- readRDS(file = paste(model.dir, '/results_Test.set_uni_DLM_non-H.RDS', sep=''))

#we can get the RMSE from the DLM on test set
et.all <- results.test.set[, grep("^[et_]", names(results.test.set), value=TRUE)]
et.all <- na.omit(et.all)
RMSE <- sqrt(mean(et.all^2))

## remove the standardization and the log transformation so that the unit of RMSE is the same as the unit of your observation
###first remove the standardization and the log transformation from the ft values
ft <- results.test.set[, grep("^[ft_]", names(results.test.set), value=TRUE)]
ft_no_s <- (ft * SD.name) + Mean.name
ft_no_s_no_l <- exp(1)^ft_no_s - 1/20000
###now subtract it to the observed values of mortality (with no standardization and the log transformation)
et.all <- Test.set$mortality.rel.20 - ft_no_s_no_l
###get RMSE
et.all <- na.omit(et.all)
RMSE <- sqrt(mean(et.all^2))


### Plot mortality of all production cycles to pdf ----

# Get the DLM parameters
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Univariate non-hierarchical DLM"
DLM_parameters <- readRDS(file = paste(model.dir, '/Uni_non-H_mu0andC0_not_improved_by_EM.RDS', sep=''))

# Get the results
model.dir <- "H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Univariate non-hierarchical DLM"
results.test.set <- readRDS(file = paste(model.dir, '/results_Test.set_uni_DLM_non-H.RDS', sep=''))

#get the right data
D <- Test.set
Standarized.factors <- DLM_parameters$Standarized.factors
for(name in relevant.names){
  Mean.name <- Standarized.factors[name, "Means"]
  SD.name <- Standarized.factors[name, "SDs"]
  D[,name] <- (D[,name] - Mean.name)/SD.name
}


#plot to pdf - using base R
setwd("H:/PhD KU/DECIDE/Scottish Salmon data/R studio/R.Data/Univariate non-hierarchical DLM/plots")
pdf("log0.mortality.rel.20_all_+6.graphs.pdf")

#specify to save plots in 3x1 grid
par(mfrow = c(3,1))

for(ID in unique(D$nseq)){ 
  
  ID.set <- subset(D, D[,identifyer] == ID)
  i <- which(results.test.set$nseq == ID)
  mt <- results.test.set[i, paste("mt_level_", relevant.names, sep="")] #only the values of mt (without the trend)
  ft <- results.test.set[i, paste("ft_", relevant.names, sep="")]
  Ct <- results.test.set[i, paste("Ct_", relevant.names, sep="")]       #only the variance of the levels
  Qt <- results.test.set[i, paste("Qt_", relevant.names, sep="")]       #variance - univariate case only 1 value
  mts <- results.test.set[i, paste("mts_", relevant.names, sep="")]     #only the values of mts (without the trend)
  Cts <- results.test.set[i, paste("Cts_", relevant.names, sep="")]     #only the variance of the levels
  
  min <- min(ft-1.96*sqrt(Qt), na.rm=T)
  max <- max(ft+1.96*sqrt(Qt), na.rm=T)
  
  #mt
  plot(ID.set[,relevant.names], type="l", xlab="month", ylab="log.mortality", 
       main=c(ID), ylim=c(min,max), lwd=2, 
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1) #observations
  
  low.limit.mt <- mt-1.96*sqrt(Ct)
  high.limit.mt <- mt+1.96*sqrt(Ct)
  lines(mt, col='red', lwd=2)                        #filtered mean
  lines(low.limit.mt, col='red', lty = "dashed")     #filtered variance
  lines(high.limit.mt, col='red', lty = "dashed")    #filtered variance
  
  legend(x = "bottomright",                          # Position
         legend = c("Obs", "mt", "95% CI"),          # Legend texts
         col = c("black", "red", "red"),             # Line colors
         lwd = 1,                                    # Line thickness
         lty = c(1,1,2),
         cex=0.7)
  
  #ft
  plot(ID.set[,relevant.names], type="l", xlab="month", ylab="log.mortality", 
       main=c(ID), ylim=c(min,max), lwd=2,
       cex.lab=0.8, cex.axis=0.8, cex.main=0.8, cex.sub=0.8) #observations
  
  low.limit.ft <- ft-1.96*sqrt(Qt)
  high.limit.ft <- ft+1.96*sqrt(Qt)
  lines(ft, col='blue', lwd=2)                       #forecasts
  lines(low.limit.ft, col='blue', lty = "dashed")    #forecast variance
  lines(high.limit.ft, col='blue', lty = "dashed")   #forecast variance
  
  legend(x = "bottomright",                          # Position
         legend = c("Obs", "ft", "95% CI"),          # Legend texts
         col = c("black", "blue", "blue"),           # Line colors
         lwd = 1,                                    # Line thickness
         lty = c(1,1,2),
         cex=0.7) 
  
  #mts (smoother)
  plot(ID.set[,relevant.names], type="l", xlab="month", ylab="log.mortality", 
       main=c(ID), ylim=c(min,max), lwd=2,
       cex.lab=1, cex.axis=1, cex.main=1, cex.sub=1) #observations
  
  low.limit.mts <- mts-1.96*sqrt(Cts)
  high.limit.mts <- mts+1.96*sqrt(Cts)
  lines(mts, col='darkgreen', lwd=2)                     #smooth mean
  lines(low.limit.mts, col='darkgreen', lty = "dashed")  #smooth variance
  lines(high.limit.mts, col='darkgreen', lty = "dashed") #smooth variance
  
  legend(x = "bottomright",                              # Position
         legend = c("Obs", "mts", "95% CI"),             # Legend texts
         col = c("black", "darkgreen", "darkgreen"),     # Line colors
         lwd = 1,                                        # Line thickness
         lty = c(1,1,2),
         cex=0.7) 
  
  
  # - remove standardization from mortality
  Standarized.factors <- DLM_parameters$Standarized.factors
  Mean.name <- Standarized.factors[relevant.names, "Means"]
  SD.name <- Standarized.factors[relevant.names, "SDs"]
  
  Obs_no_s <- subset(Test.set[,relevant.names], Test.set[,identifyer] == ID)
  mt_no_s <- (mt * SD.name) + Mean.name
  ft_no_s <- (ft * SD.name) + Mean.name
  mts_no_s <- (mts * SD.name) + Mean.name
  
  ## - For the CI's we have to remove the standardization on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt_no_s <- (low.limit.mt * SD.name) + Mean.name
  high.limit.mt_no_s <- (high.limit.mt * SD.name) + Mean.name
  low.limit.ft_no_s <- (low.limit.ft * SD.name) + Mean.name
  high.limit.ft_no_s <- (high.limit.ft * SD.name) + Mean.name
  low.limit.mts_no_s <- (low.limit.mts * SD.name) + Mean.name
  high.limit.mts_no_s <- (high.limit.mts * SD.name) + Mean.name
  
  # - remove log transformation from mortality (solve in https://www.wolframalpha.com/)
  Obs <- subset(Test.set[,"mortality.rel.20"], Test.set[,identifyer] == ID)
  mt <- exp(1)^mt_no_s - 1/20000
  ft <- exp(1)^ft_no_s - 1/20000
  mts <- exp(1)^mts_no_s - 1/20000
  
  ## - For the CI's we have again to remove log transformation on the all CI calculation, 
  ## - not individually on Ct, Qt and Cts
  low.limit.mt <- exp(1)^low.limit.mt_no_s - 1/20000
  high.limit.mt <- exp(1)^high.limit.mt_no_s - 1/20000
  low.limit.ft <- exp(1)^low.limit.ft_no_s - 1/20000
  high.limit.ft <- exp(1)^high.limit.ft_no_s - 1/20000
  low.limit.mts <- exp(1)^low.limit.mts_no_s - 1/20000
  high.limit.mts <- exp(1)^high.limit.mts_no_s - 1/20000
  
  # - plot with the real values
  min <- min(c(low.limit.ft, Obs), na.rm=T)
  max <- max(c(high.limit.ft, Obs), na.rm=T)
  
  #mt
  plot(Obs, type="l", xlab="Months", ylab="Mortality",
       main=c(ID, unique(ID.set[,"site"])), ylim=c(min,max), lwd=2)
  lines(mt, col='red', lwd=2)                         #filtered mean
  lines(low.limit.mt, col='red', lty = "dashed")      #filtered variance
  lines(high.limit.mt, col='red', lty = "dashed")     #filtered variance
  
  legend(x = "bottomright",                           # Position
         legend = c("Obs", "mt", "95% CI"),           # Legend texts
         col = c("black", "red", "red"),              # Line colors
         lwd = 1,                                     # Line thickness
         lty = c(1,1,2),
         cex=0.7)
  
  #ft
  plot(Obs, type="l", xlab="Months", ylab="Mortality", 
       main=c(ID, unique(ID.set[,"site"])), ylim=c(min,max), lwd=2)
  lines(ft, col='blue', lwd=2)                        #forecasts
  lines(low.limit.ft, col='blue', lty = "dashed")     #forecast variance
  lines(high.limit.ft, col='blue', lty = "dashed")    #forecast variance
  
  legend(x = "bottomright",                           # Position
         legend = c("Obs", "ft", "95% CI"),           # Legend texts
         col = c("black", "blue", "blue"),            # Line colors
         lwd = 1,                                     # Line thickness
         lty = c(1,1,2),
         cex=0.7)
  
  #mts
  plot(Obs, type="l", xlab="Months", ylab="Mortality", 
       main=c(ID, unique(ID.set[,"site"])), ylim=c(min,max), lwd=2)
  
  lines(mts, col='green', lwd=2)                     #smooth mean
  lines(low.limit.mts, col='green', lty = "dashed")  #smooth variance
  lines(high.limit.mts, col='green', lty = "dashed") #smooth variance
  
  legend(x = "bottomright",                          # Position
         legend = c("Obs", "mts", "95% CI"),         # Legend texts
         col = c("black", "green", "green"),         # Line colors
         lwd = 1,                                    # Line thickness
         lty = c(1,1,2),
         cex=0.7) 
}


#turn off PDF plotting
dev.off() 