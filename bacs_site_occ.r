# first row contains variable names, comma is separator 
# assign the variable id to row names
# note the / instead of \ on mswindows systems 

##############################################################################
#load R Libraries
library(beanplot)
# library(MASS) #model comparison - stepAIC
# library(leaps) #another model comparison library

#set working directory
setwd("C:/data/@projects/bacs_surveys/analysis/heat_map/r/")

#import data from exel (Access export)
habq <- read.table("BACSHabData.csv", header=TRUE, sep=",")

#create a subset with surveyed locations only
habq_survey <- habq[ which(habq$NOSURV==0),]
habq_survey$OBS <- factor(habq_survey$OBS)

#############################################################################
#re-factor habitat variables to be 2 factors instead of 3
habq_2factor <- habq_survey
habq_2factor$BA <- factor(ifelse(habq_2factor$BA<=1,0,1))
habq_2factor$FIRE <- factor(ifelse(habq_2factor$FIRE<1,0,1))
habq_2factor$WIRE <- factor(ifelse(habq_2factor$WIRE<=1,0,1))
habq_2factor$BARE <- factor(ifelse(habq_2factor$BARE<1,0,1))
habq_2factor$FORB <- factor(ifelse(habq_2factor$FORB<1,0,1))
habq_2factor$WOOD <- factor(ifelse(habq_2factor$WOOD<1,0,1))
habq_2factor$MIDHARD <- factor(ifelse(habq_2factor$MIDHARD<=1,0,1))

#create factor permutations, export grid to csv
#f2<- c("0","1")
#write.csv(expand.grid(f2, f2, f2, f2, f2, f2, f2), file="HabQuality2Factor.csv")


#############################################################################
#run full model - export coeff to csv
fit_full <- glm(BACSBIN~BA+FIRE+WIRE+FORB+WOOD+BARE+MIDHARD, data=habq_2factor, family=binomial())
summary(fit_full) #display results
write.csv(coefficients(summary(fit_full)), file="habq_2f_full_coeff.csv")


# 
#habitat characteristics where BACS Found
habq_bacs <- subset(habq_survey, TOTBACS>0,select=c(TOTBACS, BACS_BIN, BA, FIRE, WIRE, FORB, WOOD, BARE, MIDHARD, HabQIndex))
habq_nobacs <- subset(habq_survey, TOTBACS==0,select=c(TOTBACS, BACS_BIN, BA, FIRE, WIRE, FORB, WOOD, BARE, MIDHARD, HabQIndex))

pairs(~TOTBACS+HabQIndex+jitter(BA)+jitter(FIRE)+jitter(WIRE)+jitter(BARE), data=habq_survey, main="Habitat Variables & BACS Detections")

hist(habq_bacs$HabQIndex)
hist(habq_nobacs$HabQIndex)

#############################################################################
#############################################################################
#############################################################################
#begin analysis of sites using Site-Occupancy

library("R2WinBUGS") #to run winBUGS. Requires 'coda' package
library(reshape)
library(coda)
bugs.dir = "C:/WinBUGS14/"
bacs <- na.omit(subset(habq_2factor, ,select=c(PASS, ACTIVE, JULIAN, HABQI, BTIME, OBS)))	#only those sites surveyed, and with habitat quality index

y<-as.matrix(bacs[,1:2])		#detection data for passive and active listening periods
y[y>0]<-1						#convert to binary (detection/non-detection)
active<- y						#indicator for passive vs. active 
active[,1]<-0
active[,2]<-1

habqi<-bacs$HABQI				#hab quality index
mean.habqi<-mean(habqi, na.rm=TRUE)
sd.habqi<-sd(habqi)

#qqplot(habqi, rnorm(nrow(habqi),mean = mean(habqi), sd=sd(habqi))) #assess normality of habqi distribution
#beanplot(BACSBIN~HABQI,data=habq_survey)		#beanplot (boxplot alternative) of habitat quality index vs. BACS detection

#TODO - Add in time of day effect including quadratic
#TODO - add in date quadratic effect
#TODO - add in observer effect

dates<-bacs$JULIAN				#julian day of the year

#standardize dates
mean.date<-mean(dates,na.rm=TRUE)
sd.date<-sd(dates[!is.na(dates)])
DATES<-(dates-mean.date)/sd.date	#standardize
DATES[is.na(DATES)]<-0				#impute zeroes(means)

#standardize times
times <- bacs$BTIME
mean.times<-mean(times)
sd.times<-sd(times)
TIMES<-(times-mean.times)/sd.times
TIMES[is.na(TIMES)]<-0

#factorize observers
obs<-factor(bacs$OBS)



# Specify model in BUGS language
# this is a simple model accounting for date in detection probability and habitat quality index for occupancy probability
sink("model.txt")
cat("
model {

# Priors
alpha.psi ~ dnorm(0, 0.01)	#intercept occupancy
beta.psi ~ dnorm(0, 0.01)	#hab quailty index occupancy
alpha.p ~ dnorm(0, 0.01)	#intercept detection
beta1.p ~ dnorm(0, 0.01)	#date effect detection
beta2.p ~ dnorm(0, 0.01)	#active effect detection

# Likelihood
# Ecological model for the partially observed true state
for (i in 1:R) { #sites
   z[i] ~ dbern(psi[i])                # True occurrence z at site i
   psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
   lpsi.lim[i] <- min(999, max(-999, lpsi[i]))
   lpsi[i] <- alpha.psi + beta.psi * habqi[i]

   # Observation model for the observations
   for (j in 1:T) { #detections
      y[i,j] ~ dbern(mu.p[i,j])	# Detection-nondetection at i and j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
      lp[i,j] <- alpha.p + beta1.p * DATES[i] + beta2.p * active[i,j]
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])                             # Number of occupied sites
mean.p <- exp(alpha.p) / (1 + exp(alpha.p))    # Sort of average detection
}
",fill = TRUE)
sink()


# Bundle data
win.data <- list(y = y, R = nrow(y), T = ncol(y), habqi = habqi, DATES = DATES, active=active)

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)	# Good starting values crucial
inits <- function(){list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3), beta.psi=runif(1,-3,3), beta1.p=runif(1,-3,3),beta2.p=runif(1,-3,3))}

# Parameters monitored
params <- c("alpha.psi", "beta.psi", "mean.p", "occ.fs", "alpha.p", "beta1.p", "beta2.p")

# MCMC settings
ni <- 30000
nt <- 10
nb <- 20000
nc <- 3

# Call WinBUGS from R (BRT < 1 min)
out <- bugs(win.data, inits, params, "model.txt", n.chains = nc, 
n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)

write.csv(out$summary, file="bacs_site_occ.csv")


# Posterior distribution of the number of occupied sites in actual sample
hist(out$sims.list$occ.fs, nclass = 30, col = "gray", main = "", xlab = "Number of occupied sites (occ.fs)")
abline(v = 10, lwd = 2) # The observed number
 
 
#this estimates how many visits to virtually guarantee detection of spp (if present) - this needs modification due to active/passive sample sessions
# Pstar <- array(NA, dim = c(out$n.sims, 10))
# x <- cbind(rep(1, 3000), rep(2, 3000), rep(3, 3000), rep(4, 3000), rep(5, 3000), rep(6, 3000), rep(7, 3000), rep(8, 3000), rep(9, 3000), rep(10, 3000)) 
# for (i in 1:out$n.sims) {
   # for (j in 1:10){
      # Pstar[i,j] <- 1 - (1 - out$sims.list$mean.p[i])^j
      # } #j
   # } #i

# boxplot(Pstar ~ x, col = "gray", las = 1, ylab = "Pstar", xlab = "Number of surveys", outline = FALSE)
# abline(h = 0.95, lty = 2, lwd = 2)
 

#Predict effect of Hab quality index with uncertainty - THIS NEEDS MODIFICATION
mcmc.sample <- out$n.sims

original.habqi.pred <- seq(min(habqi), max(habqi), by=1)	#generate a sequence of values that matches possible habqi values
#date.pred <- (original.date.pred - mean.date)/sd.date		#example of normalization above
habqi.pred <- original.habqi.pred										#not normalized above
psi.pred.habqi <- plogis(out$mean$alpha.psi + out$mean$beta.psi * habqi.pred )

array.psi.pred.habqi <- array(NA, dim = c(length(habqi.pred), mcmc.sample))

for (i in 1:mcmc.sample){
   array.psi.pred.habqi[,i] <- plogis(out$sims.list$alpha.psi[i] + out$sims.list$beta.psi[i] * habqi.pred)
   }

# Plot for a subsample of MCMC draws
sub.set <- sort(sample(1:mcmc.sample, size = 200))

par(mfrow = c(2, 1))
plot(original.habqi.pred, psi.pred.habqi, main = "", ylab = "Occupancy probability", xlab = "Habitat Quality Index", ylim = c(0, 1),xlim=c(1,17), type = "l", lwd = 3, frame.plot = FALSE)
for (i in sub.set){
   lines(original.habqi.pred, array.psi.pred.habqi[,i], type = "l", lwd = 1, col = "gray")
   }
lines(original.habqi.pred, psi.pred.habqi, type = "l", lwd = 3, col = "blue")

#active vs passive detection probability
par(mfrow = c(2, 1))
hist(plogis(out$sims.list$alpha.p), nclass = 40, col = "gray", main = "Passive Listening", xlab = "Detection probability", xlim = c(0, 1))
hist(plogis(out$sims.list$alpha.p+ out$sims.list$beta2.p), nclass = 40, col = "gray", main = "Active Playback", xlab = "Detection probability", xlim = c(0, 1))


##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################
#method for summarizing by groups (like pivot table)
#install plyr
#library(plyr) #load library
#habq_survey$QQUAD<- factor(habq_survey$QQUAD) #create new field, convert to factor
#write.csv(ddply(habq_survey, .(QQUAD), summarize, freq=length(QQUAD)),file="C:/data/@projects/bacs_surveys/analysis/heat_map/r/QQUADCounts.csv") #count the num of points per qquad



###################################################################################
###################################################################################
#OLD CODE!

#############################################################################
#export data frame
# write.csv(habq_survey, file="C:/data/@projects/bacs_surveys/analysis/heat_map/HabQualityRExport.csv")

#view header of data
#head(habq)

#scatterplot matrix
#pairs(~TOTBACS+Main_HABQ+NOSURVEY+INSUFHAB+HabQIndex,data=habq,main="Habitat Quality Scatterplot Matrix")

#boxplot comparing habitat quality measures
#reorder hab levels to display in boxplot properly
#habq$Main_HABQ <- factor(habq$Main_HABQ, levels = c("No Habitat","Poor", "Fair", "Good", "Excellent"))
#boxplot(HabQIndex~Main_HABQ,data=habq, main="Hab Quality Index vs. Hab Quality Eval", xlab="Hab Quality Eval", ylab="HabQIndex")

#add column with binary BACS presence
#habq["BACS_BIN"] <- ifelse(habq$TOTBACS > 0, 1,0)

#alternate method - pulls out subset, and returns only specified columns
#habq_survey <- subset(habq, NOSURVEY==0 | INSUFHAB==0, select=c(POINTID,QQUAD,TOTBACS, HabQIndex))

#boxplot
#boxplot(HabQIndex~BACS_BIN, data=habq_survey, main="BACS Presence by Habitat Quality Index", xLab="BACS Presence", yLab="Habitat Quality Index")

#logistic regression
#fit <- glm(BACS_BIN~HabQIndex, data=habq_survey, family=binomial())
#summary(fit) #display results
#export coeff to csv
#write.csv(coefficients(summary(fit)), file="C:/data/@projects/bacs_surveys/analysis/heat_map/r/habq_index_coeff.csv")

#boxplot(HabQIndex~BACS_BIN, data=habq_survey, main="Habitat Quality Index at BACS Sites", xlab="BACS Presence", ylab="Habitat Quality Index")


#confint(fit) #95% CI for coefficients
#predict(fit, type="response") #predicted values
#residuals(fit, type="deviance") #residuals

#poission regression on TOTBACS
#fitp <- glm(TOTBACS~HabQIndex, data=habq_survey, family=poisson())
#summary(fitp)

#fit_win <- glm(BACS_BIN~FBA+FFIRE+FWIRE+FBARE, data=habq_factors, family=binomial())

#pfit_full <- glm(TOTBACS~BA+FIRE+WIRE+FORB+WOOD+BARE+MIDHARD, data=habq_2factor, family=poisson())
#summary(pfit_full) #display results


#controversial model comparison - stepwise variable selection based on AIC values
#library(MASS)
#step <- stepAIC(fit_full, direction="both")

#AIC values for given models
#AIC(fit_full, fit_fire, fit_wire, fit_fire_wire_ba, fit_fire_wire_bare, fit_fire_wire_bare_midhard)


#controversial model comparison - stepwise variable selection based on AIC values
#step <- stepAIC(fit_full, direction="both")
#fit_win <- glm(BACS_BIN~BA+FIRE+WIRE+BARE, data=habq_survey, family=binomial())


#alternative method comparing all possible regressions
# #explanation: http://www.stat.columbia.edu/~martin/W2024/R10.pdf
# leaps = regsubsets(BACS_BIN~BA+FIRE+WIRE+FORB+WOOD+BARE+MIDHARD, data=habq_survey, nbest=15)
# plot(leaps, scale="adjr2")
# plot(leaps, scale="aic")


# #alternative 2: stepwise comparing r-squared values
# null = glm(BACS_BIN~1, data=habq_survey, family=binomial())
# step(null, scope=list(lower=null, upper=fit_full), direction="forward") #forward
# step(fit_full, data=habq_factors,direction="backward") #backward 
