# first row contains variable names, comma is separator 
# assign the variable id to row names
# note the / instead of \ on mswindows systems 

##############################################################################
#load R Libraries
library(beanplot)
# library(stringr) #string manipulation
# library(MASS) #model comparison - stepAIC
# library(leaps) #another model comparison library

#set working directory
setwd("C:/data/@projects/bacs_surveys/analysis/bacs_statistics/")

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
#only have to do this once...
#f2<- c("0","1")
#write.csv(expand.grid(f2, f2, f2, f2, f2, f2, f2), file="HabQuality2Factor.csv")


#############################################################################
#GLM ANALYSIS
#full model - export coeff to csv
#fit_full <- glm(BACSBIN~BA+FIRE+WIRE+FORB+WOOD+BARE+MIDHARD, data=habq_2factor, family=binomial())
#summary(fit_full) #display results
#write.csv(coefficients(summary(fit_full)), file="habq_2f_full_coeff.csv")


# 
#habitat characteristics where BACS Found
# habq_bacs <- subset(habq_survey, TOTBACS>0,select=c(TOTBACS, BACS_BIN, BA, FIRE, WIRE, FORB, WOOD, BARE, MIDHARD, HabQIndex))
# habq_nobacs <- subset(habq_survey, TOTBACS==0,select=c(TOTBACS, BACS_BIN, BA, FIRE, WIRE, FORB, WOOD, BARE, MIDHARD, HabQIndex))

# pairs(~TOTBACS+HabQIndex+jitter(BA)+jitter(FIRE)+jitter(WIRE)+jitter(BARE), data=habq_survey, main="Habitat Variables & BACS Detections")

# hist(habq_bacs$HabQIndex)
# hist(habq_nobacs$HabQIndex)

#############################################################################
#############################################################################
#############################################################################
#begin analysis of sites using Site-Occupancy

library("R2WinBUGS") #to run winBUGS. Requires 'coda' package
library(reshape)
library(coda)
bugs.dir = "C:/WinBUGS14/"
bacs <- na.omit(subset(habq_2factor, ,select=c(PASS, ACTIVE, JULIAN, BA, FIRE, WIRE, FORB, WOOD, BARE, MIDHARD, HABQI, HABQ2, TMIDMIN, OBS,POINTID)))	#only those sites surveyed, and with habitat quality index

#function for converting continuous to scaled variable
#input: vector output:vector of converted values
stand <- function(x) {
	mean <- mean(x, na.rm=TRUE)
	#mean <- mean(x)
	sd <- sd(x[!is.na(x)])
	#sd <- sd(x)
	vals <- (x-mean)/sd
	result <- list(vals=vals, mean=mean, sd=sd)
	return(result)
}

# destand <- function(x, m, s) {
	# #reverse standardization back to original values
	# vals<- (x-m)/s
	# return(vals)
# }

y<-as.matrix(bacs[,1:2])		#detection data for passive and active listening periods
y[y>0]<-1						#convert to binary (detection/non-detection)
active<- y						#indicator for passive vs. active 
active[,1]<-0
active[,2]<-1

# habqi<-bacs$HABQI				#hab quality index
# mean.habqi<-mean(habqi, na.rm=TRUE)
# sd.habqi<-sd(habqi)

#qqplot(habqi, rnorm(nrow(habqi),mean = mean(habqi), sd=sd(habqi))) #assess normality of habqi distribution
#beanplot(BACSBIN~HABQI,data=habq_survey)		#beanplot (boxplot alternative) of habitat quality index vs. BACS detection

#TODO - Add in time of day effect including quadratic
#TODO - add in observer effect

#standardize dates
dates<-stand(bacs$JULIAN)				#julian day of the year
DATES <- dates$vals

#standardize times
times<- stand(bacs$TMIDMIN)
TIMES <- times$vals

# #factorize observers

# obs.jpc<-match(bacs$OBS,"JPC",nomatch=0)
# obs.ajl<-match(bacs$OBS,"AJL",nomatch=0)
# obs.bem<-match(bacs$OBS,"BEM",nomatch=0)
# obs.bs<-match(bacs$OBS,"BS",nomatch=0)
# obs.edg<-match(bacs$OBS,"EDG",nomatch=0)
# obs.gc<-match(bacs$OBS,"GC",nomatch=0)
# obs.jfm<-match(bacs$OBS,"JFM",nomatch=0)
# obs.kkj<-match(bacs$OBS,"KKJ",nomatch=0)
# obs.man<-match(bacs$OBS,"MAN",nomatch=0)
# obs.nh<-match(bacs$OBS,"NH",nomatch=0)
# obs.ska<-match(bacs$OBS,"SKA",nomatch=0)

# #ska never detected bird, combine with jfm
# obs.jfm<-obs.jfm + obs.ska


# habitat quality index
habqi<-bacs$HABQI

#Habitat variables
BA<-as.integer(bacs$BA)
FIRE<-as.integer(bacs$FIRE)
WIRE<-as.integer(bacs$WIRE)
WOOD<-as.integer(bacs$WOOD)
FORB<-as.integer(bacs$FORB)
BARE<-as.integer(bacs$BARE)
MIDHARD<-as.integer(bacs$MIDHARD)

# Specify model in BUGS language
# this is a simple model accounting for date and time in detection probability and habitat quality parameters for occupancy probability
sink("model.txt")
cat("
model {

# Priors
alpha.psi ~ dnorm(0, 0.01)	#intercept occupancy
beta1.psi ~ dnorm(0, 0.01)	#BA index occupancy
beta2.psi ~ dnorm(0, 0.01)	#FIRE index occupancy
beta3.psi ~ dnorm(0, 0.01)	#WIRE index occupancy
beta4.psi ~ dnorm(0, 0.01)	#WOOD index occupancy
beta5.psi ~ dnorm(0, 0.01)	#FORB index occupancy
beta6.psi ~ dnorm(0, 0.01)	#BARE index occupancy
beta7.psi ~ dnorm(0, 0.01)	#MIDHARD index occupancy
# betahabqi.psi ~ dnorm(0, 0.01)	#BARE index occupancy
alpha.p ~ dnorm(0, 0.01)	#intercept detection
beta1.p ~ dnorm(0, 0.01)	#active effect detection
beta2.p ~ dnorm(0, 0.01)	#date effect detection
# beta6.p ~ dnorm(0, 0.01)	#date^2 effect detection
beta3.p ~ dnorm(0, 0.01)	#time effect detection
beta5.p ~ dnorm(0, 0.01)	#time^2 effect detection
# beta4.p ~ dnorm(0, 0.01)	#obs effect detection

#individual observer betas (factors)
# beta.ajl.p ~ dnorm(0, 0.01) #obs effect factor
# beta.bem.p ~ dnorm(0, 0.01) #obs effect factor
# beta.bs.p ~ dnorm(0, 0.01) #obs effect factor
# beta.edg.p ~ dnorm(0, 0.01) #obs effect factor
# beta.gc.p ~ dnorm(0, 0.01) #obs effect factor
# beta.jfm.p ~ dnorm(0, 0.01) #obs effect factor
# beta.jpc.p ~ dnorm(0, 0.01) #obs effect factor
# beta.kkj.p ~ dnorm(0, 0.01) #obs effect factor
# beta.man.p ~ dnorm(0, 0.01) #obs effect factor
# beta.nh.p ~ dnorm(0, 0.01) #obs effect factor
# beta.ska.p ~ dnorm(0, 0.01) #obs effect factor


# Likelihood
# Ecological model for the partially observed true state
for (i in 1:R) { #sites
   z[i] ~ dbern(psi[i])                # True occurrence z at site i
   psi[i] <- 1 / (1 + exp(-lpsi.lim[i]))
   lpsi.lim[i] <- min(999, max(-999, lpsi[i]))
   lpsi[i] <- alpha.psi + beta1.psi * BA[i] + beta2.psi * FIRE[i] + beta3.psi * WIRE[i] + beta4.psi * WOOD[i] + beta5.psi * FORB[i] + beta6.psi * BARE[i] + beta7.psi * MIDHARD[i]
   # lpsi[i] <- alpha.psi + beta1.psi * BA[i] + beta2.psi * FIRE[i] + beta3.psi * WIRE[i] + beta4.psi * WOOD[i] + beta5.psi * FORB[i] + beta6.psi * BARE[i] + beta7.psi * MIDHARD[i] + beta.ajl.p * obs.ajl[i] + beta.bem.p * obs.bem[i] + beta.bs.p * obs.bs[i] + beta.edg.p * obs.edg[i] + beta.gc.p * obs.gc[i] + beta.jfm.p * obs.jfm[i] + beta.jpc.p * obs.jpc[i] + beta.kkj.p * obs.kkj[i] + beta.man.p * obs.man[i] + beta.nh.p * obs.nh[i]

   # Observation model for the observations
   for (j in 1:T) { #detections
      y[i,j] ~ dbern(mu.p[i,j])	# Detection-nondetection at i and j
      mu.p[i,j] <- z[i] * p[i,j]
      p[i,j] <- 1 / (1 + exp(-lp.lim[i,j]))
      lp.lim[i,j] <- min(999, max(-999, lp[i,j]))
      lp[i,j] <- alpha.p + beta1.p * active[i,j] + beta2.p * DATES[i] + beta3.p * TIMES[i] + beta5.p * pow(TIMES[i], 2)
      # lp[i,j] <- alpha.p + beta1.p * active[i,j] + beta2.p * DATES[i] + beta3.p * TIMES[i] + beta5.p * pow(TIMES[i], 2) + beta.ajl.p * obs.ajl[i] + beta.bem.p * obs.bem[i] + beta.bs.p * obs.bs[i] + beta.edg.p * obs.edg[i] + beta.gc.p * obs.gc[i] + beta.jfm.p * obs.jfm[i] + beta.jpc.p * obs.jpc[i] + beta.kkj.p * obs.kkj[i] + beta.man.p * obs.man[i] + beta.nh.p * obs.nh[i] + beta.ska.p * obs.ska[i]
      # lp[i,j] <- alpha.p + beta1.p * active[i,j] + beta2.p * DATES[i] + beta3.p * TIMES[i] + beta5.p * pow(TIMES[i], 2)
      } #j
   } #i

# Derived quantities
occ.fs <- sum(z[])                             # Number of occupied sites
mean.p <- exp(alpha.p) / (1 + exp(alpha.p))    # Sort of average detection
}
",fill = TRUE)
sink()


# Bundle data
win.data <- list(y = y, R = nrow(y), T = ncol(y), BA=BA, FIRE=FIRE, WIRE=WIRE, WOOD=WOOD, FORB=FORB, BARE=BARE, MIDHARD=MIDHARD, TIMES=TIMES, DATES=DATES, active=active)
#all obs factors
# win.data <- list(y = y, R = nrow(y), T = ncol(y), BA=BA, FIRE=FIRE, WIRE=WIRE, WOOD=WOOD, FORB=FORB, BARE=BARE, MIDHARD=MIDHARD, TIMES=TIMES, DATES=DATES, active=active, obs.ajl=obs.ajl, obs.bem=obs.bem, obs.bs=obs.bs, obs.edg=obs.edg, obs.gc=obs.gc, obs.jfm=obs.jfm, obs.jpc=obs.jpc, obs.kkj=obs.kkj, obs.man=obs.man, obs.nh=obs.nh, obs.ska=obs.ska)
# win.data <- list(y = y, R = nrow(y), T = ncol(y), BA=BA, FIRE=FIRE, WIRE=WIRE, WOOD=WOOD, FORB=FORB, BARE=BARE, MIDHARD=MIDHARD, TIMES=TIMES, DATES=DATES, active=active, obs.ajl=obs.ajl, obs.bem=obs.bem, obs.bs=obs.bs, obs.edg=obs.edg, obs.gc=obs.gc, obs.jfm=obs.jfm, obs.jpc=obs.jpc, obs.kkj=obs.kkj, obs.man=obs.man, obs.nh=obs.nh)

# Initial values
zst <- apply(y, 1, max, na.rm = TRUE)	# Good starting values crucial
inits <- function(){
	list(z = zst, alpha.psi=runif(1, -3, 3), alpha.p = runif(1, -3, 3), beta1.p=runif(1,-3,3), beta2.p=runif(1,-3,3), beta3.p=runif(1,-3,3), beta5.p=runif(1,-3,3))
	}

# Parameters monitored

params <- c("alpha.psi", "beta1.psi", "beta2.psi", "beta3.psi", "beta4.psi", "beta5.psi", "beta6.psi", "beta7.psi", "mean.p", "occ.fs", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta5.p")

#all obs factorized
# params <- c("alpha.psi", "beta1.psi", "beta2.psi", "beta3.psi", "beta4.psi", "beta5.psi", "beta6.psi", "beta7.psi", "mean.p", "occ.fs", "alpha.p", "beta1.p", "beta2.p", "beta3.p", "beta5.p", "beta.ajl.p","beta.bem.p","beta.bs.p","beta.edg.p","beta.gc.p","beta.jfm.p","beta.jpc.p","beta.kkj.p","beta.man.p","beta.nh.p")

# MCMC settings
ni <- 50000
nt <- 10
nb <- 30000
nc <- 3

# Call WinBUGS from R (BRT < 1 min)
out <- bugs(win.data, inits, params, "model.txt", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = TRUE, bugs.directory = bugs.dir, working.directory = getwd())

# Summarize posteriors
print(out, dig = 2)
write.csv(out$summary, file="bacs_site_occ_habs.csv")


# Posterior distribution of the number of occupied sites in actual sample
par(mfrow=c(1,1))
hist(out$sims.list$occ.fs[2500:3000], nclass = 30, col = "gray", main = "", xlab = "Number of occupied sites (occ.fs)")
abline(v = mean(out$sims.list$occ.fs[2500:3000]), lwd = 2) # The observed number
 
 
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
 
#OCCUPANCY
##########################################################
#Predict effect of Hab quality index with uncertainty

##########################################################
#graph habitat variables - beanplot?
#occupancy probability
ba.pred <-plogis(out$sims.list$alpha.psi + out$sims.list$beta1.psi)
fire.pred <-plogis(out$sims.list$alpha.psi + out$sims.list$beta2.psi)
wire.pred <-plogis(out$sims.list$alpha.psi + out$sims.list$beta3.psi)
wood.pred <-plogis(out$sims.list$alpha.psi + out$sims.list$beta4.psi)
forb.pred <-plogis(out$sims.list$alpha.psi + out$sims.list$beta5.psi)
bare.pred <-plogis(out$sims.list$alpha.psi + out$sims.list$beta6.psi)
midhard.pred <-plogis(out$sims.list$alpha.psi + out$sims.list$beta7.psi)
ba.fire.wire.pred <-plogis(out$sims.list$alpha.psi + out$sims.list$beta1.psi + out$sims.list$beta2.psi + out$sims.list$beta3.psi)
names<=c("Stand Density","Fire", "Wood", "Forb", "Bare Ground", "Mid Hardwood")

beanplot( ba.pred, fire.pred, wire.pred, wood.pred, forb.pred, bare.pred, midhard.pred, main = "Predicted Occupancy Probability", names=c("Stand Density","Fire", "Wood", "Forb", "Bare Ground", "Mid Hardwood"))
beanplot( ba.fire.wire.pred, ba.pred, fire.pred, wire.pred, wood.pred, forb.pred, bare.pred, midhard.pred, main = "Predicted Occupancy Probability", method="jitter")
beanplot( ba.fire.wire.pred, main = "Predicted Occupancy Probability", method="jitter")
# boxplot( ba.pred, fire.pred, wire.pred, wood.pred, main = "Predicted Occupancy Probability", names=c("Stand Density","Fire", "Wire", "Wood"))
# boxplot( ba.pred, fire.pred, wire.pred, wood.pred, main = "Predicted Occupancy Probability")


#DETECTABILITY
##########################################################
#graph active vs passive detection probability
par(mfrow = c(2, 1))
#plogis is the "inverse logit" function
hist(plogis(out$sims.list$alpha.p), nclass = 40, col = "gray", main = "Passive Listening", xlab = "Detection probability", xlim = c(0, 1))
abline(v=mean(plogis(out$sims.list$alpha.p)), col="blue")
hist(plogis(out$sims.list$alpha.p+ out$sims.list$beta1.p), nclass = 40, col = "gray", main = "Active Playback", xlab = "Detection probability", xlim = c(0, 1))
abline(v=mean(plogis(out$sims.list$alpha.p+ out$sims.list$beta1.p)), col="blue")



##########################################################
# Predict effect of time of day  and date with uncertainty
mcmc.sample <- out$n.sims

original.date.pred <- seq(min(bacs$JULIAN), max(bacs$JULIAN), length.out = 30)
original.time.pred <- seq(min(bacs$TMIDMIN), max(bacs$TMIDMIN), length.out = 30)
date.pred <- (original.date.pred - dates$mean)/dates$sd
time.pred <- (original.time.pred - times$mean)/times$sd
p.pred.date <- plogis(out$mean$alpha.p + out$mean$beta2.p * date.pred)
p.pred.time <- plogis(out$mean$alpha.p + out$mean$beta3.p * time.pred + out$mean$beta5.p * time.pred^2 )

array.p.pred.time <- array.p.pred.date <- array(NA, dim = c(length(time.pred), mcmc.sample))
for (i in 1:mcmc.sample){
   array.p.pred.date[,i] <- plogis(out$sims.list$alpha.p[i] + out$sims.list$beta2.p[i] * date.pred)
   array.p.pred.time[,i] <- plogis(out$sims.list$alpha.p[i] + out$sims.list$beta3.p[i] * time.pred + out$sims.list$beta5.p[i] * time.pred^2)
   }

# Plot for a subsample of MCMC draws
sub.set <- sort(sample(1:mcmc.sample, size = 200))

par(mfrow = c(2, 1))
plot(original.date.pred, p.pred.date, main = "", ylab = "Detection probability", xlab = "Julian Day", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)
for (i in sub.set){
   lines(original.date.pred, array.p.pred.date[,i], type = "l", lwd = 1, col = "gray")
   }
lines(original.date.pred, p.pred.date, type = "l", lwd = 3, col = "blue")

plot(original.time.pred, p.pred.time, main = "", ylab = "Detection probability", xlab = "Time of day (mins after midnight)", ylim = c(0, 1), type = "l", lwd = 3, frame.plot = FALSE)
for (i in sub.set){
   lines(original.time.pred, array.p.pred.time[,i], type = "l", lwd = 1, col = "gray")
   }
lines(original.time.pred, p.pred.time, type = "l", lwd = 3, col = "blue")

