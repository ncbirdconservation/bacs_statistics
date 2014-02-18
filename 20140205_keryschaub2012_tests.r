####################################################################
#Test code from Kery & Schaub 2012
#

######################chapter 13.2##################################

#set up simulated dataset
nreps <- 1000 								#no of replicates
estimates <- array(NA, dim = c(nreps, 2)) # Array to contain the estimates
R <- 250                                  # No. sites

for (i in 1:nreps) {
   cat(i, "\n")  ;  flush.console()
   x <- runif(R, 0, 10) # choose covariate values
   state<-rbinom(n = R, size = 1, prob = plogis(-3 + 1 * x)) # Occ. state
   obs <- rbinom(n = R, size = 1, prob = 0.6) * state # Observations
   fm <- glm(obs~x, family = binomial)
   estimates[i,] <- fm$coef
   }

par(mfrow = c(3, 1))
hist(estimates[,1], col = "gray", nclass = 50, main = "", xlab = "Intercept estimates", las = 1, ylab = "", freq = FALSE)
abline(v = -3, col = "red", lwd = 3)	# Truth
hist(estimates[,2], col = "gray", nclass = 50, main = "", xlab = "Slope estimates", xlim = c(0,1), las = 1, ylab = "", freq = FALSE)
abline(v = 1, col = "red", lwd = 3)		# Truth

plot(1:10, plogis(estimates[1,1] + estimates[1,2] * (1:10)), col = "gray", lwd = 1, ylab = "Occupancy probability", xlab = "Covariate value", type = "l", ylim = c(0, 1), frame.plot = FALSE, las = 1)
samp <- sample(1:nreps, 1000)
for (i in samp){
   lines(1:10, plogis(estimates[i,1] + estimates[i,2] * (1:10)), col = "gray", lwd = 1, type = "l")
   }
lines(1:10, plogis(-3 + 1 * (1:10)), col = "red", lwd = 3, type = "l")
