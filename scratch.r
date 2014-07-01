dnorm(0) == 1/sqrt(2*pi)
dnorm(1) == exp(-1/2)/sqrt(2*pi)
dnorm(1) == 1/sqrt(2*pi*exp(1))
## Using "log = TRUE" for an extended range :
plot(function(x) dnorm(0,0.01), -60, 50, main = "Normal density")
curve(dnorm(0,0.01), add = TRUE, col = "red", lwd = 2)



dist.norm.test<-dnorm(1000,0,0.01)
dist.unif.test<-runif(1000,0,0.01)

hist(dist.norm.test)