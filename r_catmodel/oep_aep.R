# How to build a simple risk model for one location in R, Dag Lohmann, Boot Camp 2009

# create 2000 events with random hazard from a normal distribution
set.seed(123)
hazard <- abs(rnorm(2000,sd=0.3))
# every hazard event has a rate of 1/1000 years, 2 events on average per year
rates  <- rep(2/length(hazard),length(hazard))

# create a parameterized vulnerabilty function, compute the damage ratio as a function of the hazard
damage.ratio <- function(h) {
	      damage.ratio <- h**8
	      damage.ratio <- ifelse(damage.ratio>1,1,damage.ratio)
}

# create exposure (e.g. 1 Million $)
exposure <- 1e6

# compute the loss for each event and average annual loss
losses <- exposure * damage.ratio(hazard)
average.annual.loss <- sum(rates*losses)

# compute oep curve, assume events are Poisson distributed
ordered.losses <- losses[order(losses)]
ordered.rates <- rates[order(losses)]
cum.rates <- cumsum(ordered.rates[length(rates):1])
oep <- 1-exp(-cum.rates)
oep <- oep[length(oep):1]
max.loss <- max(losses)

# plot oep curve on a linear-log plot
plot(ordered.losses,log10(oep),xlab="Loss [$]",ylab="EP [years]",ylim=log10(c(0.001,1)),axes=F,xlim=(c(0,max.loss)),
type="p",lwd=4)
box()
axis(1,at=pretty(c(0,max.loss),5),labels=pretty(c(0,max.loss),5))
axis(2,at=log10(c(0.001,0.002,0.01,0.02,0.1,0.2,1)),labels=c("1000","500","100","50","10","5","1"),las=1)
abline(h=log10(c(0.001,0.002,0.01,0.02,0.1,0.2,1)),col="gray",lty=2)
abline(v=pretty(c(0,max.loss),5),col="gray",lty=2)
text(0.6*max.loss,log10(0.5),paste('Average annual loss =',round(average.annual.loss,1)),cex=1.2)

# compute aep through simulation, average 10 curves for smoother results. The events are assumed to be 
# Poisson distributed. The easiest way to model that is to just draw uniformly distributed random numbers
# for the arrival times. The waiting time between events becomes then exponentially distributed.
years.to.sim <- as.integer(length(hazard)/sum(rates))
sum.aep <- rep(0,years.to.sim)
for (j in c(1:10)) {
event.vector <- rep(0,years.to.sim)
for (i in c(1:length(losses))) {
    noe <- as.integer(rates[i]*years.to.sim)
    u <- runif(noe)
    u2 <- sort(u)
    arrival <- as.integer(u2*years.to.sim+1)
    event.vector[arrival] <- event.vector[arrival] + losses[i]
}
aep1 <- sort(event.vector,decreasing=T)
rates.aep <- c(1:years.to.sim) / years.to.sim
sum.aep <- sum.aep + aep1/10
}

# plot aep curve
lines(sum.aep,log10(rates.aep),type="l",lwd=4,col='red',lty=1)
