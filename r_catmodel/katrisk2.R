# Dag Lohmann, KatRisk LLC, June 2012
# katrisk.R explains the basics of risk modeling
# The code assumes some abitrary hazard and vulnerability functions that
# would normally be peril dependent (e.g. as functions of peak wind gust
# or flood depth)
# The code contains all features of a risk model including secondary 
# uncertainty

# Here are the steps to run the code
# 1) install R
# 2) install.packages("copula","evd","fields","maps")
# type source("katrisk.R") in the R window - or use RStudio to go through 
# the code

# This is what the code does:
# 1) create hazard events
# 2) create exposure
# 3) create vulnerability
# 4) build an event loss table
# 5) define and run financial model
# 6) compute OEP and AEP curves for expected losses
# 7) Include secondary uncertainty as beta distributed around the mean loss
# 8) build an event loss table for sampled losses
# 9) run financial model on sampled losses
# 10) compute OEP and AEP curves for sampled losses
# 11) compute other analytical quantities
#     - loss cost
#     - probable maximum loss at 250 year RP
#     - contribution of return period losses to AAL
#     - EP uncertainty

# load libraries, set seed for repeatability, set plot parameter
library("copula")
library("evd")
library('fields')
library('maps')
set.seed(1)
par(mfrow=c(2,2),mgp=c(2.3,0.7,0),oma=c(0,0,2,0))

# 1) create hazard
#    This is normally done with detailed vendor specific peril models. Model
#    output can be stored by either geographical location (e.g. lat/lon for
#    1 degree grid cells), or by administrative boundaries such as postcodes.
#    The important role of geocoding is then to align the exposure to each of
#    these.

# number of grids in the hazard field (this would be equivalent to producing
# e.g. quadkeys for n grids and giving them keys
# n.grids is then just the number of different possible geolocations
# n.events is the number of events, each with a certain rate. These could be
# non-uniform of course
# rho.par is the correlation between locations. This would normally be calculated
# by a hazard model
# produce hazard on grid that looks as if model was built for the UK - but don't care whether losses are over land or ocean

x.coord <- c(-6:2)
y.coord <- c(50:59)
dim.hazard <- c(length(x.coord),length(y.coord))
n.grids  <- dim.hazard[1]*dim.hazard[2]

n.events <- 10000
rate.per.event <- 1/10000
event.rates <- rep(rate.per.event,n.events)
rho.par <- 0.5

# create n.events events of some random hazard with an elliptic copula and
# GPD marginal distribtions
myCop.norm <- ellipCopula(family="normal",dim=n.grids,dispstr="ex",param=rho.par)
gpd.par <- list(loc=1,scale=8,shape=0.5)
myMvd <- mvdc(copula=myCop.norm,margins=rep("gpd",n.grids),paramMargins=rep(list(gpd.par),n.grids))
#hazard <- rmvdc(myMvd,n.events)
hazard <- rMvdc(n.events,myMvd)

# Plot first 4 hazard fields
for (i in c(1:4)) {
  hazard.plot <- hazard[i,]
  dim(hazard.plot) <- dim.hazard
  hazard.plot <- ifelse(hazard.plot>100,100,hazard.plot)
  image.plot(x.coord,y.coord,hazard.plot,zlim=c(0,100),xlab="Longitude",ylab="Latitude")
  map(add=T,col='white',lwd=2)
}
title(main="First 4 Hazard Fields",cex.main=1.7,outer=T)

# 2) create exposure. I assume here that n.loc locations are randomly
#    distributed in the n.grids hazard data set. Normally geocoding would
#    take care of this and assign each location to a hazard grid

# number of locations in the exposure set
n.loc <- 20
# create exposure, e.g. $100k per subject at risk
exposure <- rep(100000,n.loc)
# give each of the locations a random location in the n.grids event grids
location.in.hazard <- as.integer(runif(n.loc)*n.grids)+1

# 3) Create vulnerability
# a parameterized function for a hazard between 0 and 100
damage.ratio <- function(hazard) {
  hazard <- ifelse(hazard<20,0,hazard)
  hazard <- ifelse(hazard>99.99,99.99,hazard)
  damage.ratio <- (hazard/100)**6
}

# 4) Create the event loss table
# this table contains by event the ground up losses for each location
# create top 10 events
gu.dr.loss.table <- damage.ratio(hazard[,location.in.hazard])
gu.event.loss.table <- exposure*gu.dr.loss.table
gu.loss.by.event <- apply(gu.event.loss.table,1,sum)
top10 <- order(gu.loss.by.event,decreasing=T)[1:10]
# plot top 4 events
for (i in c(1:4)) {
  hazard.plot <- hazard[top10[i],]
  dim(hazard.plot) <- dim.hazard
  hazard.plot <- ifelse(hazard.plot>100,100,hazard.plot)
  image.plot(x.coord,y.coord,hazard.plot,zlim=c(0,100),xlab="Longitude",ylab="Latitude")
  map(add=T,col='white',lwd=2)
}
title(main="Top 4 Events",cex.main=1.7,outer=T)

# 5) Define financial model
#    for simplicity same financial structure with deductible and limit
#    for all exposures
#    Compute the gross loss event loss table
l.deductible <- 1000
l.limit <- 80000
gr.event.loss.table<-ifelse(gu.event.loss.table<l.deductible,0,gu.event.loss.table-l.deductible)
gr.event.loss.table<-ifelse(gr.event.loss.table>l.limit,l.limit,gr.event.loss.table)
gr.loss.by.event <- apply(gr.event.loss.table,1,sum)

# compute average annual loss (AAL) for each location
gu.aal <- apply(gu.event.loss.table*event.rates,2,sum)
gr.aal <- apply(gr.event.loss.table*event.rates,2,sum)
# plot AAL for GU and GR loss
gu.aal.plot <- rep(0,n.grids)
gr.aal.plot <- rep(0,n.grids)
for (i in c(1:n.loc)) {
  gu.aal.plot[location.in.hazard[i]] <- gu.aal.plot[location.in.hazard[i]] + gu.aal[i]
  gr.aal.plot[location.in.hazard[i]] <- gr.aal.plot[location.in.hazard[i]] + gr.aal[i]
}

dim(gu.aal.plot) <- dim.hazard
dim(gr.aal.plot) <- dim.hazard
image.plot(x.coord,y.coord,gu.aal.plot,xlab="Longitude",ylab="Latitude")
map(add=T,col='white',lwd=2)
title(main="GU AAL",cex.main=1)
image.plot(x.coord,y.coord,gr.aal.plot,xlab="Longitude",ylab="Latitude")
map(add=T,col='white',lwd=2)
title(main="GR AAL",cex.main=1)

# 6) Compute OEP and AEP curve through simulation of n.years
#    this is also possible for some analytical distributions of frequency
#    distributions
#    OEP = occurrence exceedance probability, the distribution of the maximum
#          loss per year
#    AEP = aggregate exceedance probability, the distribution of annual
#          (normally, but could use other aggregation periods) losses

# assume a Poisson distribution for the number of occurrences
# model this with arrival times that years*U[0,1] distributed
# run each event at a minimum n.sim times 
# only simulate with events that produce a loss
# round the number of expected loss occurrences

n.sim <- 10
loss.events <- which(gu.loss.by.event>0)
n.years <- n.sim/min(event.rates)
gu.oep.losses <- rep(0,n.years)
gu.aep.losses <- rep(0,n.years)
gr.oep.losses <- rep(0,n.years)
gr.aep.losses <- rep(0,n.years)

for (i in loss.events) {
  n.occur <- round(n.years*event.rates[i])
  year.occur <- as.integer(runif(n.occur)*n.years)+1
  for (j in year.occur) {
	gu.oep.losses[j] <- max(gu.oep.losses[j],gu.loss.by.event[i])
	gr.oep.losses[j] <- max(gr.oep.losses[j],gr.loss.by.event[i])
        gu.aep.losses[j] <- gu.aep.losses[j]+gu.loss.by.event[i]
        gr.aep.losses[j] <- gr.aep.losses[j]+gr.loss.by.event[i]
	}
}

# order the losses and rates
ord.gu.oep.losses<- gu.oep.losses[order(gu.oep.losses,decreasing=T)]
ord.gr.oep.losses<- gr.oep.losses[order(gr.oep.losses,decreasing=T)]
ord.gu.aep.losses<- gu.aep.losses[order(gu.aep.losses,decreasing=T)]
ord.gr.aep.losses<- gr.aep.losses[order(gr.aep.losses,decreasing=T)]

rates <- c(1:n.years) / n.years

# plot the results
max.plot <- max(ord.gu.aep.losses)
y.lim <- log10(c(1/n.years,1))
x.lim <- c(0,max.plot)
plot(ord.gu.oep.losses,log10(rates),ylim=y.lim,xlim=x.lim,xlab="Loss",ylab="Return Period",las=1,type="l",lwd=2,axes=F)
lines(ord.gu.aep.losses,log10(rates),col="red",lwd=2)
lines(ord.gr.oep.losses,log10(rates),col="blue",lwd=2)
lines(ord.gr.aep.losses,log10(rates),col="green",lwd=2)

p.intervals = log10(n.years)+1
ep.plot <- 10**c(p.intervals:0)
axis(1,at=pretty(c(0,max.plot),p.intervals),labels=pretty(c(0,max.plot),p.intervals))
axis(2,at=log10(1/ep.plot),labels=as.character(ep.plot),las=1)
abline(h=log10(1/ep.plot),col="lightgray",lwd=2,lty=2)
abline(v=pretty(c(0,max.plot),p.intervals),lwd=2,col="lightgray",lty=2)
title(main="AAL and EP Results",cex.main=1.7,outer=T)


# 7) Include secondary uncertainty as beta distribution around the mean loss for each loaction with a correlation
# that is not dependent on location or distance. That can be changed in the copula.
# define functions for beta distribution, a and b parameter, and mdr-cv relationship
a.func <- function(mdr,cv)(1-mdr)/cv**2-mdr
b.func <- function(mdr,a) a*(1-mdr)/mdr
mdrcv <- function(mdr,p=1)(sqrt((1-mdr)/(p+mdr)))

# define n.sample as the number of samples taken
# easy to just use n.sim
n.sample <- n.sim
# define the correlation of uncertainty
rho.uncertainty <- 0.1
# calculate the number of sampled events and provision results matrix and new event rates
n.sample.events <- length(loss.events)*n.sample
sampled.damage.ratio <- matrix(0,n.sample.events,n.loc)
s.event.rates <- rep(event.rates[loss.events]/n.sample,each=n.sample)
counter <- 0
for (i in loss.events) {
	loss.locs <- which(gu.dr.loss.table[i,] > 0)
	n.l <- length(loss.locs)
	beta.cv <- mdrcv(gu.dr.loss.table[i,loss.locs])
	beta.a <- a.func(gu.dr.loss.table[i,loss.locs],beta.cv)
	beta.b <- b.func(gu.dr.loss.table[i,loss.locs],beta.a)
	if (n.l>1) {
		myUncertCop.norm <- ellipCopula(family="normal",dim=n.l,dispstr="ex",param=rho.uncertainty)
		beta.par <- list()
		for (k in c(1:n.l)) beta.par[[k]] <- list(shape1=beta.a[k],shape2=beta.b[k])
		myUncertMvd <- mvdc(copula=myUncertCop.norm,margins=rep("beta",n.l),paramMargins=beta.par)
		#mdr.sample <- rmvdc(myUncertMvd,n.sample)
                mdr.sample <- rMvdc(n.sample,myUncertMvd)
	}
	else {
		mdr.sample <- rbeta(n=n.sample,shape1=beta.a[1],shape2=beta.b[1])
	}
	sampled.damage.ratio[(n.sample*counter)+c(1:n.sample),loss.locs] <- mdr.sample			
	counter <- counter + 1
	}

# 8) Create the event loss table for sampled events
# this table contains by event the ground up losses for each location
s.gu.event.loss.table <- exposure*sampled.damage.ratio
s.gu.loss.by.event <- apply(s.gu.event.loss.table,1,sum)

# 9) Use same financial model
s.gr.event.loss.table<-ifelse(s.gu.event.loss.table<l.deductible,0,s.gu.event.loss.table-l.deductible)
s.gr.event.loss.table<-ifelse(s.gr.event.loss.table>l.limit,l.limit,s.gr.event.loss.table)
s.gr.loss.by.event <- apply(s.gr.event.loss.table,1,sum)

# compute average annual loss (AAL) for each location
# not plotted this time
s.gu.aal <- apply(s.gu.event.loss.table*s.event.rates,2,sum)
s.gr.aal <- apply(s.gr.event.loss.table*s.event.rates,2,sum)

# 10) Compute OEP and AEP curve through simulation as before

s.gu.oep.losses <- rep(0,n.years)
s.gu.aep.losses <- rep(0,n.years)
s.gr.oep.losses <- rep(0,n.years)
s.gr.aep.losses <- rep(0,n.years)

for (i in c(1:n.sample.events)) {
  n.occur <- round(n.years*s.event.rates[i])
  year.occur <- as.integer(runif(n.occur)*n.years)+1
  for (j in year.occur) {
	s.gu.oep.losses[j] <- max(s.gu.oep.losses[j],s.gu.loss.by.event[i])
	s.gr.oep.losses[j] <- max(s.gr.oep.losses[j],s.gr.loss.by.event[i])
        s.gu.aep.losses[j] <- s.gu.aep.losses[j]+s.gu.loss.by.event[i]
        s.gr.aep.losses[j] <- s.gr.aep.losses[j]+s.gr.loss.by.event[i]
	}
}

# order the losses and rates
s.ord.gu.oep.losses<- s.gu.oep.losses[order(s.gu.oep.losses,decreasing=T)]
s.ord.gr.oep.losses<- s.gr.oep.losses[order(s.gr.oep.losses,decreasing=T)]
s.ord.gu.aep.losses<- s.gu.aep.losses[order(s.gu.aep.losses,decreasing=T)]
s.ord.gr.aep.losses<- s.gr.aep.losses[order(s.gr.aep.losses,decreasing=T)]

# Plot the results
max.plot <- max(s.ord.gu.aep.losses,na.rm=T)
y.lim <- log10(c(1/n.years,1))
x.lim <- c(0,max.plot)
plot(s.ord.gu.oep.losses,log10(rates),ylim=y.lim,xlim=x.lim,xlab="Loss",ylab="Return Period",las=1,type="l",lwd=2,axes=F)
lines(s.ord.gu.aep.losses,log10(rates),col="red",lwd=2)
lines(s.ord.gr.oep.losses,log10(rates),col="blue",lwd=2)
lines(s.ord.gr.aep.losses,log10(rates),col="green",lwd=2)

p.intervals = log10(n.years)+1
ep.plot <- 10**c(p.intervals:0)
axis(1,at=pretty(c(0,max.plot),p.intervals),labels=pretty(c(0,max.plot),p.intervals))
axis(2,at=log10(1/ep.plot),labels=as.character(ep.plot),las=1)
abline(h=log10(1/ep.plot),col="lightgray",lwd=2,lty=2)
abline(v=pretty(c(0,max.plot),p.intervals),col="lightgray",lwd=2,lty=2)

# 11) Define and calculate other analytical measures
#     - loss cost
#     - probable maximum loss at 250 year RP
#     - contribution of return period losses to AAL
#
# all key figures and analytical quantities can also be derived from the
# simulated results. Here we focus on the original results (expected loss)
# These computations can be done for any loss perspective for any
# aggregation level

# Loss Cost is defined here as the normalized loss per 1000 exposure
# It can be defined for other quantities / currencies as well
loss.cost.gu.aal <- 1000* gu.aal / exposure
loss.cost.gr.aal <- 1000* gr.aal / exposure
cat(paste("loss cost GU / GR", loss.cost.gu.aal, loss.cost.gr.aal),fill=T)

# Probable maximum loss at 250 year RP
# first find the event that is the 250 year RP event
# define this on GR loss
loss.250 <- ord.gr.oep.losses[min(which(rates >= 1/250))]
# find events within 5% of that loss
pml.events <- which(abs(gr.loss.by.event - loss.250) < loss.250*0.05)
# plot 4 of those events
for (i in c(1:4)) {
  hazard.plot <- hazard[pml.events[i],]
  dim(hazard.plot) <- dim.hazard
  hazard.plot <- ifelse(hazard.plot>100,100,hazard.plot)
  image.plot(x.coord,y.coord,hazard.plot,zlim=c(0,100),xlab="Longitude",ylab="Latitude")
  map(add=T,col='white',lwd=2)
  title(main=paste("Loss = ",round(gr.loss.by.event[pml.events[i]])),cex.main=1)
}
title(main="4 PML Events",cex.main=1.7,outer=T)

# Contribution to AAL by return period for AEP
# this can also be computed for the OEP from the events themselves
# It's easy to do this straight from simulation data
#
# We split up the return periods into some intervals (rp.intervals)
# rp.index contains from the ordered AEP vector the loss for the
# corresponding return period
#
# the plot shows the cummulative sum of the contribution to the AAL
# by return period
#
rp.intervals <- c(1,5,10,20,50,100,200,500,1000,n.years)
rp.contribution <- rep(0,length(rp.intervals)-1)
rp.index <- n.years / rp.intervals

for (i in c(1:(length(rp.intervals)-1))) {
  rp.contribution[i] <- sum(ord.gr.aep.losses[c(rp.index[i]:rp.index[i+1])])
}
rp.contribution <- rp.contribution / (sum(gr.aal)*n.years)
par(mfrow=c(2,1))
barplot(cumsum(rp.contribution),names.arg=rp.intervals[2:length(rp.intervals)], cex.names=0.8)
title(main="CDF to AAL by return period",cex.main=1.7)

# EP uncertainty
# let's assume we have 20 years of observed loss data and want to place
# these on the EP curve (can be added to plot)
# at the same time we can plot the uncertainty from the modeled losses
# by splitting up the simulations into 20 year chunks
gr.oep.losses.20 <- gr.oep.losses
# create from modeled losses n.years/20 time series of losses
dim(gr.oep.losses.20) <- c(20,n.years/20)
# aal.20 contains the n.years/20 different modeled AALs
aal.20 <- apply(gr.oep.losses.20,2,mean)
gr.oep.losses.order.20 <- apply(gr.oep.losses.20,2,sort)

# Box Plots show the following (after Climate Blog):
# The rectangle shows the interquartile range (IQR); it goes from the first
# quartile (the 25th percentile) to the third quartile (the 75th percentile).
# The whiskers go from the minimum value to the maximum value unless the distance
# from the minimum value to the first quartile is more than 1.5 times the IQR.
# In that case the whisker extends out to the smallest value within 1.5 times
# the IQR from the first quartile. A similar rule is used for values larger
# than 1.5 times IQR from the third quartile. A special symbol shows the values,
# called outliers, which are smaller or larger than the
# whiskersboxplot(t(gr.oep.losses.order.20))
boxplot(t(gr.oep.losses.order.20),names=round(20/c(20:1),2),xlab="Return Period",ylab="Losses",cex=1)
title(main="Simulated EP Uncertainty",cex.main=1.7)

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

