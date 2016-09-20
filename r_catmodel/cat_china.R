rm(list = ls())
library('fields')

library('copula')
library('evd')

x.coord <- 73:136
y.coord <- 3:54

dim.hazard <- c(length(x.coord),length(y.coord))
n.grids  <- dim.hazard[1]*dim.hazard[2]

n.events <- 10000
rate.per.event <- 1/10000
event.rates <- rep(rate.per.event,n.events)
rho.par <- 0.5

myCop.norm <- ellipCopula(family="normal",dim=n.grids,dispstr="ex",param=rho.par)
gpd.par <- list(loc=1,scale=8,shape=0.5)
myMvd <- mvdc(copula=myCop.norm,margins=rep("gpd",n.grids),paramMargins=rep(list(gpd.par),n.grids))
#hazard <- rmvdc(myMvd,n.events)
hazard <- rMvdc(n.events,myMvd)

hazard.plot <- hazard[1,]
dim(hazard.plot) <- dim.hazard
hazard.plot <- ifelse(hazard.plot>100,100,hazard.plot)
image.plot(x.coord,y.coord,hazard.plot,zlim=c(0,100),xlab="Longitude",ylab="Latitude")
map(add=T,col='white',lwd=2)


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
