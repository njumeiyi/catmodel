# =======================================================================
#
# Workshop 4: Principal Component and Maximum Covariance Analyses
#
# Analysis of Climate and Weather Data
# Christoph Frei and Michael Keller
# 
# =======================================================================
#
# -----------------------------------------------------------------------
# Problem 1: Discover the North Atlantic Oscillation
# -----------------------------------------------------------------------

# Undertake a Principal Component Analysis of the winter-time sea level 
# pressure (SLP) variations over the North Atlantic and Europe. Identify 
# and interpret the leading modes of circulation variability. 

# Load library
library(ACWD)
library(pcaXcca)

# Dataset slp.eu.era40 contains monthly fields of sea level pressure. 
# The dataset was derived from a Reanalysis of observations using a 
# global numerical data assimilation system. (The ECMWF ERA40). Make 
# yourself familiar with the dataset (resolution, temporal extent, etc.)
help(slp.eu.era40)
data(slp.eu.era40)

# Extract the winter months (Nov-Mar) from the dataset and convert the 
# 3D array into a standard data matrix with grid points in columns and 
# times in rows. 
help(field2dmat)
slp.mat <- field2dmat(slp.eu.era40,months=c(1,2,3,11,12))

# view the dataset
slp.mat[1:10,1:10]
attributes(slp.mat)

# What is the dimension of the entire phasespace? How many samples 
# are there? (Number of points in the data cloud.) What is the dimension 
# of the variance-covariance matrix?
dim(slp.mat)

# Do we need to transform the data because of asymmetry / non-normality?
# Test a few example grid points!
qqnorm(slp.mat[,345])
hist(slp.mat[,345])
# Transformation is not necessary. The distribution is reasonably 
# symetric and qq-plot shows reasonable normality.

# Perform the Principal Component Analysis. What is the meaning  
# of the various elements of the resulting data object. 
help(pca)
hh <- pca(slp.mat,cor=FALSE)
names(hh)

# Make a scree plot of the result (standard deviation of the principal 
# component score as a function of the number of the component.) 
# How much of the total variance do the first 4 components explain?
help(scree.plot)
# Hint: Numbers of explained variance fraction can be directly read 
# from list elements EVF and CEVF of the object returned from pca
scree.plot(hh,num.pcs=20,main="PCA: SLP Europe/Atlantic, NDJFM")
hh$EVF       # (44.5%,18.3%,14.3%,8.7%)
hh$CEVF      # (85.7%)

# Make a plot of the first PC Loading. Interpret the circulation feature 
# of this dominant mode: How is the circulation during months with a 
# positive (negative) score of the first PC. For this purpose also make 
# a plot of the mean SLP pattern and the SLP pattern associated with one 
# standard deviation of the PC score (i.e. mean + sdev*loading). 
# (If you whish, you can also examine the loading of higher order 
# Principal Components.)
help(plot.pca.vec)
pc.num <- 1
X11(height=6,width=8)
plot.pca.vec(hh,what="loading",num=pc.num,plot.fun=plot.geo.field,
             leg.space=0.15,xmaplim=c(-180,180),ymaplim=c(0,90))
title(main=paste("PC Loading",pc.num))

X11(height=6,width=8)
plot.pca.vec(hh,what="center",
             breaks=seq(998,1024,by=2),plot.fun=plot.geo.field,
             leg.space=0.15,xmaplim=c(-180,180),ymaplim=c(0,90))
title(main="mean SLP (NDJFM)")

X11(height=6,width=8) 
plot.pca.vec(hh,what="center+loading",num=pc.num,amplitude=+1,
                   breaks=seq(998,1024,by=2),plot.fun=plot.geo.field,
                   leg.space=0.15,xmaplim=c(-180,180),ymaplim=c(0,90))
title(main=paste("Center + PC Loading",pc.num))

X11(height=6,width=8)     # open a new graphic window for comparison
plot.pca.vec(hh,what="center+loading",num=pc.num,amplitude=-1,
                   breaks=seq(998,1024,by=2),plot.fun=plot.geo.field,
                   leg.space=0.15,xmaplim=c(-180,180),ymaplim=c(0,90))
title(main=paste("Center - PC Loading",pc.num))

dev.off()       # close the graphic windows.
dev.off()
dev.off()
dev.off()

# What is the correlation between the first PC-score and the pressure 
# difference between Gibraltar and Reykjavik (the NAO-index). 
# Why is there negative correlation.
help(nao.index) 
help(extr.nao)
nao <- extr.nao(months=c("Jan", "Feb", "Mar", "Nov", "Dec"),
                start="1957-11",end="2002-03")
score1 <- hh$scores[,1]
cor(score1,nao)

# Produce a time series plot of the first Principal Component Score 
# and compare it to the time series of the North Atlantic Oscillation 
# Index.
score1.normalized <- score1/sd(score1)
str2tim <- function(str) {
       yy <- as.numeric(substr(str,1,4))
       mm <- substr(str,6,7)
       im <- switch(mm,"01"=0.1,"02"=0.3,"03"=0.5,"11"=0.7,"12"=0.9)
       return(yy+im)
}
tim <- sapply(row.names(hh$dat),FUN=str2tim)

X11(height=6,width=10)
plot(x=tim,y=-score1.normalized,type="l",col="blue",lwd=2,ylab="",xlab="")
lines(x=tim,y=(nao-mean(nao))/sd(nao),col="red",lwd=2)
legend(x=1962,y=2.5,cex=1.3,
       legend=c("minus PC score 1","NAO-Index"),
       lwd=2,col=c("blue","red"),bty="n")
dev.off()

# Feb. 1990 (Vivian) and Dec. 1999 (Anatol, Lothar, Martin) were months 
# with heavy stroms over northern Europe. Is this plausible from the 
# results of the PCA?

# Where should we truncate the PC spectrum. Consider the truncation 
# rules based on the Eigenvalue Spectrum of random noise and the 
# Preisendorfer Rule-N. How much of the total variance would the 
# truncated representation of the dataset explain? Is it feasible to 
# truncate at the suggested points according to the North et al. rule 
# of thumb. Which truncation would you choose?
help(scree.plot)
scree.plot(hh,num.pcs=50,log=TRUE,add.text=FALSE)    # be patient here !!!
lamdas <- (hh$sdev)^2
hh$CEVF
# Eigenvalue Spectrum (L=19, there is a jump from 19 to 20, 98.8%)
# Rule-N (L=7, 93.4%)
(lamdas[19]-lamdas[20])/lamdas[20] >= sqrt(2/225)
(lamdas[18]-lamdas[19])/lamdas[19] >= sqrt(2/225)
(lamdas[7]-lamdas[8])/lamdas[8] >= sqrt(2/225)
# both truncations (19, 7) would satisfy the North et al. rule of thumb,
# but a truncation at 18 would not be feasible.

# For a model aiming at simplicity the truncation at L=7 is better. 
# If a model is desired that is capable of resolving fine scale detail 
# of the slp fields we could justify L=19 but there is only a little 
# improvement in CEVF.

dev.off()
q()

# -----------------------------------------------------------------------
# Problem 2: Relationship between winter-time temperature
#            and the sea level pressure in Europe.
# -----------------------------------------------------------------------

# Investigate the relationships between the sea level pressure (SLP) 
# over the Euro-Atlantic region and temperatures over Europe in January. 
# For this purpose use the technique of Maximum Covariance Analysis. 

# Start a new R session and load library
library(ACWD)
library(pcaXcca)

# Make yourself familiar with the two multivariate datasets: 
# slp.eu.era40 for sea level pressure and t2m.eu.era40 for temps. 
# Both datasets are from the ECMWF Reanalysis (ERA40), i.e. pseudo-
# observations.
help(t2m.eu.era40)
help(slp.eu.era40)
data(t2m.eu.era40)
data(slp.eu.era40)

# Extract Januaries from the two datasets and establish the data 
# matrices needed for MCA
slp.mat <- field2dmat(slp.eu.era40,months=c(1))
t2m.mat <- field2dmat(t2m.eu.era40,months=c(1))

# Calculate the MCA  
dat.mca <- mca(left=slp.mat,right=t2m.mat)

# Make plots of the singular vectors and make some physical 
# interpretations. To gain some more insight in the associated 
# circulation also make graphs where the sea level pressure singular 
# vectors are added to the mean field. Get insight into the first 3 
# coupled modes.
ww <- 8                # preparation for graphics
asp.left <- 0.7        # preparation for graphics
asp.right <- 0.9         # preparation for graphics
nr <- 1           # number of the singular vector to plot
X11(width=ww,height=asp.left*ww)         # window for slp singular vector
plot.mca.vec(dat.mca,what="left.singular.vector",num=nr,
             plot.fun="plot.geo.field",
             nlevels=10,leg.space=0.15)
title(main=paste("Sea Level Pressure, Singular Vector ",nr,sep=""))
X11(width=ww,height=asp.right*ww)        # window for t2m singular vector
plot.mca.vec(dat.mca,what="right.singular.vector",num=nr,
             plot.fun="plot.geo.field",
             nlevels=10,leg.space=0.15,ptype="image")
title(main=paste("2m Temperature, Singular Vector ",nr,sep=""))
dev.off()
dev.off()

sig <- 1
sigstr <- ifelse(sig==1,"+","-")
X11(width=ww,height=asp.left*ww)      # window for mean +/- singular vec
plot.mca.vec(dat.mca,what="left.center+singular.vector",
             num=nr,amplitude=sig,
             plot.fun="plot.geo.field",breaks=seq(992,1024,by=2),
             leg.space=0.15)
title(main=paste("Sea Level Pressure, Mean",sigstr,
                 "Singular Vector ",nr,sep=""))
dev.off()

# How large are the correlations between the left and right coefficients 
# for the first 3 modes. 
vars.u <- apply(dat.mca$coeff.left,FUN=var,MARGIN=c(2))
vars.v <- apply(dat.mca$coeff.right,FUN=var,MARGIN=c(2))
cors <- dat.mca$singular.values/(sqrt(vars.u)*sqrt(vars.v))

# How large are the squared covariance fractions of the first 3 modes.
help(plot.evf.scf)
dat.mca$SCF              # mode 1: 92%, mode 2: 6.1%, mode 3: 1.3% 
X11(width=6,height=8)
plot.evf.scf(dat.mca,what="SCF",percentage=TRUE,num=10,
         main="squared covariance fraction",
         ylab="squared slp-temp covariance fraction (%)")
dev.off()

