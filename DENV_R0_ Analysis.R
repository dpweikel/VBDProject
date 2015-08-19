### This file will take the individual parameter samples for all the traits and use
### them to build the DENV R0 model.

## Remember to set your working directory so the code can access the data and 
## necessary supplementary code.
setwd("~/GitHub/VBDProject")

# Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')

# This file contains tools for analysis and visualization.
source("mcmc_utils.R") 

# This file contains the derivatives for the functions.
source("temp_functions.R") 
source("temp_deriv_functions.R") 

## Loading in samples from the Rsave (DENV_ParameterFits.Rsave), which have the samples
## for a, b, c, MDR, EFD, e2a, and PDR created in the DENV_IndividualParameterFitCode.R.

load("DENV_ParameterFits.Rsave")

## Also loading the mu sample made in a separate script.
load("MuFits.Rsave")
mu.samps <- mu.DENV.samps

## Next we set up the temperatures over which we will be evaluating
## R0, etc, as well as the thinning interval for the samples.

# Temperature sequence
temp = seq(5,45,by=0.1)  
t<-length(temp)

# Length of samples
n = dim(a.samps)[1]      

# Thinned samples
thinned<-seq(1, n, by=5) 

# The number of of thinned samples
lthin<-length(thinned)  

# Creating a small constant to keep denominators from being zero.
ec<-0.000001           

## Creating the function encoding the value of R0 as a function of the parameters
myR0<-function(a, b, c, PDR, MDR, EFD, e2a, mu){
  bc = (b*c)
  ((a^2*bc*(EFD*e2a*MDR/(mu + ec)^2)*exp((-mu/(PDR+ec)) + ec))/(mu + ec))^0.5
}

## The following code runs through the samples and calculates the
## posterior trajectories (i.e. as a function of temperatures) of each
## component and of R0, and sets up a matrix to store them (each
## column is a trajectory).

R0<-matrix(NA,t,lthin)
a<-b<-c<-PDR<-MDR<-EFD<-e2a<-mu<-matrix(NA,t,lthin)

for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\n")
  # calculate parameter trajectories
  i<-thinned[j]
  a[,j] = briere(temp, a.samps[i,3], a.samps[i,2], a.samps[i,1])
  PDR[,j] = briere(temp, PDR.samps[i,3], PDR.samps[i,2], PDR.samps[i,1])
  MDR[,j] = briere(temp, MDR.samps[i,3], MDR.samps[i,2], MDR.samps[i,1])
  EFD[,j] = briere(temp, EFD.samps[i,3], EFD.samps[i,2], EFD.samps[i,1])
  e2a[,j] = quad.2.trunc(temp, e2a.samps[i,1], e2a.samps[i,2], e2a.samps[i,3])
  b[,j] = briere(temp, b.samps[i,3], b.samps[i,2], b.samps[i,1])
  c[,j] = briere(temp, c.samps[i,3], c.samps[i,2], c.samps[i,1])
  mu[,j] = quad(temp, mu.samps[i,1], mu.samps[i,2], mu.samps[i,3])
  
  # Calculate Ro equation
  R0[,j]<-myR0(a[,j], b[,j], c[,j], PDR[,j], MDR[,j], EFD[,j], e2a[,j], mu[,j])
  
}

## Next we calculate the posterior mean trajectory of each component of
## R0 and R0 itself. These will be used as part of the uncertainty
## analysis.

a.M<-rowMeans(a)
b.M<-rowMeans(b)
c.M<-rowMeans(c)
PDR.M<-rowMeans(PDR)
MDR.M<-rowMeans(MDR)
EFD.M<-rowMeans(EFD)
e2a.M<-rowMeans(e2a)
mu.M<-rowMeans(mu)
R0.M<-rowMeans(R0)

# Build matrices to hold results

R0.a<-R0.b<-R0.c<-R0.EFD<-R0.e2a<-R0.MDR<-R0.mu<-R0.PDR<-matrix(NA,t,lthin)

## For uncertainty analysis: calculate posterior samples for R0 with
## all but a single component fixed the posterior mean.

for (j in 1:lthin){
  if(j%%100==0) cat("iteration ", j, "\n")
  # Calculating derivative trajectories
  i<-thinned[j]
  ## Calculating R0 with most components set to their means
  R0.a[,j] = myR0(a[,j], b.M, c.M, PDR.M, MDR.M, EFD.M, e2a.M, mu.M)
  R0.b[,j] = myR0(a.M, b[,j], c.M, PDR.M, MDR.M, EFD.M, e2a.M, mu.M)
  R0.c[,j] = myR0(a.M, b.M, c[,j], PDR.M, MDR.M, EFD.M, e2a.M, mu.M)
  R0.EFD[,j] = myR0(a.M, b.M, c.M, PDR.M, MDR.M, EFD[,j], e2a.M, mu.M)
  R0.e2a[,j] = myR0(a.M, b.M, c.M, PDR.M, MDR.M, EFD.M, e2a[,j], mu.M)
  R0.MDR[,j] = myR0(a.M, b.M, c.M, PDR.M, MDR[,j], EFD.M, e2a.M, mu.M)
  R0.mu[,j] =myR0(a.M, b.M, c.M, PDR.M, MDR.M, EFD.M, e2a.M, mu[,j])
  R0.PDR[,j] = myR0(a.M, b.M, c.M, PDR[,j], MDR.M, EFD.M, e2a.M, mu.M)
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed.

R0.q<-  apply(R0, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(R0, 1, FUN=quantile, probs=0.025, na.rm=F)

a.q<-  apply(R0.a, 1, FUN=quantile, probs=0.925)- apply(R0.a, 1, FUN=quantile, probs=0.025)
b.q<- apply(R0.b, 1, FUN=quantile, probs=0.925)- apply(R0.b, 1, FUN=quantile, probs=0.025)
c.q<- apply(R0.c, 1, FUN=quantile, probs=0.925)- apply(R0.c, 1, FUN=quantile, probs=0.025)
EFD.q<- apply(R0.EFD, 1, FUN=quantile, probs=0.925)- apply(R0.EFD, 1, FUN=quantile, probs=0.025)
e2a.q<-apply(R0.e2a, 1, FUN=quantile, probs=0.925)- apply(R0.e2a, 1, FUN=quantile, probs=0.025)
MDR.q<-  apply(R0.MDR, 1, FUN=quantile, probs=0.925)- apply(R0.MDR, 1, FUN=quantile, probs=0.025)
mu.q <-  apply(R0.mu, 1, FUN=quantile, probs=0.925)- apply(R0.mu, 1, FUN=quantile, probs=0.025)
PDR.q<- apply(R0.PDR, 1, FUN=quantile, probs=0.925)- apply(R0.PDR, 1, FUN=quantile, probs=0.025)

## Next plot relative width of quantiles 

# Creating a small constant to keep denominators from being zero.
ec<-0.1 

par(mfrow=c(1,1), bty="n")

plot(temp, a.q/(R0.q +ec), col=2, type="l", ylim=c(0,1), lwd=2,
     xlab="Temperature (C)", ylab="Relative width of quantiles", xlim=c(12,37))
lines(temp, b.q/(R0.q +ec), col=3, lwd=2)
lines(temp, c.q/(R0.q +ec), col=3, lwd=2)
lines(temp, EFD.q/(R0.q +ec), col=4, lwd=2)
lines(temp, e2a.q/(R0.q +ec), col=5, lwd=2)
lines(temp, MDR.q/(R0.q +ec), col=6, lwd=2)
lines(temp, mu.q/(R0.q +ec), col=7, lwd=3)
lines(temp, PDR.q/(R0.q +ec), col=8, lwd=3)

# Adding a legend to the plot.

leg.text<-c("a", "b", "c", "EFD", "e2a", "MDR", "mu", "PDR")
leg.col<-seq(2, 8, by=1)
legend(30, 0.95,  leg.text, col=leg.col, lwd=c(2,2,2,2,2,3,3))


# Calculate the distribution of the lower and upper limits of R0 
# and peak R0.

R0.min<-R0.max<-R0.peak<-rep(NA, length(thinned))

# Plotting the PEAK R0 distribution.
for(i in 1:length(thinned)){
  ww<-which(R0[,i]==max(R0[,i]))
  R0.peak[i]<-temp[ww[1]]
}

# Plotting the MINIMUM R0 distribution.
for(i in 1:length(thinned)){
  ww<-which(R0[,i]>0)
  R0.min[i]<-temp[ww[1]-1]
}

# Plotting the MAXIMUM R0 distribution.

for(i in 1:length(thinned)){
  ww<-which(R0[,i]>0)
  lw<-length(ww)
  R0.max[i]<-temp[ww[lw]+1]
}

# Plotting the mean R0 with it's quantiles, all scaled by max mean R0.

par(mfrow=c(1,1), bty="n")
R0.scale<-max(R0.M)
R0.q<-temp.sim.quants(R0, length(temp))##/R0.scale
plot(temp, R0.M/R0.scale, type="l", col=1, lwd=3, xlim=c(10, 40), ylim=c(0, 2.0),
     ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
add.sim.lines(temp, sim.data=NULL, q=R0.q/R0.scale, mycol=2)

hist(R0.min, xlab="Temp of min R0", freq=TRUE, main="")
hist(R0.peak, xlab="Temp of peak R0", freq=TRUE, main="")
hist(R0.max, xlab="Temp of max R0", freq=TRUE, main="")

## Plotting the data against the thermal responses fits.

# Subsetting the data.

a.data = subset(data.all, trait.name=="gcd")
a.data$trait <- 1/a.data$trait
b.data = subset(data.all, trait.name=="b")
c.data = subset(data.all, trait.name=="c")
EFD.data = subset(data.all, trait.name=="EFD")
e2a.data = subset(data.all, trait.name=="pEA")
MDR.data = subset(data.all, trait.name=="MDR")
PDR.data = subset(data.all, trait.name=="PDR")
p.days = subset(data.all, trait.name=="p/days")
p.days$trait = exp(-(p.days$trait)*(849/45) + ec)
p.almost = subset(data.all, trait.name=="p")
mu.data = rbind(p.days, p.almost)
mu.data$trait = -log(mu.data$trait + ec)

# Creating the mean responses.

a.M<-rowMeans(a)
b.M<-rowMeans(b)
c.M<-rowMeans(c)
MDR.M<-rowMeans(MDR)
PDR.M<-rowMeans(PDR)
EFD.M<-rowMeans(EFD)
e2a.M<-rowMeans(e2a)
mu.M<-rowMeans(mu)

# Getting the HPD intervals.

a.int = HPDinterval(mcmc(t(a)))
b.int = HPDinterval(mcmc(t(b)))
c.int = HPDinterval(mcmc(t(c)))
EFD.int = HPDinterval(mcmc(t(EFD)))
e2a.int = HPDinterval(mcmc(t(e2a)))
MDR.int = HPDinterval(mcmc(t(MDR)))
PDR.int = HPDinterval(mcmc(t(PDR)))
mu.int = HPDinterval(mcmc(t(mu)))

# Setting up the plot area.

par(mfrow=c(3,2))

# Code that plots each trait individually.

plot(trait~T, data=a.data, main="a", xlim=c(10,40), ylim=c(0,0.6))
lines(temp, a.M)
lines(temp, a.int[,1], lty=2, col=2)
lines(temp, a.int[,2], lty=2, col=2)

plot(trait~T, data=b.data, main="b", xlim=c(10,40), ylim=c(0,1))
lines(temp, b.M)
lines(temp, b.int[,1], lty=2, col=2)
lines(temp, b.int[,2], lty=2, col=2)

plot(trait~T, data=c.data, main="c", xlim=c(10,40), ylim=c(0,1))
lines(temp, c.M)
lines(temp, c.int[,1], lty=2, col=2)
lines(temp, c.int[,2], lty=2, col=2)

plot(trait~T, data=EFD.data, main="EFD", xlim=c(10,40), ylim=c(0,15))
lines(temp, EFD.M)
lines(temp, EFD.int[,1], lty=2, col=2)
lines(temp, EFD.int[,2], lty=2, col=2)

plot(trait~T, data=e2a.data, main="e2a", xlim=c(10,40), ylim=c(0,1))
lines(temp, e2a.M)
lines(temp, e2a.int[,1], lty=2, col=2)
lines(temp, e2a.int[,2], lty=2, col=2)

plot(trait~T, data=MDR.data, main="MDR", xlim=c(10,40), ylim=c(0,0.2))
lines(temp, MDR.M)
lines(temp, MDR.int[,1], lty=2, col=2)
lines(temp, MDR.int[,2], lty=2, col=2)

plot(trait~T, data=PDR.data, main="PDR", xlim=c(10,40), ylim=c(0,0.15))
lines(temp, PDR.M)
lines(temp, PDR.int[,1], lty=2, col=2)
lines(temp, PDR.int[,2], lty=2, col=2)

plot(trait~T, data=mu.data, main="mu", xlim=c(10,40), ylim=c(0,1))
lines(temp, mu.M)
lines(temp, mu.int[,1], lty=2, col=2)
lines(temp, mu.int[,2], lty=2, col=2)




