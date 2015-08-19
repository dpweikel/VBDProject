### This file will take resulting analysis from the CHIKV model and run 
### through a bit of sensitvity analysis. 

## Remember to set your working directory so the code can access the data and 
## necessary supplementary code.
setwd("~/GitHub/VBDProject")

# Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')

# This file contains tools for analysis and visualization.
source("mcmc_utils.R") 

# This file contains the derivatives for the functions.
source("temp_deriv_functions.R") 

## Loading in samples from the Rsave (CHIKV_ParameterFits.Rsave), which have the samples
## for a, MDR, TFD, e2a, p created in the CHIKV_IndividualParameterFitCode.R.

load("CHIKV_ParameterFits.Rsave")

## Next we set up the temperatures over which we will be evaluating
## R0, etc, as well as the thinning interval for the samples.

temp = seq(5,45,by=0.1)  ###temperature sequence
t<-length(temp)

n = dim(a.samps)[1]    	### length of samples

thinned<-seq(1, n, by=5)### thinned samples

lthin<-length(thinned)  ### number of of thinned samples

ec<-0.000001            ### small constant used to keep denominators
### from being numerically zero

## Creating the function encoding the value of R0 as a function of the parameters
myR0<-function(a, b, c, PDR, MDR, TFD, e2a, p){
  mu = -log(p+ec)
  EFD = TFD*(1/(a+ec))
  bc = (b*c)
  ((a^2*bc*(EFD*e2a*MDR/(mu)^2)*exp(-mu/(PDR+ec)))/(mu))^0.5
}

## The following code runs through the samples and calculates the
## posterior trajectories (i.e. as a function of temperatures) of each
## component and of R0, and sets up a matrix to store them (each
## column is a trajectory).

R0<-matrix(NA,t,lthin)
a<-b<-c<-PDR<-MDR<-TFD<-e2a<-bc<-p<-matrix(NA,t,lthin)

for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\n")
  # calculate parameter trajectories
  i<-thinned[j]
  a[,j] = briere(temp, a.samps[i,3], a.samps[i,2], a.samps[i,1])
  PDR[,j] = linear(temp, 1, 0)
  MDR[,j] = briere(temp, MDR.samps[i,3], MDR.samps[i,2], MDR.samps[i,1])
  TFD[,j] = briere(temp, TFD.samps[i,3], TFD.samps[i,2], TFD.samps[i,1])
  e2a[,j] = quad.2.trunc(temp, e2a.samps[i,1], e2a.samps[i,2], e2a.samps[i,3])
  b[,j] = briere(temp, b.samps[i,3], b.samps[i,2], b.samps[i,1])
  c[,j] = briere(temp, c.samps[i,3], c.samps[i,2], c.samps[i,1])
  p[,j] = quad.2.trunc(temp, p.samps[i,1], p.samps[i,2],p.samps[i,3])
  
  # Calculate Ro equation
  R0[,j]<-myR0(a[,j], b[,j], c[,j], PDR[,j], MDR[,j], TFD[,j], e2a[,j], p[,j])
  
}

## Next we calculate the posterior mean trajectory of each component of
## R0 and R0 itself. These will be used as part of the uncertainty
## analysis.

a.M<-rowMeans(a)
b.M<-rowMeans(b)
c.M<-rowMeans(c)
PDR.M<-rowMeans(PDR)
MDR.M<-rowMeans(MDR)
TFD.M<-rowMeans(TFD)
e2a.M<-rowMeans(e2a)
p.M<-rowMeans(p)
R0.M<-rowMeans(R0)

## Build matrices to hold derivatives, etc.

dR0.dT= dR0.da = dR0.db = dR0.dc= dR0.dTFD= dR0.de2a= dR0.dMDR= dR0.dp= dR0.dPDR= dR0.R0dT<-matrix(NA,t,lthin)

for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\n")
  ## calculate derivative trajectories
  i<-thinned[j]
  da.dT = d_briere(temp,a.samps[i,3],a.samps[i,2],a.samps[i,1])  
  db.dT = d_briere(temp,b.samps[i,3], b.samps[i,2], b.samps[i,1])  
  dc.dT = d_briere(temp,c.samps[i,3], c.samps[i,2], c.samps[i,1]) 
  dPDR.dT = d_linear(temp, 1, 0)
  dMDR.dT = d_briere(temp,MDR.samps[i,3],MDR.samps[i,2],MDR.samps[i,1])
  dTFD.dT = d_briere(temp,TFD.samps[i,3],TFD.samps[i,2],TFD.samps[i,1])
  de2a.dT = d_quad.2.trunc(temp,e2a.samps[i,1],e2a.samps[i,2],e2a.samps[i,3])
  dp.dT = d_quad.2.trunc(temp,p.samps[i,1],p.samps[i,2],p.samps[i,3])
  
  ## Calculate derivative equations
  ec<-0.000001
  dR0.da[,j] = 1/2*(myR0(a[,j], b.M, c.M, PDR.M, MDR.M, TFD.M, e2a.M, p.M)*2/(a[,j]+ec) * da.dT)
  dR0.db[,j] = 1/2*(myR0(a.M, b[,j], c.M, PDR.M, MDR.M, TFD.M, e2a.M, p.M)/(b[,j]+ec) * db.dT)
  dR0.dc[,j] = 1/2*(myR0(a.M, b.M, c[,j], PDR.M, MDR.M, TFD.M, e2a.M, p.M)/(c[,j]+ec) * dc.dT)
  dR0.dTFD[,j] = 1/2*(myR0(a.M, b.M, c.M, PDR.M, MDR.M, TFD[,j], e2a.M, p.M)/((TFD[,j]/(a.M+ec) + ec)*(dTFD.dT/(a.M+ec))))
  dR0.de2a[,j] = 1/2*(myR0(a.M, b.M, c.M, PDR.M, MDR.M, TFD.M, e2a[,j], p.M)/(e2a[,j]+ec) * de2a.dT)
  dR0.dMDR[,j] = 1/2*(myR0(a.M, b.M, c.M, PDR.M, MDR[,j], TFD.M, e2a.M, p.M)/(MDR[,j]+ec) * dMDR.dT)
  dR0.dp[,j] = 1/2*(myR0(a.M, b.M, c.M, PDR.M, MDR.M, TFD.M, e2a.M, p[,j])*(-3/(-log(p[,j]+ec)+ec)
                                                                            -1/(PDR.M+ec)) *(-1/p[,j])* dp.dT)
  dR0.dPDR[,j] = 1/2*(myR0(a.M, b.M, c.M, PDR[,j], MDR.M, TFD.M, e2a.M, p.M)*(-log(p.M + ec))/(PDR[,j]+ec)^2 * dPDR.dT)
  dR0.dT[,j] =  dR0.da[,j] + dR0.db[,j] + dR0.dc[,j] + dR0.dTFD[,j] + dR0.de2a[,j] + dR0.dMDR[,j] + dR0.dp[,j] + dR0.dPDR[,j]
}

## Relative sensitivities in R0 (different from mordecai paper, just
## normalized by max R0 for each curve) with regard to temperature,
## broken down into individual parameters' contributions:

R0.Med<-apply(R0, 1, FUN=median, na.rm=FALSE)

dR0.R0da=dR0.da/max(R0.M)
dR0.R0db=dR0.db/max(R0.M)
dR0.R0dc=dR0.dc/max(R0.M)
dR0.R0dTFD=dR0.dTFD/max(R0.M)
dR0.R0de2a=dR0.de2a/max(R0.M)
dR0.R0dMDR=dR0.dMDR/max(R0.M)
dR0.R0dp=dR0.dp/max(R0.M)
dR0.R0dPDR=dR0.dPDR/max(R0.M)
dR0.R0dT = dR0.dT/max(R0.M)

## Width of the quantiles of these normalized sensitivities
dR0.q<-  apply(dR0.R0dT, 1, FUN=quantile, probs=0.925, na.rm=F) - apply(dR0.R0dT, 1, FUN=quantile, probs=0.025, na.rm=F)
dR0da.q<-  apply(dR0.R0da, 1, FUN=quantile, probs=0.925)- apply(dR0.R0da, 1, FUN=quantile, probs=0.025)
dR0db.q<- apply(dR0.R0db, 1, FUN=quantile, probs=0.925)- apply(dR0.R0db, 1, FUN=quantile, probs=0.025)
dR0dc.q<- apply(dR0.R0dc, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dc, 1, FUN=quantile, probs=0.025)
dR0dTFD.q<- apply(dR0.R0.TFD, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dTFD, 1, FUN=quantile, probs=0.025)
dR0de2a.q<-apply(dR0.R0de2a, 1, FUN=quantile, probs=0.925)- apply(dR0.R0de2a, 1, FUN=quantile, probs=0.025)
dR0dMDR.q<-  apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.025)
dR0dp.q <-  apply(dR0.R0dp, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dp, 1, FUN=quantile, probs=0.025)
dR0dPDR.q<- apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.025)

plot(temp, dR0.q, type="l", xlim=c(15,35), lwd=2, xlab="Temperature (C)", ylab="width of HPD interval of dR0/dT")
lines(temp, dR0da.q, col=2, lwd=2)
lines(temp, dR0db.q, col=3, lwd=2)
lines(temp, dR0dc.q, col=4, lwd=2)
lines(temp, dR0dTFD.q, col=5, lwd=2)
lines(temp, dR0de2a.q, col=6, lwd=2)
lines(temp, dR0dMDR.q, col=7, lwd=2)
lines(temp, dR0dp.q, col=8, lwd=2)
lines(temp, dR0dPDR.q, col=9, lwd=2)

leg.text<-c("R0", "a", "b", "c", "TFD", "e2a", "MDR", "p", "PDR")
leg.col<-seq(1, 8, by=1)
legend(16.3, 0.7,  leg.text, col=leg.col, lwd=c(2,2,2,2,2,3,3))


## Now plot the relative width of the quantiles of the
## sensitivities. We again include a small offset to keep the
## denominator from being numerically zero. Figure 3(b)
ec=10^(-6)
plot(temp, dR0da.q/(dR0.q +ec), col=2, type="l", ylim=c(0,1), lwd=2, xlab="Temperature (C)", ylab="Relative width of HPD intervals", xlim=c(12, 37))
lines(temp, dR0db.q/(dR0.q +ec), col=3, lwd=2)
lines(temp, dR0dc.q/(dR0.q +ec), col=3, lwd=2)
lines(temp, dR0dTFD.q/(dR0.q +ec), col=4, lwd=2)
lines(temp, dR0de2a.q/(dR0.q +ec), col=5, lwd=2)
lines(temp, dR0dMDR.q/(dR0.q +ec), col=6, lwd=2)
lines(temp, dR0dp.q/(dR0.q +ec), col=7, lwd=3)
lines(temp, dR0dPDR.q/(dR0.q +ec), col=8, lwd=3)
## add legend
leg.text<-c("a", "b", "c", "EFD", "e2a", "MDR", "p", "PDR")
leg.col<-seq(2, 8, by=1)
legend(16.3, 1,  leg.text, col=leg.col, lwd=c(2,2,2,2,2,3,3))



















