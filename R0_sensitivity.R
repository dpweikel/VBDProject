################################################################################
################################################################################
################################################################################
### This file contains code to calculate R0 and perform the
### uncertainty analysis described in Sections 2.2.3, and 3.3 of the
### text


source("mcmc_utils.R") ### tools for analysis and visualization
source("temp_deriv_functions.R") ### this file contains the
                                 ### derivatives for the functions


## First, we identify the files where the posterior samples for each
## component of R0 are saved.

## with informative priors
posts<-c("a_posterior1.Rsave", 
         "PDR_post_prior1.Rsave", 
         "MDR_posterior1.Rsave", 
         "EFD_posterior1.Rsave", 
         "e2a_posterior1.Rsave",
         ##"bc_posterior1_briere.Rsave",  
         "bc_posterior1_quad.Rsave",  
         "mu_posterior1.Rsave")


## We will be using the quadratic response for vector competence, bc,
## and we set this preference here
bcfunc<-"q"

## Next we need to load the samples for each component. Here, we load
## the results for the posterior distributions with informative priors
samps<-samps.b<-samps.q<-NULL

for(i in 1:7){
  load(posts[i])
  if(i==1) a.samps <-samps
  if(i==2) PDR.samps  <-samps
  if(i==3) MDR.samps  <-samps
  if(i==4) EFD.samps <-samps
  if(i==5) e2a.samps <-samps
  if(i==6){
    if(bcfunc=="b") bc.samps <-samps.b
    else bc.samps <-samps.q
  }
  if(i==7) mu.samps <-samps
}


## Next we set up the temperatures over which we will be evaluating
## R0, etc, as well as the thinning interval for the samples.

temp = seq(5,45,by=0.1)	###temperature sequence
t<-length(temp)

n = dim(a.samps)[1]    	### length of samples

thinned<-seq(1, n, by=5)### thinned samples

lthin<-length(thinned)  ### number of of thinned samples

ec<-0.000001            ### small constant used to keep denominators
                        ### from being numerically zero


## function encoding the value of R0 as a function of the parameters
myR0<-function(a, PDR, MDR, EFD, e2a, bc, mu){
  ((a^2*bc*(EFD*e2a*MDR/mu^2)*exp(-mu/(PDR+ec)))/(mu+ec))^0.5
}

## The following code runs through the samples and calculates the
## posterior trajectories (i.e. as a function of temperatures) of each
## component and of R0, and sets up a matrix to store them (each
## column is a trajectory)

R0<-matrix(NA,t,lthin)
a<-PDR<-MDR<-EFD<-e2a<-bc<-mu<-matrix(NA,t,lthin)
for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\n")
  ## calculate parameter trajectories
  i<-thinned[j]
  a[,j] = briere(temp,a.samps[i,3],a.samps[i,2],a.samps[i,1])
  PDR[,j] = briere(temp,PDR.samps[i,3],PDR.samps[i,2],PDR.samps[i,1])
  MDR[,j] = briere(temp,MDR.samps[i,3],MDR.samps[i,2],MDR.samps[i,1])
  EFD[,j] = quad.2(temp,EFD.samps[i,1],EFD.samps[i,2],EFD.samps[i,3])
  e2a[,j] = quad.2.trunc(temp,e2a.samps[i,1],e2a.samps[i,2],e2a.samps[i,3]) 
  if(bcfunc=="b") bc[,j] = briere.trunc(temp,bc.samps[i,3],bc.samps[i,2],bc.samps[i,1])
  else bc[,j] = quad.2.trunc(temp,bc.samps[i,1],bc.samps[i,2],bc.samps[i,3])
  mu[,j] = quad(temp,mu.samps[i,1], -mu.samps[i,2],mu.samps[i,3])

  ## Calculate Ro and m equations
  ##m = EFD[,j]*e2a[,j]*MDR[,j]/mu[,j]^2
  ##R0[,j] = ((a[,j]^2*bc[,j]*m*exp(-mu[,j]/PDR[,j]))/(mu[,j]))^0.5
  R0[,j]<-myR0(a[,j], PDR[,j], MDR[,j], EFD[,j], e2a[,j], bc[,j], mu[,j])
  
}


## Next calculate the posterior mean trajectory of each component of
## R0 and R0 itself. These will be used as part of the sensitivity
## analysis.
a.M<-rowMeans(a)
PDR.M<-rowMeans(PDR)
MDR.M<-rowMeans(MDR)
EFD.M<-rowMeans(EFD)
e2a.M<-rowMeans(e2a)
bc.M<-rowMeans(bc)
mu.M<-rowMeans(mu)
R0.M<-rowMeans(R0)


## Build matrices to hold derivatives, etc. 

dR0.dT= dR0.da= dR0.dbc= dR0.dEFD= dR0.de2a= dR0.dMDR= dR0.dmu= dR0.dPDR= dR0.R0dT<-matrix(NA,t,lthin)

for (j in 1:lthin){
  if(j%%50==0) cat("iteration ", j, "\n")
  ## calculate derivative trajectories
  i<-thinned[j]
  da.dT = d_briere(temp,a.samps[i,3],a.samps[i,2],a.samps[i,1])  
  dPDR.dT = d_briere(temp,PDR.samps[i,3],PDR.samps[i,2],PDR.samps[i,1])
  dMDR.dT = d_briere(temp,MDR.samps[i,3],MDR.samps[i,2],MDR.samps[i,1])
  dEFD.dT = d_quad.2(temp,EFD.samps[i,1],EFD.samps[i,2],EFD.samps[i,3])
  de2a.dT = d_quad.2.trunc(temp,e2a.samps[i,1],e2a.samps[i,2],e2a.samps[i,3])
  if(bcfunc=="b") dbc.dT = d_briere.trunc(temp,bc.samps[i,3],bc.samps[i,2],bc.samps[i,1])
  else dbc.dT = d_quad.2.trunc(temp,bc.samps[i,1],bc.samps[i,2],bc.samps[i,3])
  dmu.dT = d_quad(temp,mu.samps[i,1],-mu.samps[i,2],mu.samps[i,3])
 
  ## Calculate derivative equations
  ec<-0.000001
  dR0.da[,j] = 1/2*(myR0(a[,j], PDR.M, MDR.M, EFD.M, e2a.M, bc.M, mu.M)*2/(a[,j]+ec) * da.dT)
  dR0.dbc[,j] = 1/2*(myR0(a.M, PDR.M, MDR.M, EFD.M, e2a.M, bc[,j], mu.M)/(bc[,j]+ec) * dbc.dT)
  dR0.dEFD[,j] = 1/2*(myR0(a.M, PDR.M, MDR.M, EFD[,j], e2a.M, bc.M, mu.M)/(EFD[,j]+ec) * dEFD.dT)
  dR0.de2a[,j] = 1/2*(myR0(a.M, PDR.M, MDR.M, EFD.M, e2a[,j], bc.M, mu.M)/(e2a[,j]+ec) * de2a.dT)
  dR0.dMDR[,j] = 1/2*(myR0(a.M, PDR.M, MDR[,j], EFD.M, e2a.M, bc.M, mu.M)/(MDR[,j]+ec) * dMDR.dT)
  dR0.dmu[,j] = 1/2*(myR0(a.M, PDR.M, MDR.M, EFD.M, e2a.M, bc.M, mu[,j])*(-3/(mu[,j]+ec)
           -1/(PDR.M+ec)) * dmu.dT)
  dR0.dPDR[,j] = 1/2*(myR0(a.M, PDR[,j], MDR.M, EFD.M, e2a.M, bc.M, mu.M)*mu.M/(PDR[,j]+ec)^2 * dPDR.dT)
  dR0.dT[,j] =  dR0.da[,j] + dR0.dbc[,j] + dR0.dEFD[,j] + dR0.de2a[,j] + dR0.dMDR[,j] + dR0.dmu[,j] + dR0.dPDR[,j]
}


## Relative sensitivities in R0 (different from mordecai paper, just
## normalized by max R0 for each curve) with regard to temperature,
## broken down into individual parameters' contributions:

R0.Med<-apply(R0, 1, FUN=median, na.rm=FALSE)

dR0.R0da=dR0.da/max(R0.M)
dR0.R0dbc=dR0.dbc/max(R0.M)
dR0.R0dEFD=dR0.dEFD/max(R0.M)
dR0.R0de2a=dR0.de2a/max(R0.M)
dR0.R0dMDR=dR0.dMDR/max(R0.M)
dR0.R0dmu=dR0.dmu/max(R0.M)
dR0.R0dPDR=dR0.dPDR/max(R0.M)
dR0.R0dT = dR0.dT/max(R0.M)

## Width of the quantiles of these normalized sensitivities
dR0.q<-  apply(dR0.R0dT, 1, FUN=quantile, probs=0.925, na.rm=F) - apply(dR0.R0dT, 1, FUN=quantile, probs=0.025, na.rm=F)
dR0da.q<-  apply(dR0.R0da, 1, FUN=quantile, probs=0.925)- apply(dR0.R0da, 1, FUN=quantile, probs=0.025)
dR0dbc.q<- apply(dR0.R0dbc, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dbc, 1, FUN=quantile, probs=0.025)
dR0dEFD.q<- apply(dR0.R0dEFD, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dEFD, 1, FUN=quantile, probs=0.025)
dR0de2a.q<-apply(dR0.R0de2a, 1, FUN=quantile, probs=0.925)- apply(dR0.R0de2a, 1, FUN=quantile, probs=0.025)
dR0dMDR.q<-  apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.025)
dR0dmu.q <-  apply(dR0.R0dmu, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dmu, 1, FUN=quantile, probs=0.025)
dR0dPDR.q<- apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.025)

plot(temp, dR0.q, type="l", xlim=c(15,35), lwd=2, xlab="Temperature (C)", ylab="width of HPD interval of dR0/dT")
lines(temp, dR0da.q, col=2, lwd=2)
lines(temp, dR0dbc.q, col=3, lwd=2)
lines(temp, dR0dEFD.q, col=4, lwd=2)
lines(temp, dR0de2a.q, col=5, lwd=2)
lines(temp, dR0dMDR.q, col=6, lwd=2)
lines(temp, dR0dmu.q, col=7, lwd=2)
lines(temp, dR0dPDR.q, col=8, lwd=2)

leg.text<-c("R0", "a", "bc", "EFD", "e2a", "MDR", "mu", "PDR")
leg.col<-seq(1, 8, by=1)
legend(16.3, 0.7,  leg.text, col=leg.col, lwd=c(2,2,2,2,2,3,3))


## Now plot the relative width of the quantiles of the
## sensitivities. We again include a small offset to keep the
## denominator from being numerically zero. Figure 3(b)
ec=10^(-6)
plot(temp, dR0da.q/(dR0.q +ec), col=2, type="l", ylim=c(0,1), lwd=2, xlab="Temperature (C)", ylab="Relative width of HPD intervals", xlim=c(12, 37))
lines(temp, dR0dbc.q/(dR0.q +ec), col=3, lwd=2)
lines(temp, dR0dEFD.q/(dR0.q +ec), col=4, lwd=2)
lines(temp, dR0de2a.q/(dR0.q +ec), col=5, lwd=2)
lines(temp, dR0dMDR.q/(dR0.q +ec), col=6, lwd=2)
lines(temp, dR0dmu.q/(dR0.q +ec), col=7, lwd=3)
lines(temp, dR0dPDR.q/(dR0.q +ec), col=8, lwd=3)
## add legend
leg.text<-c("a", "bc", "EFD", "e2a", "MDR", "mu", "PDR")
leg.col<-seq(2, 8, by=1)
legend(16.3, 1,  leg.text, col=leg.col, lwd=c(2,2,2,2,2,3,3))

