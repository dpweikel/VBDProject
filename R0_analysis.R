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

## with uninformative priors
uninf<-c("a_noninf.Rsave",
         "PDR_noninf.Rsave",
         "MDR_noninf.Rsave", 
         "EFD_noninf.Rsave", 
         "e2a_noninf.Rsave", 
         ##"bc_noninf_briere.Rsave",  
         "bc_noninf_quad.Rsave",  
         "mu_noninf.Rsave") 


## We will be using the quadratic response for vector competence, bc,
## and we set this preference here
bcfunc<-"q"

derivs<-TRUE

## to do the analysis for the posteriors using uninformative priors
## use this:
# posts<-uninf

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
## R0 and R0 itself. These will be used as part of the uncertainty
## analysis.
a.M<-rowMeans(a)
PDR.M<-rowMeans(PDR)
MDR.M<-rowMeans(MDR)
EFD.M<-rowMeans(EFD)
e2a.M<-rowMeans(e2a)
bc.M<-rowMeans(bc)
mu.M<-rowMeans(mu)
R0.M<-rowMeans(R0)

## Build matrices to hold results
R0.a<-R0.bc<-R0.EFD<-R0.e2a<-R0.MDR<-R0.mu<-R0.PDR<-matrix(NA,t,lthin)

## For uncertainty analysis: calculate posterior samples for R0 with
## all but a single component fixed the posterior mean
for (j in 1:lthin){
  if(j%%100==0) cat("iteration ", j, "\n")
  ## calculate derivative trajectories
  i<-thinned[j]
  ## Calculate R0 with most components set to their means
  R0.a[,j] = myR0(a[,j], PDR.M, MDR.M, EFD.M, e2a.M, bc.M, mu.M)
  R0.bc[,j] = myR0(a.M, PDR.M, MDR.M, EFD.M, e2a.M, bc[,j], mu.M)
  R0.EFD[,j] = myR0(a.M, PDR.M, MDR.M, EFD[,j], e2a.M, bc.M, mu.M)
  R0.e2a[,j] = myR0(a.M, PDR.M, MDR.M, EFD.M, e2a[,j], bc.M, mu.M)
  R0.MDR[,j] = myR0(a.M, PDR.M, MDR[,j], EFD.M, e2a.M, bc.M, mu.M)
  R0.mu[,j] =myR0(a.M, PDR.M, MDR.M, EFD.M, e2a.M, bc.M, mu[,j])
  R0.PDR[,j] = myR0(a.M, PDR[,j], MDR.M, EFD.M, e2a.M, bc.M, mu.M)
}

## Calculate the distance within the inner 95% quantile for R0 overall
## (R0.q) and for the posterior of R0 with each component held fixed
R0.q<-  apply(R0, 1, FUN=quantile, probs=0.925, na.rm=F)- apply(R0, 1, FUN=quantile, probs=0.025, na.rm=F)

a.q<-  apply(R0.a, 1, FUN=quantile, probs=0.925)- apply(R0.a, 1, FUN=quantile, probs=0.025)
bc.q<- apply(R0.bc, 1, FUN=quantile, probs=0.925)- apply(R0.bc, 1, FUN=quantile, probs=0.025)
EFD.q<- apply(R0.EFD, 1, FUN=quantile, probs=0.925)- apply(R0.EFD, 1, FUN=quantile, probs=0.025)
e2a.q<-apply(R0.e2a, 1, FUN=quantile, probs=0.925)- apply(R0.e2a, 1, FUN=quantile, probs=0.025)
MDR.q<-  apply(R0.MDR, 1, FUN=quantile, probs=0.925)- apply(R0.MDR, 1, FUN=quantile, probs=0.025)
mu.q <-  apply(R0.mu, 1, FUN=quantile, probs=0.925)- apply(R0.mu, 1, FUN=quantile, probs=0.025)
PDR.q<- apply(R0.PDR, 1, FUN=quantile, probs=0.925)- apply(R0.PDR, 1, FUN=quantile, probs=0.025)

## Next plot relative width of quantiles (Figure 3(a))
ec<-0.1 ## small constant used to keep denominators from being numerically zero

plot(temp, a.q/(R0.q +ec), col=2, type="l", ylim=c(0,1), lwd=2,
     xlab="Temperature (C)", ylab="Relative width of quantiles", xlim=c(12,37))
lines(temp, bc.q/(R0.q +ec), col=3, lwd=2)
lines(temp, EFD.q/(R0.q +ec), col=4, lwd=2)
lines(temp, e2a.q/(R0.q +ec), col=5, lwd=2)
lines(temp, MDR.q/(R0.q +ec), col=6, lwd=2)
lines(temp, mu.q/(R0.q +ec), col=7, lwd=3)
lines(temp, PDR.q/(R0.q +ec), col=8, lwd=3)
## Add legend
leg.text<-c("a", "bc", "EFD", "e2a", "MDR", "mu", "PDR")
leg.col<-seq(2, 8, by=1)
legend(12, 1,  leg.text, col=leg.col, lwd=c(2,2,2,2,2,3,3))


## Calculate the distribution of the lower and upper limits of R0 and peak R0.
R0.min<-R0.max<-R0.peak<-rep(NA, length(thinned))
## PEAK
for(i in 1:length(thinned)){
  ww<-which(R0[,i]==max(R0[,i]))
  R0.peak[i]<-temp[ww[1]]
}
## MINIMUM
for(i in 1:length(thinned)){
  ww<-which(R0[,i]>0)
  R0.min[i]<-temp[ww[1]-1]
}
## MAXIMUM
for(i in 1:length(thinned)){
  ww<-which(R0[,i]>0)
  lw<-length(ww)
  R0.max[i]<-temp[ww[lw]+1]
}

## plot mean R0 with it's quantiles, all scaled by max mean R0
par(mfrow=c(1,4), bty="n")
R0.scale<-max(R0.M)
R0.q<-temp.sim.quants(R0, length(temp))##/R0.scale
plot(temp, R0.M/R0.scale, type="l", col=1, lwd=3, xlim=c(15, 35), ylim=c(0, 3),
     ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
add.sim.lines(temp, sim.data=NULL, q=R0.q/R0.scale, mycol=2)

hist(R0.min, xlab="Temp of min R0", freq=TRUE, main="")
hist(R0.peak, xlab="Temp of peak R0", freq=TRUE, main="")
hist(R0.max, xlab="Temp of max R0", freq=TRUE, main="")

### Save all the important bits that have been calculated into a new
### data frame.

## for the posteriors calculated with informative priors
R0.posts<-list()
R0.posts$mins<-R0.min
R0.posts$max<-R0.max
R0.posts$peaks<-R0.peak
R0.posts$temp<-temp
R0.posts$q<-R0.q
R0.posts$mean<-R0.M
R0.posts$scale<-R0.scale


## for the posteriors calculated with uninformative priors
##R0.uninf<-list()
##R0.uninf$mins<-R0.min
##R0.uninf$max<-R0.max
##R0.uninf$peaks<-R0.peak
##R0.uninf$temp<-temp
##R0.uninf$q<-R0.q
##R0.uninf$mean<-R0.M
##R0.uninf$scale<-R0.scale




## Make a comparison plot to see the impact of using informative vs
## uninformative priors Figure 2

## R0 overall (Fig 2 top)
plot(temp, R0.posts$mean/R0.posts$scale,  type="l", col="red", lwd=3,
     xlim=c(15, 35), ylim=c(0, 3),
     ylab=expression(paste("relative ", R[0], sep="")), xlab="Temperature (C)")
lines(temp, R0.posts$q[1,]/R0.posts$scale, lwd=3, col="red", lty="dotted")
lines(temp, R0.posts$q[2,]/R0.posts$scale, lwd=3, col="red", lty="dotted")

lines(temp, R0.uninf$mean/R0.uninf$scale, lwd=3, col="blue", lty="dashed")
lines(temp, R0.uninf$q[1,]/R0.uninf$scale, lwd=3, col="blue", lty="dotted")
lines(temp, R0.uninf$q[2,]/R0.uninf$scale, lwd=3, col="blue", lty="dotted")


## R0 min, peak, max (Fig 2 bottom)
par(mfrow=c(1,3), bty="n")
## MINIMUM
plot(density(R0.uninf$mins),  type="l", col="blue", lwd=3,
     xlim=c(15, 25), ylim=c(0, 0.25), ylab="",
     main=expression(paste("density of lower limit of ", R[0], sep=" ")),
     xlab="Temperature (C)", lty="dashed")
lines(density(R0.posts$mins), lwd=3, col=2)
## PEAK
plot(density(R0.uninf$peaks),  type="l", col="blue", lwd=3,
     xlim=c(22, 28), ylim=c(0, 0.5), ylab="",
     main=expression(paste("density of peak ", R[0], sep=" ")),
     xlab="Temperature (C)", lty="dashed")
lines(density(R0.posts$peaks), lwd=3, col=2)
## MAXIMUM
plot(density(R0.uninf$max),  type="l", col="blue", lwd=3,
     xlim=c(24, 35), ylim=c(0, 0.4), ylab="",
     main=expression(paste("density of upper limit of ", R[0], sep=" ")),
     xlab="Temperature (C)", lty="dashed")
lines(density(R0.posts$max), lwd=3, col=2)


################################################################################
################################################################################
################################################################################
### Derivatives for the Sensitivity Analysis

a.M<-rowMeans(a)
PDR.M<-rowMeans(PDR)
MDR.M<-rowMeans(MDR)
EFD.M<-rowMeans(EFD)
e2a.M<-rowMeans(e2a)
bc.M<-rowMeans(bc)
mu.M<-rowMeans(mu)
R0.M<-rowMeans(R0)


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

### Save everything
##
save(temp, R0,dR0.da,dR0.dbc,dR0.dEFD,dR0.de2a,dR0.dMDR,dR0.dmu,dR0.dPDR,dR0.dT, #, dR0.R0da,dR0.R0dbc,dR0.R0dEFD,dR0.R0de2a,dR0.R0dMDR,dR0.R0dmu,dR0.R0dPDR,dR0.R0dT,
file = "R0_trajectories_derivs.Rsave")

##Relative change in R0 (different from mordecai paper, just
##normalized by max R0 for each curve) with regard to temperature,
##broken down into individual parameters' contributions:

R0.M<-rowMeans(R0)
R0.Med<-apply(R0, 1, FUN=median, na.rm=FALSE)

dR0.R0da=dR0.da/max(R0.M)
dR0.R0dbc=dR0.dbc/max(R0.M)
dR0.R0dEFD=dR0.dEFD/max(R0.M)
dR0.R0de2a=dR0.de2a/max(R0.M)
dR0.R0dMDR=dR0.dMDR/max(R0.M)
dR0.R0dmu=dR0.dmu/max(R0.M)
dR0.R0dPDR=dR0.dPDR/max(R0.M)
dR0.R0dT = dR0.dT/max(R0.M)

plot(temp, dR0.R0da[,1], type="l", ylim=c(-0,1), xlim=c(15, 35))
for(i in 2:100) lines(temp, dR0.R0da[,i])

plot(temp, rowMeans(dR0.R0dT, na.rm=T), type="l", col=1, ylim=c(-.5, 0.5), xlim=c(15, 35), lwd=3)
##plot(temp, apply(dR0.R0dT, 1, FUN=median, na.rm=T), type="l", col=1, ylim=c(-0.75, 0.65), xlim=c(15, 35), lwd=3)
ec<-10^(-8)
plot(temp, rowMeans(dR0.da)/(10*R0.M), col=2, lwd=1, type="l", ylim=c(-.2, 0.2), xlim=c(12, 33), xlab="Temperature (C)", ylab="")
lines(temp, rowMeans(dR0.dbc)/(10*R0.M), col=3)
lines(temp, rowMeans(dR0.dEFD)/(10*R0.M), col=4, lwd=1)
lines(temp, rowMeans(dR0.de2a)/(10*R0.M), col=5)
lines(temp, rowMeans(dR0.dMDR)/(10*R0.M), col=6)
lines(temp, rowMeans(dR0.dmu)/(10*R0.M), col=7, lwd=3)
lines(temp, rowMeans(dR0.dPDR)/(10*R0.M), col=8, lwd=3)
abline(h=0)

plot(temp, apply(dR0.da, 1, FUN=median, na.rm=T)/(10*R0.Med), col=2, lwd=1, type="l", ylim=c(-.1, 0.1), xlim=c(16, 33), xlab="Temperature (C)", ylab="")
lines(temp, apply(dR0.dbc, 1, FUN=median, na.rm=T)/(10*R0.Med), col=3)
lines(temp, apply(dR0.dEFD, 1, FUN=median, na.rm=T)/(10*R0.Med), col=4, lwd=1)
lines(temp, apply(dR0.de2a, 1, FUN=median, na.rm=T)/(10*R0.Med), col=5)
lines(temp, apply(dR0.dMDR, 1, FUN=median, na.rm=T)/(10*R0.Med), col=6)
lines(temp, apply(dR0.dmu, 1, FUN=median, na.rm=T)/(10*R0.Med), col=7, lwd=3)
lines(temp, apply(dR0.dPDR, 1, FUN=median, na.rm=T)/(10*R0.Med), col=8, lwd=3)
abline(h=0)

lines(temp, apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.025), lty="dashed", col=8, lwd=3)
lines(temp, apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.925), lty="dashed", col=8, lwd=3)

dR0.q<-  apply(dR0.R0dT, 1, FUN=quantile, probs=0.925, na.rm=F) - apply(dR0.R0dT, 1, FUN=quantile, probs=0.025, na.rm=F)
dR0da.q<-  apply(dR0.R0da, 1, FUN=quantile, probs=0.925)- apply(dR0.R0da, 1, FUN=quantile, probs=0.025)
dR0dbc.q<- apply(dR0.R0dbc, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dbc, 1, FUN=quantile, probs=0.025)
dR0dEFD.q<- apply(dR0.R0dEFD, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dEFD, 1, FUN=quantile, probs=0.025)
dR0de2a.q<-apply(dR0.R0de2a, 1, FUN=quantile, probs=0.925)- apply(dR0.R0de2a, 1, FUN=quantile, probs=0.025)
dR0dMDR.q<-  apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.025)
dR0dmu.q <-  apply(dR0.R0dmu, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dmu, 1, FUN=quantile, probs=0.025)
dR0dPDR.q<- apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.925)- apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.025)

plot(temp, dR0.q, type="l")
lines(temp, dR0da.q, col=2)
lines(temp, dR0dbc.q, col=3)
lines(temp, dR0dEFD.q, col=4)
lines(temp, dR0de2a.q, col=5)
lines(temp, dR0dMDR.q, col=6)
lines(temp, dR0dmu.q, col=7)
lines(temp, dR0dPDR.q, col=8)

## dR0.q and its components are on a much different scale (0-0.7
## instead of 0-12, for instance), and so needs a smaller ec value
ec=10^(-6)
plot(temp, dR0da.q/(dR0.q +ec), col=2, type="l", ylim=c(0,1), lwd=2, xlab="Temperature (C)", ylab="proportion", xlim=c(12, 37))
lines(temp, dR0dbc.q/(dR0.q +ec), col=3, lwd=2)
lines(temp, dR0dEFD.q/(dR0.q +ec), col=4, lwd=2)
lines(temp, dR0de2a.q/(dR0.q +ec), col=5, lwd=2)
lines(temp, dR0dMDR.q/(dR0.q +ec), col=6, lwd=2)
lines(temp, dR0dmu.q/(dR0.q +ec), col=7, lwd=3)
lines(temp, dR0dPDR.q/(dR0.q +ec), col=8, lwd=3)

leg.text<-c("a", "bc", "EFD", "e2a", "MDR", "mu", "PDR")
###"(A1) DT", expression(paste("(A2) DT (", L[m], ")")), expression(paste("(B) ", X[eff], sep="")), expression(paste("(C) ", X[eff], " (", L[m], ")", sep="")))
leg.col<-seq(2, 8, by=1)

legend(16.3, 1,  leg.text, col=leg.col, lwd=c(2,2,2,2,2,3,3))


##Relative change in R0 with regard to temperature, broken down into
##individual parameters' contributions:

dR0.R0da=dR0.da/R0
dR0.R0dbc=dR0.dbc/R0
dR0.R0dEFD=dR0.dEFD/R0
dR0.R0de2a=dR0.de2a/R0
dR0.R0dMDR=dR0.dMDR/R0
dR0.R0dmu=dR0.dmu/R0
dR0.R0dPDR=dR0.dPDR/R0
dR0.R0dT = dR0.R0da + dR0.R0dbc + dR0.R0dEFD + dR0.R0de2a + dR0.R0dMDR+ dR0.R0dmu + dR0.R0dPDR


plot(temp, apply(dR0.R0dT,  1, FUN=median, na.rm=T), type="l", col=1, ylim=c(-3.5, 3.5), xlim=c(16, 34), lwd=2)
lines(temp, apply(dR0.R0dT, 1, FUN=quantile, probs=0.025, na.rm=T), lty="dashed", lwd=2)
lines(temp, apply(dR0.R0dT, 1, FUN=quantile, probs=0.925, na.rm=T), lty="dashed", lwd=2)
      
lines(temp,  apply(dR0.R0da,  1, FUN=median, na.rm=T), col=2)
lines(temp, apply(dR0.R0da, 1, FUN=quantile, probs=0.025, na.rm=T), lty="dashed", col=2)
lines(temp, apply(dR0.R0da, 1, FUN=quantile, probs=0.925, na.rm=T), lty="dashed", col=2)

lines(temp,  apply(dR0.R0dbc,  1, FUN=median, na.rm=T), col=3)
lines(temp, apply(dR0.R0dbc, 1, FUN=quantile, probs=0.025, na.rm=T), lty="dashed", col=3)
lines(temp, apply(dR0.R0dbc, 1, FUN=quantile, probs=0.925, na.rm=T), lty="dashed", col=3)

lines(temp,  apply(dR0.R0dEFD,  1, FUN=median, na.rm=T), col=4)
lines(temp, apply(dR0.R0dEFD, 1, FUN=quantile, probs=0.025, na.rm=T), lty="dashed", col=4)
lines(temp, apply(dR0.R0dEFD, 1, FUN=quantile, probs=0.925, na.rm=T), lty="dashed", col=4)

lines(temp,  apply(dR0.R0de2a,  1, FUN=median, na.rm=T), col=5)
lines(temp, apply(dR0.R0de2a, 1, FUN=quantile, probs=0.025, na.rm=T), lty="dashed", col=5)
lines(temp, apply(dR0.R0de2a, 1, FUN=quantile, probs=0.925, na.rm=T), lty="dashed", col=5)

lines(temp,  apply(dR0.R0dMDR,  1, FUN=median, na.rm=T), col=6)
lines(temp, apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.025, na.rm=T), lty="dashed", col=6)
lines(temp, apply(dR0.R0dMDR, 1, FUN=quantile, probs=0.925, na.rm=T), lty="dashed", col=6)

lines(temp,  apply(dR0.R0dmu,  1, FUN=median, na.rm=T), col=7, lwd=2)
lines(temp, apply(dR0.R0dmu, 1, FUN=quantile, probs=0.025, na.rm=T), lty="dashed", col=7, lwd=2)
lines(temp, apply(dR0.R0dmu, 1, FUN=quantile, probs=0.925, na.rm=T), lty="dashed", col=7, lwd=2)

lines(temp,  apply(dR0.R0dPDR,  1, FUN=median, na.rm=T), col=8)
lines(temp, apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.025, na.rm=T), lty="dashed", col=8, lwd=2)
lines(temp, apply(dR0.R0dPDR, 1, FUN=quantile, probs=0.925, na.rm=T), lty="dashed", col=8, lwd=2)

