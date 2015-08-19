### MU COMP SCRIPT ~ Go Back to the Beginning (Digimon Reference)

## Load working directory
setwd("~/GitHub/VBDProject")

## Load relevant packages and materials
library(IDPmisc)
library('rjags')
source("temp_functions.R") 
source("mcmc_utils.R")

## Load data for first model fit
data.all <- read.csv("aegyptiDENVmodelTempData.csv", header=TRUE)

## Choose the response variable
data.p <- data.all[which(data.all$trait.name=="p"),]
data.days <- data.all[which(data.all$trait.name=="p/days"),]

## Manipulating the survival data so we can examine it as mu 
ec <- 0.000001
data.days$trait = exp(-(data.days$trait)*(849/45))
data <- rbind(data.days, data.p)
data$trait = -log(data$trait + ec)


## Specify the MCMC Input Parameters
n.chains <- 5
n.adapt <- 5000
n.samps <- 5000


## Set up the jags code, which contains the specifics of the quadratic model with the
## default priors. As well as a few different jags models that contain different things.

jags.quad <- jags.model('jags-quad.bug',
                            data = list('Y' = data$trait,
                                        'T' = data$T, 'N'=length(data$T)),
                            n.chains = n.chains,
                            inits=list(inter=0.6, n.slope=0.01, qd=0.1, tau=100),
                            n.adapt = n.adapt)

#jags.mu <-  jags.model('jags-mu.bug',
#                            data = list('Y' = data$trait,
#                                        'T' = data$T, 'N'=length(data$T)),
#                           n.chains = n.chains,
#                            inits=list(inter=0.6, n.slope=0.01, qd=0.1, tau=100),
#                            n.adapt = n.adapt)

## The coda.samples() functions takes n.samps new samples, and saves
## them in the coda format, which we use for visualization and
## analysis
coda.samps1 <- coda.samples(jags.quad, c('inter','n.slope', 'qd', 'sigma'), n.samps)
#coda.samps2 <- coda.samples(jags.mu, c('inter','n.slope', 'qd', 'sigma'), n.samps)

## A simple, built in command to visualize samples using coda

plot(coda.samps1,ask=T)
#plot(coda.samps2,ask=T)

## This command combines the samples from the n.chains into a format
## that we can use for further analyses
samps1<-make.pos.quad.samps(coda.samps1, nchains=n.chains, samp.lims=c(1, n.samps))
mu.DENV.samps <- samps1
#samps2<-make.pos.quad.samps(coda.samps2, nchains=n.chains, samp.lims=c(1, n.samps))

## Next we want to visualize the posterior samples of the parameters
## specifying the temp response, and compare them to the priors. We
## can look at the pair-wise joint posterior distirbution:
ipairs(samps1[,c(1:3,4)], ztransf = function(x){x[x<1] <- 1; log2(x)})
#ipairs(samps2[,c(1:3,4)], ztransf = function(x){x[x<1] <- 1; log2(x)})

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
priors1<-list()
priors1$names<-c( "inter", "n.slope", "qd","tau")
priors1$fun<-c( "gamma", "gamma", "gamma","normal")
priors1$hyper<-matrix(NA, ncol=4, nrow=3)
priors1$hyper[,1]<-c(2, 2, NA)
priors1$hyper[,2]<-c(1, 1, NA)
priors1$hyper[,3]<-c(1, 100, NA)
priors1$hyper[,4]<-c(1000, 1/500, NA)

#priors2<-list()
#priors2$names<-c( "inter", "n.slope", "qd","tau")
#priors2$fun<-c( "normal", "normal", "gamma","normal")
#priors2$hyper<-matrix(NA, ncol=4, nrow=3)
#priors2$hyper[,1]<-c(2.3, 1/0.3, NA)
#priors2$hyper[,2]<-c(0.21, 1/0.02, NA)
#priors2$hyper[,3]<-c(2, 2, NA)
#priors2$hyper[,4]<-c(1000, 1/500, NA)

## the plot.hists command can be used to plot histograms of posterior
## samples of parameters with overlying priors
plot.hists(samps1[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)
#plot.hists(samps2[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors2)

## Next we want to use the parameter samples to get posterior samples
## of the temperature responses themselves
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="quad.pos", samps1, Temps, thinned=seq(1,25000, by=5))
#out2<-make.sims.temp.resp(sim="quad.pos", samps2, Temps, thinned=seq(1,25000, by=5))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))
#q2<-temp.sim.quants(out2$fits, length(Temps))

## We can then plot the data with the fits/quantiles
mycol<-1 
plot(data$T, data$trait, xlim=c(0,45), ylim=c(0,3),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab="Mu (Survival Rate) - Quad w/o Informed Priors",
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

#mycol<-1 
#plot(data$T, data$trait, xlim=c(0,45), ylim=c(0,3),
#     pch=(mycol+20),
#     xlab="Temperature (C)",
#     ylab="Mu (Survival Rate) - Quad w/ Priors",
#     col=mycol, cex=1.5)
#
#add.sim.lines(Temps, sim.data=out2$fits, q=q2, mycol=8)

## Second Set of Data which will be used
data.all <- read.csv("albopictusCHIKVmodelTempData.csv", header=TRUE)
data <- data.all[which(data.all$trait.name=="p.succ"),]

## Manipulating the data so it can be used in the analysis
data$trait = data$trait/data$trait2
ec <- 0.000001
data$trait = -log(data$trait + ec)


## Specify the MCMC Input Parameters
n.chains<-5
n.adapt<-5000
n.samps<-5000

## Set up the jags code, which contains the specifics of the quadratic model with the
## default priors. As well as a few different jags models that contain different things.

jags.quad <- jags.model('jags-quad-trunc.bug',
                        data = list('Y' = data$trait,
                                    'T' = data$T, 'N'=length(data$T)),
                        n.chains = n.chains,
                        inits=list(inter=0.6, n.slope=0.01, qd=0.1, tau=100),
                        n.adapt = n.adapt)

#jags.mu <-  jags.model('jags-mu.bug',
#                       data = list('Y' = data$trait,
#                                   'T' = data$T, 'N'=length(data$T)),
#                       n.chains = n.chains,
#                       inits=list(inter=0.6, n.slope=0.01, qd=0.1, tau=100),
#                       n.adapt = n.adapt)

## The coda.samples() functions takes n.samps new samples, and saves
## them in the coda format, which we use for visualization and
## analysis
coda.samps1 <- coda.samples(jags.quad, c('inter','n.slope', 'qd', 'sigma'), n.samps)
#coda.samps2 <- coda.samples(jags.mu, c('inter','n.slope', 'qd', 'sigma'), n.samps)

## A simple, built in command to visualize samples using coda

plot(coda.samps1,ask=T)
#plot(coda.samps2,ask=T)

## This command combines the samples from the n.chains into a format
## that we can use for further analyses
samps1<-make.pos.quad.samps(coda.samps1, nchains=n.chains, samp.lims=c(1, n.samps))
mu.CHIKV.samps <- samps1
#samps2<-make.pos.quad.samps(coda.samps2, nchains=n.chains, samp.lims=c(1, n.samps))

## Next we want to visualize the posterior samples of the parameters
## specifying the temp response, and compare them to the priors. We
## can look at the pair-wise joint posterior distirbution:
ipairs(samps1[,c(1:3,4)], ztransf = function(x){x[x<1] <- 1; log2(x)})
#ipairs(samps2[,c(1:3,4)], ztransf = function(x){x[x<1] <- 1; log2(x)})

## This is how the priors should be specified in order to plot them
## with the histograms of samples: Default priors
priors1<-list()
priors1$names<-c( "inter", "n.slope", "qd","tau")
priors1$fun<-c( "gamma", "gamma", "gamma","normal")
priors1$hyper<-matrix(NA, ncol=4, nrow=3)
priors1$hyper[,1]<-c(1, 10, NA)
priors1$hyper[,2]<-c(1, 10, NA)
priors1$hyper[,3]<-c(1, 100, NA)
priors1$hyper[,4]<-c(1000, 1/500, NA)

#priors2<-list()
#priors2$names<-c( "inter", "n.slope", "qd","tau")
#priors2$fun<-c( "normal", "normal", "gamma","normal")
#priors2$hyper<-matrix(NA, ncol=4, nrow=3)
#priors2$hyper[,1]<-c(2.3, 1/0.3, NA)
#priors2$hyper[,2]<-c(0.21, 1/0.02, NA)
#priors2$hyper[,3]<-c(2, 2, NA)
#priors2$hyper[,4]<-c(1000, 1/500, NA)

## the plot.hists command can be used to plot histograms of posterior
## samples of parameters with overlying priors
plot.hists(samps1[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors1)
#plot.hists(samps2[,c(1:3,4)], my.par=c(2,2), n.hists=4, priors=priors2)

## Next we want to use the parameter samples to get posterior samples
## of the temperature responses themselves
Temps<-seq(0,50, by=0.1)
out1<-make.sims.temp.resp(sim="quad.pos", samps1, Temps, thinned=seq(1,25000, by=5))
#out2<-make.sims.temp.resp(sim="quad.pos", samps2, Temps, thinned=seq(1,25000, by=5))

## and then we calculate the 95% inner quantile/HPD
q1<-temp.sim.quants(out1$fits, length(Temps))
#q2<-temp.sim.quants(out2$fits, length(Temps))

## We can then plot the data with the fits/quantiles
par(mfrow=c(1,1))
mycol<-1 
plot(data$T, data$trait, xlim=c(0,45), ylim=c(0,3),
     pch=(mycol+20),
     xlab="Temperature (C)",
     ylab="Mu (Survival Rate) - Quad w/o Informed Priors",
     col=mycol, cex=1.5)

add.sim.lines(Temps, sim.data=out1$fits, q=q1, mycol=8)

#mycol<-1 
#plot(data$T, data$trait, xlim=c(0,45), ylim=c(0,3),
#     pch=(mycol+20),
#     xlab="Temperature (C)",
#     ylab="Mu (Survival Rate) - Quad w/ Priors",
#     col=mycol, cex=1.5)

#add.sim.lines(Temps, sim.data=out2$fits, q=q2, mycol=8)

## Saving the mu samps for later use since the original attempts 
## don't really work out.
save(mu.DENV.samps, mu.CHIKV.samps,
     file = "MuFits.Rsave")
