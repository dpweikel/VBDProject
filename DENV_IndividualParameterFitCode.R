### This file will process the temperature data and create the individual parameter 
### samples for all the traits that are included in the DENV R0 model. 

## Remember to set your working directory so the code can access the data and 
## necessary supplementary code.
setwd("~/Desktop/Summer '15 /Models/DENV Model")

## Loading the required packages and supplementary code for sampling and analysis.
library(IDPmisc)
library('rjags')

# This file contains tools for analysis and visualization.
source("mcmc_utils.R") 

# This file contains the derivatives for the functions.
source("temp_deriv_functions.R") 

## Loading the data; my data will be called aegyptiDENVmodelTempData.csv 
data.all <- read.csv("aegyptiDENVmodelTempData.csv", header=TRUE)

## Now the code will choose all the temperature sensitive traits and fit either a
## Briere or Quadratic model to it using MCMC sampling.

## The first trait the biting rate, a.

# This chooses the trait, GCD, which is used to model the biting rate, a.

data <- data.all[ which(data.all$trait.name=="gcd"),]

# Transforming the GCD into a, the biting rate

data$trait = (1/(data$trait))

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

# Specifing the parameters that control the MCMC (these will be used throughout the code). 

n.chains <- 5
n.adapt <- 5000
n.samps <- 5000

# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
a.samps <- samps

## The second trait to analyze is EFD, or the number of eggs laid
## per female per day.

data <- data.all[which(data.all$trait.name=="EFD"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
EFD.samps <- samps

## The next trait we'll be choosing is the vector competence, b*c.
## Unfortunately, due to the data that was collected it will be necessary
## to decompose it into its two parts, b and c.

# Choose b, the probability a human will be bitten and infected by
# an infectious mosquito (ie. transmission).

data <- data.all[ which(data.all$trait.name=="b"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
b.samps <- samps

# Next, choose c, the probability that a mosquito becomes infected
# after biting an infectious human (ie. infection). 

data <- data.all[ which(data.all$trait.name=="c"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

# Given the data we've chosen to use the Briere fucntion. Jag-briere.bug contains the 
# specifics of the Briere model with the default priors, which is then used to create 
# an MCMC sample using the data.

jags <- jags.model('jags-briere.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'= length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c', 'Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that will be used for further analyses.

samps <- make.briere.samps(coda.samps, nchains = n.chains, samp.lims = c(1, n.samps))
samps$tau <- 1/samps$sigma
c.samps <- samps

## The next trait that is fitted is MDR, the mean development time for a mosquito. 

data <- data.all[ which(data.all$trait.name=="MDR"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

# Given the data the Briere fuction is chosen. Jag-briere.bug contains the specifics of 
# the Briere model with the default priors.

jags <- jags.model('jags-briere.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N' = length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt) 

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c','Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps <- make.briere.samps(coda.samps, nchains=n.chains, samp.lims=c(1, n.samps))
samps$tau <- 1/samps$sigma
MDR.samps <-  samps

## The next parameter is pEA, the probability a mosquito will survive from hatching to 
## maturation. 

data <- data.all[which(data.all$trait.name=="pEA"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait~T, data=data)

# Given the data the Negative Quadratic function is chosen. Jags-quad-neg.bug contains
# the specifics of the Negative Quadratic model with the default priors. 

jags <- jags.model('jags-quad-neg.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'=length(data$T)),
                   n.chains = n.chains,
                   inits=list(T0=5, Tm=33, n.qd=0.005), n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis

coda.samps <- coda.samples(jags, c('T0','Tm', 'qd'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 
plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps.q <- make.quad.samps(coda.samps, nchains=n.chains, 
                           samp.lims=c(1, n.samps), sig=FALSE)
samps.q$n.qd <- samps.q$qd
e2a.samps <- samps.q

## The next parameter is p, the daily probability that an adult mosquito
## survives. 

data.p <- data.all[which(data.all$trait.name=="p"),]

# Taking the mortality rate and changing it into the survival probability
data.days <- data.all[which(data.all$trait.name=="p/days"),]

data.days$trait = exp(-(data.days$trait)*(849/45))

# Merginging two sets together
data <- rbind(data.days, data.p)

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data.p)
plot(trait ~ T, data = data.days)
plot(trait ~ T, data = data)

# Given the data the Negative Quadratic function is chosen. Jags-quad-neg.bug contains
# the specifics of the Negative Quadratic model with the default priors. 

jags <- jags.model('jags-quad-neg.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N'=length(data$T)),
                   n.chains = n.chains,
                   inits=list(T0=5, Tm=33, n.qd=0.005), n.adapt = n.adapt)

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis

coda.samps <- coda.samples(jags, c('T0','Tm', 'qd'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 
plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps.q <- make.quad.samps(coda.samps, nchains=n.chains, 
                           samp.lims=c(1, n.samps), sig=FALSE)
samps.q$n.qd <- samps.q$qd
p.samps <- samps.q

## The final parameter to fit is PDR, the parasite/virus development rate.

data <- data.all[which(data.all$trait.name=="PDR"),]

# Plot the data to see which fucntion, Briere or Quadratic, is a more suitable fit for 
# the data.

plot(trait ~ T, data = data)

# Given the data the Briere fuction is chosen. Jag-briere.bug contains the specifics of 
# the Briere model with the default priors.

jags <- jags.model('jags-briere.bug',
                   data = list('Y' = data$trait, 'T' = data$T, 'N' = length(data$T)),
                   n.chains = n.chains, inits = list(Tm = 31, T0 = 5, c = 0.00007),
                   n.adapt = n.adapt) 

# The coda.samples() function takes n.samps new samples, and saves
# them in the coda format, which we use for visualization and
# analysis.

coda.samps <- coda.samples(jags, c('c','Tm', 'T0', 'sigma'), n.samps)

# These plots are useful to asses model convergence and general diagnosticl information. 

plot(coda.samps)

# This command combines the samples from the n.chains into a format
# that we can use for further analyses

samps <- make.briere.samps(coda.samps, nchains=n.chains, samp.lims=c(1, n.samps))
samps$tau <- 1/samps$sigma
PDR.samps <-  samps

## This code is just save the MCMC samples for further analysis in the R0 model.
save(a.samps, b.samps, c.samps, MDR.samps, EFD.samps, 
     e2a.samps, p.samps, PDR.samps,
     file = "DENV_ParameterFits.Rsave")