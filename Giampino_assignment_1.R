
library(jagsUI)
library(rjags)
library(dplyr)

setwd("E:\\PhD\\Bayesian I\\Project")

forest <- read.csv("3_forestfires.csv")
summary(forest) # summary of dataset

forest <- forest %>%
  select(-X, -Y) # not consider the geographic coordinates

par(mfrow=c(1,1))
plot(density(forest$area), xlim=c(0,150)) # plot of the area
plot(density(log(forest$area)), main="Logarithm of area") # area in log

sum(forest$area == exp(log(forest$area)))

forest$area[forest$area == 0] <- 0.01 # 0 is transformed in 0.01 to apply log
forest$area <- log(forest$area) # log transformation

plot(density(forest$area))

x <- model.matrix(~.-area, forest) # covariates
y <- forest$area
n <- nrow(forest)
p <- ncol(x)

jags.m <- jagsUI::jags(model.file="multivariate_model.bugs",
             inits=NULL,
             data = list('n' = n,
                         'P' = p,
                         'x' = x, 
                         'y' = y),
             n.chains = 2,
             n.iter = 5000,
             n.burnin = 2000,
             n.adapt = 1000,
             parameters.to.save="beta")

jags.m
par(mfrow=c(1,1))
plot(jags.m) # traceplots and densities

library(plotMCMC)
# other way to obtain similar plot as above
plotTrace(jags.m$samples[[1]])
plotAuto(jags.m$samples[[1]])
plotCumu(jags.m$samples[[1]])
plotDens(jags.m$samples[[1]])

# Diagnostic check:
library(coda)
jags.mcmc <- as.mcmc(jags.m$samples[1])
geweke.diag(jags.mcmc)
geweke.plot(jags.mcmc)


jags.m2 <- jags.model( file="multivariate_model.bugs",
                      inits=NULL,
                      data = list('n' = n,
                                  'P' = p,
                                  'x' = x,
                                  'y' = y),
                      n.chains = 2,
                      n.adapt=1000)

samps <- coda.samples( jags.m2, "beta", n.iter=4000 )
summary(samps)
summary(window(samps, start=1001))  # Burn in of 1000 Start at 1001

plot(samps)
gelman.plot(samps)
gelman.diag(samps)



