
############################ GIAMPINO 790347 - FOREST FIRE ############################

rm(list=ls())

library(jagsUI)
library(rjags)
library(dplyr)
library(coda)
library(BAS)
library(dplyr)
library(tidyr)

setwd("E:\\PhD\\Bayesian I\\Project")

# Load data ----
forest <- read.csv("3_forestfires.csv")
summary(forest) # summary of dataset

# Coordinates
s <- data.frame(forest$X, forest$Y, row.names=row.names(forest))
plot(s,axes=FALSE,xlab="",ylab="",main="Monitor locations")

forest <- forest %>%
  select(-X, -Y) # not consider the geographic coordinates

# Plot area
par(mfrow=c(1,1))
plot(density(forest$area), xlim=c(0,150)) # plot of the area
plot(density(log(forest$area)), main="Logarithm of area") # area in log
forest$area <- log(forest$area+1) # log transformation
area <- forest$area

# Dummy
forest <- fastDummies::dummy_cols(forest, remove_first_dummy = TRUE)[-c(1,2)]

forest <- forest[-9]
forest <- cbind(forest, area)

# Multivariate outlier
library(car)
mod <- lm(area~., data=forest)
outlierTest(mod) # only at 5% not 1%

# Full model ----
x <- model.matrix(~.-area, forest) # covariates
y <- forest$area
n <- nrow(forest)
p <- ncol(x)

set.seed(2005)
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

jags.m # DMC + month_dec DIC = 1830
summary(jags.m)

# Diagnostic checks
par(mfrow=c(1,1))
plot(jags.m) # traceplots and densities

jags.mcmc <- as.mcmc(jags.m$samples[1])
geweke.diag(jags.mcmc)
geweke.plot(jags.mcmc)

heidel.diag(jags.mcmc)  # assess stationarity

# model with the library rjags using coda for diagnostic checks
set.seed(2020)
jags.m2 <- jags.model( file="multivariate_model.bugs",
                       inits=NULL,
                       data = list('n' = n,
                                   'P' = p,
                                   'x' = x,
                                   'y' = y),
                       n.chains = 2,
                       n.adapt=1000)

samps <- coda.samples( jags.m2, "beta", n.iter=4000 )
samps
summary(samps)
summary(window(samps, start=1001))  # Burn in of 1000 Start at 1001

plot(samps)
gelman.plot(samps)

# Model selection and model averaging ----

col_names <- c("intercept",names(forest)[-26])
variable_selection <- data.frame(names = col_names, 
                                 mean = jags.m$mean$beta, 
                                 overlap0 = jags.m$overlap0)

variable_selection # DMC, month_dec

# variable selection with DIC:

# !!! long cycle
variables <- c(1:25)

models <- list()
max_variables <- 25
n_model <- 100
count <- 1
for(k in 1:max_variables){
  for(j in 1:n_model){
    models[[count]] <- sample(variables, k)
    count <- count+1
  }
}

dic2 <- c()
dic.deviance2 <- c()
dic.penalty2 <- c()
dic.variable.length2 <- c()
count <- 1

for(model in models){
  cat("\n", count, "\n")
  count <- count + 1
  dati <- data.frame(forest[, model])
  x <- model.matrix(~.-area, data=dati)
  y <- forest$area
  n <- nrow(forest)
  p <- ncol(x)
  
  jags <- jags.model( file="multivariate_model.bugs",
                      inits=NULL,
                      data = list('n' = n,
                                  'P' = p,
                                  'x' = x,
                                  'y' = y),
                      n.chains = 2,
                      n.adapt=1000)
  
  dic.mod1 <- dic.samples(jags, 1500, "pD")
  
  
  dic2 <- c(dic2, sum(dic.mod1$deviance)+sum(dic.mod1$penalty))
  dic.deviance2 <- c(dic.deviance2, sum(dic.mod1$deviance))
  dic.penalty2 <- c(dic.penalty2, sum(dic.mod1$penalty))
  dic.variable.length2 <- c(dic.variable.length2, length(model))
  
}

dic.df2 <- data.frame(dic=dic2, dic.deviance=dic.deviance2, dic.penalty=dic.penalty2, dic.variable.length=dic.variable.length2)

write.csv(dic.df2, "dic.generation_random.csv")
dic.df2 <- read.csv("dic.generation_random.csv")[-1]

dic2.summary <- dic.df2 %>% 
  group_by(dic.variable.length) %>% 
  summarise(dic=min(dic), dic.deviance=min(dic.deviance), dic.penalty=min(dic.penalty))

plot(dic2.summary$dic, ylim = c(0, max(dic2.summary$dic)))

par(mfrow=c(1,1))
boxplot(dic~dic.variable.length,data=dic.df2, 
        main="Relationship between DIC and model complexity", 
        sub="100 random models for each number of variables",
        xlab="Number of variables", 
        ylab="DIC")

boxplot(dic.penalty~dic.variable.length,data=dic.df2, 
        main="Relationship between DIC penality and model complexity", 
        sub="100 random models for each number of variables",
        xlab="Number of variables", 
        ylab="DIC penality")


# DIC - stepwise ----

selected <- c()
results <- data.frame()
cols <- names(forest)
model <-  '
model{
    # prior
	c2 <- n
	# prior means
    for (j in 1:P){ mu.beta[j] <- 0.0 }
    # calculation of xtx
    for (i in 1:P){ for (j in 1:P){ 
       inverse.V[i,j] <- inprod( x[,i] , x[,j] )
    }}
    for(i in 1:P){ for (j in 1:P){
      prior.T[i,j] <- inverse.V[i,j] * tau /c2
    }}
    
    # likelihood
    
		for (i in 1:n){
		
		    y[i] ~ dnorm( mu[i], tau ) # stochastic componenent
			mu[i] <- inprod( beta[], x[i,] )
		}
		
	# prior distributions
    beta[1:P] ~ dmnorm( mu.beta[], prior.T[,] )
    tau    ~ dgamma( 0.01, 0.01 )  
    s2    <- 1/tau 
 	s <-sqrt(s2)
}
'


dic_function <- function(col){
  df <- forest %>% 
    select(area, col)
  x <- as.matrix(cbind(1, df[, -26]))
  y <- df$area
  n <- nrow(df)
  p <- ncol(x)
  model_data = list(n = n, y = y, x = x, P = p)
  
  model_parameters =  c("beta")
  
  jags <- jagsUI::jags(model.file=textConnection(model),
                       inits=NULL,
                       data = model_data,
                       n.chains = 2,
                       n.iter = 4000,
                       n.burnin = 500,
                       n.adapt = 1000,
                       parameters.to.save=model_parameters)
  
  print(jags)
  res <- list(col = col[length(col)], dic = jags$DIC, pd = jags$pD, 
            dev_mu = jags$mean$deviance, dev_sigma = jags$sd$deviance)
  return(res)
}

for(i in 1:26){
  
  for(col in cols[!(cols %in% selected)]){
    res <- dic_function(c(selected, col))
    res$i <- i
    results <- rbind(results, data.frame(res))
  }
  
  id_min_dic <- which.min(results[results$i == i, ]$dic)
  selected <- c(selected, as.character(results[results$i == i, ]$col[id_min_dic]))
  
}

selected

write.csv(results, "stepwise_results.csv")

boxplot(results$dic~results$i, main="Stepwise DIC", ylab="DIC", xlab="Number of variable")

# DIC - distributional assumption ----

# tau~Gamma(0.01, 0.01)

jags.m$DIC # 1830

# A) tau~Gamma(0.1, 0.1) DIC =  1830
# B) tau~Gamma( 1, 1)  DIC = 1830
# C) tau~Gamma(10, 10) DIC = 1830
# D) tau~Gamma(0.01, 1000) DIC = 2320

# A) unit-information empirical prior DIC = 1830.5 
# B) unit-information empirical prior DIC = 1830.7
# C) unit-information empirical prior DIC = 1830.3


# check structural components with DIC - Random effects model ----

burn   <- 1000
iters  <- 5000
chains <- 2 
data   <- list(y=y,x=x,n=n,P=p)

# Random effects:
set.seed(2020)
model.re  <- jagsUI::jags("multivariate_model_random_eff.bugs",
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

model.re # selected variables are the same as in the model with fixed effect

model.re$DIC # 11546

jags.m$DIC # 1830



# BAS Variable selection ----

subset <- c(1:25)
set.seed(2020)
# modelprior = beta.binomial(1,1) (default)
bas_bic <- bas.lm(area~., 
                  data=forest[, subset], 
                  prior = "BIC")

bas_gprior <- bas.lm(area~., 
                     prior="g-prior",
                     data=forest[, subset], 
                     alpha=n
)

bas_hyperg <- bas.lm(area~., 
                     prior="hyper-g",
                     data=forest[, subset], 
                     alpha=3
)

col_names <- c("intercept",names(forest)[-26])
variable_selection <- data.frame(names = col_names, 
                                 bic = round(bas_bic$probne0,3), 
                                 gprior = round(bas_gprior$probne0,3), 
                                 hyperg = round(bas_hyperg$probne0,3))

variable_selection

par(mfrow=c(1,2))
plot(bas_bic, which=4, ask=F, main="test")
image(bas_bic, rotate=F)
summary(bas_bic, 3)

# posterior inclusion > 0.2
bas_bic_selected <- which(bas_bic$probne0>0.2)
bas_gprior_selected <- which(bas_gprior$probne0>0.2)
bas_hyperg_selected <- which(bas_hyperg$probne0>0.5)
# intercept + december [+ DMC + temp + september] (hyper-g)

# modelprior = uniform
bas_bic_u <- bas.lm(area~., 
                  data=forest[, subset], 
                  modelprior = uniform(),
                  prior="BIC")

bas_gprior_u <- bas.lm(area~., 
                     data=forest[, subset], 
                     prior = "g-prior",
                     modelprior = uniform())

bas_hyperg_u <- bas.lm(area~., 
                     data=forest[, subset],
                     prior="hyper-g",
                     alpha=3,
                     modelprior = uniform())


col_names <- c("intercept",names(forest)[-26])
variable_selection_u <- data.frame(names = col_names, 
                                 bic = round(bas_bic_u$probne0,3), 
                                 gprior = round(bas_gprior_u$probne0,3),
                                 hyperg = round(bas_hyperg_u$probne0,3))

variable_selection_u

par(mfrow=c(1,2))
plot(bas_bic_u, which=4, ask=F, main="test")
image(bas_bic_u, rotate=F)
summary(bas_bic_u, 3)

# posterior inclusion > 0.2
bas_bic_selected_u <- which(bas_bic_u$probne0>0.2) 
# temp + month_dec + month_sep
bas_gprior_selected_u <- which(bas_gprior_u$probne0>0.2) 
# temp + month_dec + month_sep
bas_hyperg_selected_u <- which(bas_hyperg_u$probne0>0.5) 
# temp + month_dec + rain + month_sep

# HPM = highest probability model (or map)
# MPM = median probability model of Barbieri & Berger

col_names[which(coef(bas_bic, estimator = "HPM")$probne0>0.2)] #intercept+month_dec
col_names[which(coef(bas_bic, estimator = "MPM")$probne0>0.2)] #intercept+month_dec
col_names[which(coef(bas_bic_u, estimator = "HPM")$probne0>0.2)] #intecept+month_dec+month_sep+temp
col_names[which(fitted(bas_bic_u, estimator = "MPM")$probne0>0.2)] #none

col_names[which(coef(bas_gprior, estimator = "HPM")$probne0>0.2)] #intercept+month_dec
col_names[which(coef(bas_gprior, estimator = "MPM")$probne0>0.2)] #intercept+month_dec
col_names[which(coef(bas_gprior_u, estimator = "HPM")$probne0>0.2)] #intecept+month_dec+month_sep+temp
col_names[which(fitted(bas_gprior_u, estimator = "MPM")$probne0>0.2)] #none


col_names[which(coef(bas_hyperg, estimator = "HPM")$probne0>0.5)] 
#intercept+month_dec
col_names[which(coef(bas_hyperg, estimator = "MPM")$probne0>0.5)] #none
col_names[which(coef(bas_hyperg_u, estimator = "HPM")$probne0>0.5)] 
#intercept+DMC+temp+wind+month_dec+month_sep
col_names[which(coef(bas_hyperg_u, estimator = "MPM")$probne0>0.5)] #none


# BUGS variable selection ----

x <- model.matrix(~.-area, forest) # covariates
y <- forest$area
n <- nrow(forest)
p <- ncol(x)
beta.full <- jags.m$mean$beta
sd.full.beta <- jags.m$sd$beta
sd.full <- jags.m$mean$deviance/n


# Hyper-g - uniform
set.seed(2020)
jags.hyper.u <- jagsUI::jags(model.file="mod_hyperg_unif.bugs",
                       inits=NULL,
                       data = list('n' = n,
                                   'P' = p,
                                   'x' = x, 
                                   'y' = y,
                                   'prop.mean.beta'=rep(0, 26),
                                   'prop.sd.beta'=rep(1,26)),
                       n.chains = 2,
                       n.iter = 5000,
                       n.burnin = 2000,
                       n.adapt = 1000,
                       parameters.to.save="beta")

jags.hyper.u # month_dec, DIC = 1818

# Hyper-g - beta.binomial
set.seed(2020)
jags.hyper.bb <- jagsUI::jags(model.file="mod_hyperg_bb.bugs",
                             inits=NULL,
                             data = list('n' = n,
                                         'P' = p,
                                         'x' = x, 
                                         'y' = y,
                                         'prop.mean.beta'=rep(0, 26),
                                         'prop.sd.beta'=rep(1, 26)),
                             n.chains = 2,
                             n.iter = 5000,
                             n.burnin = 2000,
                             n.adapt = 1000,
                             parameters.to.save="beta")

jags.hyper.bb # DIC = 1823

# Unif-information empirical - uniform

beta.full <- jags.m$mean$beta
sd.full.beta <- jags.m$sd$beta
sd.full <- jags.m$mean$deviance/n

set.seed(2020)
jags.empir.u <- jagsUI::jags(model.file="mod_emp_unif.bugs",
                             inits=NULL,
                             data = list('n' = n,
                                         'P' = p,
                                         'x' = x, 
                                         'y' = y,
                                         'beta.full'=beta.full,
                                         'sd.beta.full'=sd.full.beta,
                                         'sd.full'=sd.full),
                             n.chains = 2,
                             n.iter = 5000,
                             n.burnin = 2000,
                             n.adapt = 1000,
                             parameters.to.save="beta")

jags.empir.u # DIC = 1890
col_names[which(jags.empir.u$overlap0$beta==FALSE)]
# "DMC"       "DC"        "temp"      "wind"      "month_dec" "day_sat" 

# Unif-information empirical - beta.binomial
set.seed(2020)
jags.empir.bb <- jagsUI::jags(model.file="mod_emp_bb.bugs",
                             inits=NULL,
                             data = list('n' = n,
                                         'P' = p,
                                         'x' = x, 
                                         'y' = y,
                                         'beta.full'=beta.full,
                                         'sd.beta.full'=sd.full.beta,
                                         'sd.full'=sd.full),
                             n.chains = 2,
                             n.iter = 5000,
                             n.burnin = 2000,
                             n.adapt = 1000,
                             parameters.to.save="beta")

jags.empir.bb # DIC = 1888
col_names[which(jags.empir.bb$overlap0$beta==FALSE)]
# "DMC"   "DC"   "temp" "wind" "month_dec" "month_sep" "day_sat"

# g-prior - uniform
set.seed(2020)
jags.g.u <- jagsUI::jags(model.file="mod_g_unif.bugs",
                             inits=NULL,
                             data = list('n' = n,
                                         'P' = p,
                                         'x' = x, 
                                         'y' = y,
                                         'prop.mean.beta'=rep(0,26),
                                         'prop.sd.beta'=rep(1,26)),
                             n.chains = 2,
                             n.iter = 5000,
                             n.burnin = 2000,
                             n.adapt = 1000,
                             parameters.to.save=c("beta"))

jags.g.u # DIC = 2072
col_names[which(jags.g.u$overlap0$gb==FALSE)]
# "DMC"       "month_dec"

# g-prior - beta.binomial
set.seed(2020)
jags.g.bb <- jagsUI::jags(model.file="mod_g_bb.bugs",
                         inits=NULL,
                         data = list('n' = n,
                                     'P' = p,
                                     'x' = x, 
                                     'y' = y,
                                     'prop.mean.beta'=beta.full,
                                     'prop.sd.beta'=sd.full.beta),
                         n.chains = 2,
                         n.iter = 5000,
                         n.burnin = 2000,
                         n.adapt = 1000,
                         parameters.to.save=c("beta", "g"))

jags.g.bb # DIC = 2069
col_names[which(jags.g.bb$overlap0$gb==FALSE)]
# "DMC" "month_dec"


# MAP, MPM:

unitunif.df <- as.data.frame(as.matrix(jags.g.bb$samples[1][,27:52]))
m.index <- sapply(1:26, function(p) 2^p)
model.id <- sapply(1:3000, function(i) 1 + as.matrix(unitunif.df)[i,] %*% m.index)
MAP <- as.data.frame(sort(table(model.id)/3000, decreasing = TRUE)); MAP[1,] # intercept
MPM <- MAP[which(MAP$Freq>0.5),]; MPM # intercept


# Other variable selection: BMA ----

library(BMA)
set.seed(2020)
bma <- bic.glm(area~., data=forest, glm.family="gaussian")
summary(bma) #intercept+month_dec (+month_sep+temp+DMC)


# Final model ----

# Intercept + month_dec + temp + month_sep + DMC

# All possible combination of these variables:
combination <- list()

for(i in 1:length(subsetf)){
  
  combination[[i]] <- (combn(subsetf, i))
  
}

combination
DIC <- c()

for(j in 1:length(combination)){
  
  sub <- as.data.frame(combination[j])
  
  K <- ncol(sub)
  
  for(k in K){
    subsetk <- sub[,k]
    dati <- as.data.frame(forest[,c(subsetk, 26)])
    x <- model.matrix(~.-area, dati) # covariates
    y <- forest$area
    n <- nrow(dati)
    p <- ncol(x)
    
    set.seed(2005)
    jags.f <- jagsUI::jags(model.file="multivariate_model.bugs",
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
    DIC[j] <- jags.f$DIC
    print(jags.f)
    print(DIC)
  }
  
}

# Final model:

subsetf <- c(5,10,19) # temp+month_dec+month_sep
dati <- forest[,subsetf]
x <- model.matrix(~.-area, dati) # covariates
y <- forest$area
n <- nrow(forest)
p <- ncol(x)

set.seed(2005)
jags.f <- jagsUI::jags(model.file="multivariate_model.bugs",
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

jags.f # DMC + month_dec DIC = 1803
summary(jags.f)

# Diagnostic checks
par(mfrow=c(1,1))
plot(jags.f) # traceplots and densities

jags.mcmc <- as.mcmc(jags.f$samples[1])
geweke.diag(jags.mcmc)

geweke.plot(jags.mcmc)

heidel.diag(jags.mcmc)  # assess stationarity

acf(jags.mcmc)[1]

# BAS
final_bas <- bas.lm(area~., dati)
summary(final_bas)
coef(final_bas, estimator = "MPM")
coef(final_bas, estimator = "HPM")
coef(final_bas, estimator = "BMA")

final_bma <- bic.glm(area~., data=dati, glm.family="gaussian")
summary(final_bma)

plot(jags.f$mean$beta, pch=19, ylab=expression(beta), xlab="Coefficients", ylim=c(0.0, 3), xaxt='n')
arrows(1:4, jags.f$q2.5$beta, 1:4, jags.f$q97.5$beta, angle = 90, length = 0.03)
text(c(1.07,2,3,3.9), c(1.2, 0.5, 0.5, 1),
     label=c("Intercept", "Temperature", "December", "September"))

x11()
plot(final_bas)


# Random effects:
set.seed(2020)
model.re  <- jagsUI::jags("multivariate_model_random_eff.bugs",
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

model.re # selected variables are the same as in the model with fixed effect


# Prediction for evaluation:

set.seed(2020)
sample <- sample(1:nrow(forest), 100)

forest_train <- forest[-sample,]
forest_test <- forest[sample,]
dati_test <- data.frame(forest_test[,c(10, 26)])
x_test <- model.matrix(~.-area, data=dati_test)
n_test <- nrow(x_test)
dati_train <- data.frame(forest_train[,c(subset, 26)])
x <- model.matrix(~.-area, data=dati_train) # covariates
y <- dati_train$area
n <- nrow(forest_train)
p <- ncol(x)

jags <- jagsUI::jags(model.file="multivariate_model.bugs",
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

jags
samples <- jags$samples[[1]] #save the sample from first chain
beta <- apply(samples, 2, mean)[-length(beta)]
n_test <- nrow(forest_test)
jags <- jags.model(file="multivariate_model_pred.bugs",
                   inits=NULL,
                   data = list('n_test' = n_test,
                               'x_test' = x_test, 
                               'beta' = beta),
                   n.chains = 5,
                   n.adapt = 1000)
update(jags, 1000)
samples <- jags.samples(jags, c('mu'), 4000)

posterior_means <- apply(samples$mu, 1, mean)
posterior_median <- apply(samples$mu, 1, median)

pred <- posterior_means
index <- which(pred>1000)
real <- forest_test$area

MSE_1 <- sqrt(mean((real[-index]-pred[-index])^2)) 
MSE_2 <- sqrt(mean((real[index]-pred[index])^2)) 

