
rm(list=ls())

library(jagsUI)
library(rjags)
library(dplyr)
library(coda)
library(BAS)

setwd("E:\\PhD\\Bayesian I\\Project")

forest <- read.csv("3_forestfires.csv")
summary(forest) # summary of dataset

forest <- forest %>%
  select(-X, -Y) # not consider the geographic coordinates
par(mfrow=c(1,1))
plot(density(forest$area), xlim=c(0,150)) # plot of the area
plot(density(log(forest$area)), main="Logarithm of area") # area in log

forest$area[forest$area == 0] <- 0.01 # 0 is transformed in 0.01 to apply log
forest$area <- log(forest$area) # log transformation
area <- forest$area
forest <- fastDummies::dummy_cols(forest, remove_first_dummy = TRUE)[-c(1,2)]

forest <- forest[-9]
forest <- cbind(forest, area)

n <- nrow(forest)

# mod <- lm(area~., data=forest)
# summary(mod)


# BAS Variable selection ----

subset <- c(1:25)
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

bas_hypergn <- bas.lm(area~., 
                      prior="hyper-g-n",
                      data=forest[, subset], 
                      alpha=3
)

col_names <- c("intercept",names(forest)[-26])
variable_selection <- data.frame(names = col_names, 
                                 bic = bas_bic$probne0, 
                                 gprior = bas_gprior$probne0, 
                                 hyperg = bas_hyperg$probne0,
                                 hypergn = bas_hypergn$probne0)

variable_selection

par(mfrow=c(1,2))
plot(bas_bic, which=4, ask=F, main="test")
image(bas_bic, rotate=F)
summary(bas_bic, 3)

# posterior inclusion > 0.2
bas_bic_selected <- which(bas_bic$probne0>0.2)
bas_gprior_selected <- which(bas_gprior$probne0>0.2)
bas_hyperg_selected <- which(bas_hyperg$probne0>0.2)
bas_hypergn_selected <- which(bas_hypergn$probne0>0.2)
# intercept + december [+ temp + march + september] (hyper-g)

# jags variable selection ----

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

beta_selection <- which(jags.m$overlap0$beta==FALSE & round(jags.m$mean$beta, 2) != 0)

beta_selection # 11, december

beta <- data.frame(name=col_names[beta_selection], 
                   mean=jags.m$mean$beta[beta_selection])

beta

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



# Predictions y ----

set.seed(2020)
sample <- sample(1:nrow(forest), 100)

forest_train <- forest[-sample,]
forest_test <- forest[sample,]

predictor <- function(subset, graph_id=c(1:4)){
  subset <- subset[-1]-1 # remove the intercept and take the previous value
  dati_test <- data.frame(forest_test[,c(subset, 26)])
  x_test <- model.matrix(~.-area, data=dati_test)
  n_test <- nrow(x_test)
  dati_train <- data.frame(forest_train[,c(subset, 26)])
  x <- model.matrix(~.-area, data=dati_train) # covariates
  y <- dati_train$area
  n <- nrow(forest_train)
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
  samples <- jags.m$samples[[1]] #save the sample from first chain
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
  real <- forest_test$area
  
  MEAE <- median(abs((real-pred)/real)) #median absolute error
  par(mfrow=c(1,1))
  plot(1:100, abs((real-pred)/real), main="Absolute percentage deviation", 
       xlab="Observations", ylab="Absolute per error")
  abline(h=mean(abs((real-pred)/real)), lwd=2)
  abline(h=median(abs((real-pred)/real)), lty=2, lwd=2)
  legend("topright", 
         c("Mean", "Median"),
         lty=c(1,2), 
         lwd=c(2, 2))
  
  return(list(MEAE=MEAE))
}

predictor_lm <- function(subset){
  subset <- subset[-1]-1 # remove the intercept and take the previous value
  dati_test <- data.frame(forest_test[,c(subset, 26)])
  dati_train <- data.frame(forest_train[,c(subset, 26)])
  lm <- lm(area~., data=dati_train)
  fitted <- predict(lm, dati_test)
  fitted <- exp(fitted)
  real <- dati_test$area
  MEAE <- median(abs((real-fitted)/real))#median absolute percentage deviation
  return(list(MEAE=MEAE))
}

predictor_lm(bas_bic_selected)
predictor(bas_bic_selected)

predictor_lm(bas_gprior_selected)
predictor(bas_gprior_selected)

predictor_lm(bas_hyperg_selected)
predictor(bas_hyperg_selected)

predictor_lm(bas_hypergn_selected)
predictor(bas_hypergn_selected)


# Prediction x ----

# I insert some NA in the variable temperature (temp) and use the subset of 
# covariate of hyper-g-n prior.

subset <- bas_bic_selected
subset <- subset[-1]-1 # remove the intercept and take the previous value
dati <- data.frame(forest[,c(subset, 26)])
## Delete missing data (3 in each)
set.seed(2020)
temp.miss = sample(n)[1:3]
dati$month_dec[temp.miss]= NA
x <- as.matrix(cbind(1, dati[,-ncol(dati)])) # covariates
x[temp.miss,]
y <- dati$area
n <- nrow(forest)
p <- ncol(x)

mu.x <- sum(forest$month_dec)/517



jags.mx <- jagsUI::jags(model.file="mtv_mod_x_miss.bugs",
                       inits=NULL,
                       data = list('n' = n,
                                   'P' = p,
                                   'x' = x,
                                   'mu.x' = mu.x,
                                   #'tau.x' = tau.x,
                                   'y' = y),
                       n.chains = 2,
                       n.iter = 5000,
                       n.burnin = 2000,
                       n.adapt = 1000,
                       parameters.to.save=c("beta", "x"))
jags.mx
temp.miss # 412 236 321
round(jags.mx$mean$x[412,2],2) # 0.02
round(jags.mx$mean$x[236,2],2) # 0.08
round(jags.mx$mean$x[321,2],2) # 0.05
forest$month_dec[temp.miss] # all 0

# model without missing
subset <- bas_bic_selected
subset <- subset[-1]-1 # remove the intercept and take the previous value
dati <- data.frame(forest[,c(subset, 26)])
x <- model.matrix(~.-area, dati) # covariates
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





