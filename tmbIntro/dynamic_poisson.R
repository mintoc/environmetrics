##---------------------------------------------
## Simulate and fit Poisson counts with a varying rate
## CM: Tue Apr 24 2018
## AJ: modified to include a comparison to a Bayesian fit using JAGS
##     JAGS must be installed separately to the rjags package in 
##     order for this code to run. see http://mcmc-jags.sourceforge.net
##---------------------------------------------

##------------
## SIMULATION
##------------

set.seed(1012)
## initial log rate for abundance, say
lnlambda <- 2

sigma <- 0.1

T <- 100

## simulate stochastic log rate 
for(i in 2:T){
    lnlambda[i] <- rnorm(1, mean = lnlambda[i-1], sd = sigma)
}

lambda <- exp(lnlambda)


## abundance
n <- rpois(T, lambda = lambda)

## set some to zero 
## n[30:40] <- NA
plot(n)
lines(lambda)
lines(lowess(lambda, f = 1/10), col = "grey", lty = 1)

## variable for whether n present or not
npres <- as.numeric(!is.na(n))

# set some data as missing values
npres[30:40] <- 0
n[npres==0] <- NA

##------------
## ESTIMATION
##------------
## load TMB
library(TMB)

## compile
compile("dynamic_poisson.cpp")

## load the function
dyn.load(dynlib("dynamic_poisson"))

obj <- MakeADFun(
    ## data
    data = list(
        y = n,
        ypres = npres
    ),
    ## starting parameter values
    parameters = 
        list(lnsigma = 0.1,
             lnlambda = rep(mean(2), length(n))),
    ## specify what is random
    random = "lnlambda",
    DLL = "dynamic_poisson",
    hessian = TRUE)

## run the optimisation
opt <- nlminb(objective = obj$fn, gradient = obj$gr, start = obj$par)

## calculate standard errors
rep <- sdreport(obj)
rep.summ <- summary(rep)

## extract lnlambda estimates
lnlambda.hat <- rep.summ[rownames(rep.summ) == "lnlambda", ]

## multiplier matrix for getting approximate 95% CIs
mult <- matrix(c(1, 0, 1, 2, 1, -2), nrow = 2)

lambda.hat <- exp(lnlambda.hat %*% mult)





# ------------------------------------------------------------------------------
# try a JAGS model. Quick and dirty with only one chain an no convergence 
# checking!

library(rjags)

# specify the jags model as a text string
modelstring ='
model {

# Specify the likelihood for the data, i.e. the 
# response variable whcih in this case is a set of  
# d13C values for 9 geese.
for(i in 1:N) { 
  y[i] ~ dpois(lambda[i])

  # the rate parameter lambda is the exponent of its log value
  lambda[i] <- exp(log_lambda[i])
}


# from time point 2 onwards, log(lambda) is an AR process of the 
# previous value with normal error defined by precision tau_lambda
for(i in 2:N){
  log_lambda[i] ~ dnorm(log_lambda[i-1], tau_lambda)
}

# prior on the first time point log_lambda[1]
log_lambda[1] ~ dnorm(0, 10^-6)

# prior on error around lambda is specified on standard deviation using
# a uniform distribution which is converted to precision tau
sigma_lambda ~ dunif(0, 100)
tau_lambda <- 1 / (sigma_lambda * sigma_lambda)

}
' # end of string

# specify the data to be passed to the jags model we defined
data = list(y = n,
            N = length(n)
            ) # end of data specification

# generate the model
model <- jags.model(textConnection(modelstring), data = data)

# generate the samples for the posterior and record only the lambda values
# in this example (we might also want to recover sigma_lambda but i havent
# included it here) 
output <- coda.samples(model = model,
                      variable.names = c("lambda"),
                      n.iter = 10000)

## Take a look 
plot(n, pch = 1, col = "blue", bty = "l", 
     xlab = "Time", ylab = "Abundance", 
     ylim = range(0, lambda, 1.2 * lambda.hat))
title("Dynamic Poisson time series")
## add in 95% prediction intervals
# lines(qpois(0.025, lambda = lambda.hat[,1]), lty = 2, col = "darkblue")
# lines(qpois(0.975, lambda = lambda.hat[,1]), lty = 2, col = "darkblue")
## approximate 95% CIs
polygon(c(1:T, rev(1:T)), c(lambda.hat[,2], rev(lambda.hat[,3])), 
        col = "#0000FF40", border = NA)
lines(lambda.hat[,1], col = "darkblue", lty = 2)
# lines(lambda, col = "red", lty = 2)

# add the lines for the mean and 95% credible interval
lines(colMeans(output[[1]]), type = "l")
lines(apply(output[[1]],2, function(x){quantile(x, 0.025)}))
lines(apply(output[[1]],2, function(x){quantile(x, 0.975)}))

legend("topright", lty = c(NA, 2, NA, 1), pch = c(1, NA, 15, NA),
       col = c("blue", "darkblue", "#0000FF40", "black"),
       legend = c("Count",  "lambda hat", "95% confidence interval",
                  "JAGS estimate"), 
       bty = "n")

# -----------------------------------------------------------------------------
# plot all this on a log scale
## Take a look 
plot(log(n), pch = 1, col = "blue", bty = "l", 
     xlab = "Time", ylab = "Log Abundance", 
     ylim = range(0, log(lambda), log(1.2 * lambda.hat)))
title("Dynamic Poisson time series")
## approximate 95% CIs
polygon(c(1:T, rev(1:T)), log(c(lambda.hat[,2], rev(lambda.hat[,3]))), 
        col = "#0000FF40", border = NA)
lines(log(lambda.hat[,1]), col = "darkblue", lty = 2)
# lines(log(lambda, col = "red", lty = 2))


# add the jags fit
lines(log(colMeans(output[[1]])), type = "l")
lines(log(apply(output[[1]],2, function(x){quantile(x, 0.025)})))
lines(log(apply(output[[1]],2, function(x){quantile(x, 0.975)})))

legend("topright", lty = c(NA, 2, NA, 1), pch = c(1, NA, 15, NA),
       col = c("blue", "darkblue", "#0000FF40", "black"),
       legend = c("Count",  "lambda hat", "95% confidence interval",
                  "JAGS estimate"), 
       bty = "n")
