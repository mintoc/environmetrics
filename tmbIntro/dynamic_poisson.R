##---------------------------------------------
## Simulate and fit Poisson counts with a varying rate
## CM: Tue Apr 24 2018
##
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

plot(lambda)

## abundance
n <- rpois(T, lambda = lambda)

## set some to zero 
## n[30:40] <- NA
plot(n)

## variable for whether n present or not
npres <- as.numeric(!is.na(n))

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

## Take a look 
plot(n, pch = 1, col = "blue", bty = "l", xlab = "Time", ylab = "Abundance", ylim = range(0, lambda, 1.2 * lambda.hat))
title("Dynamic Poisson time series")
## add in 95% prediction intervals
lines(qpois(0.025, lambda = lambda.hat[,1]), lty = 2, col = "darkblue")
lines(qpois(0.975, lambda = lambda.hat[,1]), lty = 2, col = "darkblue")
## approximate 95% CIs
polygon(c(1:T, rev(1:T)), c(lambda.hat[,2], rev(lambda.hat[,3])), col = "#0000FF40", border = NA)
lines(lambda.hat[,1], col = "darkblue")
legend("topright", lty = c(NA, 1, 2, NA), pch = c(1, NA, NA, 15),
       col = c("blue", "darkblue", "darkblue", "#0000FF40"),
       legend = c("Count", "Mean rate", "95% prediction interval", "95% confidence interval"))

