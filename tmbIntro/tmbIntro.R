## ----echo=FALSE----------------------------------------------------------
library(knitr)
opts_chunk$set(size="footnotesize")
## set the working directory
##setwd('../tex'); 

## ----echo=FALSE, results='asis', fig.height = 4, fig.width = 6, fig.cap = 'The function $y = 2 + 3x -x^2$ with the value of $x$ that maximises the function shown as a vertical line.', fig.pos = "h", fig.align = "center"----
curve(2 + 3 * x - x^2, from = 0, to = 3, ylab = "y", xlab = "x", 
      cex.axis = 1.5, cex.lab = 2, bty = "l")
abline(v = 3/2, lty = 2)

## ----warning = FALSE-----------------------------------------------------
fx <- function(x){
  y <- 2 + 3 * x -x^2
  ## return -y as optim automatically does 
  ## minimization
  return(-y)
}

nlminb(objective = fx, start = 0)


## ------------------------------------------------------------------------
dydx <- function(x){
  dy <- 3 - 2 * x
  return(-dy)
}

nlminb(objective = fx, gradient = dydx, start = 0)


## ------------------------------------------------------------------------
## generate some data - 10,000 observations
set.seed(101)
y <- rnorm(1e4, mean = 8)

## log-likelihood function
fmu <- function(mu){
  ll <- dnorm(y, mean = mu, sd = 2, log = TRUE)
  return(-sum(ll))
}


## ----fig.height = 5, fig.width = 6, fig.cap = "Log-likelihood of $\\mathbf{y}$ over $\\mu$ with the maximum likelihood value of $\\mu$ shown as a vertical grey line."----
fmu2 <- function(mu){
  ll <- dnorm(y, mean = mu, sd = 2, log = TRUE)
  return(sum(ll))
}

h <- Vectorize(fmu2)
curve(h, from = 5, to = 12, ylab = "Log likelihood", xlab = expression(mu), 
      cex.axis = 1.5, cex.lab = 1.5, bty = "l", mfrow = c(5, 10, 1, 1))
abline(v = mean(y), col = "grey")

## ------------------------------------------------------------------------
nlminb(objective = fmu, start = 5)

## ------------------------------------------------------------------------

dfmu <- function(mu){
  dll <- sum(2 * y / 8 - 2 * mu /8)
  return(-dll)
}

nlminb(objective = fmu, gradient = dfmu, start = 5)

## ------------------------------------------------------------------------
## load the TMB library
library(TMB)


## ----eval = FALSE, warning = FALSE---------------------------------------
## 
## ## step takes a while (on my machine anyway)
compile("unknown_mean.cpp")
## 

## ------------------------------------------------------------------------
## load the function
dyn.load(dynlib("unknown_mean"))

## ------------------------------------------------------------------------
obj <- MakeADFun(
         data = list(y = y), 
         parameters = list(mu = 5),
         DLL = "unknown_mean",
         silent = FALSE)

## compare the log-likelihoods
obj$fn(5); fmu(5)
## compare the gradients
obj$gr(5); dfmu(5)
## quite impressive!

## ------------------------------------------------------------------------

(opt <- nlminb(objective = obj$fn, gradient = obj$gr, start = obj$par))


## ------------------------------------------------------------------------

(rep <- sdreport(obj))


## ----eval = FALSE--------------------------------------------------------
compile("random_intercepts.cpp")
## 

## ------------------------------------------------------------------------

## load the function
dyn.load(dynlib("random_intercepts"))


## ------------------------------------------------------------------------

ngps <- 20

gp.means <- rnorm(ngps, mean = 8, sd = 2)

y <- c(sapply(gp.means, 
              FUN = function(x){rnorm(10, mean = x, sd = 0.5)})
       )
gps <- rep(0:(ngps - 1), each = 10)
library(ggplot2)

dat <- data.frame(gps, y)


## ------------------------------------------------------------------------

obj <- MakeADFun(
    data = list(y = y, gps = gps, ngps = ngps), 
    parameters = 
        list(mu = 5, logsigma = 0.1, logtau = 0.1, u = rep(mean(y), ngps)),
    random = "u",
    DLL = "random_intercepts",
    hessian = TRUE,
    silent = TRUE)

opt <- nlminb(objective = obj$fn, gradient = obj$gr, start = obj$par)

rep <- sdreport(obj)
rep.summ <- summary(rep)
gp.means <- rep.summ[rownames(rep.summ) == "gp_means",]
row.names(gp.means) <- NULL
pred.dat <- as.data.frame(cbind(gps = 0:(ngps-1), gp.means))
names(pred.dat)[3] <- "SE"
## compare to lme4 in R
library(lme4)
lme.fit <- lmer(y ~ 1 + (1|gps), REML = FALSE)

## TMB values
c(opt$par[1], exp(opt$par[2:3])); -opt$objective
## lme4 values
lme.fit


## ----fig.height = 4, fig.width = 6, fig.cap = "Random intercept data (grey points) and TMB-estimated mixed effect means and approximate 95\\% confidence intervals."----

theme_set(theme_bw())
ggplot(dat, aes(x = gps, y = y)) + geom_point(col = "slategrey", pch = 1) +
  xlab("Group") +
  geom_point(data = pred.dat, aes(x = gps, y = Estimate), col = "purple", size = 2) +
  geom_errorbar(
    data = pred.dat, aes(x = gps, y = Estimate, 
      ymin = Estimate - 2 * SE, ymax = Estimate + 2 * SE), 
    width = 0.4, col = "purple")
