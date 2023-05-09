##############################################################
# Alden Sears
# Community assembly under variable resource interactions
# Rosenzweig-MacArthur version
# R version 4.2.1
##############################################################

# Starting stuff

library(tidyverse)    # design and data structure
library(dplyr)        # data manipulation
library(brms)         # bayesian modeling interfacing w/ Stan
library(tidybayes)    # data manip and management in bayes
library(cmdstanr)     # R interface to CmdStan
library(posterior)    # tools for posterior dists
library(rstan)        # R interface to Stan

# Building a general R-McA model

## PUT IN CHECKS LATER ##

RMcA <- function(zN, zP, parameters, nsteps, resp, proc_error, obs_error){
  RMcaRes <- tibble(N = c(1:nsteps),
                    P = c(1:nsteps),
                    dN = c(1:nsteps),
                    dP = c(1:nsteps),
                    Time = c(1:nsteps))
  e <- parameters["e"]
  a <- parameters["a"]
  h <- parameters["h"]
  u <- parameters["u"]
  r <- parameters["r"]
  K <- parameters["K"]
  if(resp == "HII"){
    f <- expression((a * N)/(1 + a * h * N))
  } else if(resp == "I"){
    f <- expression(1/h * (1 - exp(-a * h * N)))
  } else if(resp == "Ht"){
    f <- expression(1/h * tan(a * h * N))
  }
  RMcaRes$P[1] <- rnorm(1, mean = zP, sd = obs["P"] * zP)
  RMcaRes$N[1] <- rnorm(1, mean = zN, sd = obs["N"] * zN)
  RMcaRes$dP[1] <- 0
  RMcaRes$dN[1] <- 0
  for(i in 2:nsteps){
    P <- zP
    N <- zN
    zP <- rlnorm(1, log((e * eval(f) - u) * P + P), sdlog = proc["P"])
    zN <- rlnorm(1, log(r * N * (1 - N/K) - eval(f) * P + N), sdlog = proc["N"])
    RMcaRes$P[i] <- rnorm(1, mean = zP, sd = obs["P"] * zP)
    RMcaRes$N[i] <- rnorm(1, mean = zN, sd = obs["N"] * zN)
    RMcaRes$dP[i] <- RMcaRes$P[i] - RMcaRes$P[i-1]
    RMcaRes$dN[i] <- RMcaRes$N[i] - RMcaRes$N[i-1]
  }
  return(RMcaRes)
}

# Testing it out

set.seed(666)
pars <- c("e" = 0.3, "a" = 0.01, "h" = 0.7, "u" = 0.28, "r" = 0.5, "K" = 500)
proc <- c("P" = 0.02, "N" = 0.05)
obs <- c("P" = 0.05, "N" = 0.1)

RMcADat <- RMcA(50, 10, parameters = pars, nsteps = 1000, resp = "HII", proc_error = proc, obs_error = obs)

plot(N~Time,dat=RMcADat,type="l",ylim=c(0,max(N)),lwd=1)
lines(RMcADat$P,col="red",lwd=1)














