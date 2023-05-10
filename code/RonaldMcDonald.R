##############################################################
# Alden Sears
# Community assembly under variable resource interactions
# Rosenzweig-MacArthur version
# R version 4.2.1
##############################################################

## Starting stuff

library(tidyverse)    # design and data structure
library(dplyr)        # data manipulation
library(brms)         # bayesian modeling interfacing w/ Stan
library(tidybayes)    # data manip and management in bayes
library(cmdstanr)     # R interface to CmdStan
library(posterior)    # tools for posterior dists
library(rstan)        # R interface to Stan

##############################################################
# Building a 1-spec R-McA model
##############################################################

## skipping checks for now, just need prax

Ronald <- function(N, P, parameters, nsteps, resp, proc_error, obs_error){
  
  RMcaRes <- tibble(Time = c(1:nsteps),    # number of time steps/generations
                    N = c(1:nsteps),       # starting prey density
                    P = c(1:nsteps),       # starting predator density
                    zN = c(1:nsteps),      # latent prey density
                    zP = c(1:nsteps),      # latent predator density
                    dN = c(1:nsteps),      # change in prey density
                    dP = c(1:nsteps))      # change in predator density
  
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
  
  RMcaRes$P[1] <- P
  RMcaRes$N[1] <- N
  RMcaRes$dP[1] <- 0
  RMcaRes$dN[1] <- 0
  RMcaRes$zP[1] <- P
  RMcaRes$zN[1] <- N
  zP <- P
  zN <- N
  
  for(i in 2:nsteps){
    P <- zP
    N <- zN
    zP <- rlnorm(1, log((e * eval(f) - u) * P + P), sdlog = proc["P"])
    zN <- rlnorm(1, log(r * N * (1 - N/K) - eval(f) * P + N), sdlog = proc["N"])
    RMcaRes$zP[i] <- zP
    RMcaRes$zN[i] <- zN
    RMcaRes$P[i] <- rnorm(1, mean = zP, sd = obs["P"] * zP)
    RMcaRes$N[i] <- rnorm(1, mean = zN, sd = obs["N"] * zN)
    RMcaRes$dP[i] <- RMcaRes$P[i] - RMcaRes$P[i-1]
    RMcaRes$dN[i] <- RMcaRes$N[i] - RMcaRes$N[i-1]
  }
  
  return(RMcaRes)
}

## Testing it out

set.seed(666)
pars <- c("e" = 0.3, "a" = 0.01, "h" = 0.7, "u" = 0.28, "r" = 0.5, "K" = 500)
proc <- c("P" = 0.02, "N" = 0.05)
obs <- c("P" = 0.05, "N" = 0.1)

RonDat <- Ronald(50, 10, parameters = pars, nsteps = 1000, resp = "HII", proc_error = proc, obs_error = obs)

plot(N ~ Time, data = RonDat, type = "l", ylim = c(0,max(N)), lwd = 1)
lines(RonDat$P, col = "red", lwd = 1)

##############################################################
# Building a multipredator R-McA model
##############################################################

McDonald <- function(N,               # starting density of prey
                     Ps,              # starting densities of predators
                     params,          # "general" parameters (a, r, K, u)
                     Pparams,         # predator-specific parameters (e, h)
                     nsteps,          # number of time steps/generations
                     resp,            # functional response
                     proc_error,      # process error sds (lognormal)
                     obs_error){      # observation error sds (normal)
  
  res <- tibble(Time = c(1:nsteps), N = c(1:nsteps), zN = c(1:nsteps), dN = c(1:nsteps))
  for (i in 1:length(Ps)){
    in_between <- paste("P", i, sep = "")
    zin_between <- paste("zP", i, sep = "")
    din_between <- paste("dP", i, sep = "")
    eval(res <- add_column(res, in_between = c(1:nsteps), zin_between = c(1:nsteps),
                           din_between = c(1:nsteps)))
    colnames(res)[(2+3*i):(4+3*i)] <- c(in_between, zin_between, din_between)
  }
  
}
























