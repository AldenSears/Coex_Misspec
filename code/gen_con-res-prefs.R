## Simulate community of consumer resource preferences
library(tidyverse)
library(brms)
library(mvtnorm)
library(utilities)
library(ggtern)

# set number of consumers & resources
nC <- 10000
nR <- 3

# set community preference concentration
comm_theta <- 0.2

set.seed(2392)
# draw consumer preferences
comm_param <- rdirichlet(n = nC, alpha = rep(1 * comm_theta, nR))


# backtransform to percent preferences
# rows are consumers, columns are resources
comm_perc <- comm_param %>%
  as_tibble() %>%
  rename_with(.cols = everything(), ~str_replace(., "V", "R")) %>%
  rownames_to_column(var = "C")

# plot (works for nR = 3)
ggtern(comm_perc, aes(x = R1, y = R2, z = R3))+
  geom_point(size = 2, alpha = 0.1)
