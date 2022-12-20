##################################################-
## Fitting competition models in R ----
## W.K. Petry
##################################################-
## Preliminaries ----
##################################################-
library(tidyverse)
library(nlstools)
library(brms)
library(tidybayes)
library(modelr)

theme_set(theme_bw(base_size = 18)+
            theme(panel.grid = element_blank(),
                  plot.margin = margin(10, 25, 10, 10),
                  axis.text = element_text(color = "black"),
                  axis.ticks = element_line(size = 0.5)))

##################################################-
## Specify model form ----
##################################################-
# Beverton-Holt model
bh_form <- as.formula(y ~ lami / (1 + aii * Ni + aij * Nj))
bh_norm <- function(Ni, Nj, lami, aii, aij, sd){
  out <- rnorm(n = 1, mean = lami / (1 + aii * Ni + aij * Nj), sd = sd)
  return(out)
}

##################################################-
## Simulate data ----
##################################################-
# set parameters
lami <- 2000
aii <- 2
aij <- 1.2
sd <- 5
n <- 50

set.seed(23423)
bh_dat <- tibble(Ni = sample(0:100, size = n, replace = TRUE),
                 Nj = sample(0:100, size = n, replace = TRUE)) %>%
  mutate(y = map2_dbl(.x = Ni, .y = Nj,
                      .f = ~bh_norm(Ni = .x, Nj = .y, lami = lami, aii = aii, aij = aij, sd = sd)))
bh_dat

plot(Nj ~ Ni, data = bh_dat)

##################################################-
## Fit with nls (frequentist) ----
##################################################-
bh_nls <- nls(bh_form, data = bh_dat, start = list(lami = 500, aii = 0.5, aij = 0.5))
summary(bh_nls)

plot(profile(bh_nls))
overview(bh_nls)
summary(nlsJack(bh_nls))

##################################################-
## Fit with brms ----
##################################################-
# check the priors you need
bh_form_brms <- bf(bh_form,
                   lami ~ 1,
                   aii ~ 1,
                   aij ~ 1,
                   nl = TRUE)

get_prior(bh_form_brms, family = gaussian(), data = bh_dat)

hist(abs(rnorm(1e6, 500, 10000)), breaks = 1000)

# set priors
bh_priors <- prior(normal(500, 10000), lb = 0, nlpar = "lami")+
  prior(uniform(0, 3), lb = 0, ub = 3, nlpar = "aii")+
  prior(uniform(0, 3), lb = 0, ub = 3, nlpar = "aij")

# prior predictive check
bh_brm_prior <- brm(bh_form_brms,
                    family = gaussian(), data = bh_dat,
                    sample_prior = "only",
                    backend = "cmdstanr", cores = 6,
                    prior = bh_priors)
bh_brm_prior

# fit the model to the data
bh_brm_post <- brm(bf(bh_form,
                      lami ~ 1,
                      aii ~ 1,
                      aij ~ 1,
                      nl = TRUE),
                   family = gaussian(), data = bh_dat,
                   backend = "cmdstanr", cores = 6,
                   prior = bh_priors)
bh_brm_post

# draw samples from the two models
samps_prior <- bh_dat %>%
  data_grid(Ni = 0:100,
            Nj = 0:100) %>%
  filter(Nj == 0) %>%
  add_epred_draws(bh_brm_prior, ndraws = 50)
samps_post <- bh_dat %>%
  data_grid(Ni = 0:100,
            Nj = 0:100) %>%
  filter(Nj == 0) %>%
  add_epred_draws(bh_brm_post, ndraws = 50)

# compare prior & posterior distributions
ggplot(data = NULL, aes(x = Ni, y = .epred, group = .draw))+
  geom_line(data = samps_prior, color = "green", alpha = 0.75)+
  geom_line(data = samps_post, color = "blue", alpha = 0.75)+
  scale_y_log10()
