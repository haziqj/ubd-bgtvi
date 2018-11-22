## ---- prelim ----
library(iprobit)
library(tidyverse)
library(directlabels)
# library(gganimate)
library(animation)
# library(reshape2)
# library(directlabels)
# library(rstan)
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# stan2coda <- function(fit) {
#   mcmc.list(lapply(1:ncol(fit), function(x) mcmc(as.array(fit)[,x,])))
# }
# library(ggmcmc)
# library(coda)


# Function to determine even numbers
is_even <- function(x) x %% 2 == 0
