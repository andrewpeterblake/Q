## BVAR: dummy priors
## author: Andrew Blake
## date: 31 MArch 2023

## Disclaimer: The Bank of England does not accept any liability for misleading 
## or inaccurate information or omissions in the information provided. The 
## subject matter reflects the views of the individual presenter and not the 
## wider Bank of England or its Policy Committees.

##### Libraries

library(tidyverse)
library(lubridate)

library(progress)
library(quantmod)     # For data access

library(RColorBrewer) # Colors

source("BVARdumFUNCs.R")

############################
# Data
############################

# Create a data frame named Y with:
# Date in column 1 (in date format)
# Rest of variables in next N columns

# Below retrieves data from FRED and puts in required format
series <- c('A191RO1Q156NBEA', 
            'CPALTT01USQ661S',
            'FEDFUNDS')       # FRED codes: US GDP growth, CPI, R

Growth  <- broom::tidy(getSymbols(series[1], src='FRED', auto.assign=FALSE)) %>% 
  mutate(series="Growth") 
CPI     <- broom::tidy(getSymbols(series[2], src='FRED', auto.assign=FALSE)) %>% 
  mutate(series="CPI") 
FedFunds <- broom::tidy(getSymbols(series[3], src='FRED', auto.assign=FALSE)) %>% 
  mutate(series="FedFunds") 

Y0 <- bind_rows(Growth, CPI, FedFunds) %>%
  pivot_wider(names_from=series, values_from=value) %>% 
  mutate(Inflation=100*(CPI/lag(CPI,4)-1)) %>%
  select(Date=index, Growth, Inflation, FedFunds) %>% 
  drop_na() # Drop missing obs to balance dataset

Y0 %>% 
  pivot_longer(cols = -Date) %>% 
  ggplot() +
  geom_line(aes(x=Date, y=value, color=name), show.legend = FALSE) +
  facet_wrap(~ name) +
  theme_minimal() +
  labs(x="", y="")

Y <- Y0[,-4]

#########
# Options
#########

nf <- 12 # Max forecast horizon
nb <- 21 # No. back periods plotted in graphs
l  <- 2  # Number of lags in VAR

# specify parameters of the Minnesota-type prior
tau    <- .1   # controls prior on own 1st lags (1 makes wibbly)
d      <- 1    # decay for higher lags
lambda <- 1    # prior for the constant
gamma  <- 1    # sum of coefficients unit roots
delta  <- 1    # cointegration prior

# Gibbs control
reps <- 20000 # total numbers of Gibbs iterations
burn <- 10000 # number of burn-in iterations

#### Example 1: BVAR(2) with $\tau=0.1$

pce   <- list()

Yplus <- augmentData(Y, l, tau, d, lambda, gamma, delta)
out   <- Gibbs_estimate(Yplus[[1]], Yplus[[2]], reps, burn, 1, nf)

# String with control info
controls <- paste0("Lag length ", l, ": tau=", tau, ", d=", d, 
                   ", lambda=", lambda, ", gamma=", gamma, ", delta=", delta)

# Plots
fan_chart(Y, out[[3]], controls, nb)
p <- coeff_plot(Y, l, out[[1]], out[[2]], 333, controls)
pce[[1]] <- p[[1]]

#### Example 2: BVAR(2) with $\tau=1$

tau   <- 1
Yplus <- augmentData(Y, l, tau, d, lambda, gamma, delta)
out   <- Gibbs_estimate(Yplus[[1]], Yplus[[2]], reps, burn, 1, nf)

# String with control info
controls <- paste0("Lag length ", l, ": tau=", tau, ", d=", d, 
                   ", lambda=", lambda, ", gamma=", gamma, ", delta=", delta)

# Plots
fan_chart(Y, out[[3]], controls, nb)
p <- coeff_plot(Y, l, out[[1]], out[[2]], 333, controls)
pce[[2]] <- p[[1]]

#### Example 3: BVAR(6) with $\tau=0.1$

l   <- 6
tau <- 0.1

Yplus <- augmentData(Y, l, tau, d, lambda, gamma, delta)
out   <- Gibbs_estimate(Yplus[[1]], Yplus[[2]], reps, burn, 1, nf)

# String with control info
controls <- paste0("Lag length ", l, ": tau=", tau, ", d=", d, 
                   ", lambda=", lambda, ", gamma=", gamma, ", delta=", delta)

# Plots
fan_chart(Y, out[[3]], controls, nb)
p <- coeff_plot(Y, l, out[[1]], out[[2]], 333, controls)
pce[[3]] <- p[[1]]

#### Example 4: BVAR(6) with $\tau=1$

l   <- 6
tau <- 1

Yplus <- augmentData(Y, l, tau, d, lambda, gamma, delta)
out   <- Gibbs_estimate(Yplus[[1]], Yplus[[2]], reps, burn, 1, nf)

# String with control info
controls <- paste0("Lag length ", l, ": tau=", tau, ", d=", d, 
                   ", lambda=", lambda, ", gamma=", gamma, ", delta=", delta)

# Plots
fan_chart(Y, out[[3]], controls, nb)
p <- coeff_plot(Y, l, out[[1]], out[[2]], 333, controls)
pce[[4]] <- p[[1]]

#### Coefficient estimates
  
###### Plot estimated param densities #####
pce[[1]]
pce[[2]]
pce[[3]]
pce[[4]]

## Tri-variate BVAR
  
Y   <- Y0

#### Example 5: BVAR(4) with $\tau=.1$, 3 variables

l   <- 4
tau <- .1
Yplus <- augmentData(Y, l, tau, d, lambda, gamma, delta)
out   <- Gibbs_estimate(Yplus[[1]], Yplus[[2]], reps, burn, 1, nf)

# String with control info
controls <- paste0("Lag length ", l, ": tau=", tau, ", d=", d, 
                   ", lambda=", lambda, ", gamma=", gamma, ", delta=", delta)

# Coeffs
p <- coeff_plot(Y, l, out[[1]], out[[2]], 0, controls)
p[[1]]

# Fan chart
fan_chart(Y, out[[3]], controls, nb)

#### Example 6: BVAR(6) with $\tau=1$, 3 variables

l   <- 6
tau <- 1
Yplus <- augmentData(Y, l, tau, d, lambda, gamma, delta)
out   <- Gibbs_estimate(Yplus[[1]], Yplus[[2]], reps, burn, 1, nf)

# String with control info
controls <- paste0("Lag length ", l, ": tau=", tau, ", d=", d, 
                   ", lambda=", lambda, ", gamma=", gamma, ", delta=", delta)

# Coeffs
p <- coeff_plot(Y, l, out[[1]], out[[2]], 0, controls)
p[[1]]
# Fan chart
fan_chart(Y, out[[3]], controls, nb)
