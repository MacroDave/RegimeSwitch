library(ggthemes); library(stringr);library(reshape2);library(rstan); library(dplyr); library(ggplot2); library(Quandl);
Quandl.api_key("BrfAiD2VDaxQEw2zgEfQ")

GNP <- Quandl("FRED/GNP") %>%
  mutate(DATE = as.Date(DATE)) %>% 
  arrange(DATE) %>% 
  mutate(VALUE = 100*log(VALUE/lag(VALUE, 1))) %>% 
  filter(!is.na(VALUE)) %>%
  rename(GNP = VALUE)

Y <- subset(GNP, DATE>"1952-04-01" & DATE < "1985-01-01")
Y <- Y[,-1]

compiled_model <- stan_model("C:/Users/wb398198/Documents/RegimeSwitch/regime_stan.stan")
regime_test <- sampling(compiled_model,data = list(T = 133, S = 2, L = 4, Y= Y), iter =400, chains=1, cores = 4)

library(reshape2); library(stringr); library(readr)

# Print estimated parameters from the model
print(regime_test, pars = c("constant_high", "constant_low", "your_markov_transition_matrix", "sigma", "beta"))
 