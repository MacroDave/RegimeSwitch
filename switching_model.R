#setwd("~/Downloads/")
library(rstan); library(ggplot2); library(dplyr); library(reshape2);
options(mc.cores = parallel::detectCores())

GNP <- read.csv("gnp.csv")
Y <- GNP[,-1]


# Compile model ------------------------------------------------


compiled_model <- stan_model("switching_model.stan")

# Run, plot model ---------------------------------------------------------

regime_op <- optimizing(compiled_model,data = list(T = 135, S = 2, L = 4, Y= Y))
optimized_estimates <- exp(regime_op$par[grepl("Xi", names(regime_op$par))])

regime_test <- sampling(compiled_model,data = list(T = 135, S = 2, L = 4, Y= Y), iter =2000, chains=4, cores = 4, control = list(adapt_delta = 0.95))
regime_test

print(regime_test, pars = c("sigma", "beta", "alpha", "p11", "p22"))

outcome_data <- regime_test %>% 
  as.data.frame() %>% 
  select(contains("Xi")) %>%
  melt() %>% 
  group_by(variable) %>% 
  summarise(lower = quantile(exp(value), 0.1), 
            median = median(exp(value)),
            upper = quantile(exp(value), 0.9)) %>% 
  mutate(Date = as.Date(GNP$Date, format = "%m/%d/%Y"),
         optim_estimates = optimized_estimates)


(outcome_data %>% 
  ggplot(aes(x = Date, y= optim_estimates) ) + geom_line() +
  xlim(first(outcome_data$Date), last(outcome_data$Date))) %>% nberShade()

pp <- outcome_data %>% 
  ggplot(aes(x = Date, y= median, group = 1) ) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "orange", alpha = 0.3) +
  geom_line() +
  xlim(first(outcome_data$Date), last(outcome_data$Date)) +
  #ggthemes::theme_economist() +
  ylim(0, 1) +
  ylab("Probability of being in slow state")

  library(tis)
  nberShade(pp)


# End ---------------------------------------------------------------------
#shinystan::launch_shinystan(regime_test)
