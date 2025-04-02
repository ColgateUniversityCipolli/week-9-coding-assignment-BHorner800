################################################################################
# HW9 Code
# Ben Horner
# MATH 240 - SPRING 2025
################################################################################
library(tidyverse)
library(patchwork)
library(nleqslv)


#The following is pulled from lecture16 code
################################################################################
# Precipitation in Madison County
################################################################################
dat.precip <- read_csv(file = "agacis.csv")

#####################################
# Clean Data
#####################################
dat.precip.long <- dat.precip |>    
  dplyr::select(-Annual) |>                   # Remove annual column 
  pivot_longer(cols = c(Jan, Feb, Mar, Apr,   # pivot the column data into one col
                        May, Jun, Jul, Aug, 
                        Sep, Oct, Nov, Dec), 
               values_to = "Precipitation",   # store the values in Precipitation
               names_to = "Month") |>         # store the months in Month
  mutate(Precipitation = case_when(Precipitation == "M" ~ NA_character_,
                                   TRUE                 ~ Precipitation))|>
  mutate(Precipitation = as.numeric(Precipitation))

#####################################
# Summarize the data
#####################################
ggplot(data=dat.precip.long) +
  geom_histogram(aes(x=Precipitation, y=after_stat(density)),
                 breaks=seq(0, 15, 1),
                 color="grey")+
  geom_hline(yintercept = 0)+
  theme_bw() +
  xlab("Precipitation (Inches)")
ylab("Density")

library(e1071)
dat.precip.long |>
  summarize(
    mean = mean(Precipitation, na.rm=T),
    sd = sd(Precipitation, na.rm=T),
    min = min(Precipitation, na.rm=T),
    max = max(Precipitation, na.rm=T),
    skew = skewness(Precipitation, na.rm=T),
    kurt = kurtosis(Precipitation, na.rm=T)
  )

#####################################
# MLE Gamma
#####################################
llgamma <- function(data, par, neg=F){
  alpha <- par[1]
  beta <- par[2]
  
  loglik <- sum(log(dgamma(x=data, shape=alpha, rate=beta)), na.rm = T)
  
  return(ifelse(neg, -loglik, loglik))
}

(mles <- optim(par = c(1,1),
               fn = llgamma,
               data=dat.precip.long$Precipitation,
               neg=T))
gamma.alpha.mle <- mles$par[1]
gamma.alpha.mle <- mles$par[2]

#####################################
# MLE Log-Normal
#####################################
lllognorm <- function(data, par, neg=F){
  mu <- par[1]
  sigma <- par[2]
  
  # Ensure sigma > 0 (standard deviation must be positive)
  if(sigma <= 0) {
    return(Inf)
  }
  
  loglik <- sum(log(dlnorm(data, meanlog = mu, sdlog = sigma)), na.rm = T)
  
  return(ifelse(neg, -loglik, loglik))
}

(mles <- optim(par = c(0, 1),           # Initial guesses for mu and sigma
               fn = lllognorm,          # Log-likelihood function
               data = dat.precip.long$Precipitation,  # Data for fitting
               neg = T))                # To minimize the negative log-likelihood

lognorm.mu.mle <- mles$par[1]
lognorm.sigma.mle <- mles$par[2]

