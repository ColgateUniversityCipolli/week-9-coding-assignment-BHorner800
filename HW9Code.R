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
# Maximum Likelihood (Weibull)
#####################################
llweibull <- function(par, data, neg=F){
  # a <- par[1]
  # sigma <- par[2]
  a <- exp(par[1]) # go from (-inf,inf) to (0,inf)
  sigma <- exp(par[2]) # go from (-inf,inf) to (0,inf)
  
  ll <- sum(log(dweibull(x=data, shape=a, scale=sigma)), na.rm=T)
  
  return(ifelse(neg, -ll, ll))
}

MLEs <- optim(fn = llweibull,
              par = c(1,1),
              data = dat.precip.long$Precipitation,
              neg=T)

(MLEs$par <- exp(MLEs$par)) # transform
weibull.loglike <- -MLEs$value

#####################################
# MLE Gamma
#####################################
llgamma <- function(data, par, neg=F){
  alpha <- par[1]
  beta <- par[2]
  
  loglik <- sum(log(dgamma(x=data, shape=alpha, rate=beta)), na.rm = T)
  
  return(ifelse(neg, -loglik, loglik))
}

(gamma.mles <- optim(par = c(1,1),
               fn = llgamma,
               data=dat.precip.long$Precipitation,
               neg=T))
gamma.alpha.mle <- gamma.mles$par[1]
gamma.beta.mle <- gamma.mles$par[2]

gamma.loglike <- -gamma.mles$value

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

(lognorm.mles <- optim(par = c(0, 1),           # Initial guesses for mu and sigma
               fn = lllognorm,          # Log-likelihood function
               data = dat.precip.long$Precipitation,  # Data for fitting
               neg = T))                # To minimize the negative log-likelihood

lognorm.mu.mle <- lognorm.mles$par[1]
lognorm.sigma.mle <- lognorm.mles$par[2]

lognorm.loglike <-  -lognorm.mles$value  # Log-likelihood for the Weibull distribution


#####################################
# Better fit according to likelihood ratio
#####################################
Q1 <- exp(weibull.loglike - gamma.loglike)
Q2 <- exp(weibull.loglike - lognorm.loglike)
Q3 <- exp(gamma.loglike - lognorm.loglike)

Q1
Q2
Q3
