# Load libraries
library(tidyverse)
library(lmtest)

# Load data
load("logisticdata.Rdata")

# Fit model
M <- glm(isFaculty ~ year*gender + ta + ah + advisornumadvised, data = logisticdata, family = binomial())

# Likelihood test
lrtest(M)

# Summary of model
summary(M)