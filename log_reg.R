#logistic regression analysis in 
#Temporal Dynamics of Faculty Hiring in Mathematics

#load libraries
library(tidyverse)
library(lmtest)


# Get list of top schools
load("logisticdata.Rdata")

M <- glm(isFaculty ~ year*gender + ta + ah + advisornumadvised, data = logisticdata, family = binomial())

#Likelihood Ratio Test
lrtest(M)

#Summary 
summary(M)
