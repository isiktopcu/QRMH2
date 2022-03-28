library(foreign)
library(ggplot2)
library(MASS)
library(texreg)
library(pscl)
library(lattice) #zero-inflated model

setwd("~/Desktop")
set.seed(0)
load("/Users/isiktopcu/Downloads/TurkeyPKKvotes.rda")
head(dataNewC)
names(dataNewC)

#probability density plots for the potential outcomes
par(mfrow=c(3,3))
vars <- c("attack_p22","nkill_n_p22","nwound_n_p22")
for (var in vars) plot(density(dataNewC[,var],na.rm=TRUE, main=var))
#zero-inflated and positively skewed which might suggest a very low lambda value for Poisson 
#(the bigger the lambda, the more symmetric the distribution)
#zero percentages of the outcomes 
vars <- c("attack_p22","nkill_n_p22","nwound_n_p22")
for (var in vars) print(100*sum(dataNewC[,var] == 0)/nrow(dataNewC))
#[1] 90.70796
#[1] 94.69027
#[1] 93.36283

#to use Poisson, we must check for the assumption that the mean and the variance are equal. 
vars <- c("attack_p22","nkill_n_p22","nwound_n_p22")
for (var in vars) print(mean(dataNewC[,var]))
for (var in vars) print(var(dataNewC[,var]))
#means and variance values of each distribution are completely different. We can not use Poisson because of it's mean = variance assumption) 

#overdispersion

E1 <- resid(pois, type = "pearson")
N1  <- nrow(dataNewC)
p1  <- length(coef(pois))   
sum(E1^2) / (N1 - p1)
#[1] 4.142692

E2 <- resid(atNB, type = "pearson")
N2  <- nrow(dataNewC)
p2  <- length(coef(atNB))  +1 #'+1' is for variance parameter in NB
sum(E2^2) / (N2 - p2)
#[1] 1.549871
#Looks like our model still produces overdisperion even though we used a negative binomial model. 
#Seems better with the zero-inflated negative binomial 
E3 <- resid(atZI.log, type = "pearson")
N3  <- nrow(dataNewC)
p3  <- length(coef(atZI.log)) +1  # '+1' is due to theta  
sum(E3^2) / (N3 - p3)
#[1] 0.9261557


#poisson 
pois <- glm(attack_p22 ~ margin + population + urbanization_rate + as.numeric(border),
            family="poisson", data = dataNewC)
summary(pois)
#standard errors are underestimated, very unrealistic low p values, inflated Type 1 error, false positives. 

#negative binomial for the three outcomes, we'll model on "attack_p22"
atNB <- glm.nb(attack_p22 ~ margin + population + urbanization_rate + as.numeric(border),  data = dataNewC)
summary(atNB)
texreg(atNB)  

killNB <- glm.nb(nkill_n_p22 ~ margin + population + urbanization_rate + as.numeric(border) , data = dataNewC )
summary(killNB)
texreg(killNB)
rm(killNB)

wNB <- glm.nb(nwound_n_p22 ~ margin + population + urbanization_rate + as.numeric(border) , data = dataNewC )
summary(wNB)
texreg(wNB)  
rm(wNB)

#let's give the zero inflated NG to see if it improves
#i've tried to scale and performed the log(1+bc) scale the explanatory variables,
#which sets the mean of the predictor to 0 and the SD to 1, which generally improves performance.
#because R wasn't able to compute the hessian otherwise, the SE were NA.
atZI.s <- zeroinfl(attack_p22 ~ 
                     scale(margin) +
                     scale(population) + 
                     scale(urbanization_rate)+
                     scale(as.numeric(border)),
                     dist = 'negbin',
                     data = dataNewC)
summary(atZI.s)
texreg(atZI.s) #not significant so we tried scaling differently

atZI.log <- zeroinfl(attack_p22 ~ 
                       log10(1+ margin) +
                       log10(1+ population) + 
                       log10(1 + urbanization_rate)+ 
                       log10(1+ as.numeric(border)),
                 dist = 'negbin',
                 data = dataNewC)
summary(atZI.log)
texreg(atZI.log)

#model fit 
null <- glm.nb(attack_p22 ~ population + urbanization_rate + as.numeric(border),  data = dataNewC) #explanatory variable margin excluded
comp <- model.comparison(atNB, null)
comp
#$statistics
#aic     bic bayes.factor
#atNB 812.797 841.638     6643.683
#null 835.207 859.241        0.000

#$predicted_differences
#0%    25%    50%    75%   100% 
#0.000  0.055  0.068  0.109 41.367 
#aic and bic favor the null model.


#likelihood-ratio test 

summary(atNB)
#urbanization_rate seems to be the highsest p value variable. can we exclude it?
simplermodel <- glm(attack_p22 ~ margin + population  + as.numeric(border),  data = dataNewC)
#install.packages("car")
library(lmtest)
lrtest(atNB, simplermodel)
#Likelihood ratio test
#Model 1: attack_p22 ~ margin + population + urbanization_rate + as.numeric(border)
#Model 2: attack_p22 ~ margin + population + as.numeric(border)
#Df  LogLik Df  Chisq Pr(>Chisq)    
#1   6  -400.4                         
#2   5 -1546.3 -1 2291.8  < 2.2e-16 ***  we keep the first model. 
  ---

#Matching 
#install.packages("MatchIt")
library(MatchIt)

match.out = matchit(I(margin > mean(margin)) ~ 
                      population + urbanization_rate + as.numeric(border),
                    data = dataNewC,
                    method = 'nearest', distance = 'mahalanobis',
                    replace = T)
matches = as.numeric(match.out$match.matrix)
matches
matches2 = as.numeric(row.names(match.out$match.matrix))
matches2
mean(dataNewC[matches2, 'attack_p22'] - dataNewC[matches, 'attack_p22'])
dataNewC[matches2[1], 'attack_p22'] - dataNewC[matches[1], 'attack_p22']


#anova 
#install.packages("stargazer")
library(stargazer)
stargazer(atNB, null)


