#This file containt the R code for the chapter 
#"Comparing Bayesian and frequentist models of language variation: 
#The case of help + (to) Infinitive"
#by Natalia Levshina (c)
#April 9 2019


mydata_big <- read.table("mydata_big.txt", header = T)
mydata_big$Helper <- relevel(mydata_big$Helper, ref = "Inanim")  
mydata_small <- read.table("mydata_small.txt", header = T)
mydata_small$Helper <- relevel(mydata_small$Helper, ref = "Inanim") 
mydata_small$Stress <- relevel(mydata_small$Stress, ref = "Good") 

#Section 5.1. ML model
library(lme4)
glmer_big <- glmer(Response ~ Year_new + Horror*Distance_log + MorphForm*Helpee + Helper + (1|Verb), data = mydata_big, family = binomial, glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(glmer_big)
#parsimonious model
glmer_big1 <- glmer(Response ~  Year_new + Horror*Distance_log + MorphForm*Helpee + (1|Verb), data = mydata_big, family = binomial, glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 100000)))
summary(glmer_big1)

#goodness of fit
library(MuMIn)
r.squaredGLMM(glmer_big)
glmer_predict <- predict(glmer_big, type = "response")
library(Hmisc)
somers2(glmer_predict, as.numeric(mydata_big$Response) - 1)

#Section 5.2. Bayesian models
library(brms)

brm_flat <- brm(Response ~ Year_new + Horror*Distance_log + MorphForm*Helpee + Helper + (1|Verb), data = mydata_big, family = bernoulli)
summary(brm_flat)

priors_student <- prior(student_t(3, 0, 1), class = b)
brm_student <- brm(Response ~ Year_new + Horror*Distance_log + MorphForm*Helpee + Helper + (1|Verb), data = mydata_big, family = bernoulli, prior = priors_student)
summary(brm_student)

priors_cauchy <- prior(cauchy(0, 2.5), class = b)
brm_cauchy <- brm(Response ~ Year_new + Horror*Distance_log + MorphForm*Helpee + Helper + (1|Verb), data = mydata_big, family = bernoulli, prior = priors_cauchy)
summary(brm_cauchy)

priors_norm <- prior(normal(0, 2.5), class = b)
brm_norm <- brm(Response ~ Year_new + Horror*Distance_log + MorphForm*Helpee + Helper + (1|Verb), data = mydata_big, family = bernoulli, prior = priors_norm)
summary(brm_norm)

#posterior probabilities of P (H|Data)
post_cauchy <- posterior_samples(brm_cauchy)
colnames(post_cauchy)[1:14]
mean(post_cauchy[, 1] > 0)
#etc.

###computing the C-index in a Bayesian way
fitted_brm <- extract_draws(brm_cauchy)
mymatrix_fixed <- fitted_brm$dpars$mu$fe$X
mymatrix_random <- fitted_brm$dpars$mu$re$Z[[1]]

mydraws_fixed <- fitted_brm$dpars$mu$fe$b
mydraws_random <- fitted_brm$dpars$mu$re$r$Verb

response <- fitted_brm$data$Y

C_draws <- numeric(4000)
for (i in 1:4000){
  predicted_i <- mymatrix_fixed %*% mydraws_fixed[i,] + mymatrix_random %*% mydraws_random[i,]
  prob_i <- exp(predicted_i)/(1 + exp(predicted_i))
  C_i <- somers2(predicted_i, response)[1]
  C_draws[i] <- C_i
}

mean(C_draws)
#0.8372801
quantile(C_draws, 0.025)
#2.5% 
#0.8197108 
quantile(C_draws, 0.975)
#97.5% 
#0.8548962 

#Section 6
glmer_small <- glmer(Response ~ Year_new + Horror*Distance_log + MorphForm*Helpee + Helper + Stress + (1|Verb), data = mydata_small, family = binomial, control = glmerControl(optimizer="bobyqa", optCtrl = list(maxfun = 1000000)))
summary(glmer_small)
r.squaredGLMM(glmer_small)
glmer_predict1 <- predict(glmer_small, type = "response")
somers2(glmer_predict1, as.numeric(mydata_small$Response) - 1)

priors_small <- c(prior(normal(-1, 0.2), class = "Intercept"),
                  prior(normal(-0.06, 0.01), class = b, coef = "Year_new"),
                  prior(normal(-5.7, 1.4), class = b, coef = "HorrorYes"), 
                  prior(normal(0.9, 0.4), class = b, coef = "Distance_log"),
                  prior(normal(0.6, 0.2),  class = b, coef = "MorphFormhelped"),
                  prior(normal(3.2, 0.3),  class = b, coef = "MorphFormhelping"),
                  prior(normal(1.2, 0.2),  class = b, coef = "MorphFormhelps"),
                  prior(normal(-1, 0.4), class = b, coef = "HelpeeYes"),
                  prior(normal(-0.3, 0.2), class = b, coef = "HelperAnim"),
                  prior(normal(-0.2, 0.4), class = b, coef = "HelperMissing"),
                  prior(normal(3.4, 1.3), class = b, coef = "HorrorYes:Distance_log"),
                  prior(normal(0.7, 0.4), class = b, coef = "MorphFormhelped:HelpeeYes"),
                  prior(normal(-1.8, 0.4), class = b, coef = "MorphFormhelping:HelpeeYes"),
                  prior(normal(0.2, 0.4), class = b, coef = "MorphFormhelps:HelpeeYes"),
                  prior(cauchy(0, 2.5), class = b, coef = "StressClash"),
                  prior(cauchy(0, 2.5), class = b, coef = "StressLapse"),
                  prior(normal(0, 0.5), class = sd, group = "Verb"))

brm_small <- brm(Response ~ Year_new + Horror*Distance_log + MorphForm*Helpee + Helper + Stress + (1|Verb), data = mydata_small, family = bernoulli, prior = priors_small)
summary(brm_small)

posteriors_small <- posterior_samples(brm_small)
mean(posteriors_small[, 1] > 0) #etc.
