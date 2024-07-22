## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
pnorm(9,5,2)-pnorm(1,5,2)
qnorm(c(0.025,0.975),5,2)

## -----------------------------------------------------------------------------
time <- c(6,6,6,6,7,9,10,10,11,13,16,17,19,20,22,23,25,32,32,34,35,
          1,1,2,2,3,4,4,5,5,8,8,8,8,11,11,12,12,15,17,22,23)

event <- c(0,1,1,1,1,0,0,1,0,1,1,0,0,0,1,1,0,0,0,0,0,
           1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

treat <- c(rep(0,21),rep(1,21))

## -----------------------------------------------------------------------------
library(survival)

y <- Surv(time,event)

## -----------------------------------------------------------------------------
cph_fit <- coxph(y ~ treat)

summary(cph_fit)

## -----------------------------------------------------------------------------
confint(cph_fit)
exp(confint(cph_fit))

## -----------------------------------------------------------------------------
coxph(y ~ treat,ties="efron")

coxph(y ~ treat,ties="breslow")

coxph(y ~ treat,ties="exact")

## -----------------------------------------------------------------------------
N <- sum(time) # total number of person-time observations
x <- matrix(0,N,3)  # columns id, time, and treatment group
colnames(x) <- c("id","t","treat")
y <- matrix(0,N,1) # only one column since only one competing risk

## -----------------------------------------------------------------------------
next_row <- 1  # next available row in the person-time matrices
for (i in seq_along(time)) {
  rows <- next_row:(next_row + time[i] - 1)  # person-time rows to fill
  x[rows,1] <- rep(i,time[i]) # subject id number is constant for each person
  x[rows,2] <- seq_len(time[i])  # discrete time is integers 1,2,...,time[i]
  x[rows,3] <- rep(treat[i],time[i])  # treatment group is constant
  y[rows,1] <- c(rep(0,time[i] - 1),event[i])  # outcome is 0's until study time
  next_row <- next_row + time[i]  # increment the next available row pointer
}

## -----------------------------------------------------------------------------
original_data <- data.frame(id=seq_along(time),time,event,treat)
head(original_data,2)

## -----------------------------------------------------------------------------
expanded_data <- data.frame(x,y)
head(expanded_data,12)

## -----------------------------------------------------------------------------
linear_fit <- glm(y ~ t + treat,family=binomial,data=expanded_data)
summary(linear_fit)

## -----------------------------------------------------------------------------
beta2_hat <- coef(linear_fit)["treat"]
exp(beta2_hat)

## -----------------------------------------------------------------------------
quadratic_fit <- glm(y ~ poly(t,2) + treat,family=binomial,data=expanded_data)
exp(coef(quadratic_fit)["treat"])
cubic_fit <- glm(y ~ poly(t,3) + treat,family=binomial,data=expanded_data)
exp(coef(cubic_fit)["treat"])

## -----------------------------------------------------------------------------
expanded_data$t_cat <- cut(expanded_data$t,c(0,10,20,35))

## -----------------------------------------------------------------------------
with(expanded_data,table(t,t_cat))

## -----------------------------------------------------------------------------
step_fit <- glm(y ~ t_cat + treat,family=binomial,data=expanded_data)
summary(step_fit)

## -----------------------------------------------------------------------------
exp(coef(step_fit)["treat"])

## -----------------------------------------------------------------------------
confint(step_fit)
exp(confint(step_fit)["treat",])

## ----eval=FALSE---------------------------------------------------------------
#  brea_mcmc(x, y, priors = NULL, S = 1000, B = 100, n = NULL, K = NULL, store_re = FALSE)

## -----------------------------------------------------------------------------
x_brea <- matrix(0,N,2)
x_brea[,1] <- cut(expanded_data$t,c(0,10,20,35),labels=FALSE) # grouped time t
x_brea[,2] <- expanded_data$treat + 1 # treatment group

## -----------------------------------------------------------------------------
library(brea)
set.seed(1234)
fit <- brea_mcmc(x_brea,y)

## -----------------------------------------------------------------------------
str(fit)

## -----------------------------------------------------------------------------
b_treatment <- fit$b_m_s[[2]][1,1,]
b_control <- fit$b_m_s[[2]][1,2,]
d <- b_control - b_treatment # sampled values of treatment effect on logit scale

## -----------------------------------------------------------------------------
mean(d) # posterior mean point estimate
median(d) # posterior median point estimate
sd(d) # posterior standard deviation (standard error)
summary(step_fit) # identical frequentist model fit with glm() earlier

## ----fig.width=6,fig.height=3-------------------------------------------------
par(cex=0.66,mgp=c(1.75,0.5,0),mai=c(0.4,0.4,0.1,0.1))
plot(d,type="l",xlab="MCMC Interation Number s",
     ylab="Sampled Value of Treatment Effect Parameter")

## -----------------------------------------------------------------------------
library(coda)
effectiveSize(d)

## -----------------------------------------------------------------------------
set.seed(1234)
fit_10k <- brea_mcmc(x_brea,y,S=10000,B=1000)
d <- fit_10k$b_m_s[[2]][1,2,] - fit_10k$b_m_s[[2]][1,1,]
effectiveSize(d)
mean(d)
median(d)
sd(d)
exp(median(d))

## -----------------------------------------------------------------------------
x_brea[,1] <- cut(expanded_data$t,seq(0,36,3),labels=FALSE)

## -----------------------------------------------------------------------------
set.seed(1234)
priors <- list(list("gmrf",3,0.01),list("cat",4))
fit_gmrf <- brea_mcmc(x_brea,y,priors,S=10000,B=1000)
d_gmrf <- fit_gmrf$b_m_s[[2]][1,2,] - fit_gmrf$b_m_s[[2]][1,1,]
effectiveSize(d_gmrf)

## -----------------------------------------------------------------------------
mean(d_gmrf)
median(d_gmrf)
sd(d_gmrf)

## -----------------------------------------------------------------------------
median(exp(d_gmrf))

## -----------------------------------------------------------------------------
quantile(exp(d_gmrf),c(0.025,0.975))

## ----include=FALSE, eval=FALSE------------------------------------------------
#  summary(coxph(Surv(time,event) ~ treat,ties="efron"))
#  summary(coxph(Surv(time,event) ~ treat,ties="breslow"))
#  summary(coxph(Surv(time,event) ~ treat,ties="exact"))
#  summary(linear_fit); exp(coef(linear_fit)["treat"])
#  summary(quadratic_fit); exp(coef(quadratic_fit)["treat"])
#  summary(cubic_fit); exp(coef(cubic_fit)["treat"])
#  summary(step_fit); exp(coef(step_fit)["treat"])
#  median(d); sd(d); median(exp(d))
#  median(d_gmrf); sd(d_gmrf); median(exp(d_gmrf))

