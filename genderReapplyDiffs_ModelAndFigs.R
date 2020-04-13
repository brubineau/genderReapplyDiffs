# (c) Brian Rubineau
# 2020-04-07
# Formal model of segregating effects from gender differences in Reapplication After Rejection (RAR)

# SCRIPT CONTENTS:
# 1. Set up environment
#    1.1. Clear environment
#    1.2. Load needed libraries
#
# 2. Collect case data and case-based model parameters
#    2.1. Include the data values provided in the paper
#    2.2. Add CIs for g from case data
#
# 3. Model and Variant
#
# 4. Equilibrium Analysis
#
# 5. Analytic Solutions 
#    5.1. Modeling female share among selected applicants
#    5.2 Function with the continuous differential equation model - applications only
#
# 6. Intersubjective interpretation - Sex bias 
#    6.1. Sex bias equivalent
#    6.2. Evaluating interventions - expanded pools
#
# 7. Figures
#    7.1. 4-panel Figure 2 - Change in female share by time, r, b, g parameters
#    7.2. 3-panel Figure 3 - Change in threshold (Exec Search Case) w/change in r, b, g parameters
#
# 8. Table

########################
# 1. Set up environment
# ---------------------
#    1.1. Clear work environment
rm(list=ls())
gc()

#    1.2. Load needed libraries
#         Note: if the "deSolve" package is not already installed, de-comment and run (once only) the following line:
#         install.packages("deSolve")
library(deSolve)

#    1.3. Change to working directory
setwd("~/") # ADAPT FOR YOUR OWN COMPUTING CONTEXT

######################################################
# 2. collect case data and case-based model parameters
# ----------------------------------------------------
# 2.1. Include the data values provided in the paper

cases <- data.frame(caseID=c(   "ExecSearch",    "CrowdFund",    "Patenting"),
                    # MODEL PARAMETERS
                    p=c(                0.16,           0.29,          0.088), # female share among applicants
                    r=c(                0.96,           0.64,           0.86), # rejection rate
                    b=c(              0.6878,         0.1117,         0.5441), # central reapplication rate after rejection
                    g=c(              1.0984,         1.2596,          1.047), # gender difference in reapplication after rejection
                    # OTHER CASE DATA
                    obs.m=c(           21952, round(98131*(1-0.29),0), round(3370800*(1-0.088),0)), # count of men in each case
                    obs.f=c(            1603, round(98131*0.29,0),    round(3370800*0.088,0)), # count of women in each case
                    pr.reapply.m=c( 1-0.2924,        0.13672,         0.5555), # estimates of reapplication probabilities for men
                    pr.reapply.f=c( 1-0.3558,        0.10854,         0.5305)) # estimates of reapplication probabilities for women

# ----------------------------------------------------
# 2.2. Add CIs for g from case data
# standard deviation of binary random variables with probability p is: sqrt(n*p*(1-p))
# There is no formula for the SD of a ratio of 2 random variables even where their individuals SDs are known.
# We simulate to determine the 95% CI of the ratio of the probability that men reapply to the probability that women reapply - our "g" parameter

getGs <- lapply(as.list(1:3),function(c){
  replicate(100000, # generates 100k instances of g (ratio of proportion of men who reapply to proportion of women who reapply)
            (rbinom(n=1,  # rbinom returns number of successes (reapplications) from the same number of men at men's case-identified probability of reapplying
                    size=cases$obs.m[c],
                    prob=cases$pr.reapply.m[c])/cases$obs.m[c])/ # dividing by number of men yields proportion of men who reapplied
              (rbinom(n=1, # rbinom returns number of successes (reapplications) from the same number of women at women's case-identified probability of reapplying
                      size=cases$obs.f[c],
                      prob=cases$pr.reapply.f[c])/cases$obs.f[c])) # dividing by number of women yields proportion of women who reapplied
})

# 95% CI taken from the simulated distribution
cases$glb <- cases$g - (unlist(lapply(getGs,quantile,probs=0.975))-unlist(lapply(getGs,quantile,probs=0.025)))/2 # lower bound of 95%CI
cases$gub <- cases$g + (unlist(lapply(getGs,quantile,probs=0.975))-unlist(lapply(getGs,quantile,probs=0.025)))/2 # upper bound of 95%CI

############################
# 3. MAIN MODEL AND VARIANT
# -------------------------
# 3.1. Modeling female share among selected applicants
#      Continuous model using differential equations
#      No Turnover version

model <- function (t, state, parms) {
  with(as.list(c(state, parms)), {
    #  Formulae for Male/Men (M) and Female/Women (W) as in the Paper
    dMApp <- -MApp + ((g*b/(1+(g*b)-b)) * MRej) + ((1-pctF)) # Male applicants
    dWApp <- -WApp + ((b/((1+(g*b)-b))) * WRej) + (pctF) # Female applicants
    dMRej <- pRej*MApp - MRej # Rejected Men
    dWRej <- pRej*WApp - WRej # Rejected Women
    dMAcc <- (1-pRej)*MApp # Accepted (Selected) Men
    dWAcc <- (1-pRej)*WApp # Accepted (Selected) Women
    list(c(dMApp,dWApp,dMRej,dWRej,dMAcc,dWAcc))
  })
}

# Initial values for the continuous model - all stocks initialized to zero
ini <- c(MApp=0, WApp=0, MRej=0, WRej=0, MAcc=0, WAcc=0)

# Time range over which the system is modeled (20 rounds - this will be varied later)
times <- seq(0, 20, 0.01)

# Function numericaly solves the continuous model using the:
#   specified model parameter values;
#   provided initial conditions; 
#   designated numer of rounds; and
# returns the final female share among all accepted (selected) men and women.

selModel <- function(gParam,bParam,pParam,rParam,rounds){
  parmsOpen <- c(g=gParam,b=bParam,pctF=pParam,pRej=rParam) # specified model parametres
  outOpen <- ode(y = ini, times = seq(0,rounds,0.01), func = model, parms = parmsOpen) # initial conditions, times, and parameters
  # return the final female share among all accepted (selected) Men and Women
  return(outOpen[dim(outOpen)[1],"WAcc"]/(outOpen[dim(outOpen)[1],"WAcc"]+outOpen[dim(outOpen)[1],"MAcc"]))
}

# Add the calculated female share after 20 rounds for each case to the data file
cases$sol.exit0 <- sapply(1:3,function(c){selModel(cases$g[c],cases$b[c],cases$p[c],cases$r[c],20)})

# -------------------------
# 3.2 Continuous differential equation model - applications only
#     Maximum Turnover - Instantaneous Exit

model2 <- function (t, state, parms) {
  with(as.list(c(state, parms)), {
    #  Formulae for Male/Men (M) and Female/Women (W) as in the Paper
    dMApp <- -MApp + ((g*b/(1+(g*b)-b)) * MRej) + ((1-pctF)) # Male applicants
    dWApp <- -WApp + ((b/((1+(g*b)-b))) * WRej) + (pctF) # Female applicants
    dMRej <- pRej*MApp - MRej # Rejected Men
    dWRej <- pRej*WApp - WRej # Rejected Women
    list(c(dMApp,dWApp,dMRej,dWRej))
  })
}

# Initial values for the continuous model - stocks set to zero
ini2 <- c(MApp=0, WApp=0, MRej=0, WRej=0)

# Function numericaly solves the continuous model
# using the:
#   specified model parameter values;
#   provided initial conditions; 
#   designated number of rounds; and
# returns the final female share among applicants.

appModel <- function(gParam,bParam,pParam,rParam,rounds){
  parmsOpen <- c(g=gParam,b=bParam,pctF=pParam,pRej=rParam) # specified model parametres
  outOpen <- ode(y = ini2, times = seq(0,rounds,0.01), func = model2, parms = parmsOpen) # initial conditions, times, and parameters
  # return the final female share among all accepted (selected) Men and Women
  return(outOpen[dim(outOpen)[1],"WApp"]/(outOpen[dim(outOpen)[1],"WApp"]+outOpen[dim(outOpen)[1],"MApp"]))
}

# Add the calculated female share after 20 rounds for each case to the data file
# Note: The "8" is used to designate maximum exit as it is the rotated symbol for infinity.
cases$sol.exit8 <- sapply(1:3,function(c){appModel(cases$g[c],cases$b[c],cases$p[c],cases$r[c],20)})

##################################################
# 4. Find equilibrium values for model and variant
# ------------------------------------------------
# Equilibrium value is the outcome value 
# (i.e., female share among selected - main model, or female share among applicants - variant)
# as time goes to infinity. For our purposes, we specify the number of digits desired, 
# and take the outcome once the model produces identical values 50 rounds apart to the specified number of digits in two successive rounds.
# we start at 500 rounds

# 4.1  Modeling female share among selected applicants
#      Continuous model using differential equations
#      No Turnover version

selEquil <- function(digits,gParam,bParam,pParam,rParam){
  parmsOpen <- c(g=gParam,b=bParam,pctF=pParam,pRej=rParam) # specified model parametres
  equilNotFound <- T
  equil <- NA
  rounds <- 500
  out <- ode(y = ini, times = seq(0,rounds,0.01), func = model, parms = parmsOpen) # initial conditions, times, and parameters
  outcomes <- round(out[,"WAcc"]/(out[,"WAcc"]+out[,"MAcc"]),digits+1)
  if (sum(as.numeric(outcomes[2:(length(outcomes)-5000)]==outcomes[5002:(length(outcomes))]))>0) { # test for identical outcomes 50 rounds apart (steps are 1/100th of a round)
    equilNotFound <- F
    equil <- outcomes[length(outcomes)]
  }
  while (equilNotFound){
    newIni <- c(out[dim(out)[1],"MApp"], 
                out[dim(out)[1],"WApp"], 
                out[dim(out)[1],"MRej"], 
                out[dim(out)[1],"WRej"], 
                out[dim(out)[1],"MAcc"], 
                out[dim(out)[1],"WAcc"])
    newTimes <- seq(rounds,rounds+500,0.01)
    rounds <- rounds+500
    out <- ode(y = newIni, times = newTimes, func = model, parms = parmsOpen) # initial conditions, times, and parameters
    outcomes <- round(out[,"WAcc"]/(out[,"WAcc"]+out[,"MAcc"]),digits+1)
    if (sum(as.numeric((outcomes[2:(length(outcomes)-5000)]-outcomes[5002:(length(outcomes))])==0))>0) {
      equilNotFound <- F
      equil <- outcomes[length(outcomes)]
    }
  }
  return(round(equil,digits))
}
# Add the calculated female share at equilibrium for each case to the data file (to 5 digits)
cases$selEquil <- sapply(1:3,function(c){selEquil(5,cases$g[c],cases$b[c],cases$p[c],cases$r[c])})

# -------------------------
# 4.2 Continuous differential equation model - applications only
#     Maximum Turnover - Instantaneous Exit

appEquil <- function(digits,gParam,bParam,pParam,rParam){
  parmsOpen <- c(g=gParam,b=bParam,pctF=pParam,pRej=rParam) # specified model parametres
  equilNotFound <- T
  equil <- NA
  rounds <- 500
  out <- ode(y = ini2, times = seq(0,rounds,0.01), func = model2, parms = parmsOpen) # initial conditions, times, and parameters
  outcomes <- round(out[,"WApp"]/(out[,"WApp"]+out[,"MApp"]),digits+1)
  if (sum(as.numeric(outcomes[2:(length(outcomes)-5000)]==outcomes[5002:(length(outcomes))]))>0) { # test for identical outcomes 50 rounds apart (steps are 1/100th of a round)
    equilNotFound <- F
    equil <- outcomes[length(outcomes)]
  }
  while (equilNotFound){
    newIni <- out[dim(out)[1],]
    newTimes <- seq(rounds,rounds+500,0.01)
    rounds <- rounds+500
    out <- ode(y = newIni, times = newTimes, func = model2, parms = parmsOpen) # initial conditions, times, and parameters
    outcomes <- round(out[,"WApp"]/(out[,"WApp"]+out[,"MApp"]),digits+1)
    if (sum(as.numeric((outcomes[2:(length(outcomes)-5000)]-outcomes[5002:(length(outcomes))])==0))>0) {
      equilNotFound <- F
      equil <- outcomes[length(outcomes)]
    }
  }
  return(round(equil,digits))
}
# Add the calculated female share at equilibrium for each case to the data file (to 5 digits)
cases$appEquil <- sapply(1:3,function(c){appEquil(5,cases$g[c],cases$b[c],cases$p[c],cases$r[c])})

########################################################
# 5. ANALYTICAL SOLUTIONS TO THE MAIN MODEL AND VARIANT
# -----------------------------------------------------
# 5.1. Modeling female share among selected applicants
# Solution determined using SageCell for the parameters of the executive search case 
# Reproducable as follows using https://sagecell.sagemath.org/:

# to get the formula for women in selected stock after t rounds:
# t = var('t')
# var('p r b g')
# wr = function('wr')(t)
# wa = function('wa')(t)
# ws = function('ws')(t)
# de1 = diff(wr,t) - (r*wa) + wr == 0
# de2 = diff(wa,t) - p + wa - ((b/((1+(g*b)-b)))*wr) == 0
# de3 = diff(ws,t) - ((1-r)*wa) == 0
# assume(p==0.16)
# assume(r==0.96)
# assume(0.67 < b < 0.71)
# assume(g==1.11)
# desolve_system([de1,de2,de3], [wr,wa,ws], ics=[0,0,0,0],ivar=t)
# SAVE THE "ws" formula

# to get the formula for women in selected stock after t rounds:
# t = var('t')
# var('p r b g')
# mr = function('mr')(t)
# ma = function('ma')(t)
# ms = function('ms')(t)
# de1 = diff(mr,t) - (r*ma) + mr == 0
# de2 = diff(ma,t) - (1-p) + ma - ((b*g/((1+(g*b)-b)))*mr) == 0
# de3 = diff(ms,t) - ((1-r)*ma) == 0
# assume(p==0.16)
# assume(r==0.96)
# assume(0.67 < b < 0.71)
# assume(g==1.11)
# desolve_system([de1,de2,de3], [mr,ma,ms], ics=[0,0,0,0],ivar=t)
# SAVE THE "ms" formula

# to get the final formula for the percent of women in selected stock after t rounds:
# t = var('t')
# var('p r b g')
# ws = function('ws')(t)
# ms = function('ms')(t)
# assume(p==0.16)
# assume(r==0.96)
# assume(0.67 < b < 0.71)
# assume(g==1.11)
# ws(t) = [paste the "ws" formula from above]
# ms(t) = [paste the "ms" formula from above]
# pctf = ws/(ws+ms)
# pctf

rarPctFMinExit <- function(p,r,g,b,t){ # Analytical solution for full model
  return(
    (
      (b*g - b + 1)*p*r*t/(b*g - b*r - b + 1) - (b*g - b + 1)*p*t/(b*g - b*r - b + 1) + 
        (
          (
            (3*(b^3*g^2 + b^3 - 2*b^2 - 2*(b^3 - b^2)*g + b)*p*r^2 + 
               (b^3*g^3 - 4*b^3 - 3*(2*b^3 - b^2)*g^2 + 9*b^2 + 3*(3*b^3 - 4*b^2 + b)*g - 6*b + 1)*p*r - 
               (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 + 3*b^2 + 3*(b^3 - 2*b^2 + b)*g - 3*b + 1)*p)*
              (b*g - b + 1)/(b^2*g^2 + b^2*r^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g - b^2 + b)*r - 2*b + 1) - 
              (
                (b^3*g^2 + b^3 - 2*b^2 - 2*(b^3 - b^2)*g + b)*p*r^2 + 
                  (b^3*g^3 - 2*b^3 - (4*b^3 - 3*b^2)*g^2 + 5*b^2 + (5*b^3 - 8*b^2 + 3*b)*g - 4*b + 1)*p*r - 
                  (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 + 3*b^2 + 3*(b^3 - 2*b^2 + b)*g - 3*b + 1)*p)*
              (b*g - b + 1)/(b^2*g^2 + b^2*r^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g - b^2 + b)*r - 2*b + 1))*
            sinh(sqrt((b*g - b + 1)*b*r)*t/(b*g - b + 1))/
            sqrt((b*g - b + 1)*b*r) + 
            (
              (b^3*g^2 + b^3 - 2*b^2 - 2*(b^3 - b^2)*g + b)*p*r^2 + 
                (b^3*g^3 - 2*b^3 - (4*b^3 - 3*b^2)*g^2 + 5*b^2 + (5*b^3 - 8*b^2 + 3*b)*g - 4*b + 1)*p*r - 
                (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 + 3*b^2 + 3*(b^3 - 2*b^2 + b)*g - 3*b + 1)*p)*
            cosh(
              sqrt(
                (b*g - b + 1)*b*r)*t/(b*g - b + 1))/
            (b^2*g^2 + b^2*r^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g - b^2 + b)*r - 2*b + 1))*
        exp(-t)/(b*g - b + 1) - 
        (
          (b^2*g - b^2 + b)*p*r^2 + (b^2*g^2 + 2*b^2 - (3*b^2 - 2*b)*g - 3*b + 1)*p*r - 
            (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p)/
        (b^2*g^2 + b^2*r^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g - b^2 + b)*r - 2*b + 1))/
      
      (
        (b*g - b + 1)*p*r*t/(b*g - b*r - b + 1) + b*g*t/(b*g*r - b*g + b - 1) - 
          (b*g - b + 1)*p*t/(b*g*r - b*g + b - 1) - (b*g - b + 1)*p*t/(b*g - b*r - b + 1) - 
          (b*g - (b*g - b + 1)*p - b + 1)*r*t/(b*g*r - b*g + b - 1) - b*t/(b*g*r - b*g + b - 1) + 
          (
            (
              (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 - 
                 (b^3*g^3 - 2*(b^3 - b^2)*g^2 + (b^3 - 2*b^2 + b)*g - 
                    (b^3*g^3 - 2*(b^3 - b^2)*g^2 + (b^3 - 2*b^2 + b)*g)*p)*r^2 + 3*b^2 + 3*(b^3 - 2*b^2 + b)*g - 
                 (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 + 3*b^2 + 3*(b^3 - 2*b^2 + b)*g - 3*b + 1)*p + 
                 (b^3 + (b^3 - b^2)*g^2 - 3*b^2 - 2*(b^3 - 2*b^2 + b)*g - 
                    (b^3 + (b^3 - b^2)*g^2 - 3*b^2 - 2*(b^3 - 2*b^2 + b)*g + 3*b - 1)*p + 3*b - 1)*r - 3*b + 1)*
                (b*g - b + 1)/(b^2*g^2*r^2 + b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g^2 - (b^2 - b)*g)*r - 2*b + 1) - 
                (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 - 3*(b^3*g^3 - 2*(b^3 - b^2)*g^2 + 
                                                          (b^3 - 2*b^2 + b)*g - (b^3*g^3 - 2*(b^3 - b^2)*g^2 + 
                                                                                   (b^3 - 2*b^2 + b)*g)*p)*r^2 + 3*b^2 + 
                   3*(b^3 - 2*b^2 + b)*g - (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 + 3*b^2 + 3*(b^3 - 2*b^2 + b)*g - 3*b + 1)*p + 
                   (2*b^3*g^3 + b^3 - 3*(b^3 - b^2)*g^2 - 3*b^2 - 
                      (2*b^3*g^3 + b^3 - 3*(b^3 - b^2)*g^2 - 3*b^2 + 3*b - 1)*p + 3*b - 1)*r - 3*b + 1)*
                (b*g - b + 1)/(b^2*g^2*r^2 + b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g^2 - (b^2 - b)*g)*r - 2*b + 1))*
              sinh(
                sqrt(
                  (b*g - b + 1)*b*g*r)*t/(b*g - b + 1))/
              sqrt(
                (b*g - b + 1)*b*g*r) - (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 - 
                                          (b^3*g^3 - 2*(b^3 - b^2)*g^2 + (b^3 - 2*b^2 + b)*g - 
                                             (b^3*g^3 - 2*(b^3 - b^2)*g^2 + (b^3 - 2*b^2 + b)*g)*p)*r^2 + 3*b^2 + 
                                          3*(b^3 - 2*b^2 + b)*g - (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 + 3*b^2 + 
                                                                     3*(b^3 - 2*b^2 + b)*g - 3*b + 1)*p + 
                                          (b^3 + (b^3 - b^2)*g^2 - 3*b^2 - 2*(b^3 - 2*b^2 + b)*g - 
                                             (b^3 + (b^3 - b^2)*g^2 - 3*b^2 - 2*(b^3 - 2*b^2 + b)*g + 3*b - 1)*p + 3*b - 1)*r - 
                                          3*b + 1)*
              cosh(
                sqrt(
                  (b*g - b + 1)*b*g*r)*t/(b*g - b + 1))/
              (b^2*g^2*r^2 + b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g^2 - (b^2 - b)*g)*r - 2*b + 1))*
          exp(-t)/(b*g - b + 1) + 
          (
            (
              (3*(b^3*g^2 + b^3 - 2*b^2 - 2*(b^3 - b^2)*g + b)*p*r^2 + 
                 (b^3*g^3 - 4*b^3 - 3*(2*b^3 - b^2)*g^2 + 9*b^2 + 3*(3*b^3 - 4*b^2 + b)*g - 6*b + 1)*p*r - 
                 (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 + 3*b^2 + 3*(b^3 - 2*b^2 + b)*g - 3*b + 1)*p)*
                (b*g - b + 1)/(b^2*g^2 + b^2*r^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g - b^2 + b)*r - 2*b + 1) - 
                (
                  (b^3*g^2 + b^3 - 2*b^2 - 2*(b^3 - b^2)*g + b)*p*r^2 + 
                    (b^3*g^3 - 2*b^3 - (4*b^3 - 3*b^2)*g^2 + 5*b^2 + (5*b^3 - 8*b^2 + 3*b)*g - 4*b + 1)*p*r - 
                    (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 + 3*b^2 + 3*(b^3 - 2*b^2 + b)*g - 3*b + 1)*p)*
                (b*g - b + 1)/(b^2*g^2 + b^2*r^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g - b^2 + b)*r - 2*b + 1))*
              sinh(
                sqrt(
                  (b*g - b + 1)*b*r)*t/(b*g - b + 1))/
              sqrt(
                (b*g - b + 1)*b*r) + 
              (
                (b^3*g^2 + b^3 - 2*b^2 - 2*(b^3 - b^2)*g + b)*p*r^2 + 
                  (b^3*g^3 - 2*b^3 - (4*b^3 - 3*b^2)*g^2 + 5*b^2 + (5*b^3 - 8*b^2 + 3*b)*g - 4*b + 1)*p*r - 
                  (b^3*g^3 - b^3 - 3*(b^3 - b^2)*g^2 + 3*b^2 + 3*(b^3 - 2*b^2 + b)*g - 3*b + 1)*p)*
              cosh(
                sqrt(
                  (b*g - b + 1)*b*r)*t/(b*g - b + 1))/
              (b^2*g^2 + b^2*r^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g - b^2 + b)*r - 2*b + 1))*
          exp(-t)/(b*g - b + 1) + (b^2*g^2 - (b^2*g^2 - (b^2 - b)*g - (b^2*g^2 - (b^2 - b)*g)*p)*r^2 + 
                                     b^2 - 2*(b^2 - b)*g - (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p - 
                                     (b^2 - (b^2 - b)*g - (b^2 - (b^2 - b)*g - 2*b + 1)*p - 2*b + 1)*r - 2*b + 1)/
          (b^2*g^2*r^2 + b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g^2 - (b^2 - b)*g)*r - 2*b + 1) - 
          (
            (b^2*g - b^2 + b)*p*r^2 + (b^2*g^2 + 2*b^2 - (3*b^2 - 2*b)*g - 3*b + 1)*p*r - 
              (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p)/
          (b^2*g^2 + b^2*r^2 + b^2 - 2*(b^2 - b)*g - 2*(b^2*g - b^2 + b)*r - 2*b + 1) + 
          t/(b*g*r - b*g + b - 1))
  )
}

# Add the calculated female share after 20 rounds for each case to the data file
cases$ans.exit0 <- sapply(1:3,function(c){rarPctFMinExit(cases$p[c],cases$r[c],cases$g[c],cases$b[c],t=20)})


# -------------------------
# 5.2 Function with the continuous differential equation model - applications only
# Maximum Turnover - Instantaneous Exit
# Solution determined using SageCell for the parameters of the executive search case 
# Reproducable as follows using https://sagecell.sagemath.org/:

# to get the formula for women in selected stock after t rounds:
# t = var('t')
# var('p r b g')
# wr = function('wr')(t)
# wa = function('wa')(t)
# de1 = diff(wr,t) - (r*wa) + wr == 0
# de2 = diff(wa,t) - p + wa - ((b/((1+(g*b)-b)))*wr) == 0
# assume(p==0.16)
# assume(r==0.96)
# assume(0.67 < b < 0.71)
# assume(g==1.11)
# desolve_system([de1,de2], [wr,wa], ics=[0,0,0],ivar=t)
# SAVE THE "wa" formula

# to get the formula for women in selected stock after t rounds:
# t = var('t')
# var('p r b g')
# mr = function('mr')(t)
# ma = function('ma')(t)
# de1 = diff(mr,t) - (r*ma) + mr == 0
# de2 = diff(ma,t) - (1-p) + ma - ((b*g/((1+(g*b)-b)))*mr) == 0
# assume(p==0.16)
# assume(r==0.96)
# assume(0.67 < b < 0.71)
# assume(g==1.11)
# desolve_system([de1,de2], [mr,ma], ics=[0,0,0],ivar=t)
# SAVE THE "ma" formula

# to get the final formula for the percent of women in selected stock after t rounds:
# t = var('t')
# var('p r b g')
# wa = function('wa')(t)
# ma = function('ma')(t)
# assume(p==0.16)
# assume(r==0.96)
# assume(0.67 < b < 0.71)
# assume(g==1.11)
# wa(t) = [paste the "wa" formula from above]
# ma(t) = [paste the "ma" formula from above]
# pctf = wa/(wa+ma)
# simplify(pctf)

rarPctFMaxExit <- function(p,r,g,b,t){
  return(
    ((b*g - b + 1)*p/(b*g - b*r - b + 1) - ((b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p*cosh(sqrt((b*g - b + 1)*b*r)*t/(b*g - b + 1))/(b*g - b*r - b + 1) - ((b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*(b*g - b + 1)*p/(b*g - b*r - b + 1) - ((b^2*g - b^2 + b)*p*r + (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p)*(b*g - b + 1)/(b*g - b*r - b + 1))*sinh(sqrt((b*g - b + 1)*b*r)*t/(b*g - b + 1))/sqrt((b*g - b + 1)*b*r))*exp(-t)/(b*g - b + 1))/((b*g - b + 1)*p/(b*g - b*r - b + 1) - ((b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p*cosh(sqrt((b*g - b + 1)*b*r)*t/(b*g - b + 1))/(b*g - b*r - b + 1) - ((b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*(b*g - b + 1)*p/(b*g - b*r - b + 1) - ((b^2*g - b^2 + b)*p*r + (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p)*(b*g - b + 1)/(b*g - b*r - b + 1))*sinh(sqrt((b*g - b + 1)*b*r)*t/(b*g - b + 1))/sqrt((b*g - b + 1)*b*r))*exp(-t)/(b*g - b + 1) + (((b^2*g^2 + b^2 - 2*(b^2 - b)*g - (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p + (b^2*g^2 - (b^2 - b)*g - (b^2*g^2 - (b^2 - b)*g)*p)*r - 2*b + 1)*(b*g - b + 1)/(b*g*r - b*g + b - 1) - (b^2*g^2 + b^2 - 2*(b^2 - b)*g - (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p - 2*b + 1)*(b*g - b + 1)/(b*g*r - b*g + b - 1))*sinh(sqrt((b*g - b + 1)*b*g*r)*t/(b*g - b + 1))/sqrt((b*g - b + 1)*b*g*r) + (b^2*g^2 + b^2 - 2*(b^2 - b)*g - (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p - 2*b + 1)*cosh(sqrt((b*g - b + 1)*b*g*r)*t/(b*g - b + 1))/(b*g*r - b*g + b - 1))*exp(-t)/(b*g - b + 1) - (b*g - (b*g - b + 1)*p - b + 1)/(b*g*r - b*g + b - 1))
  )
}

# Add the calculated female share after 20 rounds for each case to the data file
cases$ans.exit8 <- sapply(1:3,function(c){rarPctFMaxExit(cases$p[c],cases$r[c],cases$g[c],cases$b[c],t=20)})


#################################################
# 6. Model interpretation
# --------------------------------------------
# 6.1. Sex bias equivalent
#      How much of a hiring bias for men over women is needed to replicate the observed segregating effects from reapplication differences?

bias <- function (t, state, parms) {
  with(as.list(c(state, parms)), {
    dMAcc <- (s*a/(1+(s*a)-a))*(1-pctF)
    dWAcc <- (a/((1+(s*a)-a)))*pctF
    list(c(dMAcc,dWAcc))
  })
}
# men numerator: s*a
# women numerator: a
# a range: 0..1
# s range: 1..inf

sbu <- function(pRej,startPctF,finalPctF){
  pSbu <- 1
  pctF <- startPctF
  while (pctF > finalPctF){
    pSbu <- pSbu + 0.0001
    parms_sbu <- c(s=pSbu,a=(1-pRej),pctF=startPctF) 
    ini_sbu <- c(MAcc=(1-startPctF),WAcc=startPctF)
    times_sbu <- seq(0,20,0.01)
    out_sbu <- ode(y=ini_sbu,time=times_sbu,func=bias,parms=parms_sbu)
    pctF <- out_sbu[2000,"WAcc"]/(out_sbu[2000,"WAcc"]+out_sbu[2000,"MAcc"])
  }
  return(pSbu)
}

# sex bias parameter s, generating the equivalent segregating effect over the same period of time for the 3 cases
cases$sbu <- sapply(1:3,function(i){sbu(cases$r[i],cases$p[i],cases$sol.exit0[i])})

# sex bias interpretation: Acceptance rate for men (with a different rate for women) generating the equivalent segregating effect over the same period of time for the 3 cases, and with the same overall acceptance rate as the case
cases$accRateMen <- sapply(1:(dim(cases)[1]),function(i){
  return((1-cases$r[i])*cases$sbu[i]/(1+((1-cases$r[i])*cases$sbu[i])-(1-cases$r[i])))
})
# sex bias interpretation: Acceptance rate for women (with a different rate for men) generating the equivalent segregating effect over the same period of time for the 3 cases, and with the same overall acceptance rate as the case
cases$accRateWom <- sapply(1:(dim(cases)[1]),function(i){
  return((1-cases$r[i])/(cases$sbu[i]*(1+((1-cases$r[i])*cases$sbu[i])-(1-cases$r[i]))))
})

# sex bias interpretation: Odds ratio of the gender-specific acceptance rates
cases$sbuOdds <- sapply(1:(dim(cases)[1]),function(i){  # Odds Ratio of sex differences in acceptance rate
  return(cases$accRateMen[i]*(1-cases$accRateWom[i])/(cases$accRateWom[i]*(1-cases$accRateMen[i])))
})

# 6.2.   Evaluating interventions - expanded pools
#        A FIXED SELECTION OPPORTUNITY CONTEXT
#        how much to get back to p? (r IS changing, because fixed # of spots)
#        find p, given r,b,g,and p_eq, to a specified number of digits
#        PLUS exp - expansion factor of increased applicant pool
#        To model a selection context that is NOT fixed, x=1, as the rejection rate is not affected by the expanded pool

# That function uses the 200th round pct female among applicants as a proxy for the equilibrium composition
# These choices are to simplify and to speed the computation of terms.
# The validity of these choices are tested below using the equilibrium analysis functions, and give identical results (albeit requiring much more time)
# Number of digits by default is 5, but may be changed
getThreshP <- function(r,b,g,ss,x,digits=5){
  p_test <- ss
  epsilon <- (10^(-1*digits))/20
  while (abs(ss - rarPctFMaxExit(p_test,((x+r-1)/x),g,b,t=200)) > epsilon) {
    if (rarPctFMaxExit(p_test,((x+r-1)/x),g,b,t=200) < ss) {
      p_test <- p_test*1.1
    }
    else {
      p_test <- p_test*0.9
    }
  }
  return(round(p_test,digits))
}

# The threshold female share of new 1st-time applicants needed to yield a long-term female share of 13.8% (equilibrium from Exec Search Case)
getThreshP(0.96,0.688,1.098,0.138,1,4)
# [1] 0.1599307
# Interpretation:
# MORE THAN 16% female new 1st-time applicants are needed to move the long-term female share above 13.8%

# The threshold female share of new 1st-time applicants needed to yield a long-term female share of 20% (hypothetical example in paper)
getThreshP(0.96,0.688,1.098,0.2,1,4)
# [1] 0.2291646
# Interpretation:
# MORE THAN 22.9% female new 1st-time applicants are needed to move the long-term female share above 20%

# This version uses the same ode solving process as above, and only 100 rounds (it is much slower), and produces identical results
getThreshP2 <- function(r,b,g,ss,x,digits=5){
  p_test <- ss
  epsilon <- (10^(-1*digits))/20
  p_test.out <- appModel(gParam=g,bParam=b,pParam=p_test,rParam=((x+r-1)/x),rounds=100)
  while (abs(ss - p_test.out) > epsilon) {
    p_test.out <- appModel(gParam=g,bParam=b,pParam=p_test,rParam=((x+r-1)/x),rounds=100)
    if (p_test.out < ss) {
      p_test <- p_test*1.1
    }
    else {
      p_test <- p_test*0.9
    }
    p_test.out <- appModel(gParam=g,bParam=b,pParam=p_test,rParam=((x+r-1)/x),rounds=100)
  }
  return(p_test)
}
# Same results as above, but much slower
# getThreshP2(0.96,0.688,1.098,0.138,1,4) # - commented out because requires a long time
# [1] 0.1599307
# getThreshP2(0.96,0.688,1.098,0.2,1,4) # - commented out because requires a long time
# [1] 0.2291646


#############
# 7. Figures
# ----------
#    7.1. 4-panel Figure 2 - Change in female share by time, r, b, g parameters

bRange <- sort(cases$b) # values of the 3 cases, changed to continuum for last panel
gRange <- (100:155)/100        # 1...1.6 by steps of 0.01
pRange <- sort(cases$p)      # values of the 3 cases
rRange <- ((1:33)*3-1)/100     # 0.02..0.98 by steps of 0.03 
tRange <- 1:200              # focus is t=20

# Plotting as a function of t - time
png(width=10,height=7,units="in",res=300,filename=paste0("Fig2",gsub("-","",Sys.Date()),".PNG"))
par(mar=c(4,5,1,1)+0.1,mfrow=c(2,2))
plot(x=tRange,y=rep(0,length(tRange)),ylim=c(0.0,0.32),col="white",
     ylab="Female Share (of Selected or Applicants)",
     xlab="time, in terms of application and reapplication rounds",cex.axis=1,cex.lab=1,yaxt="n")
yLabels <- seq(0, 0.3, 0.05)
axis(2, at=yLabels, labels=sprintf(round(100*yLabels), fmt="%2.0f%%"), cex.axis=1,cex.lab=1)
for(c in 1:3){
  lines(x=tRange,y=sapply(tRange,
                          function(timex){rarPctFMaxExit(cases$p[c],cases$r[c],cases$g[c],cases$b[c],t=timex)}),col=gray(0.15*c),lty=3)
  abline(h=cases$p[c],lty=3,col=gray(0.2),lwd=2)
  lines(x=rep(20,2),
        y=c(cases$p[c],
            rarPctFMaxExit(cases$p[c],cases$r[c],cases$g[c],cases$b[c],t=20)),
        lwd=4,col=gray(0.15*c))
  lines(x=tRange,y=sapply(tRange,
                          function(timex){rarPctFMinExit(cases$p[c],cases$r[c],cases$g[c],cases$b[c],t=timex)}),col=gray(0.15*c),lwd=2)
  arrows(20, rarPctFMinExit(cases$p[c],cases$r[c],cases$gub[c],cases$b[c],t=20), 
         20, rarPctFMinExit(cases$p[c],cases$r[c],cases$glb[c],cases$b[c],t=20), length=0.05, angle=90, code=3)
}
text(x=40, y=0.17,labels="Executive Search {p=0.16,r=0.96,b=0.69,g=1.10}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))
text(x=40,y=0.30,labels="Crowdfunding {p=0.29,r=0.64,b=0.11,g=1.26}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))
text(x=40,y=0.1,labels="Patenting {p=0.088,r=0.86,b=0.54,g=1.05}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))

# Plotting as a function of g - gender difference in RAR
par(mar=c(4,2,1,1))
plot(x=gRange,y=rep(0,length(gRange)),ylim=c(0.0,0.32),col="white",
     ylab="",
     xlab="Gender difference in reapplication after rejection parameter, g",cex.axis=1,cex.lab=1,yaxt="n")
yLabels <- seq(0, 0.3, 0.05)
axis(2, at=yLabels, labels=sprintf(round(100*yLabels), fmt="%2.0f%%"), cex.axis=1,cex.lab=1)
for(c in 1:3){
  lines(x=gRange,y=sapply(gRange,
                          function(gp){rarPctFMaxExit(cases$p[c],cases$r[c],gp,cases$b[c],t=20)}),col=gray(0.15*c),lty=3)
  abline(h=cases$p[c],lty=3,col=gray(0.2),lwd=2)
  polygon(x=c(rep(cases$glb[c],2),rep(cases$gub[c],2)),
          y=c(cases$p[c],
              rarPctFMaxExit(cases$p[c],cases$r[c],cases$glb[c],cases$b[c],t=20),
              rarPctFMaxExit(cases$p[c],cases$r[c],cases$gub[c],cases$b[c],t=20),
              cases$p[c]),
          lwd=1,col=gray(0.45 + 0.15*c))
  lines(x=rep(cases$g[c],2),
        y=c(cases$p[c],
            rarPctFMaxExit(cases$p[c],cases$r[c],cases$g[c],cases$b[c],t=20)),
        lwd=4,col=gray(0.15*c))
  lines(x=gRange,y=sapply(gRange,
                          function(gp){rarPctFMinExit(cases$p[c],cases$r[c],gp,cases$b[c],t=20)}),col=gray(0.15*c),lwd=2)
}
text(x=1.01, y=0.17,labels="Executive Search {b=0.69,p=0.16,r=0.96}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))
text(x=1.01,y=0.30,labels="Crowdfunding {b=0.11,p=0.29,r=0.64}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))
text(x=1.01,y=0.1,labels="Patenting {b=0.54,p=0.088,r=0.86}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))

# Plotting as a function of r - rejection rate
# png(width=960,height=640,filename=paste0("rParam",gsub("-","",Sys.Date()),".PNG"))
# par(mar=c(5,6,4,2)+0.1)
par(mar=c(4,5,1,1))
plot(x=rRange,y=rep(0,length(rRange)),ylim=c(0.0,0.32),col="white",ylab="Female Share (of Selected or Applicants)",xlab="Overall rejection rate parameter, r",cex.axis=1,cex.lab=1,yaxt="n")
yLabels <- seq(0, 0.3, 0.05)
axis(2, at=yLabels, labels=sprintf(round(100*yLabels), fmt="%2.0f%%"), cex.axis=1,cex.lab=1)
for(c in 1:3){
  lines(x=rRange,y=sapply(rRange,
                          function(rp){rarPctFMinExit(cases$p[c],rp,cases$g[c],cases$b[c],t=20)}),col=gray(0.15*c),lwd=2)
  lines(x=rRange,y=sapply(rRange,
                          function(rp){rarPctFMaxExit(cases$p[c],rp,cases$g[c],cases$b[c],t=20)}),col=gray(0.15*c),lty=3)
  abline(h=cases$p[c],lty=3,col=gray(0.2),lwd=2)
  lines(x=rep(cases$r[c],2),
        y=c(cases$p[c],
            rarPctFMaxExit(cases$p[c],cases$r[c],cases$g[c],cases$b[c],t=20)),
        lwd=4,col=gray(0.15*c))
  arrows(cases$r[c], rarPctFMinExit(cases$p[c],cases$r[c],cases$gub[c],cases$b[c],t=20), cases$r[c], rarPctFMinExit(cases$p[c],cases$r[c],cases$glb[c],cases$b[c],t=20), length=0.05, angle=90, code=3)
}
text(x=0.1, y=0.17,labels="Executive Search {b=0.69,p=0.16,r=0.96}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))
text(x=0.1,y=0.30,labels="Crowdfunding {b=0.11,p=0.29,r=0.64}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))
text(x=0.1,y=0.1,labels="Patenting {b=0.54,p=0.088,r=0.86}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))

# Plotting as a function of b - baseline RAR
bRange <- (1:99)/100

par(mar=c(4,2,1,1))
plot(x=bRange,y=rep(0,length(bRange)),ylim=c(0.0,0.32),col="white",
     ylab="",
     xlab="Baseline reapplication after rejection, b",cex.axis=1,cex.lab=1,yaxt="n")
yLabels <- seq(0, 0.3, 0.05)
axis(2, at=yLabels, labels=sprintf(round(100*yLabels), fmt="%2.0f%%"), cex.axis=1,cex.lab=1)
for(c in 1:3){
  lines(x=bRange,y=sapply(bRange,
                          function(bp){rarPctFMinExit(cases$p[c],cases$r[c],cases$g[c],bp,t=20)}),col=gray(0.15*c),lwd=2)
  lines(x=bRange,y=sapply(bRange,
                          function(bp){rarPctFMaxExit(cases$p[c],cases$r[c],cases$g[c],bp,t=20)}),col=gray(0.15*c),lty=3)
  abline(h=cases$p[c],lty=3,col=gray(0.2),lwd=2)
  lines(x=rep(cases$b[c],2),
        y=c(cases$p[c],
            rarPctFMaxExit(cases$p[c],cases$r[c],cases$g[c],cases$b[c],t=20)),
        lwd=4,col=gray(0.15*c))
  arrows(cases$b[c], rarPctFMinExit(cases$p[c],cases$r[c],cases$gub[c],cases$b[c],t=20), 
         cases$b[c], rarPctFMinExit(cases$p[c],cases$r[c],cases$glb[c],cases$b[c],t=20), length=0.05, angle=90, code=3)
}
text(x=0.1, y=0.17,labels="Executive Search {p=0.16,r=0.96,g=1.10}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))
text(x=0.1,y=0.30,labels="Crowdfunding {p=0.29,r=0.64,g=1.26}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))
text(x=0.1,y=0.1,labels="Patenting {p=0.088,r=0.86,g=1.05}",col=gray(0.2),cex.axis=1,cex.lab=1,cex=1,adj = c(0, NA))
dev.off()


# ----------
#    7.2. 3-panel Figure 3 - Change in threshold (Exec Search Case, steady-state=0.138) w/change in r, b, g parameters

steadyState <- 0.138
fig3 <- sapply(1+((0:100)/100),function(ex){getThreshP(cases$r[1],cases$b[1],cases$g[1],steadyState,ex)})
fig3_80r <- sapply(1+((0:100)/100),function(ex){getThreshP(0.8*cases$r[1],cases$b[1],cases$g[1],steadyState,ex)})
fig3_60r <- sapply(1+((0:100)/100),function(ex){getThreshP(0.6*cases$r[1],cases$b[1],cases$g[1],steadyState,ex)})

fig3_80g <- sapply(1+((0:100)/100),function(ex){getThreshP(cases$r[1],cases$b[1],1+((cases$g[1]-1)*0.8),steadyState,ex)})
fig3_60g <- sapply(1+((0:100)/100),function(ex){getThreshP(cases$r[1],cases$b[1],1+((cases$g[1]-1)*0.6),steadyState,ex)})

fig3_80b <- sapply(1+((0:100)/100),function(ex){getThreshP(cases$r[1],0.8*cases$b[1],cases$g[1],steadyState,ex)})
fig3_60b <- sapply(1+((0:100)/100),function(ex){getThreshP(cases$r[1],0.6*cases$b[1],cases$g[1],steadyState,ex)})


png(width=1500,height=2400,filename=paste0("Fig3_gbr_",gsub("-","",Sys.Date()),".PNG"),res=300)
par(mar=c(3,3,1,0)+0.1,mfrow=c(3,1))
# g first
# PLOTTING CHANGES in g - gender difference in reapplication
# png(width=1200,height=640,filename=paste0("xParam_g",gsub("-","",Sys.Date()),".PNG"))
# par(mar=c(5,6,4,2)+0.1)
plot(x=1+((0:100)/100),
     y=fig3,
     ylim=c(0.138,0.162),type="l",
     ylab="Threshold %Female", xlab="expanded pool relative size, x",cex.axis=0.7,cex.lab=0.7,yaxt="n",lwd=1.5)
mtext("expanded pool relative size, x",side=1,line=2,cex=0.7)
mtext("Threshold %Female",side=2,line=2,cex=0.7)
yLabels <- seq(0.138, 0.162, 0.002)
axis(2, at=yLabels, labels=sprintf(round(100*yLabels,1), fmt="%2.1f%%"), cex.axis=0.7,cex.lab=0.7)
lines(x=c(1,2),y=rep(0.138,2),lty=3,col="darkGray")
lines(x=c(1,2),y=rep(0.16,2),lty=3,col="darkGray")

lines(x=1+((0:100)/100),
      y=fig3_80g,lwd=1.5,lty=5,col="black")
lines(x=1+((0:100)/100),
      y=fig3_60g,lwd=1.5,lty=2,col="black")
text(x=1.05,y=0.1545,labels=paste0("80% of original gender difference in RAR, g=",round(1+((cases$g[1]-1)*0.8),2)),cex=0.7,adj=0)
text(x=1.4,y=0.15,labels=paste0("60% of original gender difference in RAR, g=",round(1+((cases$g[1]-1)*0.6),2)),cex=0.7,adj=0)
text(x=1,y=0.1615,labels="Rising threshold as pool expands",cex=0.7,adj=0)
text(x=1.8,y=0.139,labels="Steady state, 13.8%",cex=0.7,adj=0)
text(x=1.4,y=0.159,labels="Threshold, p=16% for steady state=13.8%, at original pool size",cex=0.7,adj=0)

# PLOTTING CHANGES in b - baseline reapplication rate
# png(width=1200,height=640,filename=paste0("xParam_b",gsub("-","",Sys.Date()),".PNG"))
# par(mar=c(5,6,4,2)+0.1)
plot(x=1+((0:100)/100),
     y=fig3,
     ylim=c(0.138,0.162),type="l",
     ylab="Threshold %F", xlab="expanded pool relative size, x",cex.axis=0.7,cex.lab=0.7,yaxt="n",lwd=1.5)
mtext("expanded pool relative size, x",side=1,line=2,cex=0.7)
mtext("Threshold %Female",side=2,line=2,cex=0.7)
yLabels <- seq(0.138, 0.162, 0.002)
axis(2, at=yLabels, labels=sprintf(round(100*yLabels,1), fmt="%2.1f%%"), cex.axis=0.7,cex.lab=0.7)
lines(x=c(1,2),y=rep(0.138,2),lty=3,col="darkGray")
lines(x=c(1,2),y=rep(0.16,2),lty=3,col="darkGray")
lines(x=1+((0:100)/100),
      y=fig3_80b,lwd=1.5,lty=5,col="black")
lines(x=1+((0:100)/100),
      y=fig3_60b,lwd=1.5,lty=2,col="black")
text(x=1.05,y=0.1525,labels=paste0("80% of original base reapplication rate, b=",round(0.8*cases$b[1],2)),cex=0.7,adj=0)
text(x=1.4,y=0.144,labels=paste0("60% of original base reapplication rate, b=",round(0.6*cases$b[1],2)),cex=0.7,adj=0)
text(x=1,y=0.1615,labels="Rising threshold as pool expands",cex=0.7,adj=0)
text(x=1.8,y=0.139,labels="Steady state, 13.8%",cex=0.7,adj=0)
text(x=1.4,y=0.159,labels="Threshold, p=16% for steady state=13.8%, at original pool size",cex=0.7,adj=0)

# PLOTTING CHANGES in r - baseline rejection rate
plot(x=1+((0:100)/100),
     y=fig3,
     ylim=c(0.138,0.162),type="l",
     ylab="Threshold %Female", xlab="expanded pool relative size, x",cex.axis=0.7,cex.lab=0.7,yaxt="n",lwd=1.5)
mtext("expanded pool relative size, x",side=1,line=2,cex=0.7)
mtext("Threshold %Female",side=2,line=2,cex=0.7)
yLabels <- seq(0.138, 0.162, 0.002)
axis(2, at=yLabels, labels=sprintf(round(100*yLabels,1), fmt="%2.1f%%"), cex.axis=0.7,cex.lab=0.7)

lines(x=c(1,2),y=rep(0.138,2),lty=3,col="darkGray")
lines(x=c(1,2),y=rep(0.16,2),lty=3,col="darkGray")

lines(x=1+((0:100)/100),
      y=fig3_80r,lwd=1.5,lty=5,col="black")
lines(x=1+((0:100)/100),
      y=fig3_60r,lwd=1.5,lty=2,col="black")
text(x=1.05,y=0.154,labels=paste0("80% of original rejection rate, r=",round(0.8*cases$r[1],2)),cex=0.7,adj=0)
text(x=1.5,y=0.1475,labels=paste0("60% of original rejection rate, r=",round(0.6*cases$r[1],2)),cex=0.7,adj=0)
text(x=1,y=0.1615,labels="Rising threshold as pool expands",cex=0.7,adj=0)
text(x=1.8,y=0.139,labels="Steady state, 13.8%",cex=0.7,adj=0)
text(x=1.4,y=0.159,labels="Threshold, p=16% for steady state=13.8%, at original pool size",cex=0.7,adj=0)

dev.off()

##########
# 8. Table
# --------

tab2 <- round(rbind(t(cases[,c("p","r","b","g")]),
                    g_ci = ((cases$gub-cases$g)+(cases$g-cases$glb))/2,
                    after20 = cases$sol.exit0,
                    after20ci = ((rarPctFMinExit(cases$p,cases$r,cases$glb,cases$b,t=20) - cases$sol.exit0)+
                                   (cases$sol.exit0-rarPctFMinExit(cases$p,cases$r,cases$gub,cases$b,t=20)))/2,
                    delta_abs = cases$sol.exit0-cases$p,
                    delta_rel = (cases$sol.exit0-cases$p)/cases$p,
                    t(cases[,c("selEquil","sbu","accRateMen","accRateWom","sbuOdds")]))[,c(1,3,2)],4)
colnames(tab2) <- cases$caseID[c(1,3,2)]
tab2
