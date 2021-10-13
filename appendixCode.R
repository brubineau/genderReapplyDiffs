# Code for Paper

rm(list=ls())
gc()
library(deSolve)

# Differential Equation Model
# equation taken from Table 1

model <- function (t, state, parms) {
  with(as.list(c(state, parms)), {
    #  Formulae for Male/Men (M) and Female/Women (W) as in the Paper
    dAM <- -AM + 1 - p + RM*b*g/(1+(b*g)-b) # Male applicants
    dAW <- -AW + p + RW * b                 # Female applicants
    dRM <- r*AM - RM                        # Rejected Men
    dRW <- r*AW - RW                        # Rejected Women
    list(c(dAM,dAW,dRM,dRW))
  })
}

# Initial values for the continuous model - all stocks set to zero
ini <- c(AM=0, AW=0, RM=0, RW=0)

#####################################################
# Numeric solution
# ----------------
# appModel Function numericaly solves the continuous model
# using the:
#   specified model parameter values;
#   provided initial conditions; 
#   designated number of rounds; and
# returns the final female share among applicants.

appModel <- function(rounds,pParam,rParam,bParam,gParam){
  parmsOpen <- c(p=pParam,r=rParam,b=bParam,g=gParam) # specified model parametres
  outOpen <- ode(y = ini,                    # initial conditions, 
                 times = seq(0,rounds,0.01), # time horizon for model
                 func = model,               # call to model
                 parms = parmsOpen)          # parameters
  # return the final women's share among all applicants
  return(outOpen[dim(outOpen)[1],"AW"]/(outOpen[dim(outOpen)[1],"AW"]+outOpen[dim(outOpen)[1],"AM"]))
}

#####################################################
# Analytic solution
# -----------------
# women in applicant pool at time t
wT <- function(t,p,r,b){
  return((b*p*r*cosh(sqrt(b*r)*t)/(b*r - 1) + b*p*r*sinh(sqrt(b*r)*t)/((b*r - 1)*sqrt(b*r)))*exp(-t) - p/(b*r - 1))  
}
# men in applicant pool at time t
mT <- function(t,p,r,b,g){
  return(((b^2*g^2 - (b^2 - b)*g - (b^2*g^2 - (b^2 - b)*g)*p)*(b*g - b + 1)*r*sinh(sqrt((b*g - b + 1)*b*g*r)*t/(b*g - b + 1))/(sqrt((b*g - b + 1)*b*g*r)*(b*g*r - b*g + b - 1)) + (b^2*g^2 - (b^2 - b)*g - (b^2*g^2 - (b^2 - b)*g)*p)*r*cosh(sqrt((b*g - b + 1)*b*g*r)*t/(b*g - b + 1))/(b*g*r - b*g + b - 1))*exp(-t)/(b*g - b + 1) - (b*g - (b*g - b + 1)*p - b + 1)/(b*g*r - b*g + b - 1))
}
# women's share among applicants at time t
wShareT <- function(t,p,r,b,g){
  return(wT(t,p,r,b)/(wT(t,p,r,b)+mT(t,p,r,b,g)))
}


eqT <- function(p,r,b,g){
  eqt <- log(((sqrt(b*r)*b*g*r - ((sqrt(b*r)*r^2 - sqrt(b*r)*r)*g^2 + sqrt(b*r)*g*r)*b^2)*cosh(sqrt(b*r)*t)*sinh(sqrt((b^2*g^2 - (b^2 - b)*g)*r)*t/(b*g - b + 1)) + ((sqrt((b^2*g^2 - (b^2 - b)*g)*r)*b^2*g*r^2 - sqrt((b^2*g^2 - (b^2 - b)*g)*r)*b*g*r)*cosh(sqrt((b^2*g^2 - (b^2 - b)*g)*r)*t/(b*g - b + 1)) + ((g^2*r^2 - g*r^2)*b^3 - (g^2*r^2 - g*r^2)*b^2)*sinh(sqrt((b^2*g^2 - (b^2 - b)*g)*r)*t/(b*g - b + 1)))*sinh(sqrt(b*r)*t))/(((sqrt((b^2*g^2 - (b^2 - b)*g)*r)*g*r - sqrt((b^2*g^2 - (b^2 - b)*g)*r)*r)*b^2 - (sqrt((b^2*g^2 - (b^2 - b)*g)*r)*g - sqrt((b^2*g^2 - (b^2 - b)*g)*r)*r - sqrt((b^2*g^2 - (b^2 - b)*g)*r))*b - sqrt((b^2*g^2 - (b^2 - b)*g)*r))*sinh(sqrt(b*r)*t) - (((sqrt(b*r)*r - sqrt(b*r))*g^2 + sqrt(b*r)*g)*b - sqrt(b*r)*g)*sinh(sqrt((b^2*g^2 - (b^2 - b)*g)*r)*t/(b*g - b + 1))))
  return(eqt)
}

#####################################################
# Illustrating equivalence of both solution methods
# -------------------------------------------------

# function compareSols takes no arguments
# it samples a random point in the model's parameter space (including time)
# and finds the women's share among applicants at that 
compareSols <- function(){
  # pick random values for model parameters within finite ranges
  pRand <- runif(1,min=0.01,max=0.49) # women's share of first-time applicants
  rRand <- runif(1,min=0.01,max=0.99) # rejection rate
  bRand <- runif(1,min=0.01,max=0.99) # baseline (women's) reapplication rate
  gRand <- runif(1,min=1,max=10)      # odds ratio 
  tRand <- runif(1,min=10,max=1000)   # time 
  numSol <- appModel(tRand,pRand,rRand,bRand,gRand) # get the outcome - women's share of applicants - via the numeric solution
  frmSol <- wShareT(tRand,pRand,rRand,bRand,gRand)  # get the same outcome from the formal analytic solution
  return(c(time=tRand,diff=as.numeric(numSol - frmSol))) # return the difference between the 2 solutions
}

set.seed(100) # to reproduce shared results, use the same random seed
              # to generate new results, set a new seed or omit the prior line
diff1k <- replicate(1000,compareSols()) # Note - this may require close to an hour to complete
mean(diff1k[2,],na.rm=T) #  2.0 * 10^-6
max(diff1k[2,],na.rm=T)  #  6.7 * 10^-4
min(diff1k[2,],na.rm=T)  # -3.7 * 10^-15
plot(density(diff1k[2,],na.rm=T)) # note almost all results are indistinguishable from zero
plot(diff1k[1,],diff1k[2,],xlab="model time value (t)",ylab="Numeric - Analytic Result",type="p")


# time to equilibrium
time2Eq <- function(maxTime){
  # pick random values for model parameters within finite ranges
  pRand <- runif(1,min=0.01,max=0.49) # women's share of first-time applicants
  rRand <- runif(1,min=0.01,max=0.99) # rejection rate
  bRand <- runif(1,min=0.01,max=0.99) # baseline (women's) reapplication rate
  gRand <- runif(1,min=1,max=10)      # gender difference odds ratio in reapplication
  tRand <- runif(1,min=10,max=maxTime)   # time 
  frmSol <- wShareT(tRand,pRand,rRand,bRand,gRand)  # get the same outcome from the formal analytic solution
return(c(rParam=rRand,bParam=bRand,gParam=gRand,pParam=pRand,time=tRand,out=frmSol,out100=wShareT(tRand+100,pRand,rRand,bRand,gRand))) # return the difference between the 2 solutions
}

time2Eq1kMax10k <- replicate(10000,time2Eq(10000))
plot(time2Eq1kMax10k[5,],time2Eq1kMax10k[6,]-time2Eq1kMax10k[7,])
time2Eq1kMax1k <- replicate(10000,time2Eq(1000))
plot(time2Eq1kMax1k[5,],time2Eq1kMax1k[6,]-time2Eq1kMax1k[7,])

wShareEq <- function(p,r,b,g){
  return(wT(500,p,r,b)/(wT(500,p,r,b)+mT(500,p,r,b,g)))
}


# Sex Bia Equivalent
getSBEWO <- function(pParam,rParam,wParam,oParam){
  gParam <- oParam/(1-wParam+(oParam*wParam))
  bParam <- wParam/(1+wParam-(gParam*wParam))
  q <- equilPctFWO(pParam,rParam,wParam,oParam)
  r <- rParam
  p <- pParam
  return(c(oddsR = round(((r*q) - q + p)*(1-q) / (q *(q - (q*r) + r - p)),4),
           accM=round((r-1)*(q-1)/(1-p),4),
           accW=round((q-(q*r))/p,4)))
}
