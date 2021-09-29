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
library(plotly)
library(lattice)
library(ContourFunctions)

#    1.3. Change to working directory
# setwd("~/") # UN-COMMENT & ADAPT FOR YOUR OWN COMPUTING CONTEXT

############################
# 3. MAIN MODEL AND VARIANT
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

##################################################
# 4. Find equilibrium values for model and variant
# ------------------------------------------------
# Equilibrium value is the outcome value 
# (i.e., female share among selected - main model, or female share among applicants - variant)
# as time goes to infinity. For our purposes, we specify the number of digits desired, 
# and take the outcome once the model produces identical values 50 rounds apart to the specified number of digits in two successive rounds.
# we start at 500 rounds

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

rarPctFMaxExit <- function(p,r,g,b,t){
  return(
    ((b*g - b + 1)*p/(b*g - b*r - b + 1) - ((b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p*cosh(sqrt((b*g - b + 1)*b*r)*t/(b*g - b + 1))/(b*g - b*r - b + 1) - ((b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*(b*g - b + 1)*p/(b*g - b*r - b + 1) - ((b^2*g - b^2 + b)*p*r + (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p)*(b*g - b + 1)/(b*g - b*r - b + 1))*sinh(sqrt((b*g - b + 1)*b*r)*t/(b*g - b + 1))/sqrt((b*g - b + 1)*b*r))*exp(-t)/(b*g - b + 1))/((b*g - b + 1)*p/(b*g - b*r - b + 1) - ((b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p*cosh(sqrt((b*g - b + 1)*b*r)*t/(b*g - b + 1))/(b*g - b*r - b + 1) - ((b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*(b*g - b + 1)*p/(b*g - b*r - b + 1) - ((b^2*g - b^2 + b)*p*r + (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p)*(b*g - b + 1)/(b*g - b*r - b + 1))*sinh(sqrt((b*g - b + 1)*b*r)*t/(b*g - b + 1))/sqrt((b*g - b + 1)*b*r))*exp(-t)/(b*g - b + 1) + (((b^2*g^2 + b^2 - 2*(b^2 - b)*g - (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p + (b^2*g^2 - (b^2 - b)*g - (b^2*g^2 - (b^2 - b)*g)*p)*r - 2*b + 1)*(b*g - b + 1)/(b*g*r - b*g + b - 1) - (b^2*g^2 + b^2 - 2*(b^2 - b)*g - (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p - 2*b + 1)*(b*g - b + 1)/(b*g*r - b*g + b - 1))*sinh(sqrt((b*g - b + 1)*b*g*r)*t/(b*g - b + 1))/sqrt((b*g - b + 1)*b*g*r) + (b^2*g^2 + b^2 - 2*(b^2 - b)*g - (b^2*g^2 + b^2 - 2*(b^2 - b)*g - 2*b + 1)*p - 2*b + 1)*cosh(sqrt((b*g - b + 1)*b*g*r)*t/(b*g - b + 1))/(b*g*r - b*g + b - 1))*exp(-t)/(b*g - b + 1) - (b*g - (b*g - b + 1)*p - b + 1)/(b*g*r - b*g + b - 1))
  )
}

appEquil2 <- function(digits,gParam,bParam,pParam,rParam){
  # define equilibrium as value when result is the same after 1000 rounds
  equil <- NA
  equilNotFound <- T
  rounds <- 1000
  outcome <- round(rarPctFMaxExit(p=pParam,r=rParam,g=gParam,b=bParam,t=rounds),digits+1)
  if (abs(outcome-round(rarPctFMaxExit(p=pParam,r=rParam,g=gParam,b=bParam,t=rounds+1000),digits+1)) < (1/(10^(digits+1)))) { # test for identical outcomes within specified rounding level
    equilNotFound <- F
    equil <- outcome
  }
  while (equilNotFound){
    rounds <- rounds+1000
    outcome <- round(rarPctFMaxExit(p=pParam,r=rParam,g=gParam,b=bParam,t=rounds),digits+1)
    if (abs(outcome-round(rarPctFMaxExit(p=pParam,r=rParam,g=gParam,b=bParam,t=rounds+1000),digits+1))<(1/(10^(digits+1)))) {
      equilNotFound <- F
      equil <- outcome
    }
  }
  return(round(equil,digits))
}

# analytic solution?
# share of return applicants
fShare <- function(t,p,r,b,g){
  k <- r*b/(1+(b*g)-b)
  return((exp((t+log(k*p))*(k-1))-p)/(k-1))
}
mShare <- function(t,p,r,b,g){
  k <- r*b/(1+(b*g)-b)
  return((exp((t+(log(g*k-g*p*k)/(g*k-1)))*(g*k-1))-(1-p))/(g*k-1))
}

equilPctF <- function(p,r,b,g){
#  fShare <- (exp(100*((r*b/(1+(b*g)-b))-1))-p)/((r*b/(1+(b*g)-b))-1)
#  mShare <- (exp(100*((r*b*g/(1+(b*g)-b))-1))-(1-p))/((r*b*g/(1+(b*g)-b))-1)
  return(fShare(1000,p,r,b,g)/(fShare(1000,p,r,b,g)+mShare(1000,p,r,b,g)))
}
# > equilPctF(0.3,0.8,0.7,1.2)
# [1] 0.2569546
# > appEquil(digits=6,gParam=1.2,bParam=0.7,pParam=0.3,rParam=0.8)
# [1] 0.256955

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

findSbuEquil <- function(startPctF,rParam,sParam){
  ini_sbu <- c(MAcc=(1-startPctF),WAcc=startPctF)
  parms_sbu <- c(s=sParam,a=(1-rParam),pctF=startPctF) 
  equilNotFound <- T
  equil <- NA
  rounds <- 500
  out <- ode(y = ini_sbu, times = seq(0,rounds,0.01), func = bias, parms = parms_sbu) # initial conditions, times, and parameters
  outcomes <- round(out[,"WAcc"]/(out[,"WAcc"]+out[,"MAcc"]),4)
  if (sum(as.numeric(outcomes[2:(length(outcomes)-5000)]==outcomes[5002:(length(outcomes))]))>0) { # test for identical outcomes 50 rounds apart (steps are 1/100th of a round)
    equilNotFound <- F
    equil <- outcomes[length(outcomes)]
  }
  while (equilNotFound){
    newIni <- out[dim(out)[1],]
    newTimes <- seq(rounds,rounds+500,0.01)
    rounds <- rounds+500
    out <- ode(y = newIni, times = newTimes, func = bias, parms = parms_sbu) # initial conditions, times, and parameters
    outcomes <- round(out[,"WApp"]/(out[,"WApp"]+out[,"MApp"]),digits+1)
    if (sum(as.numeric((outcomes[2:(length(outcomes)-5000)]-outcomes[5002:(length(outcomes))])==0))>0) {
      equilNotFound <- F
      equil <- outcomes[length(outcomes)]
    }
  }
  return(c(equil=round(equil,4),steps=rounds))
}

findS <- function(pParam,rParam,target){
  pSbu <- 1
  sFound <- F
  sEquil <- as.numeric(findSbuEquil(pParam,rParam,pSbu)[1])
  while (round(target,4) != round(sEquil,4)){
    if (target < sEquil){
      pSbu <- pSbu * 1.5
      sEquil <- as.numeric(findSbuEquil(pParam,rParam,pSbu)[1])
    }
    else {
      pSbu <- pSbu/1.25
      sEquil <- as.numeric(findSbuEquil(pParam,rParam,pSbu)[1])
    }
  }
  sFound <- T
  return(c(sParam=pSbu,equil=round(sEquil,4)))
}

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
# cases$sbu <- sapply(1:3,function(i){sbu(cases$r[i],cases$p[i],cases$sol.exit0[i])})

# sex bias interpretation: Acceptance rate for men (with a different rate for women) generating the equivalent segregating effect over the same period of time for the 3 cases, and with the same overall acceptance rate as the case
# cases$accRateMen <- sapply(1:(dim(cases)[1]),function(i){
#   return((1-cases$r[i])*cases$sbu[i]/(1+((1-cases$r[i])*cases$sbu[i])-(1-cases$r[i])))
# })

# sex bias interpretation: Acceptance rate for women (with a different rate for men) generating the equivalent segregating effect over the same period of time for the 3 cases, and with the same overall acceptance rate as the case
# cases$accRateWom <- sapply(1:(dim(cases)[1]),function(i){
#   return((1-cases$r[i])/(cases$sbu[i]*(1+((1-cases$r[i])*cases$sbu[i])-(1-cases$r[i]))))
# })

# sex bias interpretation: Odds ratio of the gender-specific acceptance rates
# cases$sbuOdds <- sapply(1:(dim(cases)[1]),function(i){  # Odds Ratio of sex differences in acceptance rate
#   return(cases$accRateMen[i]*(1-cases$accRateWom[i])/(cases$accRateWom[i]*(1-cases$accRateMen[i])))
# })

getS <- function(pParam,rParam,target){
  sParam <- pParam*(1-target)/(target*(1-pParam))
  accM <- sParam*(1-rParam)/(1+(sParam*(1-rParam))-(1-rParam))
  accW <- (1-rParam)/(1+(sParam*(1-rParam))-(1-rParam))
  return(c(sParam=round(sParam,4),
           accM=round(accM,4),
           accW=round(accW,4),
           oddsR = round((accM/(1-accM))/(accW/(1-accW)),4)))
}

getSBE <- function(pParam,rParam,bParam,gParam){
  q <- equilPctF(pParam,rParam,bParam,gParam)
  r <- rParam
  p <- pParam
  return(c(oddsR = round(((r*q) - q + p)*(1-q) / (q *(q - (q*r) + r - p)),4),
           accM=round((r-1)*(q-1)/(1-p),4),
           accW=round((q-(q*r))/p,4)))
}

# Examine Paramater space
# # First pass
# # # Non-continuous dimensions:
# # # p - percent female among applicants: 10 30 50 70 90
# # # p - percent female among applicants: 10 20 30 40 50
# # # r - overall rejection rate: 10 30 50 70 90
# # # Axis dimenations:
# # # b - overall reapplication rate: 0..1
# # # g - gender difference in reapplication 1..2

# # # Axis dimenations:
# # # b - overall reapplication rate: 0..1
# # # g - gender difference in reapplication 1..2
pRange <- (10+10*(0:4))/100
rRange <- (10+20*(0:4))/100
gRange <- seq(1,2,0.01)
bRange <- seq(0.01,0.99,0.01)

z <- expand.grid(list(pRange,rRange))
z$vName <- unlist(sapply(1:(dim(z)[1]),function(r){paste0('z.p',z$Var1[r]*100,'r',z$Var2[r]*100)}))

# Surface / heatmap plot of segregating effects (as % of p, then OR)

data = data.frame(
  x = rep(gRange, each=length(bRange)),
  y = rep(bRange, length(gRange))
)

z.mat <- matrix(NA,nrow=dim(data)[1],ncol=dim(z)[1])
colnames(z.mat) <- z$vName
data <- cbind(data,z.mat)

for(c in 3:(dim(data)[2])){
  data[,c] <- unlist(sapply(1:(dim(data)[1]),function(r){
#    appEquil(digits=4,rParam=z$Var2[c-2],pParam=z$Var1[c-2],gParam=data$x[r],bParam=data$y[r])
    equilPctF(p=z$Var1[c-2],r=z$Var2[c-2],b=data$y[r],g=data$x[r])
  }))
}
save(data,file="dataFor3D.RData")
load("dataFor3D.RData")
wireframe(z.p30r90 ~ x * y, data=data)
par(mfrow=c(1,1))
plotContour <- function(pParam,rParam){
  x <- sort(unique(data$x)) 
  y <- sort(unique(data$y)) 
  z <- as.matrix(reshape(data[,c("x","y",paste0('z.p',pParam*100,'r',rParam*100))],
                        direction="wide",idvar="x",timevar="y"))[,-1]
#  fig <- plot_ly(x=x,y=y,z=z,type="contour")
#  fig <- filled.contour(x=x,y=y,z=z,zlim=c(0,1),nlevels=50)
  fig <- contour(x=x,y=y,z=z,zlim=c(0,1),nlevels=50)
  return(fig)
}
fig <- plotContour(pParam=0.3,rParam=0.9)
fig
contourplot(z.p30r90 ~ x * y, data=data)


getDataMat <- function(pParam,rParam){
  x <- sort(unique(data$x)) 
  y <- sort(unique(data$y)) 
  z <- as.matrix(reshape(data[,c("x","y",paste0('z.p',pParam*100,'r',rParam*100))],
                         direction="wide",idvar="x",timevar="y"))[,-1]
  row.names(z) <- as.character(x)
  colnames(z) <- as.character(y)
  return(z)
}
datMat.lst <- list(unlist(lapply(as.list(pRange),function(p){
  lapply(as.list(rRange),function(r){getDataMat(p,r)})
})))


fbyf.lst <- list(unlist(lapply(as.list(pRange),function(p){
  lapply(as.list(rRange),function(r){plotContour(p,r)})
})))
try5x5 <- subplot(fbyf.lst,nrows=5)

dev.off()

par(mfrow=c(5,5),mar=c(2,2,0,0))
for(p in pRange){
  for(r in rRange){
#    form <- as.formula(paste0('z.p',p*100,'r',r*100,' ~ x * y'))
#    print(formula)
    plotContour(p,r)
  }
}
par(mfrow=c(1,1))

png(filename=paste0("fiveXfive",gsub("-","",Sys.Date()),".PNG"),width=1200,height=1200)
par(mfrow=c(5,5),mar=c(2,2,0,0))
x <- sort(unique(data$x)) 
y <- sort(unique(data$y))
for(p in pRange){
  for(r in rRange){
    z <- as.matrix(reshape(data[,c("x","y",paste0('z.p',p*100,'r',r*100))],
                           direction="wide",idvar="x",timevar="y"))[,-1]
    contour(x=x,y=y,z=z,zlim=c(0,1),
           nlevels=50, axes=F,labcex=1.25)
    axis(1,at=seq(min(gRange),max(gRange),0.1),labels=(p==max(pRange)),line=0,cex=1.5,cex.axis=1.5)
    axis(2,at=seq(min(bRange),max(bRange),0.1),labels=(r==min(rRange)),line=0,cex=1.5,cex.axis=1.5)
    text(x=1.25,y=0.2,labels=paste0('Share female\n@start: ',p),cex=1.5)
    polygon(x=c(1,1,1.2,1.2),
            y=c(min(bRange),max(bRange),max(bRange),min(bRange)),
            col=rgb(0.5,0.5,0.5,0.2))
    polygon(x=c(min(gRange),min(gRange),max(gRange),max(gRange)),
            y=c(min(bRange),max(bRange),max(bRange),min(bRange)),
            col=rgb(0,0,0,0))
  }
}
par(mfrow=c(1,1))
dev.off()
for(p in 0.3){
  for(r in 0.7){
    z <- as.matrix(reshape(data[,c("x","y",paste0('z.p',p*100,'r',r*100))],
                           direction="wide",idvar="x",timevar="y"))[,-1]
    contour(x=x,y=y,z=z,zlim=c(0,1),
            nlevels=50, axes=F, labcex=1.25)
    axis(1,at=seq(min(gRange),max(gRange),0.1),labels=(p==0.3),line=0,cex=1.5,cex.axis=1.5)
    axis(2,at=seq(min(bRange),max(bRange),0.1),labels=(r==0.7),line=0,cex=1.5,cex.axis=1.5)
    text(x=1.25,y=0.2,labels=paste0('Share female\n@start: ',p),cex=1.5)
    polygon(x=c(1,1,1.2,1.2),
            y=c(min(bRange),max(bRange),max(bRange),min(bRange)),
            col=rgb(0.5,0.5,0.5,0.2))
    polygon(x=c(min(gRange),min(gRange),max(gRange),max(gRange)),
            y=c(min(bRange),max(bRange),max(bRange),min(bRange)),
            col=rgb(0,0,0,0))
  }
}
png(filename=paste0("cf_4dim_",gsub("-","",Sys.Date()),".PNG"),width=1200,height=1200)
par(mar=c(1,1,0,0))
cf_4dim(function(x){equilPctF(x[1],x[2],x[3],x[4])},
        over=c(2,1),nover1=5,nover2=5,
        low=c(0.1,0.1,0.01,1),
        high=c(0.9,0.9,0.99,2),
        nlevels=32,same_scale=T,
        var_names = c("p","r","b","g"),
#        color.palette = gray.colors,
        color.palette = rainbow,
        # same_scale=T,
        # bar_width=1, 
        with_lines=T,
        axes=T,
#        axes=list(axis(1,at=seq(1,2,0.2)),axis(2,at=seq(0,1,0.2))),
        n=5)
axis(1)
dev.off()

png(filename=paste0("cf_4dimOR_",gsub("-","",Sys.Date()),".PNG"),width=1200,height=1200)
par(mar=c(2,2,0,0))
cf_4dim(function(x){getS(x[1],x[2],target=equilPctF(x[1],x[2],x[3],x[4]))[4]},
        over=c(2,1),nover1=5,nover2=5,
        low=c(0.1,0.1,0.01,1),
        high=c(0.9,0.9,0.99,2),
        nlevels=32,
        var_names = c("p","r","b","g"),
        color.palette = rainbow,
        bar_width=1, with_lines=T,
        axes={ axis(1, seq(1, 2, by = 0.2))
          axis(2, seq(0, 1, by = 0.2)) },
        n=5)
dev.off()

for (x in seq(-0.5,2,0.5)) {
  for (y in seq(-0.5,2,0.5)){
    text(x=x,y=y,labels=paste0("x:",x,", y:",y))
  }
}

plots25.lst <- lapply(as.list(pRange),function(p){
  lapply(as.list(rRange),function(r){
    plotContour(p,r)
  })
})
contourplot(z.p10r10 ~ x * y, data=data)
filled.contour(z.p10r10 ~ x * y, data=data)

par(mfrow=c(1,1))

fbyf

par(mfrow=c(1,1))

data.or <- data
for(c in names(data.or)[-1*1:2]){
  pParam <- as.numeric(substr(c,4,5))/100
  rParam <- as.numeric(substr(c,7,8))/100
  data.or[,c] <- as.numeric(
    unlist(
      sapply(
        data.or[,c],function(t){
          getS(pParam,rParam,t)[4]
          })))
}

png(filename=paste0("fiveXfiveOR",gsub("-","",Sys.Date()),".PNG"),width=1200,height=1200)
par(mfrow=c(5,5),mar=c(2,2,0,0))
x <- sort(unique(data.or$x)) 
y <- sort(unique(data.or$y))
for(p in pRange){
  for(r in rRange){
    z <- as.matrix(reshape(data.or[,c("x","y",paste0('z.p',p*100,'r',r*100))],
                           direction="wide",idvar="x",timevar="y"))[,-1]
    contour(x=x,y=y,z=z,zlim=c(1,10),
            nlevels=50, axes=F,labcex=1.25)
    axis(1,at=seq(min(gRange),max(gRange),0.1),labels=(p==max(pRange)),line=0,cex=1.5,cex.axis=1.5)
    axis(2,at=seq(min(bRange),max(bRange),0.1),labels=(r==min(rRange)),line=0,cex=1.5,cex.axis=1.5)
    text(x=1.25,y=0.2,labels=paste0('Share female\n@start: ',p),cex=1.5)
    polygon(x=c(1,1,1.2,1.2),
            y=c(min(bRange),max(bRange),max(bRange),min(bRange)),
            col=rgb(0.5,0.5,0.5,0.2))
    polygon(x=c(min(gRange),min(gRange),max(gRange),max(gRange)),
            y=c(min(bRange),max(bRange),max(bRange),min(bRange)),
            col=rgb(0,0,0,0))
  }
}
par(mfrow=c(1,1))
dev.off()

paramSpace <- expand.grid(list(pRange,rRange,gRange,bRange))
cf_4dim(function(paramSpace){equilPctF(paramSpace[,1],paramSpace[,2],paramSpace[,3],paramSpace[,4])})

cf_4dim(function(x) {equilPctF(x[1],x[2],x[3],x[4])},nover=c(4,4),over=c(1,2),low=c(.1,.1,1,0),high=c(.9,.9,2,1))

#f4outer <- function(xv,yv){openModel(gParam=xv,bParam=yv,pParam=0.25,rParam=0.8)}
#zmat <- outer(X=gRange,Y=bRange,FUN=f4outer)

zmat50 <- reshape(data[,c("x","y","z50")],idvar="x",timevar="y",v.names="z50",direction="wide")
zmat25 <- reshape(data[,c("x","y","z25")],idvar="x",timevar="y",v.names="z25",direction="wide")
zmat75 <- reshape(data[,c("x","y","z75")],idvar="x",timevar="y",v.names="z75",direction="wide")

plot_ly(x=gRange,y=bRange,z=as.matrix(zmat25[,-1])) %>% add_surface()
plot_ly(x=gRange,y=bRange,z=as.matrix(zmat50[,-1])) %>% add_surface()
plot_ly(x=gRange,y=bRange,z=as.matrix(zmat75[,-1])) %>% add_surface()

x = rep(gRange, each=length(bRange))
y = rep(bRange, length(gRange))
z.df <- expand.grid(gRange,bRange)
z.df$z <- equilPctF(0.3,0.8,b=z.df$Var2,g=z.df$Var1)
data3d <- xyz.coords(z.df)
filled.contour2(x=gRange,
               y=bRange,
               z=matrix(data3d$z,nrow=length(gRange),ncol=length(bRange)),
               zlim=c(0,1),
               levels=seq(0,1,0.01),
               color.palette = rainbow
               )
contour(x=gRange,
               y=bRange,
               z=matrix(data3d$z,nrow=length(gRange),ncol=length(bRange)),
               zlim=c(0,1),
               levels=seq(0,1,0.01),
        add=T)

plot_ly(x=gRange,y=bRange,z=matrix(data3d$z,nrow=length(gRange),ncol=length(bRange)),
        type="contour",
        colorscale='rainbow',
        autocontour=F,
        contours=list(start=0.01,end=0.99,size=0.01,showlabels=T),
        line=list(smoothing=0)
        )

getPctFPanel <- function(p,r){
  fig <- plot_ly(
    x=gRange,
    y=bRange,
    z=matrix(equilPctF(p,
                       r,
                       b=data3d$y,
                       g=data3d$x),nrow=length(gRange),
             ncol=length(bRange)),
    type="contour",
    colorscale=cbind(seq(0, 1, by=0.01), rainbow(101)),
    autocontour=F,
    contours=list(start=0.01,end=0.99,size=0.01,showlabels=T),
    line=list(smoothing=0),
    showscale=F
  )
  fig <- fig %>% layout(
    xaxis=list(title="g"),
    yaxis=list(title="b")
  )
  return(fig)
}

getORPanel <- function(p,r){
  fig <- plot_ly(
    x=gRange,y=bRange,
    z=matrix(apply(as.data.frame(data3d[1:2]),
                   MAR=1,
                   FUN=function(paramRow){
                     getSBE(pParam=p,
                          rParam=r,
                          bParam=paramRow[2],
                          gParam=paramRow[1])[1]
      }),nrow=length(gRange),ncol=length(bRange)),
    type="contour",
    colorscale=cbind(seq(0, 1, by=0.01), rainbow(101)),
    autocontour=F,
    contours=list(start=1,end=10,size=0.1,showlabels=T),
    line=list(smoothing=0),
    showscale=F
  )
  fig <- fig %>% layout(
    xaxis=list(title="g"),
    yaxis=list(title="b")
  )
  return(fig)
}

getSubsPanel <- function(p,r,thresh=1.1,maxg=1.3){ # contour plot highlighting substantive effect threshold (1.1) and max value of g (1.3)
  xRange <- seq(1,maxg,0.01)
  z.df <- expand.grid(xRange,bRange)
  z.df$z <- equilPctF(p=p,r=r,b=z.df$Var2,g=z.df$Var1)
  coords3d <- xyz.coords(z.df)
  fig <- plot_ly(
    x=xRange,
    y=bRange,
    z=t(matrix(apply(as.data.frame(coords3d[1:2]),
                   MAR=1,
                   FUN=function(paramRow){
                     getSBE(pParam=p,
                            rParam=r,
                            bParam=paramRow[2],
                            gParam=paramRow[1])[1]
                   }),nrow=length(xRange),ncol=length(bRange))),
    type="contour",
#    colorscale=cbind(seq(1,thresh, by=(thresh-1)), gray.colors(2)),
    autocontour=F,
    contours=list(start=1,end=thresh,size=(thresh-1),showlabels=T),
    line=list(smoothing=0),
    showscale=F
  )
    fig <- fig %>% layout(
    xaxis=list(title="g",xlim=c(1,thresh)),
    yaxis=list(title="b")
  )
  return(fig)
}

paramsPR <- expand.grid(pRange,rRange)
panelsPctF <- lapply(
  lapply(1:(dim(paramsPR)[1]),
         function(ps){
           c(paramsPR[ps,1],paramsPR[ps,2])
           }),
  function(pr){
    getPctFPanel(pr[1],pr[2])
    })
panelsOR <- lapply(
  lapply(1:(dim(paramsPR)[1]),
         function(ps){
           c(paramsPR[ps,1],paramsPR[ps,2])
           }),
  function(pr){
    getORPanel(pr[1],pr[2])
    })
gridFigPctF <- subplot(panelsPctF[1:25],nrows=5)
gridFigPctF
gridFigOR <- subplot(panelsOR[1:25],nrows=5)
gridFigOR

panelsSubs <- lapply(
  lapply(1:(dim(paramsPR)[1]),
         function(ps){
           c(paramsPR[ps,1],paramsPR[ps,2])
         }),
  function(pr){
    getSubsPanel(pr[1],pr[2])
  })
gridFigSubs <- subplot(panelsSubs[1:25],nrows=5)
gridFigSubs

panPctF3.7 <- getPctFPanel(0.3,0.7)
panPctF3.7
panOR3.7 <- getORPanel(0.3,0.7)
panOR3.7
panSubs3.7 <- getSubsPanel(0.3,0.7)
panSubs3.7


pan37 <- getPanel(0.3,0.8)
pan1.1 <- panels[[1]]
pan3.7 <- panels[[17]]

gridFig <- subplot(lapply(panels,function(p){p[1]}),plot_ly(),nrows=5)
fig <- subplot(panels[1:25],nrows=5)
fig <- fig %>% layout(coloraxis=list(colorscale='Rainbow'))
fig

# # Second pass
# # # Non-continuous dimensions:
# # # p - percent female among applicants
# # # b - overall reapplication rate: 10 30 50 70 90
# # # Axis dimenations:
# # # b - overall reapplication rate: 0..1
# # # g - gender difference in reapplication 1..2

filled.contour3 <-
  function (x = seq(0, 1, length.out = nrow(z)),
            y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes,mar, ...) 
  {
    # modification by Ian Taylor of the filled.contour function
    # to remove the key and facilitate overplotting with contour()
    # further modified by Carey McGilliard and Bridget Ferris
    # to allow multiple plots on one page
    
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    # on.exit(par(par.orig))
    # w <- (3 + mar.orig[2]) * par("csi") * 2.54
    # par(las = las)
    # mar <- mar.orig
    plot.new()
    # par(mar=mar)
    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
    if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
      stop("no proper 'z' matrix specified")
    if (!is.double(z)) 
      storage.mode(z) <- "double"
    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
                            col = col))
    if (missing(plot.axes)) {
      if (axes) {
        title(main = "", xlab = "", ylab = "")
        Axis(x, side = 1)
        Axis(y, side = 2)
      }
    }
    else plot.axes
    if (frame.plot) 
      box()
    if (missing(plot.title)) 
      title(...)
    else plot.title
    invisible()
  }

filled.legend <-
  function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                         length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
            ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
            levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
            col = color.palette(length(levels) - 1), plot.title, plot.axes, 
            key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
            axes = TRUE, frame.plot = axes, ...) 
  {
    # modification of filled.contour by Carey McGilliard and Bridget Ferris
    # designed to just plot the legend
    if (missing(z)) {
      if (!missing(x)) {
        if (is.list(x)) {
          z <- x$z
          y <- x$y
          x <- x$x
        }
        else {
          z <- x
          x <- seq.int(0, 1, length.out = nrow(z))
        }
      }
      else stop("no 'z' matrix specified")
    }
    else if (is.list(x)) {
      y <- x$y
      x <- x$x
    }
    if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
      stop("increasing 'x' and 'y' values expected")
    #  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
    #  on.exit(par(par.orig))
    #  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
    #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
    #  par(las = las)
    #  mar <- mar.orig
    #  mar[4L] <- mar[2L]
    #  mar[2L] <- 1
    #  par(mar = mar)
    # plot.new()
    plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                yaxs = "i")
    rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
    if (missing(key.axes)) {
      if (axes) 
        axis(4)
    }
    else key.axes
    box()
  }
#
#    if (!missing(key.title)) 
#        key.title
#    mar <- mar.orig
#    mar[4L] <- 1
#    par(mar = mar)
#    plot.new()
#    plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
#    if (!is.matrix(z) || nrow(z) <= 1L || ncol(z) <= 1L) 
#        stop("no proper 'z' matrix specified")
#    if (!is.double(z)) 
#        storage.mode(z) <- "double"
#    .Internal(filledcontour(as.double(x), as.double(y), z, as.double(levels), 
#        col = col))
#    if (missing(plot.axes)) {
#        if (axes) {
#            title(main = "", xlab = "", ylab = "")
#            Axis(x, side = 1)
#            Axis(y, side = 2)
#        }
#    }
#    else plot.axes
#    if (frame.plot) 
#        box()
#    if (missing(plot.title)) 
#        title(...)
#    else plot.title
#    invisible()
#}

