# biased selection to equilibrium

bias <- function (t, state, parms) {
  with(as.list(c(state, parms)), {
    dMAcc <- (s*a/(1+(s*a)-a))*(1-pctF)
    dWAcc <- (a/((1+(s*a)-a)))*pctF
    list(c(dMAcc,dWAcc))
  })
}

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
  return(round(equil,4))
  
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

getS <- function(pParam,rParam,target){
  sParam <- pParam*(1-target)/(target*(1-pParam))
  accM <- sParam*(1-rParam)/(1+(sParam*(1-rParam))-(1-rParam))
  accW <- (1-rParam)/(1+(sParam*(1-rParam))-(1-rParam))
  return(c(sParam=round(sParam,4),
           accM=round(accM,4),
           accW=round(accW,4),
           oddsR = round((accM/(1-accM))/(accW/(1-accW)),4)))
}
