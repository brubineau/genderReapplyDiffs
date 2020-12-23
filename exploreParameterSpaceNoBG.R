########################
# 1. Set up environment
# ---------------------
#    1.1. Clear work environment
rm(list=ls())
gc()

#    1.2. Load needed libraries
#         Note: if the "deSolve" package is not already installed, de-comment and run (once only) the following line:
#         install.packages("deSolve")
library(plotly)
library(ContourFunctions)

#    1.3. Change to working directory
# setwd("~/") # UN-COMMENT & ADAPT FOR YOUR OWN COMPUTING CONTEXT

############################
# 3. MAIN MODEL AND VARIANT
# -------------------------
# 3.2 Continuous differential equation model - applications only
#     Maximum Turnover - Instantaneous Exit

# need to re-work
# w is women's reapplication rate = b/(1+bg-b)
# m is men's reapplication rate = bg/(1+bg-b) = w*g
# OR is odds ratio:  (m*(1-w))/(w*(1-m))
# m = OR*w/(w-OR*w-1)

# analytic solution?
# share of return applicants
fShareWO <- function(t,p,r,w,o){
  g <- o/(1-w+(o*w))
  b <- w/(1+w-(g*w))
  k <- r*b/(1+(b*g)-b)
  return((exp((t+log(k*p))*(k-1))-p)/(k-1))
}
mShareWO <- function(t,p,r,w,o){
  g <- o/(1-w+(o*w))
  b <- w/(1+w-(g*w))
  k <- r*b/(1+(b*g)-b)
  return((exp((t+(log(g*k-g*p*k)/(g*k-1)))*(g*k-1))-(1-p))/(g*k-1))
}

equilPctFWO <- function(p,r,w,o){
  g <- o/(1-w+(o*w))
  b <- w/(1+w-(g*w))
#  fShare <- (exp(100*((r*b/(1+(b*g)-b))-1))-p)/((r*b/(1+(b*g)-b))-1)
#  mShare <- (exp(100*((r*b*g/(1+(b*g)-b))-1))-(1-p))/((r*b*g/(1+(b*g)-b))-1)
  return(fShareWO(1000,p,r,w,o)/(fShareWO(1000,p,r,w,o)+mShareWO(1000,p,r,w,o)))
}

fShareMW <- function(t,p,r,m,w){
  o <- (m*(1-w))/(w*(1-m))
  g <- o/(1-w+(o*w))
  b <- w/(1+w-(g*w))
  k <- r*b/(1+(b*g)-b)
  return((exp((t+log(k*p))*(k-1))-p)/(k-1))
}
mShareMW <- function(t,p,r,m,w){
  o <- (m*(1-w))/(w*(1-m))
  g <- o/(1-w+(o*w))
  b <- w/(1+w-(g*w))
  k <- r*b/(1+(b*g)-b)
  return((exp((t+(log(g*k-g*p*k)/(g*k-1)))*(g*k-1))-(1-p))/(g*k-1))
}

equilPctFMW <- function(p,r,m,w){
  o <- (m*(1-w))/(w*(1-m))
  g <- o/(1-w+(o*w))
  b <- w/(1+w-(g*w))
  return(fShareMW(1000,p,r,m,w)/(fShareMW(1000,p,r,m,w)+mShareMW(1000,p,r,m,w)))
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

getSBEMW <- function(pParam,rParam,mParam,wParam){
  oParam <- (mParam*(1-wParam))/(wParam*(1-mParam))
  gParam <- oParam/(1-wParam+(oParam*wParam))
  bParam <- wParam/(1+wParam-(gParam*wParam))
  q <- equilPctFMW(pParam,rParam,mParam,wParam)
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
wRange <- seq(0.01,0.99,0.01)
mRange <- seq(0.01,0.99,0.01)
oRange <- seq(1,2,0.1)

# x = rep(mRange, each=length(mRange))
# y = rep(wRange, length(wRange))
zMW.df <- expand.grid(mRange,wRange)
zMW.df$z <- equilPctFMW(0.3,0.8,w=zMW.df$Var2,m=zMW.df$Var1)
data3dMW <- xyz.coords(zMW.df)

# x = rep(oRange, each=length(oRange))
# y = rep(wRange, length(wRange))
zWO.df <- expand.grid(oRange,wRange)
zWO.df$z <- equilPctFWO(0.3,0.8,w=zWO.df$Var2,o=zWO.df$Var1)
data3dWO <- xyz.coords(zWO.df)

getPctFPanelMW <- function(p,r){
  zmat <- t(matrix(equilPctFMW(p,
                             r,
                             w=data3dMW$y,
                             m=data3dMW$x),
                   nrow=length(mRange),
                   ncol=length(wRange)))
  zmat_nas <- t(matrix(as.numeric(data3dMW$y>data3dMW$x),
                       nrow=length(mRange),ncol=length(wRange)))
  zmat[which(zmat_nas==1)]<- NA
  fig <- plot_ly(
    x=mRange,
    y=wRange,
    z=zmat,
    type="contour",
    colorscale=cbind(seq(0, 1, by=0.01), rainbow(101)),
    autocontour=F,
    contours=list(start=0.01,end=0.99,size=0.01,showlabels=T),
    line=list(smoothing=0),
    showscale=F
  )
  fig <- fig %>% layout(
    xaxis=list(title="Men's reapplication rate"),
    yaxis=list(title="Women's reapplication rate")
  )
  return(fig)
}

getORPanelMW <- function(p,r){
  zmat <- t(matrix(apply(as.data.frame(data3dMW[1:2]),
                         MAR=1,
                         FUN=function(paramRow){
                           getSBEMW(pParam=p,
                                  rParam=r,
                                  wParam=paramRow[2],
                                  mParam=paramRow[1])[1]
                         }),nrow=length(mRange),ncol=length(wRange)))
  zmat_nas <- t(matrix(as.numeric(data3dMW$y>data3dMW$x),
                       nrow=length(mRange),ncol=length(wRange)))
  zmat[which(zmat_nas==1)]<- NA
  fig <- plot_ly(
    x=mRange,y=wRange,
    z=zmat,
    type="contour",
    colorscale=cbind(seq(0, 1, by=0.01), rainbow(101)),
    autocontour=F,
    contours=list(start=1,end=10,size=0.1,showlabels=T),
    line=list(smoothing=0),
    showscale=F
  )
  fig <- fig %>% layout(
    xaxis=list(title="Men's reapplication rate"),
    yaxis=list(title="Women's reappication rate")
  )
  return(fig)
}

getPctFPanelWO <- function(p,r){
  zmat <- t(matrix(equilPctFWO(p,
                             r,
                             w=data3dWO$y,
                             o=data3dWO$x),
                   nrow=length(oRange),
                   ncol=length(wRange)))
  zmat_nas <- t(matrix(as.numeric(data3dWO$x<1),
                       nrow=length(mRange),ncol=length(wRange)))
  zmat[which(zmat_nas==1)]<- NA
  fig <- plot_ly(
    x=oRange,
    y=wRange,
    z=zmat,
    type="contour",
    colorscale=cbind(seq(0, 1, by=0.01), rainbow(101)),
    autocontour=F,
    contours=list(start=0.01,end=0.99,size=0.01,showlabels=T),
    line=list(smoothing=0),
    showscale=F
  )
  fig <- fig %>% layout(
    xaxis=list(title="Odds Ratio of men's higher reapplication rate"),
    yaxis=list(title="Women's reapplication rate")
  )
  return(fig)
}

getORPanelWO <- function(p,r){
  zmat <- t(matrix(apply(as.data.frame(data3dWO[1:2]),
                         MAR=1,
                         FUN=function(paramRow){
                           getSBEWO(pParam=p,
                                  rParam=r,
                                  wParam=paramRow[2],
                                  oParam=paramRow[1])[1]
                         }),nrow=length(oRange),ncol=length(wRange)))
  zmat_nas <- t(matrix(as.numeric(data3dWO$x<1),
                       nrow=length(oRange),ncol=length(wRange)))
  zmat[which(zmat_nas==1)]<- NA
  fig <- plot_ly(
    x=oRange,y=wRange,
    z=zmat,
    type="contour",
    colorscale=cbind(seq(0, 1, by=0.01), rainbow(101)),
    autocontour=F,
    contours=list(start=1,end=10,size=0.1,showlabels=T),
    line=list(smoothing=0),
    showscale=F
  )
  fig <- fig %>% layout(
    xaxis=list(title="Odds Ratio of men's higher reapplication rate"),
    yaxis=list(title="Women's reappication rate")
  )
  return(fig)
}

paramsPR <- expand.grid(pRange,rRange)

panelsPctFMW <- lapply(
  lapply(1:(dim(paramsPR)[1]),
         function(ps){
           c(paramsPR[ps,1],paramsPR[ps,2])
           }),
  function(pr){
    getPctFPanelMW(pr[1],pr[2])
    })

panelsORMW <- lapply(
  lapply(1:(dim(paramsPR)[1]),
         function(ps){
           c(paramsPR[ps,1],paramsPR[ps,2])
           }),
  function(pr){
    getORPanelMW(pr[1],pr[2])
    })

panelsPctFWO <- lapply(
  lapply(1:(dim(paramsPR)[1]),
         function(ps){
           c(paramsPR[ps,1],paramsPR[ps,2])
         }),
  function(pr){
    getPctFPanelWO(pr[1],pr[2])
  })

panelsORWO <- lapply(
  lapply(1:(dim(paramsPR)[1]),
         function(ps){
           c(paramsPR[ps,1],paramsPR[ps,2])
         }),
  function(pr){
    getORPanelWO(pr[1],pr[2])
  })

gridFigPctFMW <- subplot(panelsPctFMW[1:25],nrows=5)
gridFigPctFMW

gridFigORMW <- subplot(panelsORMW[1:25],nrows=5)
gridFigORMW

gridFigPctFWO <- subplot(panelsPctFWO[1:25],nrows=5)
gridFigPctFWO

gridFigORWO <- subplot(panelsORWO[1:25],nrows=5)
gridFigORWO

panPctMWF3.7 <- getPctFPanelMW(0.3,0.7)
panPctMWF3.7
panORMW3.7 <- getORPanelMW(0.3,0.7)
panORMW3.7
panPctFWO3.7 <- getPctFPanelWO(0.3,0.7)
panPctFWO3.7
panORWO3.7 <- getORPanelWO(0.3,0.7)
panORWO3.7

pan37 <- getPanel(0.3,0.8)
pan1.1 <- panels[[1]]
pan3.7 <- panels[[17]]

gridFig <- subplot(lapply(panels,function(p){p[1]}),plot_ly(),nrows=5)
fig <- subplot(panels[1:25],nrows=5)
fig <- fig %>% layout(coloraxis=list(colorscale='Rainbow'))
fig


getORcalcMW <- function(){
  zmat <- t(matrix(data3dMW$x*(1-data3dMW$y)/(data3dMW$y*(1-data3dMW$x)),
                   nrow=length(mRange),
                   ncol=length(wRange)))
  zmat_nas <- t(matrix(as.numeric(data3dMW$y>data3dMW$x),
                       nrow=length(mRange),ncol=length(wRange)))
  zmat[which(zmat_nas==1)]<- NA
  fig <- plot_ly(
    x=mRange,
    y=wRange,
    z=zmat,
    type="contour",
    colorscale=cbind(seq(0, 1, by=0.01), rainbow(101)),
    autocontour=F,
    contours=list(start=1,end=100,size=0.5,showlabels=T),
    line=list(smoothing=0),
    showscale=F
  )
  fig <- fig %>% layout(
    xaxis=list(title="Men's reapplication rate"),
    yaxis=list(title="Women's reapplication rate")
  )
  return(fig)
}
getORcalcMW()
