randSbe <- function(p=runif(1,min=0.01,max=0.49),
                    r=runif(1,min=0.01,max=0.99),
                    b=runif(1,min=0.01,max=0.99),
                    g=runif(1,min=1.01,max=1.5)){
  return(c(sub=as.numeric(sbe(p,r,b,g)[1]>=1.2),
           sbe=sbe(p,r,b,g)[1]))
}

paramEffectsP <- function(iters=1000){
  paramRange <- as.numeric(seq(0.01,0.99,0.02))
  pctSub <- sapply(paramRange,function(param){
    apply(replicate(iters,randSbe(p=param)),
          MAR=1,FUN=mean,na.rm=T)})
  return(data.frame(pRange=paramRange,
                    pctSub=pctSub[1,],
                    meanOR=pctSub[2,])) # return the percent the sampled point in parameter space
}
pEval <- paramEffectsP(10000)

paramEffectsR <- function(iters=1000){
  paramRange <- as.numeric(seq(0.01,0.99,0.02))
  pctSub <- sapply(paramRange,function(param){
    apply(replicate(iters,randSbe(r=param)),
          MAR=1,FUN=mean,na.rm=T)})
  return(data.frame(rRange=paramRange,
                    pctSub=pctSub[1,],
                    meanOR=pctSub[2,])) # return the percent the sampled point in parameter space
}
rEval <- paramEffectsR(10000)

paramEffectsB <- function(iters=1000){
  paramRange <- as.numeric(seq(0.01,0.99,0.02))
  pctSub <- sapply(paramRange,function(param){
    apply(replicate(iters,randSbe(b=param)),
          MAR=1,FUN=mean,na.rm=T)})
  return(data.frame(bRange=paramRange,
                    pctSub=pctSub[1,],
                    meanOR=pctSub[2,])) # return the percent the sampled point in parameter space
}
bEval <- paramEffectsB(10000)

paramEffectsG <- function(iters=1000){
  paramRange <- as.numeric(seq(1.01,1.49,0.02))
  pctSub <- sapply(paramRange,function(param){
    apply(replicate(iters,randSbe(g=param)),
          MAR=1,FUN=mean,na.rm=T)})
  return(data.frame(gRange=paramRange,
                    pctSub=pctSub[1,],
                    meanOR=pctSub[2,])) # return the percent the sampled point in parameter space
}
gEval <- paramEffectsG(10000)

png(filename=paste0("paramInterp_",gsub("-","",Sys.Date()),".PNG"),
    width=6.5,height=8,units="in",res=300)

par(mar=c(4.1,5.1,4.1,0.5),mfrow=c(4,2))


plot(x=pEval$pRange,
     y=pEval$meanOR,
     ylim=c(1,1.25),type="l",
     main="Effects of p parameter:\nWomen's share of first-time applicants",
     xlab="p parameter value",
     ylab="Mean ORm - sex bias equivalent\nodds ratio favoring men in selection")
abline(v=cases$p,lt=3)
text(x=cases$p-0.02,y=1.15,labels=cases$caseID,srt=90)
plot(x=pEval$pRange,
     y=pEval$pctSub,
     ylim=c(0,0.3), type="l",
     main="Effects of p parameter:\nWomen's share of first-time applicants",
     xlab="p parameter value",
     ylab="Portion of sampled parameter space\nwith substantive segregating effects")
abline(v=cases$p,lt=3)
text(x=cases$p-0.025,y=0.18,labels=cases$caseID,srt=90)


plot(x=rEval$rRange,
     y=rEval$meanOR,
     ylim=c(1,1.25),type="l",
     main="Effects of r parameter: Rejection rate",
     xlab="r parameter value",
     ylab="Mean ORm - sex bias equivalent\nodds ratio favoring men in selection")
abline(v=cases$r,lt=3)
text(x=cases$r-0.02,y=1.15,labels=cases$caseID,srt=90)
plot(x=rEval$rRange,
     y=rEval$pctSub,
     ylim=c(0,0.3),type="l",
     main="Effects of r parameter: Rejection rate",
     xlab="r parameter value",
     ylab="Portion of sampled parameter space\nwith substantive segregating effects")
abline(v=cases$r,lt=3)
text(x=cases$r-0.025,y=0.18,labels=cases$caseID,srt=90)

plot(x=bEval$bRange,
     y=bEval$meanOR,
     ylim=c(1,1.25),type="l",
     main="Effects of b parameter:\nWomen's repplication rate",
     xlab="b parameter value",
     ylab="Mean ORm - sex bias equivalent\nodds ratio favoring men in selection")
abline(v=cases$b,lt=3)
text(x=cases$b-0.02,y=1.15,labels=cases$caseID,srt=90)
plot(x=bEval$bRange,
     y=bEval$pctSub,
     ylim=c(0,0.3),type="l",
     main="Effects of b parameter:\nWomen's repplication rate",
     xlab="b parameter value",
     ylab="Portion of sampled parameter space\nwith substantive segregating effects")
abline(v=cases$b,lt=3)
text(x=cases$b-0.025,y=0.18,labels=cases$caseID,srt=90)

plot(x=gEval$gRange,
     y=gEval$meanOR,
     ylim=c(1,1.25),type="l",
     main="Effects of g parameter: Odds ratio of\nmen's greater reapplication rate",
     xlab="g parameter value",
     ylab="Mean ORm - sex bias equivalent\nodds ratio favoring men in selection")
abline(v=cases$g,lt=3)
text(x=cases$g-0.01,y=1.15,labels=cases$caseID,srt=90)
plot(x=gEval$gRange,
     y=gEval$pctSub,
     ylim=c(0,0.3),type="l",
     main="Effects of g parameter: Odds ratio of\nmen's greater reapplication rate",
     xlab="g parameter value",
     ylab="Portion of sampled parameter space\nwith substantive segregating effects")
abline(v=cases$g,lt=3)
text(x=cases$g-0.01,y=0.18,labels=cases$caseID,srt=90)

dev.off()

prop1Sbe <- function(p=runif(1,min=0.01,max=0.49),
                    r=runif(1,min=0.8,max=0.99),
                    b=runif(1,min=0.55,max=0.85),
                    g=runif(1,min=1.3,max=1.5)){
  return(as.numeric(sbe(p,r,b,g)[1]>=1.2))
}
mean(replicate(100000,prop1Sbe()))

