randSbe <- function(p=runif(1,min=0.01,max=0.49),
                    r=runif(1,min=0.01,max=0.99),
                    b=runif(1,min=0.01,max=0.99),
                    g=runif(1,min=1.01,max=1.5)){
  sbe1 <- sbe(p,r,b,g)[1]
  return(c(sub=as.numeric(sbe1>=1.2),
           sbe=sbe1))
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
                    b=runif(1,min=0.40,max=0.99),
                    g=runif(1,min=1.3,max=1.5)){
  sbe1 <- sbe(p,r,b,g)[1]
  return(c(sub=as.numeric(sbe1>=1.2),
           sbe=sbe1))
}

rowMeans(replicate(100000,prop1Sbe()))

prop1EffectsB <- function(iters=1000){
  paramRange <- as.numeric(seq(0.01,0.99,0.02))
  pctSub <- sapply(paramRange,function(param){
    rowMeans(replicate(iters,prop1Sbe(b=param)),na.rm=T)})
  return(data.frame(bRange=paramRange,
                    pctSub=pctSub[1,],
                    meanOR=pctSub[2,])) # return the percent the sampled point in parameter space
}

bEvalP1 <- prop1EffectsB(10000)

plot(x=bEvalP1$bRange,
     y=bEvalP1$pctSub,
     type="l",
     main="Effects of b parameter:\nWomen's repplication rate",
     xlab="b parameter value",
     ylab="Portion of sampled parameter space\nwith substantive segregating effects")
abline(h=0.34,lt=2)
abline(v=c(0.54,0.88),lt=2)
# abline(v=(5:9)/10,lt=2)
abline(v=cases$b,lt=3)
text(x=cases$b-0.025,y=0.18,labels=cases$caseID,srt=90)



prop1EffectsP <- function(iters=1000){
  paramRange <- as.numeric(seq(0.01,0.49,0.01))
  pctSub <- sapply(paramRange,function(param){
    rowMeans(replicate(iters,prop1Sbe(p=param)),na.rm=T)})
  return(data.frame(pRange=paramRange,
                    pctSub=pctSub[1,],
                    meanOR=pctSub[2,])) # return the percent the sampled point in parameter space
}
pEvalP1 <- prop1EffectsP(10000)

prop1EffectsR <- function(iters=1000){
  paramRange <- as.numeric(seq(0.01,0.99,0.02))
  pctSub <- sapply(paramRange,function(param){
    rowMeans(replicate(iters,prop1Sbe(r=param)),na.rm=T)})
  return(data.frame(rRange=paramRange,
                    pctSub=pctSub[1,],
                    meanOR=pctSub[2,])) # return the percent the sampled point in parameter space
}
rEvalP1 <- prop1EffectsR(10000)

prop1EffectsB <- function(iters=1000){
  paramRange <- as.numeric(seq(0.01,0.99,0.02))
  pctSub <- sapply(paramRange,function(param){
    rowMeans(replicate(iters,prop1Sbe(b=param)),na.rm=T)})
  return(data.frame(bRange=paramRange,
                    pctSub=pctSub[1,],
                    meanOR=pctSub[2,])) # return the percent the sampled point in parameter space
}
bEvalP1 <- prop1EffectsB(10000)

prop1EffectsG <- function(iters=1000){
  paramRange <- as.numeric(seq(1.01,1.49,0.01))
  pctSub <- sapply(paramRange,function(param){
    rowMeans(replicate(iters,prop1Sbe(g=param)),na.rm=T)})
  return(data.frame(gRange=paramRange,
                    pctSub=pctSub[1,],
                    meanOR=pctSub[2,])) # return the percent the sampled point in parameter space
}
gEvalP1 <- prop1EffectsG(10000)


png(filename=paste0("paramProp1_",gsub("-","",Sys.Date()),".PNG"),
    width=6.5,height=8,units="in",res=300)

par(mar=c(4.1,5.1,4.1,0.5),mfrow=c(4,2))

plot(x=pEvalP1$pRange,
     y=pEvalP1$meanOR,
     ylim=c(1,1.4),type="l",
     main="Effects of p parameter\nWhen r>0.8, 0.40<b<1, and g>1.3",
     xlab="p parameter value",
     ylab="Mean ORm - sex bias equivalent\nodds ratio favoring men in selection")
abline(v=cases$p,lt=3)
text(x=cases$p-0.02,y=1.2,labels=cases$caseID,srt=90)
plot(x=pEvalP1$pRange,
     y=pEvalP1$pctSub,
     ylim=c(0,1), type="l",
     main="Effects of p parameter\nWhen r>0.8, 0.40<b<1, and g>1.3",
     xlab="p parameter value",
     ylab="Portion of sampled parameter space\nwith substantive segregating effects")
abline(v=cases$p,lt=3)
text(x=cases$p-0.025,y=0.5,labels=cases$caseID,srt=90)

plot(x=rEvalP1$rRange,
     y=rEvalP1$meanOR,
     ylim=c(1,1.4),type="l",
     main="Effects of r parameter\nWhen 0.40<b<1 and g>1.3",
     xlab="r parameter value",
     ylab="Mean ORm - sex bias equivalent\nodds ratio favoring men in selection")
abline(v=cases$r,lt=3)
text(x=cases$r-0.02,y=1.2,labels=cases$caseID,srt=90)
plot(x=rEvalP1$rRange,
     y=rEvalP1$pctSub,
     ylim=c(0,1),type="l",
     main="Effects of r parameter\nWhen 0.40<b<1 and g>1.3",
     xlab="r parameter value",
     ylab="Portion of sampled parameter space\nwith substantive segregating effects")
abline(v=cases$r,lt=3)
text(x=cases$r-0.025,y=0.5,labels=cases$caseID,srt=90)

plot(x=bEvalP1$bRange,
     y=bEvalP1$meanOR,
     ylim=c(1,1.4),type="l",
     main="Effects of b parameter\nWhen r>0.8 and g>1.3",
     xlab="b parameter value",
     ylab="Mean ORm - sex bias equivalent\nodds ratio favoring men in selection")
abline(v=cases$b,lt=3)
text(x=cases$b-0.02,y=1.2,labels=cases$caseID,srt=90)
plot(x=bEvalP1$bRange,
     y=bEvalP1$pctSub,
     ylim=c(0,1),type="l",
     main="Effects of b parameter\nWhen r>0.8 and g>1.3",
     xlab="b parameter value",
     ylab="Portion of sampled parameter space\nwith substantive segregating effects")
abline(v=cases$b,lt=3)
text(x=cases$b-0.025,y=0.5,labels=cases$caseID,srt=90)
abline(h=0.34,lt=2,lwd=0.5,col="gray")
abline(v=c(0.54,0.88),lt=2,lwd=0.5,col="gray")

plot(x=gEvalP1$gRange,
     y=gEvalP1$meanOR,
     ylim=c(1,1.4),type="l",
     main="Effects of g parameter\nWhen r>0.8 and 0.40<b<1",
     xlab="g parameter value",
     ylab="Mean ORm - sex bias equivalent\nodds ratio favoring men in selection")
abline(v=cases$g,lt=3)
text(x=cases$g-0.01,y=1.2,labels=cases$caseID,srt=90)
plot(x=gEvalP1$gRange,
     y=gEvalP1$pctSub,
     ylim=c(0,1),type="l",
     main="Effects of g parameter\nWhen r>0.8 and 0.40<b<1",
     xlab="g parameter value",
     ylab="Portion of sampled parameter space\nwith substantive segregating effects")
abline(v=cases$g,lt=3)
text(x=cases$g-0.01,y=0.5,labels=cases$caseID,srt=90)

dev.off()
