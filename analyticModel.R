# Empirical Grounding Table

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

# analytic solution?
# share of return applicants
fShare <- function(t,p,r,b,g){
  k <- r*b/(1+(b*g)-b)
  return((exp((t+(log(k*p)/(k-1)))*(k-1))-p)/(k-1))
}
mShare <- function(t,p,r,b,g){
  k <- r*b/(1+(b*g)-b)
  return((exp((t+(log((g*k)-(g*p*k))/((g*k)-1)))*((g*k)-1))-(1-p))/((g*k)-1))
}

equilPctF <- function(p,r,b,g){
  #  fShare <- (exp(100*((r*b/(1+(b*g)-b))-1))-p)/((r*b/(1+(b*g)-b))-1)
  #  mShare <- (exp(100*((r*b*g/(1+(b*g)-b))-1))-(1-p))/((r*b*g/(1+(b*g)-b))-1)
  return(fShare(1000,p,r,b,g)/(fShare(1000,p,r,b,g)+mShare(1000,p,r,b,g)))
}

cases$sol.exit0 <- unlist(sapply(1:3,function(c){
  equilPctF(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c])}))

getSBE <- function(pParam,rParam,bParam,gParam){
  q <- equilPctF(pParam,rParam,bParam,gParam)
  r <- rParam
  p <- pParam
  return(c(oddsR = round(((r*q) - q + p)*(1-q) / (q *(q - (q*r) + r - p)),4),
           accM=round((r-1)*(q-1)/(1-p),4),
           accW=round((q-(q*r))/p,4)))
}

sbe <- sapply(1:3,function(c){
  getSBE(pParam=cases$p[c],rParam=cases$r[c],bParam=cases$b[c],gParam=cases$g[c])})

cases$accRateMen <- sbe[2,]
cases$accRateWom <- sbe[3,]
cases$sbuOddsRatio <- sbe[1,]

##########
# 8. Table
# --------

tab2 <- round(rbind(t(cases[,c("p","r","b","g")]),
                    g_ci = ((cases$gub-cases$g)+(cases$g-cases$glb))/2,
                    after20 = unlist(sapply(1:3,function(c){
                      f <- fShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=20)
                      m <- mShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=20)
                      return(f/(f+m))})),
                    after20ci = unlist(sapply(1:3,function(c){
                      fub <- fShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$gub[c],t=20)
                      mub <- mShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$gub[c],t=20)
                      flb <- fShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$glb[c],t=20)
                      mlb <- mShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$glb[c],t=20)
                      return(((flb/(flb+mlb))-(fub/(fub+mub)))/2)})),
                    atEquil = cases$sol.exit0,
                    atEquilci = unlist(sapply(1:3,function(c){
                      eub <- equilPctF(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$gub[c])
                      elb <- equilPctF(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$glb[c])
                      return((elb-eub)/2)})),
                    delta_abs = cases$sol.exit0-cases$p,
                    delta_rel = (cases$sol.exit0-cases$p)/cases$p,
                    t(cases[,c("accRateMen","accRateWom","sbuOddsRatio")]))[,c(1,3,2)],4)
colnames(tab2) <- cases$caseID[c(1,3,2)]
tab2

unlist(sapply(1:3,function(c){equilPctF(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$gub[c])}))
unlist(sapply(1:3,function(c){equilPctF(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$glb[c])}))
unlist(sapply(1:3,function(c){
  f <- fShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=100)
  m <- mShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=100)
  return(f/(f+m))}))

plot(1:100, 
     unlist(
       sapply(
         1:100,
         function(tp){
           c <- 1
           f <- fShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=tp)
           m <- mShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=tp)
           return(f/(f+m))
         })),type="l", ylim=c(0,0.5))

lines(1:100, 
     unlist(
       sapply(
         1:100,
         function(tp){
           c <- 2
           f <- fShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=tp)
           m <- mShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=tp)
           return(f/(f+m))
         })),type="l")

lines(1:100, 
      unlist(
        sapply(
          1:100,
          function(tp){
            c <- 3
            f <- fShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=tp)
            m <- mShare(p=cases$p[c],r=cases$r[c],b=cases$b[c],g=cases$g[c],t=tp)
            return(f/(f+m))
          })),type="l",col="green")
