#don't remember which one i'm using so load them both
library(mnormt)
library(MESS)

# set up code. run as many times as you have trials/hazard rates/etc
A <- 2
t.length <- 2.33
nsim <- 50000
null.rate <-.693
hz  <- .75

control.rate <- ceiling(500/(1+A))
expt.rate <- floor(500*A/(1+A))
n.control <- ceiling(control.rate*t.length)
n.expt <- ceiling(expt.rate*t.length)
t.rate <- null.rate*hz

control.entry <- seq(0, t.length, length=n.control)
expt.entry <- seq(0, t.length, length=n.expt)
ctr.v <- c(rep(0, n.control), rep(1, n.expt))

#need to change the t to run different sims
pval.t1 <- rep(NA, nsim)

set.seed(3511)

for (i in 1:nsim) {
  ctrl.time <- rexp(n.control, rate=null.rate)
  tx.time <- rexp(n.expt, rate = t.rate)
  
  ctrl.list <- ctrl.time+control.entry
  tx.list <- tx.time +expt.entry

  censc <- ifelse(ctrl.list>t.length, 0, 1)
  censt <- ifelse(tx.list>t.length,0,1)
  ttec <- ifelse(censc, ctrl.time,t.length-control.entry)
  ttet <- ifelse(censt, tx.time, t.length-expt.entry)
  test <- summary(coxph(Surv(c(ttec, ttet), c(censc, censt))~ctr.v))
  #another place to update pval.t[x] with a different number if desired
  pval.t1[i] <- test$coef[5]
  if (i %% 5000 == 0) {
    cat("Simulation", i*100/nsim,"%", "donezo!", "\n")
  } 
};beep()

#two arms
corr<- .66

vcov <- matrix(data = c(1, corr, corr, 1), nrow = 2, ncol = 2)
mn <- rep(0, 2)

#can set length but this should cover it
sig <- 1-seq(.025, .035, by=.001)
ls <- length(sig)

props <- array(NA, dim = c(ls, ls))

for (j in 1:ls) {
    for (k in 1:ls) {
      sigf <- c(qnorm(sig[j]), qnorm(sig[k]))
      props[j,k] <- 1-pmnorm(sigf, mn, vcov)
  }
};beep()

#this is for correlation .66 specifically
#you have to change the loop depending what you're summing to 
#or do it manually as in the 3 and 4 case
pwr1 <- array(NA, dim = c(2, 9))

for (i in 1:9) {
  j <- .024+(i/1000)
  k <- .059-j
  pwr1[1,i] <- mean(pval.t1<=j)
  pwr1[2,i] <- mean(pval.t2<=k)
}

disjp <- rep(NA, 9)
conjp <- rep(NA, 9)

for (i in 1:9) {
  cval <- c(qnorm(pwr1[1,i]), qnorm(pwr1[2,i]))
  conjp[i] <- pmnorm(cval, mn, vcov)
  dval <- c(qnorm(1-pwr1[1,i]), qnorm(1-pwr1[2,i]))
  disjp[i] <- 1 - pmnorm(dval, mn, vcov)
}

#three arm case
corr1 <- .58 #arm 1/arm 2 corr 
corr2 <- .41 #arm 1/arm 3 corr
corr3 <- .49 #arm 2/arm 3 corr

vcov <- matrix(data = c(1, corr1, corr2, corr1, 1, corr3, 
                        corr2, corr3, 1), nrow = 3, ncol = 3)
mn <- rep(0, 3)

#rounded up to .17 on the bfrni correction
sig <- 1-seq(.017, .03, by=.001)
ls <- length(sig)

props <- array(NA, dim = c(ls, ls, ls))

for (i in 1:ls){
  for (j in 1:ls) {
    for (k in 1:ls) {
      sigf <- c(qnorm(sig[j]), qnorm(sig[k]), qnorm(sig[i]))
      props[i, j,k] <- 1-pmnorm(sigf, mn, vcov)
    }
  }
};beep()

conjp2 <- rep(NA, 20)
disjp2 <- rep(NA, 20)

#set these according to grid search
pwr2 <- c(
  mean(pval.t1<=.024),
  mean(pval.t2<=.017),
  mean(pval.t3<=.017)
)

#mess with the index i to store results in the vector for comparison
i <- 1
cval2 <- c(qnorm(pwr2[1]), qnorm(pwr2[2]), qnorm(pwr2[3]))
conjp2[i] <- pmnorm(cval2, mn, vcov)
dval2 <- c(qnorm(1-pwr2[1]), qnorm(1-pwr2[2]),qnorm(1-pwr2[3]))
disjp2[i] <- 1 - pmnorm(dval2, mn, vcov)


#four arm case
corr1 <- .58 #arm 1/arm 2 corr
corr2 <- .41 #arm 1/arm 3 corr
corr3 <- .49 #arm 2/arm 3 corr
corr4 <- .33 #arm 3/arm 4 corr
corr5 <- .14 #arm 2/arm 4 corr
corr6 <- .09 #arm 1/arm 4 corr

vcov <- matrix(data = c(1, corr1, corr2, corr6, corr1, 1, corr3, corr5, 
                        corr2, corr3, 1, corr4, corr6, corr5, corr4, 1), nrow = 4, ncol = 4)
mn <- rep(0, 4)

#rounded up again
sig <- 1-seq(.013, .025, by=.001)
ls <- length(sig)

props <- array(NA, dim = c(ls, ls, ls, ls))

for (i in 1:ls){
  for (j in 1:ls) {
    for (k in 1:ls) {
      for (l in 1:ls) {
        sigf <- c(qnorm(sig[i]), qnorm(sig[j]), qnorm(sig[k]), qnorm(sig[l]))
        props[i, j,k,l] <- 1-pmnorm(sigf, mn, vcov)
      }
    }
  }
  print(i)
};beep()
#there is no good way to look at these

conjp3 <- rep(NA, 20)
disjp3 <- rep(NA, 20)

pwr3 <- c(
  mean(pval.t1<=.0145),
  mean(pval.t2<=.0145),
  mean(pval.t3<=.0145),
  mean(pval.t4<=.0135)
)

i <- 1
cval3 <- c(qnorm(pwr3[1]), qnorm(pwr3[2]), qnorm(pwr3[3]), qnorm(pwr3[4]))
conjp3[i] <- pmnorm(cval3, mn, vcov)
dval3 <- c(qnorm(1-pwr3[1]), qnorm(1-pwr3[2]),qnorm(1-pwr3[3]), qnorm(1-pwr3[4]))
disjp3[i] <- 1 - pmnorm(dval3, mn, vcov)
