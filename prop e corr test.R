library(survival)
library(survminer)
library(janitor)
library(beepr)

nsim <- 5000
null.hz <- .3
tx.hz <- .75*null.hz
nevents <- 20
start2 <- ceiling(nevents*.2)
ntot <- c(ceiling(nevents*.8), ceiling(nevents*.8))
arm_no <- c(2, 6)
n.arms <- c(5,6,2)
gshape <- 1
grate <- .693
grate_tx <- grate*1
rate <- 20

sh.ctrl <- rep(NA, nsim)
ctrl.evs <- rep(NA, nsim)

ts1 <- rep(NA, nsim)
ts5 <- rep(NA, nsim)

for (i in 1:nsim) {
  entry <- seq(1/rate,100, by=1/rate)
  ctrl1entry.1 <- entry[seq(1,length(entry), by=n.arms[1])]
  tx1entry.1 <- entry[seq(arm_no[1],length(entry),by=n.arms[1])]
  ctrl1evs.1 <- rgamma(length(ctrl1entry.1), gshape, rate=grate)#rexp(length(ctrl1entry.1), null.hz)
  tx1evs.1 <- rgamma(length(tx1entry.1), gshape, rate=grate_tx)#rexp(length(tx1entry.1), tx.hz)
  
  ctrl1time.1 <- ctrl1entry.1+ctrl1evs.1
  tx1time.1 <- tx1entry.1 + tx1evs.1
  
  comp1 <- c(ctrl1time.1, tx1time.1)
  censTime1.1 <- comp1[order(comp1)][start2]
  
  inds1 <- which(ctrl1entry.1<=censTime1.1)
  inds2 <- which(tx1entry.1<=censTime1.1)
  
  ctrl1evs.1 <- ctrl1evs.1[inds1]
  tx1evs.1 <- tx1evs.1[inds2]
  ctrl1time.1 <- ctrl1time.1[inds1]
  tx1time.1 <- tx1time.1[inds2]
  ctrl1entry.1 <- ctrl1entry.1[inds1]
  tx1entry.1 <- tx1entry.1[inds2]
  
  nu_st <- round_to_fraction(censTime1.1, rate, 3)
  
  entry2 <- seq(nu_st+1/rate,10+nu_st, by=1/rate)
  ctrl1entry.2 <- entry2[seq(1,length(entry2), by=n.arms[2])]
  tx1entry.2 <- entry2[seq(arm_no[1],length(entry2),by=n.arms[2])]
  ctrl5entry.1 <- entry[seq(1,length(entry), by=n.arms[2])]
  tx5entry.1 <- entry[seq(arm_no[2], length(entry), by=n.arms[2])]
  
  ctrl.sh.evs <- rgamma(length(ctrl1entry.2), gshape, rate=grate)#rexp(length(ctrl1entry.2), null.hz)
  tx1evs.2 <- rgamma(length(tx1entry.2), gshape, rate=grate_tx)#rexp(length(tx1entry.2), tx.hz)
  tx5evs.1 <- rgamma(length(tx5entry.1), gshape, rate=grate_tx)#rexp(length(tx5entry.1), tx.hz)
  
  ctrl1time.2 <- ctrl1entry.2+ctrl.sh.evs
  tx1time.2 <- tx1entry.2 + tx1evs.2
  ctrl5time.1 <- ctrl5entry.1 + ctrl.sh.evs
  tx5time.1 <- tx5entry.1 + tx5evs.1
  
  comp1 <- c(ctrl1time.1, tx1time.1, ctrl1time.2, tx1time.2)
  censTime1.2 <- comp1[order(comp1)][ntot[1]]
  
  inds1 <- which(ctrl1entry.2<censTime1.2)
  inds2 <- which(tx1entry.2<censTime1.2)
  inds3 <- which(tx5entry.1<(censTime1.2-nu_st))
  
  ctrl.sh.evs <- ctrl.sh.evs[inds1]
  tx1evs.2 <- tx1evs.2[inds2]
  tx5evs.1 <- tx5evs.1[inds3]
  ctrl1time.2 <- ctrl1time.2[inds1]
  ctrl5time.1 <- ctrl5time.1[inds1]
  tx1time.2 <- tx1time.2[inds2]
  tx5time.1 <- tx5time.1[inds3]
  ctrl1entry.2 <- ctrl1entry.2[inds1]
  ctrl5entry.1 <- ctrl5entry.1[inds1]
  tx1entry.2 <- tx1entry.2[inds2]
  tx5entry.1 <- tx5entry.1[inds3]
  
  ctrl1.evs <- c(ctrl1evs.1, ctrl.sh.evs)
  tx1.evs <- c(tx1evs.1, tx1evs.2)
  ctrl1.time <- c(ctrl1time.1, ctrl1time.2)
  tx1.time <- c(tx1time.1, tx1time.2)
  ctrl1.entry <- c(ctrl1entry.1, ctrl1entry.2)
  tx1.entry <- c(tx1entry.1, tx1entry.2)
  
  comp1 <- c(ctrl1.time, tx1.time)
  censTime1 <- comp1[order(comp1)][nevents]
  
  c1stat <- ctrl1.time <= censTime1
  t1stat <- tx1.time <= censTime1
  
  nu_st2 <- round_to_fraction((censTime1.2-nu_st), 500, 3)
  
  entry3 <- seq(nu_st2+1/rate,10, by=1/rate)
  ctrl5entry.2 <- entry3[seq(1,length(entry3), by=n.arms[3])]
  tx5entry.2 <- entry3[seq(2, length(entry3), by=n.arms[3])]
  
  ctrl5evs.2 <- rgamma(length(ctrl5entry.2), gshape, rate=grate)#rexp(length(ctrl5entry.2), null.hz)
  tx5evs.2 <- rgamma(length(tx5entry.2), gshape, rate=grate_tx)#rexp(length(tx5entry.2), tx.hz)
  
  ctrl5time.2 <- ctrl5entry.2+ctrl5evs.2
  tx5time.2 <- tx5entry.2 + tx5evs.2
  
  comp5 <- c(ctrl5time.1, tx5time.1, ctrl5time.2, tx5time.2)
  censTime2.2 <- comp5[order(comp5)][ntot[2]]
  
  inds1 <- which(ctrl5entry.2<censTime2.2)
  inds2 <- which(tx5entry.2<censTime2.2)
  
  ctrl5evs.2 <- ctrl5evs.2[inds1]
  tx5evs.2 <- tx5evs.2[inds2]
  ctrl5time.2 <- ctrl5time.2[inds1]
  tx5time.2 <- tx5time.2[inds2]
  ctrl5entry.2 <- ctrl5entry.2[inds1]
  tx5entry.2 <- tx5entry.2[inds2]
  
  ctrl5.evs <- c(ctrl.sh.evs, ctrl5evs.2)
  tx5.evs <- c(tx5evs.1, tx5evs.2)
  ctrl5.time <- c(ctrl5time.1, ctrl5time.2)
  tx5.time <- c(tx5time.1, tx5time.2)
  ctrl5.entry <- c(ctrl5entry.1, ctrl5entry.2)
  tx5.entry <- c(tx5entry.1, tx5entry.2)
  
  comp5 <- c(ctrl5.time, tx5.time)
  censTime5 <- comp5[order(comp5)][nevents]
  
  c5stat <- ctrl5.time <= censTime5
  t5stat <- tx5.time <= censTime5
  
  evTimesC1 <- ctrl1.evs * c1stat + (censTime1 - ctrl1.entry) * (1 - c1stat)
  evTimesT1 <- tx1.evs * t1stat + (censTime1 - tx1.entry) * (1 - t1stat)
  evTimesC5 <- ctrl5.evs * c5stat + (censTime5 - ctrl5.entry) * (1 - c5stat)
  evTimesT5 <- tx5.evs * t5stat + (censTime5 - tx5.entry) * (1 - t5stat)
  
  sh.ctrl[i] <- sum(ctrl1time.2 <= censTime1)
  ctrl.evs[i] <- sum(c5stat)+sum(ctrl1time.1<=censTime1)
  
  test1c <- summary(coxph(Surv(time = c(evTimesC1,evTimesT1) , event = c(c1stat, t1stat))
                          ~c(rep(0, length(evTimesC1)), rep(1, length(evTimesT1)))))
  
  test5c <- summary(coxph(Surv(time = c(evTimesC5,evTimesT5) , event = c(c5stat, t5stat))~
                            c(rep(0, length(evTimesC5)), rep(1, length(evTimesT5)))))
  
  ts1[i] <- test1c$coefficients[4]
  ts5[i] <- test5c$coefficients[4]
  
};beep()


mean(sh.ctrl)/(2*mean(ctrl.evs))
cor(ts1, ts5)

t <- c(.001, .005, .01, .05, .1, .2, .3, .4, .5, .6, .65)
c1 <- c(.4569915,.4230184,.3926285,.2928085,.2289988, 0.1522665, .1012968, .06530908, 
        .03761772, .03767777, .007556579)
c2 <- c(0.4430902, 0.4208719, 0.3967856, .3323371, 0.2657181, 0.2051999, 0.1440322, 0.1049914, 
        0.05786682, .05727152, 0.006542993)

plot(t, c1, type="l")
lines(t, c2, col=2)

