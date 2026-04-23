library(survival)
library(survminer)
library(beepr)

nsim <- 5000
null.hz <- .693
tx.hz <- .75*null.hz
nevents <- 524
ntot <- ceiling(nevents*1.5)
halfn <- ntot/2
gpc <- .693
gp1 <- gpc*.75
gp2 <- 1
t_st <- 300/524
rte <- 500

n.arms <- c(5,6,2)
#number of arms at any given time
first_push <- round(ntot*t_st)

ctrl.n.t1 <- ceiling(first_push/2)
t1_first <- floor(first_push/2)
ctrl.n.shared <- floor((ntot-first_push)/2)
t.n.shared <- ceiling((ntot-first_push)/2)
ctrl.n.t2 <- ceiling(ntot/2)-ctrl.n.shared
t2_second <- floor(ntot/2)-t.n.shared

entry <- seq(1/rte,100, by=1/rte)
ctrl1.entry <- c(entry[seq(1, n.arms[1]*ctrl.n.t1, by=n.arms[1])], 
                 entry[seq(n.arms[1]*ctrl.n.t1+1, n.arms[1]*ctrl.n.t1+n.arms[2]*ctrl.n.shared, by=n.arms[2])])
tx1.entry <- c(entry[seq(2,n.arms[1]*t1_first, by=n.arms[1])],
               entry[seq(n.arms[1]*t1_first+2,n.arms[1]*t1_first+n.arms[2]*t.n.shared, by=n.arms[2])])
ctrl5.entry <- c(entry[seq(1, n.arms[2]*ctrl.n.shared, by=n.arms[2])], 
                 entry[seq(n.arms[2]*ctrl.n.shared+1, n.arms[2]*ctrl.n.shared+n.arms[3]*ctrl.n.t2, by=n.arms[3])])
tx5.entry <- c(entry[seq(n.arms[2],n.arms[2]*t.n.shared, by=n.arms[2])],
               entry[seq(n.arms[2]*t.n.shared+2,n.arms[2]*t.n.shared+n.arms[3]*t2_second, by=n.arms[3])])

#stat <- rep(NA, nsim)
p.expt_1c <- rep(NA, nsim)
p.expt_5c <- rep(NA, nsim)
p.expt_15 <- rep(NA, nsim)

sh.ctrl <- rep(NA, nsim)
ctrl.evs <- rep(NA, nsim)

Tstat1 <- rep(NA, nsim)
Tstat5 <- rep(NA, nsim)

for (i in 1:nsim) {
  ctrl.shared.evs <- rgamma(ctrl.n.shared, shape = gp2, rate = gpc) 
  ctrl1.evs <- c(rgamma(ctrl.n.t1, shape = gp2, rate = gpc), ctrl.shared.evs) 
  ctrl5.evs <- c(ctrl.shared.evs, rgamma(ctrl.n.t2, shape = gp2, rate = gpc))
  tx1.evs <- rgamma(length(tx1.entry), shape = gp2, rate = gp1)
  tx5.evs <- rgamma(length(tx5.entry), shape = gp2, rate = gp1)
  
  ctrl1.time <- ctrl1.entry+ctrl1.evs
  ctrl5.time <- ctrl5.entry+ctrl5.evs
  tx1.time <- tx1.entry+tx1.evs
  tx5.time <- tx5.entry+tx5.evs
  
  comp1 <- c(ctrl1.time, tx1.time)
  censTime1 <- comp1[order(comp1)][nevents]
  comp5 <- c(ctrl5.time, tx5.time)
  censTime5 <- comp5[order(comp5)][nevents]
  
  c1stat <- ctrl1.time <= censTime1
  t1stat <- tx1.time <= censTime1
  c5stat <- ctrl5.time <= censTime5
  t5stat <- tx5.time <= censTime5
  
  sh1 <- (ctrl.n.t1+1):(ctrl.n.shared+ctrl.n.t1)
  
  sh.ctrl[i] <- sum(ctrl1.time[sh1]<=censTime1)
  ctrl.evs[i] <- sum(ctrl1.time[-sh1]<=censTime1)+sum(c5stat)
  
  evTimesC1 <- ctrl1.evs * c1stat + (censTime1 - ctrl1.entry) * (1 - c1stat)
  evTimesT1 <- tx1.evs * t1stat + (censTime1 - tx1.entry) * (1 - t1stat)
  evTimesC5 <- ctrl5.evs * c5stat + (censTime5 - ctrl5.entry) * (1 - c5stat)
  evTimesT5 <- tx5.evs * t5stat + (censTime5 - tx5.entry) * (1 - t5stat)
  
  test1c <- summary(coxph(Surv(time = c(evTimesC1,evTimesT1) , event = c(c1stat, t1stat))
                          ~c(rep(0, length(evTimesC1)), rep(1, length(evTimesT1)))))
  
  test5c <- summary(coxph(Surv(time = c(evTimesC5,evTimesT5) , event = c(c5stat, t5stat))~
                            c(rep(0, length(evTimesC5)), rep(1, length(evTimesT5)))))
  
  Tstat1[i] <- test1c$coef[4]
  Tstat5[i] <- test5c$coef[4]
  
};beep()


t_st
mean(sh.ctrl/ctrl.evs)*A/(A+1)
cor(Tstat1, Tstat5)

t <- c(1/524, 0.01908397,0.09541985,0.1908397, 0.2862595,0.3816794,0.4770992,
       0.5725191,0.6679389,0.7633588,0.8587786,0.9541985)
c1 <- c(0.4960042, 0.4727439,0.3838493,0.2954787,0.2234227,0.1671546,0.1205747,
        0.08372938,0.05341577,0.0307505,0.01369667,0.003089916)
c2 <- c(0.4657339, 0.4608015,0.3987601,0.3406108,0.277623,0.2222681,0.1832036,
        0.1187485,0.08677928,0.06356928,0.02448023,-0.005740668)

plot(t, c1, type="l", main="Fixed N exp corr")
lines(t, c2, col=2, lty="dashed")

