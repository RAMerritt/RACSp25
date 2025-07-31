library(survival)
library(beepr)

A <- 2
t.length <- 2.33
nsim <- 10000
arm2.start <- seq(0,2,by=.2)
null.rate <-.693
hz1 <- hz2 <- .75
t1.rate <- null.rate*hz1
t2.rate <- null.rate*hz2

crate1 <- ceiling(500/(1+A))
crate2 <- 100
txrate1 <- floor(500*A/(1+A))
txrate2 <- 200

stat1 <- list()
p.expt1 <- list()
stat2 <- list()
p.expt2 <- list()

set.seed(123)

for (j in 1:length(arm2.start)) {
  stat1[[j]] <- rep(NA, nsim)
  p.expt1[[j]] <- rep(NA, nsim)
  stat2[[j]] <- rep(NA, nsim)
  p.expt2[[j]] <- rep(NA, nsim)
  
start2 <-arm2.start[j]

control.entry <- c(seq(0, start2, length=round(crate1*start2)), 
                   seq(start2, t.length, length=round(crate2*(t.length-start2))),
                   seq(t.length, t.length+2, length=round(crate1*start2)))
expt1.entry <- c(seq(0, start2, length=round(txrate1*start2)), 
                 seq(start2, t.length, length=round(txrate2*(t.length-start2))))
expt2.entry <- c(seq(start2, t.length, length=round(txrate2*(t.length-start2))),
                 seq(t.length, t.length+start2, length=round(txrate1*start2)))
ctr.v <- c(rep(0, round(crate1*start2) + round(crate2*(t.length-start2))),
           rep(1, round(txrate1*start2)+round(txrate2*(t.length-start2))))

comp2 <- round(crate1*start2) + round(crate2*(t.length-start2))
comp1 <- which(control.entry>arm2.start[j])[[1]]
c3 <- comp1 + comp2-1

ctrl.list <- list()
tx1.list <- list()
tx2.list <- list()
ctrl.time <- list()
tx1.time <- list()
tx2.time <- list()

for (i in 1:nsim) {
  ctrl.time[[i]] <- c(rexp(round(crate1*start2), null.rate),
                      rexp(round(crate2*(t.length-start2)), rate=null.rate),
                      rexp(round(crate1*start2), rate=null.rate))
  tx1.time[[i]] <- c(rexp(round(txrate1*start2), t1.rate),
                     rexp(round(txrate2*(t.length-start2)), t1.rate))
  tx2.time[[i]] <- c(rexp(round(txrate2*(t.length-start2)), t2.rate),
                     rexp(round(txrate1*start2), t2.rate))
  
  ctrl.list[[i]] <- ctrl.time[[i]]+control.entry
  tx1.list[[i]] <- tx1.time[[i]] +expt1.entry
  tx2.list[[i]] <- tx2.time[[i]] +expt2.entry

  censc <- ifelse(ctrl.list[[i]][1:comp2]>t.length, 0, 1)
  censt <- ifelse(tx1.list[[i]]>t.length,0,1)
  ttec <- ifelse(censc, ctrl.time[[i]][1:comp2],t.length-control.entry[1:comp2])
  ttet <- ifelse(censt, tx1.time[[i]], t.length-expt1.entry)
  test <- summary(coxph(Surv(c(ttec, ttet), c(censc, censt))~ctr.v))
  stat1[[j]][i] <- test$coef[4]
  p.expt1[[j]][i] <- test$coef[5]

  censc <- ifelse(ctrl.list[[i]][comp1:c3]>t.length+start2, 0, 1)
  censt <- ifelse(tx2.list[[i]]>t.length+start2,0,1)
  ttec <- ifelse(censc, ctrl.time[[i]][comp1:c3],t.length+start2-control.entry[comp1:c3])
  ttet <- ifelse(censt, tx2.time[[i]], t.length+start2-expt2.entry)
  test <- summary(coxph(Surv(c(ttec, ttet), c(censc, censt))~ctr.v))
  stat2[[j]][i] <- test$coef[4]
  p.expt2[[j]][i] <- test$coef[5]
}
  cat("Simulation", j*10,"%", "donezo!", "\n")
};beep()

hist(p.expt1)
hist(p.expt2)
cor(stat1, stat2)

for (i in 1:11) {
  print(cor(stat1[[i]], stat2[[i]]))
}
for (i in 1:11) {
  print(mean(p.expt2[[i]]<.05&p.expt1[[i]]<.05))
  
}
