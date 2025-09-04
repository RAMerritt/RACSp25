library(survival)
library(beepr)

A <- 2
t.length <- 2.33
nsim <- 50000
arm2.start <- seq(0,2,by=.2)
null.rate <-.693
hz1 <- hz2 <- .75

control.rate <- ceiling(500/(1+A))
expt.rate <- floor(500*A/(1+A))
n.control <- ceiling(control.rate*(t.length+2))
n.fin <- ceiling(control.rate*t.length)
n.expt <- ceiling(expt.rate*t.length)
t1.rate <- null.rate*hz1
t2.rate <- null.rate*hz2

control.entry <- seq(0, (t.length+2), length=n.control)
control.entry[1]<-.00001
expt.entry <- seq(0, t.length, length=n.expt)
expt.entry[1] <-.00001
ctr.v <- c(rep(0, n.fin), rep(1, n.expt))

ctrl.list <- list()
tx1.list <- list()
tx2.list <- list()
ctrl.time <- list()
tx1.time <- list()
tx2.time <- list()

set.seed(3511)

for (i in 1:nsim) {
  ctrl.time[[i]] <- rexp(n.control, rate=null.rate)
  tx1.time[[i]] <- rexp(n.expt, rate = t1.rate)
  tx2.time[[i]] <- rexp(n.expt, rate = t2.rate)
  
  ctrl.list[[i]] <- ctrl.time[[i]]+control.entry
  tx1.list[[i]] <- tx1.time[[i]] +expt.entry
  tx2.list[[i]] <- tx2.time[[i]] +expt.entry
};beep()

stat1 <- rep(NA, nsim)
p.expt1 <- rep(NA, nsim)

for (i in 1:nsim) {
  censc <- ifelse(ctrl.list[[i]][1:n.fin]>t.length, 0, 1)
  censt <- ifelse(tx1.list[[i]]>t.length,0,1)
  ttec <- ifelse(censc, ctrl.time[[i]][1:n.fin],t.length-control.entry[1:n.fin])
  ttet <- ifelse(censt, tx1.time[[i]], t.length-expt.entry)
  test <- summary(coxph(Surv(c(ttec, ttet), c(censc, censt))~ctr.v))
  stat1[i] <- test$coef[4]
  p.expt1[i] <- test$coef[5]
  if (i %% 5000 == 0) {
    cat("Simulation", i*100/nsim,"%", "donezo!", "\n")
  } 
};beep()

stat2 <- list()
p.expt2 <- list()

for (i in 1:length(arm2.start)) {
  stat2[[i]] <- rep(NA, nsim)
  p.expt2[[i]] <- rep(NA, nsim)
  t <- which(control.entry>arm2.start[i])[[1]]
  t1 <- t+n.fin-1
  for (j in 1:nsim) {
    censc <- ifelse(ctrl.list[[j]][t:t1]>t.length+arm2.start[i], 0, 1)
    censt <- ifelse(tx2.list[[j]]>t.length,0,1)
    ttec <- ifelse(censc, ctrl.time[[j]][t:t1],t.length+arm2.start[i]-control.entry[t:t1])
    ttet <- ifelse(censt, tx2.time[[j]], t.length-expt.entry)
    test <- summary(coxph(Surv(c(ttec, ttet), c(censc, censt))~ctr.v))
    stat2[[i]][j] <- test$coef[4]
    p.expt2[[i]][j] <- test$coef[5]
  }

    cat("Simulation", i*100/11,"%", "donezo!", "\n")

};beep()

cor.vals <- rep(NA, length(arm2.start))

for (i in 1:length(arm2.start)) {
  cor.vals[i] <-cor(stat.hold, stat2[[i]])
  print(cor(stat.hold, stat2[[i]]))
}
