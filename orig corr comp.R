library(survival)
library(beepr)

A <- 1
rec.time <- 2
fu.time <- 1.5
rec.rate <- 500
nsim <- 15000
arm2.start <- seq(0,2,by=.2)
null.rate <-.693
hz1 <- hz2 <- .75

t.time <- rec.time+fu.time
control.rate <- ceiling(rec.rate/(1+A))
expt.rate <- floor(rec.rate*A/(1+A))
n.control <- ceiling(control.rate*(2*rec.time))
n.fin <- ceiling(control.rate*rec.time)
n.expt <- ceiling(expt.rate*rec.time)
t1.rate <- null.rate*hz1
t2.rate <- null.rate*hz2

control.entry <- seq(0, (2*rec.time), length=n.control)
control.entry[1]<-.00001
expt.entry <- seq(0, rec.time, length=n.expt)
expt.entry[1] <-.00001
ctr.v <- c(rep(0, n.fin), rep(1, n.expt))

ctrl.list <- list()
tx1.list <- list()
tx2.list <- list()
ctrl.time <- list()
tx1.time <- list()
tx2.time <- list()

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
  censc <- ifelse(ctrl.list[[i]][1:n.fin]>t.time, 0, 1)
  censt <- ifelse(tx1.list[[i]]>t.time,0,1)
  ttec <- ifelse(censc, ctrl.time[[i]][1:n.fin],t.time-control.entry[1:n.fin])
  ttet <- ifelse(censt, tx1.time[[i]], t.time-expt.entry)
  test <- summary(coxph(Surv(c(ttec, ttet), c(censc, censt))~ctr.v))
  stat1[i] <- test$coef[4]
  p.expt1[i] <- test$coef[5]
  if (i %% 5000 == 0) {
    cat("Simulation", i*100/nsim,"%", "donezo!", "\n")
  } 
};beep()


stat2 <- list()
p.expt2 <- list()
ctrl.ev <- list()
sh.ctrl.ev <- list()

for (i in 1:length(arm2.start)) {
  stat2[[i]] <- rep(NA, nsim)
  p.expt2[[i]] <- rep(NA, nsim)
  t <- which(control.entry>arm2.start[i])[[1]]
  t1 <- t+n.fin-1
  ctrl.ev[[i]] <- rep(NA, nsim)
  sh.ctrl.ev[[i]] <- rep(NA, nsim)
  for (j in 1:nsim) {
    censc <- ifelse(ctrl.list[[j]][t:t1]>t.time+arm2.start[i], 0, 1)
    censt <- ifelse(tx2.list[[j]]>t.time,0,1)
    ttec <- ifelse(censc, ctrl.time[[j]][t:t1],t.time+arm2.start[i]-control.entry[t:t1])
    ttet <- ifelse(censt, tx2.time[[j]], t.time-expt.entry)
    test <- summary(coxph(Surv(c(ttec, ttet), c(censc, censt))~ctr.v))
    stat2[[i]][j] <- test$coef[4]
    p.expt2[[i]][j] <- test$coef[5]
    ctrl.ev[[i]][j] <- sum(censc)+ifelse(t==1, 0, sum(ifelse(ctrl.list[[j]][1:t]>t.time, 0, 1))-1)
    sh.ctrl.ev[[i]][j] <- sum(ifelse(ctrl.list[[j]][t:n.fin]>t.time, 0, 1))
  }
  
  cat("Simulation", i*100/11,"%", "donezo!", "\n")
  
};beep()

cor.vals <- rep(NA, length(arm2.start))
rat <- rep(NA, length(arm2.start))

for (i in 1:length(arm2.start)) {
  cor.vals[i] <-cor(stat1, stat2[[i]])
  rat[i] <- mean(sh.ctrl.ev[[i]]/ctrl.ev[[i]])*A/(A+1)
}

tab <- cbind(cor.vals, rat)
head(tab,11)

