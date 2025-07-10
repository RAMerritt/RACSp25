library(survival)

control.rate <- 334
expt.rate <- 166

n.control <- 789
n.expt <- 391
nsim <- 50000
arm2.start <- seq(0,2,by=.2)

control.entry <- seq(0, 2.36, length=n.control)
control.entry[1]<-.00001
expt.entry <- seq(0, 2.36, length=n.expt)
expt.entry[1] <-.00001
exp2.start <- which(expt.entry>arm2.start[l])

control.event <- rexp(n.control, rate=0.693)
expt.event <- rexp(n.expt, rate = .51975)
expt.event2 <- rexp(n.expt, rate = .51975)

ctrl.event.time <- control.entry+control.event
expt.event.time <- expt.entry+expt.event
expt.event.time2 <- expt.entry+expt.event2

ctrl.censor <- as.integer(ctrl.event.time < 2.36)
expt.censor <- as.integer(expt.event.time < 2.36)

ctr.v.tx1 <- c(rep(0, n.control), rep(1, n.expt))

for (i in 1:length(ctrl.event.time)) {
  if(ctrl.censor[i]==0){
    ctrl.event.time[i] <- 2.36
  }
}

event.time1 <- c(ctrl.event.time, expt.event.time)
censor.time1 <- c(ctrl.censor, expt.censor)

survdiff(Surv(event.time1, censor.time1) ~ ctr.v.tx1)$chisq

ctrl.list <- list()
tx1.list <- list()
tx2.list <- list()
ctrl.time <- list()
tx1.time <- list()
tx2.time <- list()

set.seed(123)

for (i in 1:nsim) {
  ctrl.list[[i]] <- rexp(n.control, rate=0.693)
  tx1.list[[i]] <- rexp(n.expt, rate = .51975)
  tx2.list[[i]] <- rexp(n.expt, rate = .51975)
  ctrl.time[[i]] <- as.integer(ctrl.list[[i]]<2.36)
  tx1.time[[i]] <- as.integer(tx1.list[[i]] < 2.36)
  tx2.time[[i]] <- as.integer(tx2.list[[i]] < 2.36)
  for (j in 1:n.control) {
    if (ctrl.list[[i]][j]+control.entry[j]<2.36){
      ctrl.list[[i]][j] <- ctrl.list[[i]][j]+control.entry[j]
    }
    else {
      ctrl.list[[i]][j] <- 2.36
    }
  }
  for (j in 1:n.expt) {
    if (tx1.list[[i]][j] +expt.entry[j]<2.36){
      tx1.list[[i]][j] <- tx1.list[[i]][j] +expt.entry[j]
    }
    else {
      tx1.list[[i]][j] <- 2.36
    }
    if (tx2.list[[i]][j] +expt.entry[j]<2.36){
      tx2.list[[i]][j] <- tx2.list[[i]][j] +expt.entry[j]
    }
    else {
      tx2.list[[i]][j] <- 2.36
    }
  }
}

diff.chi <- rep(NA, nsim)

for (i in 1:nsim) {
  event.time <- c(ctrl.list[[i]], tx1.list[[i]])
  censor.time <- c(ctrl.time[[i]], tx1.time[[i]])
  diff.chi[i] <- survdiff(Surv(event.time, censor.time) ~ ctr.v.tx1)$chisq
}

diff.chi2 <- list()

for (i in 1:length(arm2.start)) {
  diff.chi2[[i]] <- rep(NA, nsim)
  s <- which(expt.entry>arm2.start[i])[[1]]
  t <- which(control.entry>arm2.start[i])[[1]]
  for (j in 1:nsim) {
    event.time <- c(ctrl.list[[j]][t:789], tx2.list[[j]][s:391])
    censor.time <- c(ctrl.time[[j]][t:789], tx2.time[[j]][s:391])
    c.tx <- c(rep(0, length(t:789)), rep(1, length(s:391)))
    diff.chi2[[i]][j] <- survdiff(Surv(event.time, censor.time) ~ c.tx)$chisq
  }
}

cor.vals <- rep(NA, length(arm2.start))

for (i in 1:length(arm2.start)) {
  cor.vals[i] <- cor(diff.chi, diff.chi2[[i]])
}

for (i in 1:nsim) {
  event.time <- c(ctrl.list[[i]], tx2.list[[i]])
  censor.time <- c(ctrl.time[[i]], tx2.time[[i]])
  diff.chi2[i] <- survdiff(Surv(event.time, censor.time) ~ ctr.v.tx1)$chisq
}

cor(diff.chi, diff.chi2)
