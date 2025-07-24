library(survival)
library(beepr)

control.rate <- 167
expt.rate <- 333

n.control <- 1045
n.fin <- 545
n.expt <- 545
nsim <- 50000
arm2.start <- seq(0,2,by=.2)
t.length <- 2.18

control.entry <- seq(0, (t.length+2), length=n.control)
control.entry[1]<-.00001
expt.entry <- seq(0, t.length, length=n.expt)
expt.entry[1] <-.00001
ctr.v <- c(rep(0, n.fin), rep(1, n.expt))

#testing code
# control.event <- rexp(n.control, rate=0.693)
# expt.event <- rexp(n.expt, rate = .51975)
# expt.event2 <- rexp(n.expt, rate = .51975)
# 
# ctrl.event.time <- control.entry+control.event
# expt.event.time <- expt.entry+expt.event
# expt.event.time2 <- expt.entry+expt.event2
# 
# ctrl.censor <- as.integer(ctrl.event.time < 2.36)
# expt.censor <- as.integer(expt.event.time < 2.36)
# 
# ctr.v.tx1 <- c(rep(0, n.control), rep(1, n.expt))
# 
# for (i in 1:length(ctrl.event.time)) {
#   if(ctrl.censor[i]==0){
#     ctrl.event.time[i] <- 2.36
#   }
# }
# 
# event.time1 <- c(ctrl.event.time, expt.event.time)
# censor.time1 <- c(ctrl.censor, expt.censor)
# 
# survdiff(Surv(event.time1, censor.time1) ~ ctr.v.tx1)$chisq

ctrl.list <- list()
tx1.list <- list()
tx2.list <- list()
ctrl.time <- list()
tx1.time <- list()
tx2.time <- list()

set.seed(123)

for (i in 1:nsim) {
  ctrl.list[[i]] <- rexp(n.control, rate=0.693)+control.entry
  tx1.list[[i]] <- rexp(n.expt, rate = 0.51975)+expt.entry
  tx2.list[[i]] <- rexp(n.expt, rate = 0.51975)+expt.entry
  # ctrl.list[[i]] <- rep(NA, n.control)
  # tx1.list[[i]] <- rep(NA, n.expt)
  # tx2.list[[i]] <- rep(NA, n.expt)
  # 
  # ctrl.list[[i]] <- ctrl.time[[i]]+control.entry
  # tx1.list[[i]] <- tx1.time[[i]] +expt.entry
  # tx2.list[[i]] <- tx2.time[[i]] +expt.entry
};beep()

diff.chi <- rep(NA, nsim)
p.expt1 <- rep(NA, nsim)


for (i in 1:nsim) {
  # ctrl.list.tmp <- sapply(ctrl.list[[i]][1:789], function(x) if(x>2.36) x=2.36 else x=x)
  # tx1.list.tmp <- sapply(tx1.list[[i]], 
  #                        function(x) if(x>2.36) x=2.36 else x=x)
  # event.time <- c(sapply(ctrl.list[[i]][1:789], 
  #                        function(x) if(x>2.36) x=2.36 else x=x), sapply(tx1.list[[i]], 
  #                                                                        function(x) if(x>2.36) x=2.36 else x=x))
  # censor.time <- c(ctrl.time[[i]][1:789], tx1.time[[i]])
  diff.chi[i] <- survdiff(Surv(c(sapply(ctrl.list[[i]][1:n.fin], 
                                        function(x) if(x>t.length) x=t.length else x=x), sapply(tx1.list[[i]], 
                                                                                        function(x) if(x>t.length) x=t.length else x=x)), c(sapply(ctrl.list[[i]][1:n.fin], function(x) if(x>t.length) x=0 else x=1), sapply(tx1.list[[i]], function(x) if(x>t.length) x=0 else x=1))) ~ ctr.v)$chisq
  p.expt1[i] <- survdiff(Surv(c(sapply(ctrl.list[[i]][1:n.fin], 
                function(x) if(x>t.length) x=t.length else x=x), sapply(tx1.list[[i]], 
                function(x) if(x>t.length) x=t.length else x=x)), c(sapply(ctrl.list[[i]][1:n.fin], 
                function(x) if(x>t.length) x=0 else x=1), sapply(tx1.list[[i]], function(x) if(x>t.length) x=0 else x=1))) ~ ctr.v)$pval
};beep()

diff.chi2 <- list()
p.expt2 <- list()

tx2.dat <- list()
tx2.cens <- list()

for (i in 1:nsim) {
  tx2.dat[[i]] <- sapply(tx2.list[[i]], function(x) if(x>t.length) x=t.length else x=x)
  tx2.cens[[i]] <- sapply(tx2.list[[i]], function(x) if(x>t.length) x=0 else x=1)
};beep()

for (i in 1:length(arm2.start)) {
  diff.chi2[[i]] <- rep(NA, nsim)
  p.expt2[[i]] <- rep(NA, nsim)
  t <- which(control.entry>arm2.start[i])[[1]]
  t1 <- t+n.fin-1
  for (j in 1:nsim) {
    # event.time <- c(sapply(ctrl.list[[j]][t:t1], function(x) if(x>(2.36+arm2.start[i])) x=(2.36+arm2.start[i]) else x=x), sapply(tx2.list[[j]][s:s1], function(x) if(x>(2.36+arm2.start[i])) x=(2.36+arm2.start[i]) else x=x))
    # censor.time <- c(ctrl.time[[j]][t:t1], tx2.time[[j]][s:s1])
    #c.tx <- c(rep(0, length(t:789)), rep(1, length(s:391)))
    diff.chi2[[i]][j] <- survdiff(Surv(c(sapply(ctrl.list[[j]][t:t1], function(x) if(x>(t.length+arm2.start[i])) x=t.length else x=x-arm2.start[i]), tx2.dat[[j]]), c(sapply(ctrl.list[[j]][t:t1], function(x) if(x>(t.length+arm2.start[i])) x=0 else x=1), tx2.cens[[j]])) ~ ctr.v)$chisq
    p.expt2[[i]][j] <- survdiff(Surv(c(sapply(ctrl.list[[j]][t:t1], function(x) if(x>(t.length+arm2.start[i])) x=t.length else x=x-arm2.start[i]), tx2.dat[[j]]), c(sapply(ctrl.list[[j]][t:t1], function(x) if(x>(t.length+arm2.start[i])) x=0 else x=1), tx2.cens[[j]])) ~ ctr.v)$pval
  }
};beep()

# storage <- rep(NA, 200)
# storage2 <- rep(NA, 200)
# 
# for (j in 1:200) {
#   storage[j] <- survdiff(Surv(c(sapply(ctrl.list[[j]][1:789], function(x) if(x>2.36) x=2.36 else x=x), sapply(tx2.list[[j]], function(x) if(x>2.36) x=2.36 else x=x)), c(sapply(ctrl.list[[j]][1:789], function(x) if(x>2.36) x=0 else x=1), sapply(tx2.list[[j]], function(x) if(x>2.36) x=0 else x=1))) ~ ctr.v.tx1)$chisq
#   storage2[j] <- survdiff(Surv(c(sapply(ctrl.list[[j]][1:789], function(x) if(x>2.36) x=2.36 else x=x), sapply(tx2.list[[j]], function(x) if(x>2.36) x=2.36 else x=x)), c(sapply(ctrl.list[[j]][1:789], function(x) if(x>2.36) x=0 else x=1), sapply(tx2.list[[j]], function(x) if(x>2.36) x=0 else x=1))) ~ ctr.v.tx1)$pval
# };beep()
# 
# for (j in 1:nsim) {
#   diff.chi2[[1]][j] <- survdiff(Surv(c(sapply(ctrl.list[[j]][1:n.fin], function(x) if(x>2.36) x=2.36 else x=x), tx2.dat[[j]]), c(sapply(ctrl.list[[j]][1:n.fin], function(x) if(x>2.36) x=0 else x=1), tx2.cens[[j]])) ~ ctr.v)$chisq
#   p.expt2[[1]][j] <- survdiff(Surv(c(sapply(ctrl.list[[j]][1:n.fin], function(x) if(x>2.36) x=2.36 else x=x), tx2.dat[[j]]), c(sapply(ctrl.list[[j]][1:n.fin], function(x) if(x>2.36) x=0 else x=1), tx2.cens[[j]])) ~ ctr.v)$pval
# };beep()
# 
# for (j in 1:nsim) {
#   diff.chi2[[2]][j] <- survdiff(Surv(c(sapply(ctrl.list[[j]][68:856], function(x) if(x>2.56) x=2.36 else x=x-.2), sapply(tx2.list[[j]], function(x) if(x>2.36) x=2.36 else x=x)), c(sapply(ctrl.list[[j]][68:856], function(x) if(x>2.56) x=0 else x=1), sapply(tx2.list[[j]], function(x) if(x>2.36) x=0 else x=1))) ~ ctr.v.tx1)$chisq
#   p.expt2[[2]][j] <- survdiff(Surv(c(sapply(ctrl.list[[j]][68:856], function(x) if(x>2.56) x=2.36 else x=x-.2), sapply(tx2.list[[j]], function(x) if(x>2.36) x=2.36 else x=x)), c(sapply(ctrl.list[[j]][68:856], function(x) if(x>2.56) x=0 else x=1), sapply(tx2.list[[j]], function(x) if(x>2.36) x=0 else x=1))) ~ ctr.v.tx1)$pval
# };beep()
# 
# cor(diff.chi, diff.chi2[[2]])


cor.vals <- rep(NA, length(arm2.start))

for (i in 1:length(arm2.start)) {
  cor.vals[i] <- cor(sqrt(diff.chi), sqrt(diff.chi2[[i]]))
}

null.p.t1 
null.p.t2 
stat.t1 
stat.t2 
sig.p.t1 
sig.p.t2 
sig.stat.t1 
sigstat.t2 
null2.p.t1 
null2.p.t2 
null2.stat.t1 
null2.stat.t2 
sig2.p.t1
sig2.p.t2 
sig2.stat.t1 
sig2.stat.t2 
null3.p.t1
null3.p.t2 
null3.stat.t1 
null3.stat.t2 
sig3.p.t1 <- p.expt1
sig3.p.t2 <- p.expt2
sig3.stat.t1 <- diff.chi
sig3.stat.t2 <- diff.chi2

# for (i in 1:nsim) {
#   event.time <- c(ctrl.list[[i]], tx2.list[[i]])
#   censor.time <- c(ctrl.time[[i]], tx2.time[[i]])
#   diff.chi2[i] <- survdiff(Surv(event.time, censor.time) ~ ctr.v.tx1)$chisq
# }
# 
# cor(diff.chi, diff.chi2)

for (i in 1:length(sig.p.t2)) {
  print(mean(sig.p.t2[[i]]<.05))
}

for (i in 1:11) {
  power <- rep(NA, nsim)
  for (j in 1:nsim) {
    if (sig3.p.t1[j]<.025 & sig3.p.t2[[i]][j]<.025){
      power[j] <- 1
    }
    else {
      power[j] <-0
    }
  }
  print(mean(power))
}
