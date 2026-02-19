library(survival)
library(survminer)

nsim <- 500
null.hz <- .3
tx.hz <- .75*null.hz
nevents <- 524
ntot <- ceiling(nevents*1.5)
halfn <- ntot/2

n.arms <- c(5,6,2)
#number of arms at any given time
first_push <- round(ntot*7/8)

ctrl.n.t1 <- ceiling(first_push/2)
t1_first <- floor(first_push/2)
ctrl.n.shared <- floor((ntot-first_push)/2)
t.n.shared <- ceiling((ntot-first_push)/2)
ctrl.n.t2 <- ceiling(ntot/2)-ctrl.n.shared
t2_second <- floor(ntot/2)-t.n.shared

entry <- seq(1/500,6, by=1/500)
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

slow_null_ct_t <- list()
slow_null_tx_t <- list()
fast_null_ct_t <- list()
fast_null_tx_t <- list()

CI_slow_null_ct_l <- list()
CI_slow_null_ct_u <- list()
CI_slow_null_tx_l <- list()
CI_slow_null_tx_u <- list()
CI_fast_null_ct_l <- list()
CI_fast_null_ct_u <- list()
CI_fast_null_tx_l <- list()
CI_fast_null_tx_u <- list()

t_length <- data.frame(arm1=0, arm2=0)

mdiff1 <- rep(NA, nsim)
mdiff5 <- rep(NA, nsim)

surv1T <- list()
surv5T <- list()

diff1 <- list()
diff5 <- list()
time1 <- list()
time5 <- list()

for (i in 1:nsim) {
  ctrl.shared.evs <- rexp(ctrl.n.shared, null.hz)
  ctrl1.evs <- c(rexp(ctrl.n.t1, null.hz), ctrl.shared.evs)
  ctrl5.evs <- c(ctrl.shared.evs, rexp(ctrl.n.t2, null.hz))
  tx1.evs <- rexp(length(tx1.entry), tx.hz)
  tx5.evs <- rexp(length(tx5.entry), tx.hz)
  
  ctrl1.time <- ctrl1.entry+ctrl1.evs
  ctrl5.time <- ctrl5.entry+ctrl5.evs
  tx1.time <- tx1.entry+tx1.evs
  tx5.time <- tx5.entry+tx5.evs
  
  comp1 <- c(ctrl1.time, tx1.time)
  censTime1 <- comp1[order(comp1)][nevents]
  comp5 <- c(ctrl5.time, tx5.time)
  censTime5 <- comp5[order(comp5)][nevents]
  
  c1stat <- ctrl1.time < censTime1
  t1stat <- tx1.time < censTime1
  c5stat <- ctrl5.time < censTime5
  t5stat <- tx5.time < censTime5
  
  evTimesC1 <- ctrl1.evs * c1stat + (censTime1 - ctrl1.entry) * (1 - c1stat)
  evTimesT1 <- tx1.evs * t1stat + (censTime1 - tx1.entry) * (1 - t1stat)
  evTimesC5 <- ctrl5.evs * c5stat + (censTime5 - ctrl5.entry) * (1 - c5stat)
  evTimesT5 <- tx5.evs * t5stat + (censTime5 - tx5.entry) * (1 - t5stat)
  
  t_length[i,] <- cbind(censTime1, censTime5)
  
  test1c <- summary(coxph(Surv(time = c(evTimesC1,evTimesT1) , event = c(c1stat, t1stat))
                          ~c(rep(0, length(evTimesC1)), rep(1, length(evTimesT1)))))
  
  test5c <- summary(coxph(Surv(time = c(evTimesC5,evTimesT5) , event = c(c5stat, t5stat))~
                            c(rep(0, length(evTimesC5)), rep(1, length(evTimesT5)))))
  
  test15 <- summary(coxph(Surv(time = c(evTimesT1,evTimesT5) , event = c(t1stat, t5stat))~
                                 c(rep(0, length(evTimesT1)), rep(1, length(evTimesT5)))))
  
  p.expt_1c[i] <- test1c$coef[5]
  p.expt_5c[i] <- test5c$coef[5]
  p.expt_15[i] <- test15$coef[5]
  
  surv1T[[i]] <- survfit(Surv(time = evTimesT1, event = t1stat)~1)
  surv5T[[i]] <- survfit(Surv(time = evTimesT5, event = t5stat)~1)
  
  diff1[[i]] <- surv1T[[i]]$upper-surv1T[[i]]$lower
  diff5[[i]] <- surv5T[[i]]$upper-surv5T[[i]]$lower
  
  mdiff1[i] <- max(surv1T[[i]]$upper-surv1T[[i]]$lower)
  mdiff5[i] <- max(surv5T[[i]]$upper-surv5T[[i]]$lower)
  
  time1[[i]] <- round(surv1T[[i]]$time, digits=2)
  time5[[i]] <- round(surv5T[[i]]$time, digits=2)
  
  CI_slow_null_tx_l[[i]] <- surv1T[[i]]$lower
  CI_slow_null_tx_u[[i]] <- surv1T[[i]]$upper
  CI_fast_null_tx_l[[i]] <- surv5T[[i]]$lower
  CI_fast_null_tx_u[[i]] <- surv5T[[i]]$upper
  
}

mean(p.expt_1c<.05)
mean(p.expt_5c<.05)
mean(p.expt_15<.05)

t.test(t_length$arm1, t_length$arm2, paired = T)
t.test(mdiff1, mdiff5, paired = T)

hist(t_length$arm1)
hist(t_length$arm2)

ks.test(diff1[[4]], diff5[[4]])

plot(ecdf(diff1[[18]]))
lines(ecdf(diff5[[18]]), col=2)

# hist(diff1[[8]])
# hist(diff5[[8]])

slodiff <- cbind(unlist(diff1), round(unlist(time1), digits=2))
fastdiff <- cbind(unlist(diff5), round(unlist(time5), digits=2))

slodiff <- slodiff[complete.cases(slodiff),]
fastdiff <- fastdiff[complete.cases(fastdiff),]

sdind <- unique(slodiff[,2])
fdind <- unique(fastdiff[,2])

fst_CI <- data.frame()
slo_CI <- data.frame()

for (i in fdind) {
  ind <- which(fastdiff[,2]==i)
  qt_l <- quantile(fastdiff[ind,1], .025)
  qt_u <- quantile(fastdiff[ind,1], .975)
  fst_CI <- rbind(fst_CI, data.frame(qt_l,qt_u,i))
}

for (i in sdind) {
  ind <- which(slodiff[,2]==i)
  qt_l <- quantile(slodiff[ind,1], .025)
  qt_u <- quantile(slodiff[ind,1], .975)
  slo_CI <- rbind(slo_CI, data.frame(qt_l,qt_u,i))
}

plot(sort(slo_CI$i), sort(slo_CI$qt_l),type='l', col='red')
lines(sort(slo_CI$i), sort(slo_CI$qt_u), col='red')
lines(sort(fst_CI$i), sort(fst_CI$qt_l),type='l', col='blue')
lines(sort(fst_CI$i), sort(fst_CI$qt_u), col='blue')

# slo_ct_l <- cbind(unlist(slow_null_ct_t),unlist(CI_slow_null_ct_l))
# slo_ct_u <- cbind(unlist(slow_null_ct_t),unlist(CI_slow_null_ct_u))
# fst_ct_l <- cbind(unlist(fast_null_ct_t),unlist(CI_fast_null_ct_l))
# fst_ct_u <- cbind(unlist(fast_null_ct_t),unlist(CI_fast_null_ct_u))
slo_tx_l <- cbind(unlist(time1),unlist(CI_slow_null_tx_l))
slo_tx_u <- cbind(unlist(time1),unlist(CI_slow_null_tx_u))
fst_tx_l <- cbind(unlist(time5),unlist(CI_fast_null_tx_l))
fst_tx_u <- cbind(unlist(time5),unlist(CI_fast_null_tx_u))

# slo_ct_l <- slo_ct_l[complete.cases(slo_ct_l),]
# slo_ct_u <- slo_ct_u[complete.cases(slo_ct_u),]
# fst_ct_l <- fst_ct_l[complete.cases(fst_ct_l),]
# fst_ct_u <- fst_ct_u[complete.cases(fst_ct_u),]
slo_tx_l <- slo_tx_l[complete.cases(slo_tx_l),]
slo_tx_u <- slo_tx_u[complete.cases(slo_tx_u),]
fst_tx_l <- fst_tx_l[complete.cases(fst_tx_l),]
fst_tx_u <- fst_tx_u[complete.cases(fst_tx_u),]

#fnct <- unique(fst_ct_l[,1])
fntt <- unique(fst_tx_l[,1])
# snct <- unique(slo_ct_l[,1])
 sntt <- unique(slo_tx_l[,1])

#fst_ct_ci <- data.frame()
fst_tx_ci <- data.frame()
#slo_ct_ci <- data.frame()
slo_tx_ci <- data.frame()

# for (i in fnct) {
#   ind <- which(fst_ct_l[,1]==i)
#   qt_l_l <- quantile(fst_ct_l[ind,2], .025)
#   qt_l_u <- quantile(fst_ct_l[ind,2], .975)
#   qt_u_l <- quantile(fst_ct_u[ind,2], .025)
#   qt_u_u <- quantile(fst_ct_u[ind,2], .975)
#   fst_ct_ci <- rbind(fst_ct_ci, data.frame(qt_l_l,qt_l_u, qt_u_l,qt_u_u,i))
# }

for (i in fntt) {
  ind <- which(fst_tx_l[,1]==i)
  qt_l_l <- quantile(fst_tx_l[ind,2], .025)
  qt_l_u <- quantile(fst_tx_l[ind,2], .975)
  qt_u_l <- quantile(fst_tx_u[ind,2], .025)
  qt_u_u <- quantile(fst_tx_u[ind,2], .975)
  fst_tx_ci <- rbind(fst_tx_ci, data.frame(qt_l_l,qt_l_u, qt_u_l,qt_u_u,i))
}

# for (i in snct) {
#   ind <- which(slo_ct_l[,1]==i)
#   qt_l_l <- quantile(slo_ct_l[ind,2], .025)
#   qt_l_u <- quantile(slo_ct_l[ind,2], .975)
#   qt_u_l <- quantile(slo_ct_u[ind,2], .025)
#   qt_u_u <- quantile(slo_ct_u[ind,2], .975)
#   slo_ct_ci <- rbind(slo_ct_ci, data.frame(qt_l_l,qt_l_u, qt_u_l,qt_u_u,i))
# }
# 
for (i in sntt) {
  ind <- which(slo_tx_l[,1]==i)
  qt_l_l <- quantile(slo_tx_l[ind,2], .025)
  qt_l_u <- quantile(slo_tx_l[ind,2], .975)
  qt_u_l <- quantile(slo_tx_u[ind,2], .025)
  qt_u_u <- quantile(slo_tx_u[ind,2], .975)
  slo_tx_ci <- rbind(slo_tx_ci, data.frame(qt_l_l,qt_l_u, qt_u_l,qt_u_u,i))
}

plot(sort(slo_tx_ci$i), sort(slo_tx_ci$qt_l_l, decreasing=T),type='l', col='red')
lines(sort(slo_tx_ci$i), sort(slo_tx_ci$qt_l_u, decreasing=T), col='red')
lines(sort(fst_tx_ci$i), sort(fst_tx_ci$qt_l_l, decreasing=T),type='l', col='blue')
lines(sort(fst_tx_ci$i), sort(fst_tx_ci$qt_l_u, decreasing=T), col='blue')
# lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_l_l, decreasing=T),type='l', col='purple')
# lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_l_u, decreasing=T), col='purple')
# lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_l_l, decreasing=T),type='l', col='forestgreen')
# lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_l_u, decreasing=T), col='forestgreen')

lines(sort(fst_tx_ci$i), sort(fst_tx_ci$qt_u_l, decreasing=T), col='blue')
lines(sort(fst_tx_ci$i), sort(fst_tx_ci$qt_u_u, decreasing=T), col='blue')
lines(sort(slo_tx_ci$i), sort(slo_tx_ci$qt_u_l, decreasing=T),lty='dashed', col='red')
lines(sort(slo_tx_ci$i), sort(slo_tx_ci$qt_u_u, decreasing=T),lty='dashed', col='red')
# lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_u_l, decreasing=T),lty='dashed', col='purple')
# lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_u_u, decreasing=T),lty='dashed', col='purple')
# lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_u_l, decreasing=T), col='forestgreen')
# lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_u_u, decreasing=T), col='forestgreen')
