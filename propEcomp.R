library(survival)
library(survminer)
library(janitor)
library(beepr)

nsim <- 500
null.hz <- .3
tx.hz <- 1*null.hz
nevents <- 524
start2 <- ceiling(nevents*.35)
ntot <- ceiling(nevents*.7)
arm_no <- c(2, 6)
n.arms <- c(5,6,2)
#number of arms at any given time

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
n_pat <- data.frame(arm1=0, arm2=0)

surv1T <- list()
surv5T <- list()

diff1 <- list()
diff5 <- list()
time1 <- list()
time5 <- list()

mdiff1 <- rep(NA, nsim)
mdiff5 <- rep(NA, nsim)

for (i in 1:nsim) {
  entry <- seq(1/500,6, by=1/500)
  ctrl1entry.1 <- entry[seq(1,3000, by=n.arms[1])]
  tx1entry.1 <- entry[seq(arm_no[1],3000,by=n.arms[1])]
  ctrl1evs.1 <- rexp(length(ctrl1entry.1), null.hz)
  tx1evs.1 <- rexp(length(tx1entry.1), tx.hz)
  
  ctrl1time.1 <- ctrl1entry.1+ctrl1evs.1
  tx1time.1 <- tx1entry.1 + tx1evs.1
  
  comp1 <- c(ctrl1time.1, tx1time.1)
  censTime1.1 <- comp1[order(comp1)][start2]
  
  inds1 <- which(ctrl1entry.1<censTime1.1)
  inds2 <- which(tx1entry.1<censTime1.1)
  
  ctrl1evs.1 <- ctrl1evs.1[inds1]
  tx1evs.1 <- tx1evs.1[inds2]
  ctrl1time.1 <- ctrl1time.1[inds1]
  tx1time.1 <- tx1time.1[inds2]
  ctrl1entry.1 <- ctrl1entry.1[inds1]
  tx1entry.1 <- tx1entry.1[inds2]
  
  nu_st <- round_to_fraction(censTime1.1, 500, 3)
  
  entry2 <- seq(nu_st+1/500,6+nu_st, by=1/500)
  ctrl1entry.2 <- entry2[seq(1,length(entry2), by=n.arms[2])]
  tx1entry.2 <- entry2[seq(arm_no[1],length(entry2),by=n.arms[2])]
  ctrl5entry.1 <- entry[seq(1,length(entry), by=n.arms[2])]
  tx5entry.1 <- entry[seq(arm_no[2], length(entry), by=n.arms[2])]
  #ctrl.5.st <- ctrl.sh.entry[1]
  #tx5.st <- tx5entry.1[1]
  
  ctrl.sh.evs <- rexp(length(ctrl1entry.2), null.hz)
  tx1evs.2 <- rexp(length(tx1entry.2), tx.hz)
  tx5evs.1 <- rexp(length(tx5entry.1), tx.hz)
  
  ctrl1time.2 <- ctrl1entry.2+ctrl.sh.evs
  tx1time.2 <- tx1entry.2 + tx1evs.2
  ctrl5time.1 <- ctrl5entry.1 + ctrl.sh.evs
  tx5time.1 <- tx5entry.1 + tx5evs.1
  
  comp1 <- c(ctrl1time.1, tx1time.1, ctrl1time.2, tx1time.2)
  censTime1.2 <- comp1[order(comp1)][ntot]
  
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
  
  c1stat <- ctrl1.time < censTime1
  t1stat <- tx1.time < censTime1
  
  nu_st2 <- round_to_fraction((censTime1.2-nu_st), 500, 3)
  
  entry3 <- seq(nu_st2+1/500,8, by=1/500)
  ctrl5entry.2 <- entry3[seq(1,length(entry3), by=n.arms[3])]
  tx5entry.2 <- entry3[seq(2, length(entry3), by=n.arms[3])]
  
  ctrl5evs.2 <- rexp(length(ctrl5entry.2), null.hz)
  tx5evs.2 <- rexp(length(tx5entry.2), tx.hz)
  
  ctrl5time.2 <- ctrl5entry.2+ctrl5evs.2
  tx5time.2 <- tx5entry.2 + tx5evs.2
  
  comp5 <- c(ctrl5time.1, tx5time.1, ctrl5time.2, tx5time.2)
  censTime2.2 <- comp5[order(comp5)][ntot]
  
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
  ctrl5.time <- c(ctrl5time.1, ctrl5time.2)#-ctrl.5.st
  tx5.time <- c(tx5time.1, tx5time.2)#-tx5.st
  ctrl5.entry <- c(ctrl5entry.1, ctrl5entry.2)#-nu_st
  tx5.entry <- c(tx5entry.1, tx5entry.2)#-nu_st
  
  comp5 <- c(ctrl5.time, tx5.time)
  censTime5 <- comp5[order(comp5)][nevents]
  
  c5stat <- ctrl5.time < censTime5
  t5stat <- tx5.time < censTime5
  
  evTimesC1 <- ctrl1.evs * c1stat + (censTime1 - ctrl1.entry) * (1 - c1stat)
  evTimesT1 <- tx1.evs * t1stat + (censTime1 - tx1.entry) * (1 - t1stat)
  evTimesC5 <- ctrl5.evs * c5stat + (censTime5 - ctrl5.entry) * (1 - c5stat)
  evTimesT5 <- tx5.evs * t5stat + (censTime5 - tx5.entry) * (1 - t5stat)
  
  t_length[i,] <- cbind(censTime1, censTime5)
  n_pat[i, ] <- cbind(length(comp1), length(comp5))
  
  test1c <- summary(coxph(Surv(time = c(evTimesC1,evTimesT1) , event = c(c1stat, t1stat))
                          ~c(rep(0, length(evTimesC1)), rep(1, length(evTimesT1)))))
  
  test5c <- summary(coxph(Surv(time = c(evTimesC5,evTimesT5) , event = c(c5stat, t5stat))~
                            c(rep(0, length(evTimesC5)), rep(1, length(evTimesT5)))))
  
  test15 <- summary(coxph(Surv(time = c(evTimesT1,evTimesT5) , event = c(t1stat, t5stat))~
                            c(rep(0, length(evTimesT1)), rep(1, length(evTimesT5)))))
  
  #stat[i] <- test$coef[4]
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
  
};beep()

mean(p.expt_1c<.05)
mean(p.expt_5c<.05)
mean(p.expt_15<.05)

hist(p.expt_1c)
hist(p.expt_5c)
hist(p.expt_15)

t.test(t_length$arm1, t_length$arm2, paired = T)
t.test(n_pat$arm1, n_pat$arm2, paired=T)
t.test(mdiff1, mdiff5, paired = T)

ks.test(diff1[[4]], diff5[[4]])

plot(ecdf(diff1[[18]]))
lines(ecdf(diff5[[18]]), col=2)

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
