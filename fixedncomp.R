library(survival)
library(survminer)

nsim <- 500
null.hz <- .3
tx.hz <- 1*null.hz
nevents <- 524
halfn <- nevents*3/4


entry <- seq(0,6, by=1/500)
ctrl.entry <- entry[c(seq(1,1425, by=5),seq(1426,2072, by=6),seq(2073,2644, by=2))]
tx1.entry <- entry[c(seq(2,1425, by=5),seq(1427,2072, by=6))]
tx5.entry <- entry[c(seq(1431,2072, by=6), seq(2074,2644, by=2))]

ctrl.l <- length(ctrl.entry)
five_st <- entry[1431]
ctrl.5.st <- entry[1426]

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

set.seed(2349)

for (i in 1:nsim) {
  ctrl.evs <- rexp(length(ctrl.entry), null.hz)
  tx1.evs <- rexp(length(tx1.entry), tx.hz)
  tx5.evs <- rexp(length(tx5.entry), tx.hz)
  
  ctrl.time <- ctrl.entry+ctrl.evs
  tx1.time <- tx1.entry+tx1.evs
  tx5.time <- tx5.entry+tx5.evs
  
  comp1 <- rbind(cbind(ctrl.time[1:halfn], ctrl.evs[1:halfn],ctrl.entry[1:halfn],
                       rep(0, halfn)), cbind(tx1.time, tx1.evs, tx1.entry, rep(1,halfn)))
  comp1 <- comp1[order(comp1[,1]),]
  comp1 <- cbind(comp1, c(rep(1,nevents),rep(0,nrow(comp1)-nevents)))
  test1 <- ifelse(comp1[1:nrow(comp1),1]>comp1[nevents], comp1[nevents]-comp1[1:nrow(comp1),3], 
                 comp1[1:nrow(comp1),2])
  comp1[,2] <- replace(comp1[,2],1:nrow(comp1),test1)
  
  comp5 <- rbind(cbind(ctrl.time[(ctrl.l-halfn+1):ctrl.l]-ctrl.5.st, 
                       ctrl.evs[(ctrl.l-halfn+1):ctrl.l], ctrl.entry[(ctrl.l-halfn+1):ctrl.l]-ctrl.5.st,
                       rep(0, halfn)), cbind(tx5.time-five_st, tx5.evs, tx5.entry-five_st, rep(1,halfn)))
  comp5 <- comp5[order(comp5[,1]),]
  comp5 <- cbind(comp5, c(rep(1,nevents),rep(0,nrow(comp5)-nevents)))
  test5 <- ifelse(comp5[1:nrow(comp5),1]>comp5[nevents], comp5[nevents]-comp5[1:nrow(comp5),3], 
                 comp5[1:nrow(comp5),2])
  comp5[,2] <- replace(comp5[,2],1:nrow(comp5),test5)
  
  t_length[i,] <- cbind(comp1[nevents], comp5[nevents])
  
  test1c <- summary(coxph(Surv(comp1[,2], comp1[,5])~comp1[,4]))
  test5c <- summary(coxph(Surv(comp5[,2], comp5[,5])~comp5[,4]))
  test15 <- summary(coxph(Surv(c(comp1[comp1[,4]==1,2],comp5[comp5[,4]==1,2]),
                               c(comp1[comp1[,4]==1,5],comp5[comp5[,4]==1,5]))~
                               c(rep(0,length(comp1[comp1[,4]==1,2])),
                                     rep(1,length(comp5[comp5[,4]==1,2])))))
  
  #stat[i] <- test$coef[4]
  p.expt_1c[i] <- test1c$coef[5]
  p.expt_5c[i] <- test5c$coef[5]
  p.expt_15[i] <- test15$coef[5]
  
  obj1 <- survfit(Surv(comp1[comp1[,4]==1,2], comp1[comp1[,4]==1,5])~1)
  
  obj2 <- survfit(Surv(comp1[comp1[,4]==0,2], comp1[comp1[,4]==0,5])~1)
  check <- all.equal(obj2$n,length(obj2$time),length(obj2$lower),length(obj2$upper))
  
  obj3 <- survfit(Surv(comp5[comp5[,4]==1,1], comp5[comp5[,4]==1,5])~1)
  check <- all.equal(obj3$n,length(obj3$time),length(obj3$lower),length(obj3$upper))
  
  obj4 <- survfit(Surv(comp5[comp5[,4]==0,1], comp5[comp5[,4]==0,5])~1)
  
  slow_null_ct_t[[i]] <- round(obj2$time, digits=1)
  slow_null_tx_t[[i]] <- round(obj1$time, digits=1)
  fast_null_ct_t[[i]] <- round(obj4$time, digits=1)
  fast_null_tx_t[[i]] <- round(obj3$time, digits=1)
  
  CI_slow_null_ct_l[[i]] <- obj2$lower
  CI_slow_null_ct_u[[i]] <- obj2$upper
  CI_slow_null_tx_l[[i]] <- obj1$lower
  CI_slow_null_tx_u[[i]] <- obj1$upper
  CI_fast_null_ct_l[[i]] <- obj4$lower
  CI_fast_null_ct_u[[i]] <- obj4$upper
  CI_fast_null_tx_l[[i]] <- obj3$lower
  CI_fast_null_tx_u[[i]] <- obj3$upper
  
}

mean(p.expt_1c<.05)
mean(p.expt_5c<.05)
mean(p.expt_15<.05)

t.test(t_length$arm1, t_length$arm2)

slo_ct_l <- cbind(unlist(slow_null_ct_t),unlist(CI_slow_null_ct_l))
slo_ct_u <- cbind(unlist(slow_null_ct_t),unlist(CI_slow_null_ct_u))
fst_ct_l <- cbind(unlist(fast_null_ct_t),unlist(CI_fast_null_ct_l))
fst_ct_u <- cbind(unlist(fast_null_ct_t),unlist(CI_fast_null_ct_u))
slo_tx_l <- cbind(unlist(slow_null_tx_t),unlist(CI_slow_null_tx_l))
slo_tx_u <- cbind(unlist(slow_null_tx_t),unlist(CI_slow_null_tx_u))
fst_tx_l <- cbind(unlist(fast_null_tx_t),unlist(CI_fast_null_tx_l))
fst_tx_u <- cbind(unlist(fast_null_tx_t),unlist(CI_fast_null_tx_u))

slo_ct_l <- slo_ct_l[complete.cases(slo_ct_l),]
slo_ct_u <- slo_ct_u[complete.cases(slo_ct_u),]
fst_ct_l <- fst_ct_l[complete.cases(fst_ct_l),]
fst_ct_u <- fst_ct_u[complete.cases(fst_ct_u),]
slo_tx_l <- slo_tx_l[complete.cases(slo_tx_l),]
slo_tx_u <- slo_tx_u[complete.cases(slo_tx_u),]
fst_tx_l <- fst_tx_l[complete.cases(fst_tx_l),]
fst_tx_u <- fst_tx_u[complete.cases(fst_tx_u),]

fnct <- unique(fst_ct_l[,1])
fntt <- unique(fst_tx_l[,1])
snct <- unique(slo_ct_l[,1])
sntt <- unique(slo_tx_l[,1])

fst_ct_ci <- data.frame()
fst_tx_ci <- data.frame()
slo_ct_ci <- data.frame()
slo_tx_ci <- data.frame()

for (i in fnct) {
  ind <- which(fst_ct_l[,1]==i)
  qt_l_l <- quantile(fst_ct_l[ind,2], .025)
  qt_l_u <- quantile(fst_ct_l[ind,2], .975)
  qt_u_l <- quantile(fst_ct_u[ind,2], .025)
  qt_u_u <- quantile(fst_ct_u[ind,2], .975)
  fst_ct_ci <- rbind(fst_ct_ci, data.frame(qt_l_l,qt_l_u, qt_u_l,qt_u_u,i))
}

for (i in fntt) {
  ind <- which(fst_tx_l[,1]==i)
  qt_l_l <- quantile(fst_tx_l[ind,2], .025)
  qt_l_u <- quantile(fst_tx_l[ind,2], .975)
  qt_u_l <- quantile(fst_tx_u[ind,2], .025)
  qt_u_u <- quantile(fst_tx_u[ind,2], .975)
  fst_tx_ci <- rbind(fst_tx_ci, data.frame(qt_l_l,qt_l_u, qt_u_l,qt_u_u,i))
}

for (i in snct) {
  ind <- which(slo_ct_l[,1]==i)
  qt_l_l <- quantile(slo_ct_l[ind,2], .025)
  qt_l_u <- quantile(slo_ct_l[ind,2], .975)
  qt_u_l <- quantile(slo_ct_u[ind,2], .025)
  qt_u_u <- quantile(slo_ct_u[ind,2], .975)
  slo_ct_ci <- rbind(slo_ct_ci, data.frame(qt_l_l,qt_l_u, qt_u_l,qt_u_u,i))
}

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
lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_l_l, decreasing=T),type='l', col='purple')
lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_l_u, decreasing=T), col='purple')
lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_l_l, decreasing=T),type='l', col='forestgreen')
lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_l_u, decreasing=T), col='forestgreen')

lines(sort(fst_tx_ci$i), sort(fst_tx_ci$qt_u_l, decreasing=T), col='blue')
lines(sort(fst_tx_ci$i), sort(fst_tx_ci$qt_u_u, decreasing=T), col='blue')
lines(sort(slo_tx_ci$i), sort(slo_tx_ci$qt_u_l, decreasing=T),lty='dashed', col='red')
lines(sort(slo_tx_ci$i), sort(slo_tx_ci$qt_u_u, decreasing=T),lty='dashed', col='red')
lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_u_l, decreasing=T),lty='dashed', col='purple')
lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_u_u, decreasing=T),lty='dashed', col='purple')
lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_u_l, decreasing=T), col='forestgreen')
lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_u_u, decreasing=T), col='forestgreen')

