library(survival)
library(survminer)

nsim <- 1000
null.hz <- .3
tx.hz <- 1*null.hz
nevents <- 524
boy <- ceiling(nevents*.35)
boy2 <- ceiling(nevents*.7)

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

set.seed(2349)

for (i in 1:nsim) {
  entry <- seq(1/500,6, by=1/500)
  ctrl.entry <- entry[seq(1,3000, by=5)]
  tx1.entry <- entry[seq(2,3000,by=5)]
  ctrl.1 <- rexp(length(ctrl.entry), null.hz)
  tx.1 <- rexp(length(tx1.entry), tx.hz)
  
  comp1 <- rbind(cbind(ctrl.1+ctrl.entry,ctrl.1, ctrl.entry, rep(0,length(ctrl.entry))),
                 cbind(tx.1+tx1.entry,tx.1, tx1.entry, rep(1,length(tx1.entry))))
  comp1 <- comp1[order(comp1[,1]),]
  cut <- comp1[boy,1]
  inds <- which(comp1[,3]>cut)
  comp1 <- comp1[-inds,]
  
  entry2 <- seq(cut,6, by=1/500)
  ctrl.entry <- entry2[seq(1,length(entry2), by=6)]
  tx1.entry <- entry2[seq(2,length(entry2),by=6)]
  tx5.entry <- entry2[seq(6, length(entry2), by=6)]
  ctrl.5.st <- ctrl.entry[1]
  tx5.st <- tx5.entry[1]
  
  ctrl.1 <- rexp(length(ctrl.entry), null.hz)
  tx.1 <- rexp(length(tx1.entry), tx.hz)
  tx.5 <- rexp(length(tx5.entry), tx.hz)
  
  comp1 <- rbind(comp1, 
                 rbind(cbind(ctrl.1+ctrl.entry,ctrl.1, ctrl.entry, rep(0,length(ctrl.entry))),
                       cbind(tx.1+tx1.entry,tx.1, tx1.entry, rep(1,length(tx1.entry)))))
  comp1 <- comp1[order(comp1[,1]),]
  cut <- comp1[boy2,1]
  inds <- which(comp1[,3]>cut)
  comp1 <- comp1[-inds,]
  comp1 <- cbind(comp1, c(rep(1,nevents),rep(0,nrow(comp1)-nevents)))
  test1 <- ifelse(comp1[1:nrow(comp1),1]>comp1[nevents], comp1[nevents]-comp1[1:nrow(comp1),3], 
                  comp1[1:nrow(comp1),2])
  comp1[,2] <- replace(comp1[,2],1:nrow(comp1),test1)
  
  comp5 <- rbind(cbind(ctrl.1+ctrl.entry,ctrl.1, ctrl.entry, rep(0,length(ctrl.entry))),
                 cbind(tx.5+tx5.entry,tx.5, tx5.entry, rep(1,length(tx5.entry))))
  comp5 <- comp5[order(comp5[,1]),]
  inds <- which(comp5[,3]>cut)
  comp5 <- comp5[-inds,]
  
  entry3 <- seq(cut,8, by=1/500)
  ctrl.entry <- entry3[seq(1,length(entry3), by=2)]
  tx5.entry <- entry3[seq(2, length(entry3), by=2)]
  
  ctrl.1 <- rexp(length(ctrl.entry), null.hz)
  tx.5 <- rexp(length(tx5.entry), tx.hz)
  
  comp5 <- rbind(comp5, 
                 rbind(cbind(ctrl.1+ctrl.entry,ctrl.1, ctrl.entry, rep(0,length(ctrl.entry))),
                       cbind(tx.5+tx5.entry,tx.5, tx5.entry, rep(1,length(tx5.entry)))))
  comp5 <- comp5[order(comp5[,1]),]
  cut <- comp5[boy2,1]
  inds <- which(comp5[,3]>cut)
  comp5 <- comp5[-inds,]
  comp5 <- cbind(comp5, c(rep(1,nevents),rep(0,nrow(comp5)-nevents)))
  comp5[comp5[,4]==0,1] <- comp5[comp5[,4]==0,1]-ctrl.5.st
  comp5[comp5[,4]==0,3] <- comp5[comp5[,4]==0,3]-ctrl.5.st
  comp5[comp5[,4]==1,1] <- comp5[comp5[,4]==1,1]-tx5.st
  comp5[comp5[,4]==1,3] <- comp5[comp5[,4]==1,3]-tx5.st
  test5 <- ifelse(comp5[1:nrow(comp5),1]>comp5[nevents], comp5[nevents]-comp5[1:nrow(comp5),3], 
                  comp5[1:nrow(comp5),2])
  comp5[,2] <- replace(comp5[,2],1:nrow(comp5),test5)
  
  t_length[i,] <- cbind(comp1[nevents], comp5[nevents])
  n_pat[i, ] <- cbind(nrow(comp1), nrow(comp5))
  
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

hist(p.expt_1c)
hist(p.expt_5c)
hist(p.expt_15)

t.test(t_length$arm1, t_length$arm2)
t.test(n_pat$arm1, n_pat$arm2)

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

