#insert packages here

n.arms <- 3
null.rate <-.693/12
hr <- rep(1, n.arms)
beta <- .9
alpha <- rep(1-(.95)^(1/n.arms), n.arms)
A <- 1
rate <- 700
tx.rate <- hr*null.rate


test <- mean(rexp(10000, null.rate))


sum(500*5/6,
500/4,
100,
500/3)
sum(500*5/6,
    500/4,
    100)

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


for (i in 1:10000) {
  ctrl.ev <- rexp(ctrl.n, rate=null.rate)
  tx1.ev <- rexp(tx1.n, rate=tx.rate[1])
  tx2.ev <- rexp(tx2.n, rate=tx.rate[2])
  tx3.ev <- rexp(tx3.n, rate=tx.rate[3])
  
  ctrl.entry <- seq_maker(ctrl.patients,1,length(ctrl.patients[1:11]))
  tx1.entry <- seq_maker(ctrl.patients,st[1],end[1])
  tx3.entry <- seq_maker(ctrl.patients,st[6],end[6])
  
  ctrl.time <- ctrl.ev + ctrl.entry
  tx1.time <- tx1.ev + tx1.entry
  tx3.time <- tx3.ev + tx3.entry-st[6]+1
  
  comp1 <- rbind(cbind(ctrl.time[1:length(tx1.ev)], ctrl.ev[1:length(tx1.ev)],
                       rep(0, length(tx1.ev))), cbind(tx1.time, tx1.ev, rep(1,length(tx1.ev))))
  comp1 <- comp1[order(comp1[,1]),]
  comp1 <- cbind(comp1, c(rep(1,evs[1]),rep(0,nrow(comp1)-evs[1])))
  test <- ifelse(comp1[(evs[1]+1):nrow(comp1),2]>comp1[evs[1]], comp1[evs[1]], 
                 comp1[(evs[1]+1):nrow(comp1),2])
  comp1[,1] <- replace(comp1[,1],(evs[1]+1):nrow(comp1),test)
  
  comp3 <- rbind(cbind(ctrl.time[(length(ctrl.entry)-length(tx3.ev)+1):(length(ctrl.entry))]-st[6]+1, 
                       ctrl.ev[(length(ctrl.entry)-length(tx3.ev)+1):(length(ctrl.entry))],
                       rep(0, length(tx3.ev))), cbind(tx3.time, tx3.ev, rep(1,length(tx3.ev))))
  comp3 <- comp3[order(comp3[,1]),]
  comp3 <- cbind(comp3, c(rep(1,evs[3]),rep(0,nrow(comp3)-evs[3])))
  test <- ifelse(comp3[(evs[3]+1):nrow(comp3),2]>comp3[evs[3]], comp3[evs[3]], 
                 comp3[(evs[3]+1):nrow(comp3),2])
  comp3[,1] <- replace(comp3[,1],(evs[3]+1):nrow(comp3),test)
  
  
  obj1 <- survfit(Surv(comp1[comp1[,3]==1,1], comp1[comp1[,3]==1,4])~1)
  
  obj2 <- survfit(Surv(comp1[comp1[,3]==0,1], comp1[comp1[,3]==0,4])~1)

  
  obj3 <- survfit(Surv(comp3[comp3[,3]==1,1], comp3[comp3[,3]==1,4])~1)
  
  obj4 <- survfit(Surv(comp3[comp3[,3]==0,1], comp3[comp3[,3]==0,4])~1)
  
  slow_null_ct_t[[i]] <- round(sort(obj2$time), digits=1)
  slow_null_tx_t[[i]] <- round(sort(obj1$time), digits=1)
  fast_null_ct_t[[i]] <- round(sort(obj4$time), digits=1)
  fast_null_tx_t[[i]] <- round(sort(obj3$time), digits=1)
  
  CI_slow_null_ct_l[[i]] <- sort(obj2$lower, decreasing=T)
  CI_slow_null_ct_u[[i]] <- sort(obj2$upper, decreasing=T)
  CI_slow_null_tx_l[[i]] <- sort(obj1$lower, decreasing=T)
  CI_slow_null_tx_u[[i]] <- sort(obj1$upper, decreasing=T)
  CI_fast_null_ct_l[[i]] <- sort(obj4$lower, decreasing=T)
  CI_fast_null_ct_u[[i]] <- sort(obj4$upper, decreasing=T)
  CI_fast_null_tx_l[[i]] <- sort(obj3$lower, decreasing=T)
  CI_fast_null_tx_u[[i]] <- sort(obj3$upper, decreasing=T)
}

slo_ct_l <- cbind(unlist(slow_null_ct_t),unlist(CI_slow_null_ct_l))
slo_ct_u <- cbind(unlist(slow_null_ct_t),unlist(CI_slow_null_ct_u))
fst_ct_l <- cbind(unlist(fast_null_ct_t),unlist(CI_fast_null_ct_l))
fst_ct_u <- cbind(unlist(fast_null_ct_t),unlist(CI_fast_null_ct_u))
slo_tx_l <- cbind(unlist(slow_null_tx_t),unlist(CI_slow_null_tx_l))
slo_tx_u <- cbind(unlist(slow_null_tx_t),unlist(CI_slow_null_tx_u))
fst_tx_l <- cbind(unlist(fast_null_tx_t),unlist(CI_fast_null_tx_l))
fst_tx_u <- cbind(unlist(fast_null_tx_t),unlist(CI_fast_null_tx_u))

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



plot(sort(fst_tx_ci$i), sort(fst_tx_ci$qt_u_l, decreasing=T),type='l', col='blue')
lines(sort(fst_tx_ci$i), sort(fst_tx_ci$qt_u_u, decreasing=T), col='blue')
lines(sort(slo_tx_ci$i), sort(slo_tx_ci$qt_u_l, decreasing=T),lty='dashed', col='red')
lines(sort(slo_tx_ci$i), sort(slo_tx_ci$qt_u_u, decreasing=T),lty='dashed', col='red')
lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_u_l, decreasing=T),lty='dashed', col='purple')
lines(sort(slo_ct_ci$i), sort(slo_ct_ci$qt_u_u, decreasing=T),lty='dashed', col='purple')
lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_u_l, decreasing=T), col='forestgreen')
lines(sort(fst_ct_ci$i), sort(fst_ct_ci$qt_u_u, decreasing=T), col='forestgreen')

