process <- cumsum(rexp(2328, 500))

goodup <- function(x, level=1) round(x + 5*10^(-level-1), level)

A <- 2
t.length <- 2.33
nsim <- 15000
arm2.start <- seq(0,2,by=.2)
null.rate <-.693
hz1 <- hz2 <- .75

control.rate <- ceiling(500/(1+A))
expt.rate <- floor(500*A/(1+A))
n.control <- ceiling(control.rate*(t.length+2))
n.expt <- ceiling(expt.rate*t.length)
t1.rate <- null.rate*hz1
t2.rate <- null.rate*hz2

tlength <- rep(NA, 11)
stat1 <- list()
p.expt1 <- list()
stat2 <- list()
p.expt2 <- list()


for (i in 1:length(arm2.start)) {
  stat1[[i]] <- rep(NA, nsim)
  p.expt1[[i]] <- rep(NA, nsim)
  stat2[[i]] <- rep(NA, nsim)
  p.expt2[[i]] <- rep(NA, nsim)
  
  ind1 <- which(process>arm2.start[i])[[1]]
  
  if (ind1==1) {
    ntot <- n.expt*2 + round(n.expt/A)
    tmp <- process[1:ntot]
    
    tlength[i] <- goodup(max(tmp),1)+.2
    
    ind2 <- sample(1:ntot, round(n.expt/A))
    ctrl.entry <- tmp[ind2]
    
    tmp <- tmp[-ind2]
    ind2 <- sample(1:length(tmp), n.expt)
    exp1.entry <- tmp[ind2]
    exp2.entry <- tmp[-ind2]
    
    n.st <-1
    n.fin <- length(ctrl.entry)
    
  }else {
    tmp <- process[1:ind1]
    ind2 <- sample(1:ind1, floor(length(tmp)*A/(1+A)))
    
    ctrl.entry <- tmp[-ind2]
    exp1.entry <- tmp[ind2]
  
  n.st <- length(ctrl.entry)+1
  
  ind2 <- ind1+(n.expt-length(exp1.entry))*2 + round((n.expt-length(exp1.entry))/A)
  ind3 <- ind1+1
  
  tmp <- process[ind3:ind2]
  
  tlength[i] <- goodup(max(tmp),1)+.2
  
  tmp2 <- ind2-ind3+1
  ind3 <- sample(1:tmp2, round((n.expt-length(exp1.entry))/A))
  
  ctrl.entry <- c(ctrl.entry, tmp[ind3])
  n.fin <- length(ctrl.entry)
  
  tmp <- tmp[-ind3]
  ind4 <- sample(1:length(tmp), (n.expt-length(exp1.entry)))
  
  exp1.entry <- c(exp1.entry, tmp[ind4])
  exp2.entry <- tmp[-ind4]
  
  ind3 <- ind2+(n.expt-length(exp2.entry))+ n.st-1
  ind2 <- ind2+1
  tmp <- process[ind2:ind3]
  tmp2 <- ind3-ind2+1
  ind3 <- sample(1:tmp2, n.st-1)
  
  ctrl.entry <- c(ctrl.entry, tmp[ind3])
  exp2.entry <- c(exp2.entry, tmp[-ind3])
  }
  
  ctrl.entry <- sort(ctrl.entry)
  exp1.entry <- sort(exp1.entry)
  exp2.entry <- sort(exp2.entry)
  
  ctrl.list <- list()
  tx1.list <- list()
  tx2.list <- list()
  ctrl.time <- list()
  tx1.time <- list()
  tx2.time <- list()
  
  ctr.v <- c(rep(0, n.fin), rep(1, n.expt))
  
  for (j in 1:nsim) {
    ctrl.time[[j]] <- rexp(length(ctrl.entry), rate=null.rate)
    tx1.time[[j]] <- rexp(n.expt, rate = t1.rate)
    tx2.time[[j]] <- rexp(n.expt, rate = t2.rate)
    
    ctrl.list[[j]] <- ctrl.time[[j]]+ctrl.entry
    tx1.list[[j]] <- tx1.time[[j]] +exp1.entry
    tx2.list[[j]] <- tx2.time[[j]] +exp2.entry
    
    censc <- ifelse(ctrl.list[[j]][1:n.fin]>tlength[i], 0, 1)
    censt <- ifelse(tx1.list[[j]]>tlength[i],0,1)
    ttec <- ifelse(censc, ctrl.time[[j]][1:n.fin],tlength[i]-ctrl.entry[1:n.fin])
    ttet <- ifelse(censt, tx1.time[[j]], tlength[i]-exp1.entry)
    test <- summary(coxph(Surv(c(ttec, ttet), c(censc, censt))~ctr.v))
    stat1[[i]][j] <- test$coef[4]
    p.expt1[[i]][j] <- test$coef[5]
    
    censc <- ifelse(ctrl.list[[j]][n.st:length(ctrl.entry)]>tlength[i]+arm2.start[i], 0, 1)
    censt <- ifelse(tx2.list[[j]]>tlength[i]+arm2.start[i],0,1)
    ttec <- ifelse(censc, ctrl.time[[j]][n.st:length(ctrl.entry)],tlength[i]+arm2.start[i]-ctrl.entry[n.st:length(ctrl.entry)])
    ttet <- ifelse(censt, tx2.time[[j]], tlength[i]+arm2.start[i]-exp2.entry)
    test <- summary(coxph(Surv(c(ttec, ttet), c(censc, censt))~ctr.v))
    stat2[[i]][j] <- test$coef[4]
    p.expt2[[i]][j] <- test$coef[5]
    
  }
  
  cat("Simulation", i*100/11,"%", "donezo!", "\n")
  
};beep()

for (i in 1:11) {
  print(cor(stat1[[i]], stat2[[i]]))
}

for (i in 1:11) {
  print(mean(p.expt1[[i]]<.025|p.expt2[[i]]<.025))
}
