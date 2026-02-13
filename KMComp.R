library(mnormt)
library(purrr)
library(survival)
library(gsDesign)
library(plyr)
library(tidyverse)
library(survminer)

n.arms <- 3
null.rate <-.693/12
hr <- rep(.75, n.arms)
beta <- .9
alpha <- rep(1-(.95)^(1/n.arms), n.arms)
A <- 1
rate <- 500
tx.rate <- hr*null.rate
st <- c(1, 1, 1,1, 1,6)
end <- c(9, 9,9,9,9, 11)


ctrl.n <- sum(ctrl.patients[1:11])
tx1.n <- tx2.n <- sum(ctrl.patients[st[1]:end[1]])
tx3.n <- sum(ctrl.patients[st[6]:end[6]])

ctrl.ev <- rexp(ctrl.n, rate=null.rate)
tx1.ev <- rexp(tx1.n, rate=tx.rate[1])
tx2.ev <- rexp(tx2.n, rate=tx.rate[2])
tx3.ev <- rexp(tx3.n, rate=tx.rate[3])

seq_maker <- function(entry.times, start, end){
  seql <- c()
  for (i in start:end) {
    seql <- append(seql, (seq(i-1,i-1/entry.times[i], by=1/entry.times[i])))
  }
  return(seql)
}

ctrl.entry <- seq_maker(ctrl.patients,1,length(ctrl.patients[1:11]))
tx1.entry <- tx2.entry <- seq_maker(ctrl.patients,st[1],end[1])
tx3.entry <- seq_maker(ctrl.patients,st[6],end[6])

ctrl.time <- ctrl.ev + ctrl.entry
tx1.time <- tx1.ev + tx1.entry
#tx2.time <- tx2.ev + tx2.entry
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

ggsurvplot(obj1, data = as.data.frame(cbind(comp1[comp1[,3]==1,1], comp1[comp1[,3]==1,4])),
           conf.int = T,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

obj2 <- survfit(Surv(comp1[comp1[,3]==0,1], comp1[comp1[,3]==0,4])~1)

ggsurvplot(obj2, data = as.data.frame(cbind(comp1[comp1[,3]==0,1], comp1[comp1[,3]==0,4])),
           conf.int = T,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

obj3 <- survfit(Surv(comp3[comp3[,3]==1,1], comp3[comp3[,3]==1,4])~1)

ggsurvplot(obj3, data = as.data.frame(cbind(comp3[comp3[,3]==1,1], comp3[comp3[,3]==1,4])),
           conf.int = T,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

obj4 <- survfit(Surv(comp3[comp3[,3]==0,1], comp3[comp3[,3]==0,4])~1)

ggsurvplot(obj4, data = as.data.frame(cbind(comp3[comp3[,3]==0,1], comp3[comp3[,3]==0,4])),
           conf.int = T,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())
