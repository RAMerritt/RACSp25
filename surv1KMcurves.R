library(plyr)
library(survival)
library(tidyverse)
library(survminer)
per_time <- 500
pers <- 3
null.rate <-.693
hz <- .75
t.rate <- null.rate*hz
unit_adj <- 1/365
breaks <- seq(per_time, per_time*3, per_time)
t1.time <- 3
t2.time <- 5


unif.entry.tot <- round_any(sort(c(runif(per_time), 1+runif(per_time), 
                                   2+runif(per_time))), unit_adj)
ctrl.inds <- c(seq(1, breaks[1], 3), seq(breaks[1]+1, breaks[2], 5), 
               seq(breaks[2]+1, breaks[3], 3))
t1.inds <- c(seq(2,breaks[1]+2,3), seq(breaks[1]+2, breaks[2], 5))
t2.inds <- c(seq(3,breaks[1]+3,3), seq(breaks[1]+3, breaks[2], 5))
t3.inds <- c(seq(breaks[1]+4, breaks[2], 5), seq(breaks[2]+2, breaks[3], 3))
t4.inds <- c(seq(breaks[1]+5, breaks[2], 5), seq(breaks[2]+3, breaks[3], 3))

unif.ctrl.entry <- unif.entry.tot[ctrl.inds]
unif.entry.t1 <- unif.entry.tot[t1.inds]
unif.entry.t2 <- unif.entry.tot[t2.inds]
unif.entry.t3 <- unif.entry.tot[t3.inds]
unif.entry.t4 <- unif.entry.tot[t4.inds]

n.ctrl <- length(ctrl.inds)
n.t1 <- length(t1.inds)
n.t2 <- length(t2.inds)
n.t3 <- length(t3.inds)
n.t4 <- length(t4.inds)

ctrl.events <- round_any(rexp(n.ctrl, null.rate), unit_adj)
t1.events <- round_any(rexp(n.t1, t.rate),unit_adj)
t2.events <- round_any(rexp(n.t2, t.rate),unit_adj)
t3.events <- round_any(rexp(n.t3, t.rate),unit_adj)
t4.events <- round_any(rexp(n.t4, t.rate),unit_adj)

ctrl.tot <- ctrl.events+unif.ctrl.entry
t1.tot <- t1.events+unif.entry.t1
t2.tot <- t2.events+unif.entry.t2
t3.tot <- t3.events+unif.entry.t3
t4.tot <- t4.events+unif.entry.t4

ctrl.block <- length(c(seq(1, breaks[1], 3), seq(breaks[1]+1, breaks[2], 5)))

ctrl.cens1 <- ifelse(ctrl.tot[1:ctrl.block]>t1.time, 0, 1)
ctrl.obs1 <- ifelse(ctrl.cens1, ctrl.events[1:ctrl.block], 
                    t1.time-unif.ctrl.entry[1:ctrl.block])
t1.cens <- ifelse(t1.tot > t1.time, 0, 1)
t1.obs <- ifelse(t1.cens, t1.events, t1.time-unif.entry.t1)
t2.cens <- ifelse(t2.tot > t1.time, 0, 1)
t2.obs <- ifelse(t2.cens, t2.events, t1.time-unif.entry.t2)

ctrl.cens2 <- ifelse(ctrl.tot[(length(ctrl.tot)-ctrl.block+1):(length(ctrl.tot))]>t2.time, 0, 1)
ctrl.obs2 <- ifelse(ctrl.cens2, ctrl.events[(length(ctrl.tot)-ctrl.block+1):(length(ctrl.tot))],
                    t2.time-unif.ctrl.entry[(length(ctrl.tot)-ctrl.block+1):(length(ctrl.tot))])
t3.cens <- ifelse(t3.tot > t2.time, 0, 1)
t3.obs <- ifelse(t3.cens, t3.events, t2.time-unif.entry.t3)
t4.cens <- ifelse(t4.tot > t2.time, 0, 1)
t4.obs <- ifelse(t4.cens, t4.events, t2.time-unif.entry.t4)

obj1 <- survfit(Surv(t1.obs, t1.cens)~1)
summary(obj1)

ggsurvplot(obj1, data = as.data.frame(cbind(t1.obs, t1.cens)),
           conf.int = FALSE,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

obj2 <- survfit(Surv(t2.obs, t2.cens)~1)
summary(obj2)

ggsurvplot(obj2, data = as.data.frame(cbind(t2.obs, t2.cens)),
           conf.int = FALSE,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

obj3 <- survfit(Surv(t3.obs, t3.cens)~1)
summary(obj3)

ggsurvplot(obj3, data = as.data.frame(cbind(t3.obs, t3.cens)),
           conf.int = FALSE,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

obj4 <- survfit(Surv(t4.obs, t4.cens)~1)
summary(obj4)

ggsurvplot(obj4, data = as.data.frame(cbind(t4.obs, t4.cens)),
           conf.int = FALSE,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

pois.num <- rpois(pers, per_time)
breaks <- cumsum(pois.num)
exp.entry.tot <- round_any(cumsum(unlist(map(.x=pois.num, ~rexp(n=., rate=.)))), unit_adj)

exp.ctrl.entry <- exp.entry.tot[ctrl.inds]
exp.entry.t1 <- exp.entry.tot[t1.inds]
exp.entry.t2 <- exp.entry.tot[t2.inds]
exp.entry.t3 <- exp.entry.tot[t3.inds]
#exp.entry.t4 <- exp.entry.tot[t4.inds]

n.ctrl <- length(ctrl.inds)
n.t1 <- length(t1.inds)
n.t2 <- length(t2.inds)
n.t3 <- length(t3.inds)
#n.t4 <- length(t4.inds)

ctrl.events <- round_any(rexp(n.ctrl, null.rate), unit_adj)
t1.events <- round_any(rexp(n.t1, t.rate),unit_adj)
t2.events <- round_any(rexp(n.t2, t.rate),unit_adj)
t3.events <- round_any(rexp(n.t3, t.rate),unit_adj)
#t4.events <- round_any(rexp(n.t4, t.rate),unit_adj)

ctrl.tot <- ctrl.events+exp.ctrl.entry
t1.tot <- t1.events+exp.entry.t1
t2.tot <- t2.events+exp.entry.t2
t3.tot <- t3.events+exp.entry.t3
#t4.tot <- t4.events+exp.entry.t4

ctrl.block <- length(c(seq(1, breaks[1], 3), seq(breaks[1]+1, breaks[2], 5)))

ctrl.cens1 <- ifelse(ctrl.tot[1:ctrl.block]>t1.time, 0, 1)
ctrl.obs1 <- ifelse(ctrl.cens1, ctrl.events[1:ctrl.block], 
                    t1.time-unif.ctrl.entry[1:ctrl.block])
t1.cens <- ifelse(t1.tot > t1.time, 0, 1)
t1.obs <- ifelse(t1.cens, t1.events, t1.time-unif.entry.t1)
t2.cens <- ifelse(t2.tot > t1.time, 0, 1)
t2.obs <- ifelse(t2.cens, t2.events, t1.time-unif.entry.t2)

ctrl.cens2 <- ifelse(ctrl.tot[(length(ctrl.tot)-ctrl.block+1):(length(ctrl.tot))]>t2.time, 0, 1)
ctrl.obs2 <- ifelse(ctrl.cens2, ctrl.events[(length(ctrl.tot)-ctrl.block+1):(length(ctrl.tot))],
                    t2.time-unif.ctrl.entry[(length(ctrl.tot)-ctrl.block+1):(length(ctrl.tot))])
t3.cens <- ifelse(t3.tot > t2.time, 0, 1)
t3.obs <- ifelse(t3.cens, t3.events, t2.time-unif.entry.t3)
#t4.cens <- ifelse(t4.tot > t2.time, 0, 1)
#t4.obs <- ifelse(t4.cens, t4.events, t2.time-unif.entry.t4)

obj1 <- survfit(Surv(t1.obs, t1.cens)~1)

ggsurvplot(obj1, data = as.data.frame(cbind(t1.obs, t1.cens)),
           conf.int = FALSE,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

obj2 <- survfit(Surv(t2.obs, t2.cens)~1)

ggsurvplot(obj2, data = as.data.frame(cbind(t2.obs, t2.cens)),
           conf.int = FALSE,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

obj3 <- survfit(Surv(t3.obs, t3.cens)~1)

ggsurvplot(obj3, data = as.data.frame(cbind(t3.obs, t3.cens)),
           conf.int = FALSE,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())

obj4 <- survfit(Surv(t4.obs, t4.cens)~1)

ggsurvplot(obj4, data = as.data.frame(cbind(t4.obs, t4.cens)),
           conf.int = FALSE,
           risk.table = TRUE,
           break.time.by = 5,
           ggtheme = theme_light())
