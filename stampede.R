library(tidyverse)
library(survival)
library(beepr)

one_arm <- function(nsim, n, rate) {
  map(1:nsim, ~rexp(n, rate=rate))
}

nsim <- 10000

ctrl_arm <- one_arm(10000,4754, 1)
b <- one_arm(10000,593,1) #7.25 yr #2 fu
c <- one_arm(10000,592,1) #7.25 yr #2 fu
d <- one_arm(10000,321,1) #5.25yr #4.75 fu
e <- one_arm(10000,593,1) #7.25 yr #2 fu
f <- one_arm(10000,311,1) #5.25yr #4.75 fu
g <- one_arm(10000,960,1) #2.25 yr #3.3 fu
h <- one_arm(10000,1032,1) #3.6 yr #2 fu
j <- one_arm(10000,989,1) #1.75 yr #4fu
k <- one_arm(10000,1059,1) #6.5 yr #3.6fu
l <- one_arm(10000,347,1) #6 yr #3.6fu


#med fu: 43 mo

ctrl_entry <- seq(0, 17.5, length=4754)
b_entry <- seq(0, 7.25, length=593)
c_entry <- seq(0, 7.25, length=592)
d_entry <- seq(0, 5.25, length=321)
e_entry <- seq(0, 7.25, length=593)
f_entry <- seq(0, 5.25, length=311)
g_entry <- seq(6, 8.25, length=960)
h_entry <- seq(7.25, 10.85, length=1032)
j_entry <- seq(8.75, 10.5, length=989)
k_entry <- seq(11, 17.5, length=1059)
l_entry <- seq(11.5, 17.5, length=347)

ctrl.list <- list()
b.list <- list()
c.list <- list()
d.list <- list()
e.list <- list()
f.list <- list()
g.list <- list()
h.list <- list()
j.list <- list()
k.list <- list()
l.list <- list()


for (i in 1:nsim) {
  
  ctrl.list[[i]] <- ctrl_arm[[i]]+ctrl_entry
  b.list[[i]] <- b[[i]] +b_entry
  c.list[[i]] <- c[[i]] +c_entry
  d.list[[i]] <- d[[i]] +d_entry
  e.list[[i]] <- e[[i]] +e_entry
  f.list[[i]] <- f[[i]] +f_entry
  g.list[[i]] <- g[[i]] +g_entry
  h.list[[i]] <- h[[i]] +h_entry
  j.list[[i]] <- j[[i]] +j_entry
  k.list[[i]] <- k[[i]] +k_entry
  l.list[[i]] <- l[[i]] +l_entry
};beep()

censor <- function(data, entry, fu) {
  ifelse(data>entry+fu, 0,1)
}

fu=3.6

ctrl_cens <- lapply(ctrl.list, censor, entry=ctrl_entry, fu=fu)
b_cens <- lapply(b.list, censor, entry=b_entry, fu=fu)
c_cens <- lapply(c.list, censor, entry=c_entry, fu=fu)
d_cens <- lapply(d.list, censor, entry=d_entry, fu=fu)
e_cens <- lapply(e.list, censor, entry=e_entry, fu=fu)
f_cens <- lapply(f.list, censor, entry=f_entry, fu=fu)
g_cens <- lapply(g.list, censor, entry=g_entry, fu=fu)
h_cens <- lapply(h.list, censor, entry=h_entry, fu=fu)
j_cens <- lapply(j.list, censor, entry=j_entry, fu=fu)
k_cens <- lapply(k.list, censor, entry=k_entry, fu=fu)
l_cens <- lapply(l.list, censor, entry=l_entry, fu=fu)

ctrl_be <- c(rep(0, length(ctrl_entry[1:which.min(ctrl_entry<7.25)-1])), 
             rep(1, length(b_entry)))
ctrl_c <- c(rep(0, length(ctrl_entry[1:which.min(ctrl_entry<7.25)-1])), 
            rep(1, length(c_entry)))
ctrl_d <- c(rep(0, length(ctrl_entry[1:which.min(ctrl_entry<5.25)-1])), 
            rep(1, length(d_entry)))
ctrl_f <- c(rep(0, length(ctrl_entry[1:which.min(ctrl_entry<5.25)-1])), 
            rep(1, length(f_entry)))
ctrl_g <- c(rep(0, length(ctrl_entry[which.min(ctrl_entry<6):which.min(ctrl_entry<8.25)-1])),
            rep(1, length(g_entry)))
ctrl_h <- c(rep(0, length(ctrl_entry[which.min(ctrl_entry<7.25):which.min(ctrl_entry<10.85)-1])),
            rep(1, length(h_entry)))
ctrl_j <- c(rep(0, length(ctrl_entry[which.min(ctrl_entry<8.75):which.min(ctrl_entry<10.5)-1])),
            rep(1, length(j_entry)))
ctrl_k <- c(rep(0, length(ctrl_entry[which.min(ctrl_entry<11):which.min(ctrl_entry<17.5)-1])),
            rep(1, length(k_entry)))
ctrl_l <- c(rep(0, length(ctrl_entry[which.min(ctrl_entry<11.5):which.min(ctrl_entry<17.5)-1])),
            rep(1, length(l_entry)))

getps <- function(censc, censt, ctrl.time, ctrl1, ctrl2, tx.time, fu, compsn, nsim){
  p <- rep(NA, nsim)
  for (i in 1:nsim) {
    ttec <- ifelse(censc[[i]][ctrl1:ctrl2], ctrl.time[[i]][ctrl1:ctrl2], fu)
    ttet <-ifelse(censt[[i]], tx.time[[i]], fu)
    test <- summary(coxph(Surv(c(ttec, ttet), c(censc[[i]][ctrl1:ctrl2], censt[[i]]))~compsn))
    p[i] <- test$coef[5]
  }
  return(p)
}

pb <- getps(ctrl_cens, b_cens, ctrl_arm, 1, which.min(ctrl_entry<7.25)-1, 
            b, fu, ctrl_be, nsim);beep()
pc <- getps(ctrl_cens, c_cens, ctrl_arm, 1, which.min(ctrl_entry<7.25)-1, 
            c, fu, ctrl_c, nsim);beep()
pd <- getps(ctrl_cens, d_cens, ctrl_arm, 1, which.min(ctrl_entry<5.25)-1, 
            d, fu, ctrl_d, nsim);beep()
pe <- getps(ctrl_cens, e_cens, ctrl_arm, 1, which.min(ctrl_entry<7.25)-1, 
            e, fu, ctrl_be, nsim);beep()
pf <- getps(ctrl_cens, f_cens, ctrl_arm, 1, which.min(ctrl_entry<5.25)-1, 
            f, fu, ctrl_f, nsim);beep()
pg <- getps(ctrl_cens, g_cens, ctrl_arm,which.min(ctrl_entry<6),which.min(ctrl_entry<8.25)-1, 
            g, fu, ctrl_g, nsim);beep()
ph <- getps(ctrl_cens, h_cens, ctrl_arm,ctrl1=which.min(ctrl_entry<7.25),ctrl2=which.min(ctrl_entry<10.85)-1, 
            h, fu, ctrl_h, nsim);beep()
pj <- getps(ctrl_cens, j_cens, ctrl_arm,which.min(ctrl_entry<8.75),which.min(ctrl_entry<10.5)-1, 
            j, fu, ctrl_j, nsim);beep()
pk <- getps(ctrl_cens, k_cens, ctrl_arm,which.min(ctrl_entry<11),which.min(ctrl_entry<17.5)-1, 
            k, fu, ctrl_k, nsim);beep()
pl <- getps(ctrl_cens, l_cens, ctrl_arm,which.min(ctrl_entry<11.5),which.min(ctrl_entry<17.5)-1, 
            l, fu, ctrl_l, nsim);beep()

rej_b <- ifelse(pb<.05, 1, 0)
rej_c <- ifelse(pc<.05, 1, 0)
rej_d <- ifelse(pd<.05, 1, 0)
rej_e <- ifelse(pe<.05, 1, 0)
rej_f <- ifelse(pf<.05, 1, 0)
rej_g <- ifelse(pg<.05, 1, 0)
rej_h <- ifelse(ph<.05, 1, 0)
rej_j <- ifelse(pj<.05, 1, 0)
rej_k <- ifelse(pk<.05, 1, 0)
rej_l <- ifelse(pl<.05, 1, 0)

tot <- rej_b+rej_c+rej_d+rej_e+rej_f+rej_g+rej_h+rej_j+rej_k+rej_l
flat <- ifelse(tot>=1, 1, 0)
mean(flat)
