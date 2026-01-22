library(mnormt)
library(purrr)
library(survival)
library(gsDesign)

n.arms <- 6
arm.start <- c(0,0,0,0,0,6)
null.rate <-.693/12
hr <- rep(.75, n.arms)
beta <- .9
alpha <- rep(1-(.95)^(1/n.arms), n.arms)
A <- 1
rate <- 500
prop <- .7

tx.rate <- hr*null.rate

expt_events <- function(n_in, t=1, haz){
  if(n_in==0){
    return(0)
  }else{}
  u <- seq(0, t-1/n_in, 1/n_in)
  expt <- n_in - sum(exp(haz*(u-t)))
  return(expt)
}

mort <- function(input, hz, t){
  return((1-exp(-hz*t))*input)
}

curr.arms <- arm.start==0
act <- 1 + sum(curr.arms)
evs <- ceiling(nEvents(hr=hr, alpha=alpha, beta=1-beta, ratio=A))

ctrl.events <- list()
tx.events <- list()
ctrl.patients <- 0

for (i in 1:n.arms) {
  tx.events[[i]] <- data.frame(patients = 0,
                                 cum_pt = 0,
                                 events = 0,
                                 cum_events = 0)
  ctrl.events[[i]] <- data.frame(at_risk=0)
}

t<-1

while (sum(curr.arms)>0) {
  pat <- floor(rate/sum(curr.arms,1))
  ctrl.patients[t] <- pat
  
  for (j in which(curr.arms==T)) {
    tx.events[[j]][t,1] <- pat
    tx.events[[j]][t,2] <- sum(tx.events[[j]][1:t,1], na.rm = T)
    if(t==1){
      tx.events[[j]][1,3] <-round(expt_events(pat,1,tx.rate[j]))
    }else{
      cumpt <- ifelse(is.na(tx.events[[j]][t-1,2]), 0,tx.events[[j]][t-1,2])
      cumev <- ifelse(is.na(tx.events[[j]][t-1,4]), 0,tx.events[[j]][t-1,4])
      tx.events[[j]][t,3] <-round(expt_events(pat,1,tx.rate[j]))+
        round(mort(cumpt-cumev,tx.rate[j],1))
    }
    tx.events[[j]][t,4] <- sum(tx.events[[j]][1:t,3], na.rm = T)
  }
  
  
  for (i in 1:n.arms) {
    if (t!=1){
      ctrl.events[[i]][,ncol(ctrl.events[[i]])+1] <- NA
    }
    for (j in 1:t) {
        if (j<ncol(ctrl.events[[i]])){
          ctrl.events[[i]][j,t-j+1] <- ctrl.events[[i]][j,t-j]-
            round(mort(ctrl.events[[i]][j,t-j],null.rate,t-j+1))
        }else{
          if (curr.arms[i]==F){
            ctrl.events[[i]][j,t-j+1] <- 0
          }else{
          ctrl.events[[i]][j,t-j+1] <- pat-round(expt_events(pat,1,null.rate))
        }
      }
    } 
  }
  
  for (k in which(curr.arms==T)) {
    on <- which(ctrl.events[[k]][,1]!=0)
    if (length(on)==0){
      c.evs <- 0
    }else if (length(on)==1){
      c.evs <- ctrl.patients[on]-min(ctrl.events[[k]][on,], na.rm = T)
    }else {
      c.evs <- sum(ctrl.patients[on]-apply(ctrl.events[[k]][on,], 1, FUN = min, na.rm = TRUE))
    }
    n.evs <- c.evs + max(tx.events[[k]]$cum_events, na.rm = T)
    if (n.evs > prop*evs[k]) {
      curr.arms[k] <-F
    }
  }
  
  for (k in which(curr.arms==F)) {
    if(arm.start[k]==t){
      curr.arms[k]<-T
    }
  }
  
  t<-t+1
}

(cpat <- sum(ctrl.patients))
for (i in 1:n.arms) {
  print(max(tx.events[[i]]$cum_pt, na.rm = T))
}

