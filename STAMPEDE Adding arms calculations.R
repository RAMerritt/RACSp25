A <- 0.5
start2 <- 0.5
expt.length <- 2.36
control.length <- expt.length + start2

patients.per.unit.time <- 500
control.rate <- round(1 / (A + 1) * patients.per.unit.time)
expt.rate <- round(A / (A + 1) * patients.per.unit.time)

n.control <- round(control.length * control.rate)
n.expt1 <- round(expt.length * expt.rate)
n.expt2 <- round(expt.length * expt.rate)

# Entry times treated as deterministic here--should they be
# random?


control.entry <- seq(0, control.length, length = n.control)
expt1.entry <- seq(0, expt.length, length = n.expt1)
expt2.entry <- seq(start2, start2 + expt.length, length = n.expt2)

# Treatment events should be similarly generated, but with
# a rate corresponding to whatever value of the hazard ratio
# we want to explore

h1 <- h2 <- 0.75
lam0 <- 0.693
lam1 <- lam0 * h1
lam2 <- lam0 * h2

nsim <- 10000
stat1 <- stat2 <- pval1 <- pval2 <- rep(0, nsim)
e0 <- rep(0, nsim)
for (i in 1:nsim) {
  control.event <- rexp(n.control, rate = lam0)
  expt1.event <- rexp(n.expt1, rate = lam1)
  expt2.event <- rexp(n.expt2, rate = lam2)
  
  # Calendar time of events
  control.event.cal <- control.entry + control.event
  expt1.event.cal <- expt1.entry + expt1.event
  expt2.event.cal <- expt2.entry + expt2.event
  
  # Separate control into patients who accrued by end of first
  # study, and patients who accrued by end of second study
  control1.event.cal <- control.event.cal[control.entry <= expt.length]
  control2.event.cal <- control.event.cal[control.entry >= start2]
  
  control.event1 <- control.event[control.entry <= expt.length]
  control.event2 <- control.event[control.entry >= start2]
  
  control.entry1 <- control.entry[control.entry <= expt.length]
  control.entry2 <- control.entry[control.entry >= start2]
  
  e0[i] <- sum(control.event.cal < 2.36)
  
  # Censored if the calendar event time is greater than
  # expt.length
  control1.death <- control1.event.cal <= expt.length
  control2.death <- control2.event.cal <= (expt.length + start2)
  expt1.death <- expt1.event.cal <= expt.length
  expt2.death <- expt2.event.cal <= (expt.length + start2)
  
  # Observation times:
  # expt.length - entry time if censored
  # event.time if not censored
  control1.obsTime <- ifelse(control1.death, control.event1, expt.length - control.entry1)
  control2.obsTime <- ifelse(control2.death,
                             control.event2,
                             expt.length + start2 - control.entry2)
  expt1.obsTime <- ifelse(expt1.death, expt1.event, expt.length - expt1.entry)
  expt2.obsTime <- ifelse(expt2.death, expt2.event, expt.length + start2 - expt2.entry)
  
  
  event.time1 <- c(control1.obsTime, expt1.obsTime)
  delta.1 <- c(control1.death, expt1.death)
  grp1 <- rep(c(0, 1), c(length(control1.obsTime), length(expt1.obsTime)))
  
  event.time2 <- c(control2.obsTime, expt2.obsTime)
  delta.2 <- c(control2.death, expt2.death)
  grp2 <- rep(c(0, 1), c(length(control2.obsTime), length(expt2.obsTime)))
  
  rslt1 <- summary(coxph(Surv(event.time1, delta.1) ~ grp1))
  rslt2 <- summary(coxph(Surv(event.time2, delta.2) ~ grp2))
  
  stat1[i] <- rslt1$coef[4]
  stat2[i] <- rslt2$coef[4]
  
  pval1[i] <- rslt1$coef[5]
  pval2[i] <- rslt2$coef[5]
  
  if (i %% 500 == 0) {
    cat("Simulation", i, "donzo!", "\n")
  }
}

mean(pval1 < 0.05)
mean(pval2 < 0.05)

# Disjunctive Power
mean(pval1 < 0.05 | pval2 < 0.05)
# Conjunctive Power
mean(pval1 < 0.05 & pval2 < 0.05)
# FWER
mean(pval1 < 0.025 | pval2 < 0.025)
# Test Stat Correlation
cor(stat1, stat2)


# They seem to be using alpha 0.025 for FWER calculations but 
# alpha = 0.05 for Power

#######

# Confirm that on average we get 401 events by 
# the end of the study
sum(control.event.cal < 2.36)

