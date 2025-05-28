control.rate <- 333
expt.rate <- 167

n.control <- 789
n.expt <- 374

control.entry <- seq(0, 2.36, length=n.control)
expt.entry <- seq(0, 2.36, length=n.expt)

control.event <- rexp(n.control, rate=0.693)
control.event.cal <- control.entry + control.event

sum(control.event.cal < 2.36)


n.control2 = round(control.rate*0.36)

control.entry.overlap <- seq(2, 2.36, length=n.control2)
control.event.overlap <- rexp(n.control2, rate=0.693)
control.event.overlap.cal <- control.entry.overlap + control.event.overlap

sum(control.event.overlap.cal < 2.36)
