version = "V7N500"

n = 500
p = 10
npc = 2*(p-1)
beta0=c(0.7,0.7,0.7,-0.5,-0.5,-0.5,0.3,0.3,0.3,0)
sigmaZ = 1

Nsetup = 1
setup = c("H")
isetup = 1


m = 500

Nrep = 500
Nperturb = 1000
Nlam = 100
min.lam = 1e-4
Nfold = 5
repid = formatC(1:Nrep,width = 1+floor(log(Nrep,10)), flag = 0)

source("source/sim-V7-run.R")