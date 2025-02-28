rm(list=objects())
version = "pcaPO-N500"

source("source/utility.R")
source("source/sim_gen_rev.R")
n = 500
p = 10
npc = 2*(p-1)
beta0=c(0.7,0.7,0.7,-0.5,-0.5,-0.5,0.3,0.3,0.3,0)
sigmaZ = 1

Nsetup = 1
setup = c("S")
isetup = 1

misT = noise.loglinear.S
link = expit
dlink = dexpit
inv.link = logit

m = 500

Nrep = 500
Nperturb = 1000
Nlam = 100
min.lam = 1e-4
Nfold = 5
repid = formatC(1:Nrep,width = 1+floor(log(Nrep,10)), flag = 0)

source("source/sim-adaPCA-run.R")
rm(list=objects())