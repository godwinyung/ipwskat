args <- (commandArgs(TRUE))
eval(parse(text=args[[1]])) #taskID1
eval(parse(text=args[[2]])) #taskID2
eval(parse(text=args[[3]])) #njobs = total number of jobs to be submitted to the cluster
source("ipwskat.r")
source("sim-cts.R")

#####
### PACKAGES
#####

library(SKAT)
library(geepack)

#####
### INPUT
#####

input          <- "/n/home01/godwin/Research/ipwskat/data/cosi/"
regions        <- read.table(paste(input,"10Kregions.txt",sep=""))
SNPInfo        <- read.table(paste(input,"out.pos-1",sep=""),header=T)

#####
### SETTING
#####

c.Y = c(0,log(5)/2,log(2.5)/2)
c.D = c(0,log(5)/2,log(2.5)/2)
beta.Y = c(0,log(2)/2,log(2))
kappa = 0.1
link = "logit"
MAF.threshold = 0.03
prop.Y = c(0,0.20,0.50)
prop.pos.Y = 0.50
prop.D = c(0,0.20,0.50)
prop.pos.D = 0.50
n1 = 1000 # c(100,250,500,1000)
nhap  = 10000  # number of haplotypes
# njobs = 100
interval = 1000/njobs
R.min = nhap/njobs*(taskID1-1)+interval*(taskID2-1)+1
R.max = nhap/njobs*(taskID1-1)+interval*(taskID2-1)+interval
S=1

sink(paste("settings/sim",taskID1,"_",taskID2,".txt",sep=""))
cat("Simulation for continuous (v1).R\n\n")
cat("Settings:\n")
cat("n1 =", n1, "\n")
cat("R.min =", R.min, "\n")
cat("R.max =", R.max, "\n")
cat("S =", S, "\n")
sink()

#####
### MAIN
#####

sim.cts(c.Y = c.Y[1], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y, prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=seq(0,1,0.1), R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, filename=paste("../results/t1e",j,k,"_",taskID1,"_",taskID2,".RData",sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[1], beta.Y = beta.Y[2], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y, prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=seq(0,1,0.1), R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, filename=paste("../results/t1e",j,k,"_",taskID1,"_",taskID2,".RData",sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[1], beta.Y = beta.Y[3], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y, prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=seq(0,1,0.1), R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, filename=paste("../results/t1e",j,k,"_",taskID1,"_",taskID2,".RData",sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[2], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y, prop.D=prop.D[2], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=seq(0,1,0.1), R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, filename=paste("../results/t1e",j,k,"_",taskID1,"_",taskID2,".RData",sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[2], beta.Y = beta.Y[2], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y, prop.D=prop.D[2], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=seq(0,1,0.1), R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, filename=paste("../results/t1e",j,k,"_",taskID1,"_",taskID2,".RData",sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[2], beta.Y = beta.Y[3], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y, prop.D=prop.D[2], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=seq(0,1,0.1), R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, filename=paste("../results/t1e",j,k,"_",taskID1,"_",taskID2,".RData",sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[3], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y, prop.D=prop.D[3], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=seq(0,1,0.1), R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, filename=paste("../results/t1e",j,k,"_",taskID1,"_",taskID2,".RData",sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[3], beta.Y = beta.Y[2], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y, prop.D=prop.D[3], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=seq(0,1,0.1), R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, filename=paste("../results/t1e",j,k,"_",taskID1,"_",taskID2,".RData",sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[3], beta.Y = beta.Y[3], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y, prop.D=prop.D[3], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=seq(0,1,0.1), R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, filename=paste("../results/t1e",j,k,"_",taskID1,"_",taskID2,".RData",sep=""))

for (r in R.min:R.max) {
     fname = paste(input,"10Khaplotypes/haplotypes",r,".txt",sep="")
     if (file.exists(fname)) {
          file.remove(fname)
          cat(paste("",r),file="output/done.txt",append=T)
     }
}

