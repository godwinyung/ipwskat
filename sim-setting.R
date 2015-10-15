args <- (commandArgs(TRUE))
eval(parse(text=args[[1]])) #taskID
eval(parse(text=args[[2]])) #interval = number of haplotypes analyzed per job in the cluster

#####
### PACKAGES
#####

library(SKAT)
library(geepack)

#####
### INPUT AND SOURCES
#####

file.regions   <- "/n/regal/xlin/godwinyung/ipwskat/data/regions.txt"
file.pos       <- "/n/regal/xlin/godwinyung/ipwskat/data/out.pos-1"
folder.hap     <- "/n/regal/xlin/godwinyung/ipwskat/data/haplotypes"
source("/n/regal/xlin/godwinyung/ipwskat/code/ipwskat.r")
source("/n/regal/xlin/godwinyung/ipwskat/code/sim-cts.r")

#####
### SETTING
#####

MAF.threshold  = 0.03

# primary model (D)
n1        = 1000 # c(100,250,500,1000)
kappa     = 0.10
link      = "logit"
beta.Y    = c(0,log(2)/2,log(2))
c.D       = c(0,log(5)/2,log(2.5)/2)
prop.D    = c(0,0.20,0.50)
prop.pos.D= 0.50

# secondary model (Y)
c.Y       = c(0,log(2),log(2)/2)
prop.Y    = c(0,0.20,0.50)
prop.pos.Y= c(0.5,0.8,1.0)
r.corr    = seq(0,1,0.1)

R.min = interval*(taskID-1)+1
R.max = interval*taskID
S=1

#####
### MAIN
#####

# type I error

# sim.cts(c.Y = c.Y[1], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("t1e",1,1,"_",taskID,sep=""))
# 
# sim.cts(c.Y = c.Y[1], c.D = c.D[1], beta.Y = beta.Y[2], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("t1e",1,2,"_",taskID,sep=""))
# 
# sim.cts(c.Y = c.Y[1], c.D = c.D[1], beta.Y = beta.Y[3], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("t1e",1,3,"_",taskID,sep=""))
# 
# sim.cts(c.Y = c.Y[1], c.D = c.D[2], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[2], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("t1e",2,1,"_",taskID,sep=""))
# 
# sim.cts(c.Y = c.Y[1], c.D = c.D[2], beta.Y = beta.Y[2], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[2], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("t1e",2,2,"_",taskID,sep=""))
# 
# sim.cts(c.Y = c.Y[1], c.D = c.D[2], beta.Y = beta.Y[3], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[2], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("t1e",2,3,"_",taskID,sep=""))
# 
# sim.cts(c.Y = c.Y[1], c.D = c.D[3], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[3], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("t1e",3,1,"_",taskID,sep=""))
# 
# sim.cts(c.Y = c.Y[1], c.D = c.D[3], beta.Y = beta.Y[2], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[3], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("t1e",3,2,"_",taskID,sep=""))
# 
# sim.cts(c.Y = c.Y[1], c.D = c.D[3], beta.Y = beta.Y[3], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[3], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("t1e",3,3,"_",taskID,sep=""))

# power

sim.cts(c.Y = c.Y[1], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("power",1,1,"_",taskID,sep=""))

sim.cts(c.Y = c.Y[2], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[2], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("power",1,2,"_",taskID,sep=""))

sim.cts(c.Y = c.Y[3], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[3], prop.pos.Y = prop.pos.Y[1], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("power",1,3,"_",taskID,sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[2], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("power",2,1,"_",taskID,sep=""))

sim.cts(c.Y = c.Y[2], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[2], prop.pos.Y = prop.pos.Y[2], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("power",2,2,"_",taskID,sep=""))

sim.cts(c.Y = c.Y[3], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[3], prop.pos.Y = prop.pos.Y[2], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("power",2,3,"_",taskID,sep=""))

sim.cts(c.Y = c.Y[1], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[1], prop.pos.Y = prop.pos.Y[3], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("power",3,1,"_",taskID,sep=""))

sim.cts(c.Y = c.Y[2], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[2], prop.pos.Y = prop.pos.Y[3], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("power",3,2,"_",taskID,sep=""))

sim.cts(c.Y = c.Y[3], c.D = c.D[1], beta.Y = beta.Y[1], kappa = kappa, link = link, MAF.threshold = MAF.threshold, prop.Y = prop.Y[3], prop.pos.Y = prop.pos.Y[3], prop.D=prop.D[1], prop.pos.D = prop.pos.D, n1=n1, n.resampling=10000, r.corr=r.corr, R.min=R.min, R.max=R.max, S=S, track=T, sr.mode=T, file.regions=file.regions, file.pos=file.pos, folder.hap=folder.hap, file.out=paste("power",3,3,"_",taskID,sep=""))