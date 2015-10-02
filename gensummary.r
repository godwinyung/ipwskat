###
# FUNCTIONS 
###

combine <- function(kmin,kmax,suffix) {
  # initialize
  out <- local(get(load(paste(suffix,kmin,".RData",sep=""))))
  
  # concatenate results
  for (k in (kmin+1):kmax) {
    fname <- paste(suffix,k,".RData",sep="")
    if (file.exists(fname)) {
      temp <- local(get(load(fname)))
      out$p.values <- rbind(out$p.values,temp$p.values)
      out$characteristics <- rbind(out$characteristics, temp$characteristics)
    }
  }
  
  # remove missing results
  has.na <- which(rowSums(is.na(out$p.values))>0)
  if (length(has.na)>0) {
    out$p.values <- out$p.values[-has.na,]
    out$characteristics <- out$characteristics[-has.na,]
  }
  
  return(out)
}

calc.t1e <- function(t1e, alpha=0.05) {
  out <- c()
  for (i in 1:length(t1e)) {
    for (j in 1:length(t1e[[i]])) {
      out <- cbind(out, colMeans(t1e[[i]][[j]]$p.values<alpha))
    }
  }
  out <- data.frame(out)
  return(out)
}

###
# MAIN
###

# combine scenario-specific results and save
t1e <- list()
for (i in 1:3) {
  t1e[[i]] <- list()
  for (j in 1:3) {
    print(paste(i,j),quote=F)
    print("Combining ...", quote=F)
    temp <- combine(1,10000,paste("t1e",i,j,"_",sep=""))
    print("Saving ...", quote=F)
    save(temp,file=paste("t1e",i,j,"_all.RData",sep=""))
    t1e[[i]][[j]] <- temp
  }
}

# load scenario-specific results
t1e <- list()
for (i in 1:3) {
  t1e[[i]] <- list()
  for (j in 1:3) {
    t1e[[i]][[j]] <- local(get(load(paste("t1e",i,j,"_all.RData",sep=""))))
  }
}

t1e1 <- calc.t1e(t1e,0.05)
t1e2 <- calc.t1e(t1e,0.01)

round(t1e1[c("full","naive","ctrl","case","adj","ropt.asy1","ropt.adj1","ropt.asy2","ropt.adj2"),],3)
round(t1e2[c("full","naive","ctrl","case","adj","ropt.asy1","ropt.adj1","ropt.asy2","ropt.adj2"),],3)