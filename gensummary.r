###
# FUNCTIONS 
###

calc.power <- function(pval, alpha=0.05) {
  out <- c()
  for (i in 1:length(pval)) {
    for (j in 1:length(pval[[i]])) {
      out <- cbind(out, colMeans(pval[[i]][[j]]<alpha))
    }
  }
  out <- data.frame(out)
  return(out)
}

###
# MAIN
###

# combine scenario-specific results
pval <- list()
for (i in 1:3) {
  pval[[i]] <- list()
  for (j in 1:3) {
    pval[[i]][[j]] <- read.table(paste("power",i,j,".pvalues",sep=""), header=T)[,3:74]
  }
}

power1 <- calc.power(pval,0.05)
power2 <- calc.power(pval,0.01)

round(power1[c("r1.ctrl","r1.case","r1.naive","r1.joint","r1.ipw","r1.aipw"),],3)
round(power1[c("r11.ctrl","r11.case","r11.naive","r11.joint","r11.ipw","r11.aipw"),],3)
round(power1[c("ropt.ctrl","ropt.case","ropt.naive","ropt.joint","ropt.ipw","ropt.aipw"),],3)


