sim.cts <- function(c.Y = 0, c.D = 0, beta.Y = 0, kappa = 0.1, link = "logit", MAF.threshold = 0.03, prop.Y = 0.05, prop.pos.Y = 0.50, prop.D=0.05, prop.pos.D = 0.50, n1=1000, n.resampling=1000, r.corr=seq(0,1,0.1), R.min=1, R.max=100, S=10, track=F, sr.mode=F, file.regions=NA, file.pos=NA, folder.hap=NA, file.out=NA) {
     
     n.corr    <- length(r.corr)
     regions   <- read.table(file.regions)
     SNPInfo   <- read.table(file.pos,header=T)
     
     # define link function
     if (link == "logit") {
          gD        <- logit
          gD.inv    <- expit
          lambda.D   <- 1
     } else {
          gD        <- qnorm
          gD.inv    <- pnorm
          lambda.D   <- sqrt(3)/pi
     }
     
     # set initial baseline
     beta.0 <- gD(kappa) - lambda.D*(0.5^2+beta.Y*0.5^2)
     
     # calculate required sample size for simulations
     N <- calc.N(kappa,1e-4,n1) 
     
     # initialize matrix of results
     file.out.rdata <- paste("results/RData/",file.out,".RData",sep="")
     if (sr.mode & file.exists(file.out.rdata)) {
          results   <- local(get(load(file.out.rdata)))
          ssim      <- which(rowSums(is.na(results$p.values))==0) # successful simulations
          if (length(ssim)>0) {
               lastsim   <- results$p.values[max(ssim),] # last successful simulation  
               r = lastsim$r-R.min+1
               s = lastsim$s+1
               if (s > S) {
                    r = r+1
                    s = 1
               }
          } else {
               r = 1
               s = 1
          }
     } else {
          p.values        <- data.frame(matrix(NA,nrow=(R.max-R.min+1)*S,ncol=2+6*(1+n.corr)))
          characteristics <- data.frame(matrix(0,nrow=(R.max-R.min+1)*S,ncol=13))
          
          names(p.values) <- c("r", "s", paste("r",paste(c(seq(n.corr),"opt"),".ctrl",sep=""),sep=""), paste("r",paste(c(seq(n.corr),"opt"),".case",sep=""),sep=""), paste("r",paste(c(seq(n.corr),"opt"),".naive",sep=""),sep=""), paste("r",paste(c(seq(n.corr),"opt"),".joint",sep=""),sep=""), paste("r",paste(c(seq(n.corr),"opt"),".ipw",sep=""),sep=""), paste("r",paste(c(seq(n.corr),"opt"),".aipw",sep=""),sep=""))
          names(characteristics) <- c("kappa", "n.variants", "n.variants.03", "n.variants.01", "n.variants.001",  "n.ovariants", "n.ovariants.03", "n.ovariants.01", "n.ovariants.001", "n.cvariants", "n.cvariants.03", "n.cvariants.01", "n.cvariants.001")
          
          p.values$r <- rep(seq(R.min,R.max),each=S)
          p.values$s <- rep(seq(1,S),R.max-R.min+1)
          results    <- list(p.values=p.values,characteristics=characteristics)
          r = 1
          s = 1
     }
     
     # perform simulations
     while(r <= (R.max-R.min+1)) {
          
          results$p.values[((r-1)*S+1):(r*S),]$r <- R.min+r-1
          
          # haplotype matrix
          fname     <- paste(folder.hap,"/haplotypes",R.min+r-1,".txt",sep="")
          haplotype <- data.matrix(read.table(fname))
          haplotype <- subset(haplotype,rowSums(is.na(haplotype))==0)
          n.hap     <- dim(haplotype)[1]
          n.variants<- dim(haplotype)[2]
          
          # statistics of variants in subregion
          info  <- SNPInfo[regions[R.min+r-1,1]:regions[R.min+r-1,2],]
          FREQ  <- info$FREQ1 # allele frequencies
          cFREQ <- unlist(lapply(FREQ,function(x) min(x,1-x))) # minor allele frequencies
          results$characteristics[((r-1)*S+1):(r*S),]$n.variants       <- n.variants
          results$characteristics[((r-1)*S+1):(r*S),]$n.variants.03    <- sum(cFREQ < .03)
          results$characteristics[((r-1)*S+1):(r*S),]$n.variants.01    <- sum(cFREQ < .01)
          results$characteristics[((r-1)*S+1):(r*S),]$n.variants.001   <- sum(cFREQ < .001)
          
          # rare variants in subregion
          rvariants      <- which(cFREQ < MAF.threshold)   # index wrt variants in subregion
          n.rvariants    <- length(rvariants)              # number of rare variants
#           info           <- temp[rvariants,]                               # information on rare variants
          
          # generate population genotype matrix G
          hap1 <- sample(n.hap, N, replace=TRUE)
          hap2 <- sample(n.hap, N, replace=TRUE)
          G    <- haplotype[hap1,] + haplotype[hap2,]
          G[,FREQ>0.5] <- 2-G[,FREQ>0.5]
#           G    <- vG[,rvariants]
          
          # randomly select a set of rare variants to be causally associated with D
          cvariants <- sample(rvariants, n.rvariants*prop.D)
          beta.G <- matrix(0, nrow=n.variants, ncol=1)
          if (c.D != 0) { # randomly select a subset of the causal variants to have positive association with D
               positive       <- sample(cvariants, length(cvariants)*prop.pos.D)
               pcvariants     <- cvariants[cvariants %in% positive]  # causal variant w/ positive association (T/F)
               ncvariants     <- cvariants[!cvariants %in% positive] # causal variant w/ negative association (T/F)
               beta.G[pcvariants,] <- -1*c.D*log10(cFREQ[pcvariants])/2
               beta.G[ncvariants,] <- c.D*log10(cFREQ[ncvariants])/2
          }
          
          # randomly select a set of rare variants to be causally associated with Y
          cvariants <- sample(rvariants, n.rvariants*prop.Y)
          alpha.G <- matrix(0, nrow=n.variants, ncol=1)
          if (c.Y != 0) { # randomly select a subset of the causal variants to have positive association with Y
               positive       <- sample(cvariants, length(cvariants)*prop.pos.Y)          
               pcvariants     <- cvariants[cvariants %in% positive]  # causal variant w/ positive association (T/F)
               ncvariants     <- cvariants[!cvariants %in% positive] # causal variant w/ negative association (T/F)
               alpha.G[pcvariants,] <- -1*c.Y*log10(cFREQ[pcvariants])/2
               alpha.G[ncvariants,] <- c.Y*log10(cFREQ[ncvariants])/2
          }
          
          while (s <= S) {
               
               # Generate population X,Y,D
               X <- matrix(c(rnorm(N), rbinom(N, 1, 0.5)), ncol=2)
               Y <- X%*%matrix(0.5,nrow=2) + G%*%alpha.G + rnorm(N)
               D <- rbinom(N, 1, gD.inv(beta.0 + lambda.D*(X%*%matrix(0.5,nrow=2) + beta.Y*Y + G%*%beta.G)))
               kappa.hat <- mean(D)
               while (abs(kappa-kappa.hat) > kappa*5e-3) {
                    beta.0 <- beta.0 + kappa - kappa.hat
                    D <- rbinom(N, 1, gD.inv(beta.0 + lambda.D*(X%*%matrix(0.5,nrow=2) + beta.Y*Y + G%*%beta.G)))
                    kappa.hat <- mean(D)
               }
               
               if (sum(D) >= n1) {
                    
                    # track progress
                    if (track) {
                         print(paste(R.min+r-1,s), quote=F)
                    }
                    
                    results$p.values[(r-1)*S+s,]$s <- s
                    
                    # randomly sample n1 cases and n0=n1 controls
                    cases     <- sample(which(D==1), n1)
                    controls  <- sample(which(D==0), n1)
                    
                    results$characteristics[(r-1)*S+s,]$kappa            <- mean(D)                                           # disease prevalence
                    temp <- colSums(G[c(cases,controls),])>0
                    results$characteristics[(r-1)*S+s,]$n.ovariants      <- sum(temp)
                    results$characteristics[(r-1)*S+s,]$n.ovariants.03   <- sum(temp & (cFREQ < 0.03))
                    results$characteristics[(r-1)*S+s,]$n.ovariants.01   <- sum(temp & (cFREQ < 0.01))
                    results$characteristics[(r-1)*S+s,]$n.ovariants.001  <- sum(temp & (cFREQ < 0.001))
                    temp <- colSums(G[c(cases,controls),cvariants])>0
                    results$characteristics[(r-1)*S+s,]$n.cvariants      <- sum(temp)
                    results$characteristics[(r-1)*S+s,]$n.cvariants.03   <- sum(temp & (cFREQ[cvariants] < 0.03))
                    results$characteristics[(r-1)*S+s,]$n.cvariants.01   <- sum(temp & (cFREQ[cvariants] < 0.01))
                    results$characteristics[(r-1)*S+s,]$n.cvariants.001  <- sum(temp & (cFREQ[cvariants] < 0.001))
                    
                    # control-only SKAT-O
                    Y1   <- Y[controls] 
                    X1   <- X[controls, ]
                    ovariants <- which( colSums( G[controls, ] ) > 0 )
                    G1   <- G[controls, ovariants]
                    MAF  <- cFREQ[ovariants]
                    obj  <- SKAT_Null_Model(Y1~X1)
                    out  <- try(SKAT(G1, obj, weights=dbeta(MAF,1,25), r.corr=r.corr))
                    if ( class(out) != "try-error" ) {
                      results$p.values[(r-1)*S+s,3:(3+n.corr)] <- c(out$param$p.val.each, out$p.value)
                    }
                    
                    # case-only SKAT-O
                    Y1   <- Y[cases] 
                    X1   <- X[cases, ]
                    ovariants <- which( colSums( G[cases, ] ) > 0 )
                    G1   <- G[cases, ovariants]
                    MAF  <- cFREQ[ovariants]
                    obj  <- SKAT_Null_Model(Y1~X1)
                    out  <- try(SKAT(G1, obj, weights=dbeta(MAF,1,25), r.corr=r.corr))
                    if ( class(out) != "try-error" ) {
                      results$p.values[(r-1)*S+s,(4+n.corr):(4+2*n.corr)] <- c(out$param$p.val.each, out$p.value)
                    }
                    
                    # naive SKAT-O
                    Y1   <- Y[c(cases,controls)] 
                    X1   <- X[c(cases,controls), ]
                    ovariants <- which( colSums( G[c(cases,controls), ] ) > 0 )
                    G1   <- G[c(cases,controls), ovariants]
                    MAF  <- cFREQ[ovariants]
                    obj  <- SKAT_Null_Model(Y1~X1)
                    out  <- try(SKAT(G1, obj, weights=dbeta(MAF,1,25), r.corr=r.corr))
                    if ( class(out) != "try-error" ) {
                      results$p.values[(r-1)*S+s,(5+2*n.corr):(5+3*n.corr)] <- c(out$param$p.val.each, out$p.value)
                    }
                    
                    # joint SKAT-O
                    X1   <- cbind(X,D)[c(cases,controls), ]
                    obj  <- SKAT_Null_Model(Y1~X1)
                    out  <- try(SKAT(G1, obj, weights=dbeta(MAF,1,25), r.corr=r.corr))
                    if ( class(out) != "try-error" ) {
                      results$p.values[(r-1)*S+s,(6+3*n.corr):(6+4*n.corr)] <- c(out$param$p.val.each, out$p.value)
                    }
                    
                    # IPW SKAT-O
                    X1   <- X[c(cases,controls), ]
                    D1   <- D[c(cases,controls)]
                    weights.ipw <- D1*2*mean(D) + (1-D1)*2*mean(1-D)
                    obj <- IPWSKAT_Null_Model(Y1~X1,weights.ipw=weights.ipw,ccstatus=D1,n.resampling=n.resampling)
                    opt.asy <- try(IPWSKAT(G1, obj, weights.snp=dbeta(MAF,1,25), r.corr=r.corr, adjustment=F))
                    opt.adj <- try(IPWSKAT(G1, obj, weights.snp=dbeta(MAF,1,25), r.corr=r.corr, adjustment=T))
                    if (inherits(opt.asy,"try-error") | inherits(opt.adj,"try-error")) {
                         s = s - 1
                    } else {
                         results$p.values[(r-1)*S+s,(7+4*n.corr):(8+6*n.corr)]  <- c(opt.asy$param$p.val.each, opt.asy$p.value, opt.adj$param$p.val.each, opt.adj$p.value)
                    }
                    
                    # save progress
                    if (sr.mode) {
                         save(results, file=file.out.rdata)
                    }
                    
                    s = s+1
                    
               }
          }
          
          r = r+1
          s = 1
          
     }
     
     save(results, file=file.out.rdata)
     write.table(results$p.values, file=paste("results/pvalues/",file.out,".pvalues",sep=""), quote=F, row.names=F)
     write.table(results$characteristics, file=paste("results/char/",file.out,".char",sep=""), quote=F, row.names=F)
     
}
