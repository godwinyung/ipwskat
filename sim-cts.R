sim.cts <- function(c.Y = 0, c.D = 0, beta.Y = 0, kappa = 0.1, link = "logit", MAF.threshold = 0.03, prop.Y = 0.05, prop.pos.Y = 0.50, prop.D=0.05, prop.pos.D = 0.50, n1=1000, n.resampling=1000, r.corr=seq(0,1,0.1), R.min=1, R.max=100, S=10, track=F, sr.mode=F, file.regions=NA, file.pos=NA, folder.hap=NA, file.out=NA) {
     
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
          p.values        <- data.frame(matrix(NA,nrow=(R.max-R.min+1)*S,ncol=7+2*(1+length(r.corr))))
          characteristics <- data.frame(matrix(0,nrow=(R.max-R.min+1)*S,ncol=11))
          
          names(p.values) <- c("r", "s", "full", "naive", "ctrl", "case", "adj", paste("r",paste(c(seq(length(r.corr)),"opt"),".asy1",sep=""),sep=""), paste("r",paste(c(seq(length(r.corr)),"opt"),".adj1",sep=""),sep=""))
          names(characteristics) <- c("kappa", "n.variants", "n.variants.03", "n.variants.01", "n.variants.001", "n.ovariants.cc", "n.ovariants.case", "n.ovariants.ctrl", "p.cc", "p.case", "p.ctrl")
          
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
          
          # variants in subregion
          temp <- SNPInfo[regions[R.min+r-1,1]:regions[R.min+r-1,2],]
          results$characteristics[((r-1)*S+1):(r*S),]$n.variants       <- dim(temp)[1]
          results$characteristics[((r-1)*S+1):(r*S),]$n.variants.03    <- sum(temp$FREQ1 < .03)
          results$characteristics[((r-1)*S+1):(r*S),]$n.variants.01    <- sum(temp$FREQ1 < .01)
          results$characteristics[((r-1)*S+1):(r*S),]$n.variants.001   <- sum(temp$FREQ1 < .001)
          
          # rare variants in subregion
          rvariants      <- which(temp$FREQ1 < MAF.threshold)              # index wrt variants in subregion
          n.rvariants    <- length(rvariants)                              # number of rare variants
          info           <- temp[rvariants,]                               # information on rare variants
          
          # Generate population genotype matrix vG (variants) and G (rare variants)
          hap1 <- sample(n.hap, N, replace=TRUE)
          hap2 <- sample(n.hap, N, replace=TRUE)
          vG   <- haplotype[hap1,] + haplotype[hap2,]
          G    <- vG[,rvariants]
          
          cvariants      <- sample(n.rvariants, n.rvariants*prop.D)        # randomly select a set of rare variants to be causally associated with D
          if (c.D != 0) { # randomly select a subset of the causal variants to have positive association with D
               positive       <- sample(cvariants, length(cvariants)*prop.pos.D)          # index wrt rvariants
               pcvariants     <- seq(1:n.rvariants) %in% positive                         # causal variant w/ positive association (T/F)
               ncvariants     <- seq(1:n.rvariants) %in% setdiff(cvariants, positive)     # causal variant w/ negative association (T/F)
               beta.G         <- matrix((-1*pcvariants + ncvariants) * c.D * log10(info$FREQ1), ncol=1)
          } else {
               beta.G         <- matrix(0, nrow=n.rvariants, ncol=1)
          }
          
          cvariants      <- sample(n.rvariants, n.rvariants*prop.Y)        # randomly select a set of rare variants to be causally associated with Y
          if (c.Y != 0) { # randomly select a subset of the causal variants to have positive association with Y
               positive       <- sample(cvariants, length(cvariants)*prop.pos.Y)          # index wrt rvariants
               pcvariants     <- seq(1:n.rvariants) %in% positive                         # causal variant w/ positive association (T/F)
               ncvariants     <- seq(1:n.rvariants) %in% setdiff(cvariants, positive)     # causal variant w/ negative association (T/F)
               alpha.G        <- matrix((-1*pcvariants + ncvariants) * c.Y * log10(info$FREQ1), ncol=1)
          } else {
               alpha.G        <- matrix(0, nrow=n.rvariants, ncol=1)
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
                    
                    results$characteristics[(r-1)*S+s,]$kappa            <- mean(D)                                             # disease prevalence
                    results$characteristics[(r-1)*S+s,]$n.ovariants.cc   <- sum(colSums(vG[c(cases,controls),])>0)              # number of variants observed by cases and ctrls
                    results$characteristics[(r-1)*S+s,]$n.ovariants.case <- sum(colSums(vG[cases,])>0)                          
                    results$characteristics[(r-1)*S+s,]$n.ovariants.ctrl <- sum(colSums(vG[controls,])>0)
                    results$characteristics[(r-1)*S+s,]$p.cc             <- sum(colSums(G[c(cases,controls),cvariants])>0)      # number of causal variants observed by cases and ctrls
                    results$characteristics[(r-1)*S+s,]$p.case           <- sum(colSums(G[cases,cvariants])>0)
                    results$characteristics[(r-1)*S+s,]$p.ctrl           <- sum(colSums(G[controls,cvariants])>0)
                    
                    MAF  <- info$FREQ1
                    
                    # analysis 1 : full cohort SKAT
                    Y1 <- Y
                    ovariants <- which(colSums(G) > 0)
                    G1 <- G[,ovariants]
                    MAF1 <- MAF[ovariants]
                    obj  <- SKAT_Null_Model(Y~X)
                    results$p.values[(r-1)*S+s,]$full <- SKAT(G1, obj, weights=dbeta(MAF1,1,25))$p.value
                    
                    # analysis 3 : control-only SKAT
                    Y1   <- Y[controls] 
                    X1   <- X[controls, ]
                    ovariants <- which( colSums( G[controls, ] ) > 0 )
                    G1   <- G[controls, ovariants]
                    MAF1 <- MAF[ovariants]
                    obj  <- SKAT_Null_Model(Y1~X1)
                    results$p.values[(r-1)*S+s,]$ctrl <- SKAT(G1, obj, weights=dbeta(MAF1,1,25))$p.value
                    
                    # analysis 4 : case-only SKAT
                    Y1   <- Y[cases] 
                    X1   <- X[cases, ]
                    ovariants <- which( colSums( G[cases, ] ) > 0 )
                    G1   <- G[cases, ovariants]
                    MAF1 <- MAF[ovariants]
                    obj  <- SKAT_Null_Model(Y1~X1)
                    results$p.values[(r-1)*S+s,]$case <- SKAT(G1, obj, weights=dbeta(MAF1,1,25))$p.value
                    
                    # analysis 2 : naive SKAT
                    Y1   <- Y[c(cases,controls)] 
                    X1   <- X[c(cases,controls), ]
                    ovariants <- which( colSums( G[c(cases,controls), ] ) > 0 )
                    G1   <- G[c(cases,controls), ovariants]
                    MAF1 <- MAF[ovariants]
                    obj  <- SKAT_Null_Model(Y1~X1)
                    results$p.values[(r-1)*S+s,]$naive <- SKAT(G1, obj, weights=dbeta(MAF1,1,25))$p.value
                    
                    # analysis 5 : adjusted SKAT
                    X1   <- cbind(X,D)[c(cases,controls), ]
                    obj  <- SKAT_Null_Model(Y1~X1)
                    results$p.values[(r-1)*S+s,]$adj <- SKAT(G1, obj, weights=dbeta(MAF1,1,25))$p.value
                    
                    # analysis 6 : ipw SKAT
                    X1   <- X[c(cases,controls), ]
                    D1   <- D[c(cases,controls)]
                    mu.D <- glm(D1~X1+Y1,family="binomial")$fitted.values
                    weights.snp <- dbeta(MAF1,1,25)
                    weights.ipw <- D1*2*mean(D) + (1-D1)*2*mean(1-D)
                    obj <- IPWSKAT_Null_Model(Y1~X1,weights.ipw=weights.ipw,ccstatus=D1,n.resampling=n.resampling)
                    opt.asy1 <- try(IPWSKAT(G1, obj, weights.snp.beta=dbeta(MAF1,1,25), r.corr=r.corr, adjustment=F))
                    opt.adj1 <- try(IPWSKAT(G1, obj, weights.snp.beta=dbeta(MAF1,1,25), r.corr=r.corr, adjustment=T))
                    if (inherits(opt.asy1,"try-error") | inherits(opt.adj1,"try-error")) {
                         s = s - 1
                    } else {
                         results$p.values[(r-1)*S+s,8:(7+2*(1+length(r.corr)))]  <- c(opt.asy1$param$p.val.each, opt.asy1$p.value, opt.adj1$param$p.val.each, opt.adj1$p.value)
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
