#############
# LIBRARIES #
#############

library(SKAT)
library(geepack)

#############
# FUNCTIONS #
#############

logit <- function(x) log(x/(1-x))

expit <- function(x) exp(x)/(1 + exp(x))

# Determine sample size needed to ensure that a sufficient number of cases is available for sampling during simulation
calc.N <- function(kappa=NULL, alpha=NULL, n1=NULL) {
     a <- kappa
     b <- qnorm(alpha)*sqrt(kappa*(1-kappa))
     c <- -n1
     return( ceiling( (-b+sqrt(b^2-4*a*c))^2/4/a^2 ) )
}


IPWSKAT_Null_Model <- function(formula, data=NULL, out_type="C", weights.ipw=NULL, ccstatus=NULL, prev=NULL, n.resampling=0) {
     
     # find individuals with missing phenotype or covariates
     obj1 <- model.frame(formula, na.action = na.omit, data)          
     obj2 <- model.frame(formula, na.action = na.pass, data)
     n1   <- dim(obj1)[1]
     n2   <- dim(obj2)[1]
     if(n2 - n1 > 0){
          MSG <- sprintf("%d samples have either missing phenotype or missing covariates. They are excluded from the analysis!", n2 - n1)
          warning(MSG,call.=FALSE)
     }
     id_include <- SKAT:::SKAT_Null_Model_Get_Includes(obj1,obj2)     # individuals to include in subsequent analysis
     
     # define the IPW weights for each individual
     if (is.null(ccstatus)) {
          stop("Please provide the case-control status for each individual.")
     } else if (is.null(prev) & is.null(weights.ipw)) {  # weight every individual equally
          weights.ipw <- rep(1,n1)
     } else if (length(weights.ipw)!=n2) {
          if (length(ccstatus)==n2 & all(unique(ccstatus) %in% c(0,1)) & !is.null(prev)) {     # use case-control status and disease prevalence to calculate weights
               weights.ipw <- 2*ccstatus*prev + 2*(1-ccstatus)*(1-prev)
               weights.ipw <- weights.ipw[id_include]
          } else {
               stop("Provide either (1) weights for each individual or (2) the disease prevalence.")
          }
     }
     ccstatus = ccstatus[id_include]
     
     # fit the null model
     id  <- seq(n2)
     data2 <- data.frame(cbind(obj1,weights.ipw,id))
     mu  <- NULL
     res <- NULL
     s2  <- NULL
     if (out_type=="C") { # weighted linear regression
          mod  <- geeglm(formula, data=data2, weights=weights.ipw, id=id)
          mu   <- mod$fitted.values
          res  <- matrix(mod$resid, ncol=1)
          s2   <- summary(mod)$dispersion$Estimate
     } else { # weighted logistic regression
          mod  <- geeglm(formula, data=data2, weights=weights.ipw, id=id, family="binomial")
          mu   <- mod$fitted.values
          res  <- matrix(mod$y - mu, ncol=1)
     }
     X1  <- model.matrix(formula,data=data)
     
     # bootstrap resampling
     res.moments <- NULL
     if (n.resampling>0) {
          mu1 = mu[ccstatus==1]
          mu0 = mu[ccstatus==0]
          Y1  = mod$y[ccstatus==1]
          Y0  = mod$y[ccstatus==0]
          res.moments <- matrix(0,nrow=n2,ncol=n.resampling)
          if (out_type=="C") {
               hat.matrix <- X1 %*% solve(t(X1)%*%(weights.ipw*X1)) %*% t(weights.ipw*X1)
               res1 = res[ccstatus==1]
               res0 = res[ccstatus==0]
               size1 = sum(ccstatus)
               size0 = sum(1-ccstatus)
               for (r in 1:n.resampling) {
                    Y.boot <- rep(0,n1)
                    Y.boot[ccstatus==1] <- mu1 + sample(res1,size1)
                    Y.boot[ccstatus==0] <- mu0 + sample(res0,size0)
                    res.moments[,r] <- Y.boot - hat.matrix %*% Y.boot
               }
          } else {
               for (r in 1:n.resampling) {
                    Y.boot <- rep(0,n1)
                    Y.boot[ccstatus==1][resample_bernoulli(Y1,mu1)] <- 1
                    Y.boot[ccstatus==0][resample_bernoulli(Y0,mu0)] <- 1
#                     data3 <- data.frame(cbind(Y.boot,X1,weights.ipw,id))
#                     mod.boot <- geeglm(Y.boot~-1+X1, data=data3, weights=weights.ipw, id=id, family="binomial")
#                     res.moments[,r] <- Y.boot - mod.boot$fitted.values
                    res.moments[,r] <- Y.boot - mu
               }
          }
     }
     
     return(list(res=res, X1=X1, res.moments=res.moments, weights.ipw=weights.ipw, out_type=out_type, n.resampling=n.resampling, id_include=id_include, s2=s2, mu=mu))
     
}


resample_bernoulli <- function(Y,mu) {
     
     n = length(mu)
     n.case = sum(Y)
     res.out1<-rbinom(n,1,mu)
     res.out2<-rbinom(n,1,mu)
     
     id_case1<-which(res.out1==1)
     id_case2<-which(res.out2==1)
     
     id_c1<-intersect(id_case1,id_case2)
     id_c2<-union(setdiff(id_case1,id_case2),setdiff(id_case2,id_case1))
     if(n.case <= length(id_c1)){
          id_case<-sample(id_c1,n.case)
     }else if (n.case > length(id_c1) && n.case <= length(id_c1)+length(id_c2)){
          id_c3<-sample(id_c2,n.case - length(id_c1))
          id_case<-c(id_c1,id_c3)
     }else {					
          id_case3<-union(id_c1,id_c2)
          id_c4<-setdiff(1:n,id_case3)
          n.needed<-n.case - length(id_case3)
          
          id_c5<-sample(id_c4,n.needed,prob=mu[id_c4])
          id_case<-union(id_case3,id_c5)
     }
     
     return(id_case)
     
}


IPWSKAT <- function(Z, obj, weights.snp=NULL, weights.snp.beta=c(1,25), r.corr=0, is_check_genotype=TRUE, missing_cutoff=0.15, adjustment=TRUE) {
     
     # check genotype matrix Z
     n <- dim(Z)[1]
     m <- dim(Z)[2]
     out.z <- SKAT:::SKAT_MAIN_Check_Z(Z, n, obj$id_include, SetID=NULL, weights.snp, weights.snp.beta, "fixed", is_check_genotype, is_dosage=FALSE, missing_cutoff) 
     
     # no polymorphic SNP
     if(out.z$return ==1){
          out.z$param$n.marker <- m
          return(out.z)
     }
     
     # edit r.corr
     if(length(r.corr) > 1 && dim(out.z$Z.test)[2] <= 1){
          r.corr=0
     }
     r.corr[which(r.corr>0.999)]=0.999  # for computational purposes, set upper limit of r.corr to 0.999
     
     Z    <- t(t(out.z$Z.test)*out.z$weights) # Z = G %*% W.snp
     X    <- obj$X1
     if (obj$out_type == "C") {    # V = diag{a(psi)*v(mu)*[g'(mu)^2]}
          V.inv <- 1
     } else {
          V.inv <- c(obj$mu*(1-obj$mu))
     }
     D    <- c(obj$weights.ipw*abs(obj$res))
     DPZ  <- D*Z - D*X %*% solve(t(X)%*%(obj$weights.ipw*V.inv*X)) %*% t(obj$weights.ipw*V.inv*X) %*% Z
     
     out.Q<- IPWSKAT_Get_Q(obj$weights.ipw*Z, obj$res, r.corr, obj$res.moments, obj$out_type)
     Q.r  <- out.Q$Q.r
     Q.sim<- out.Q$Q.sim
     
     if (adjustment & obj$n.resampling>0) {   # small sample size adjustment
          out  <- IPWSKAT_Get_Pvalue(Q.r, DPZ, obj, r.corr, Q.sim)
     } else if (adjustment & obj$n.resampling==0) {
          stop("Please refit null model, specifying number of resamples.")   
     } else {
          out  <- SKAT:::SKAT_Optimal_Get_Pvalue(Q.r, DPZ, r.corr, method="davies")
     }
     
     minp <- min(out$p.val.each)
     rho_est <- r.corr[which(out$p.val.each==minp)]
     
     param <- list(p.val.each=out$p.val.each, q.val.each=Q.r, rho=r.corr, minp=minp, rho_est=rho_est)
     
     return(list(p.value=out$p.value, param=param))
     
}


IPWSKAT_Get_Q <- function(Z1, res, r.corr, res.moments=NULL, out_type) {
     
     # Z1 = W.ipw %*% G %*% W.snp
     n.r <- length(r.corr)
     p.m <- dim(Z1)[2]
     
     temp<- t(res) %*% Z1            # t(Y-mu) %*% W.ipw %*% G %*% W.snp
     Q.r <- matrix(0,ncol=n.r)
     for (i in 1:n.r) {
          r.temp    <- r.corr[i]
          Q.skat    <- (1-r.temp) * rowSums(temp^2)
          Q.burden  <- r.temp * p.m^2 * rowMeans(temp)^2 
          Q.r[i]    <- Q.skat + Q.burden
     }
     
     Q.sim <- NULL
     if (!is.null(res.moments)) {
          temp <- t(res.moments) %*% Z1
          n.moments <- dim(res.moments)[2]
          Q.sim <- matrix(0,nrow=n.moments,ncol=n.r)
          for (i in 1:n.r) {
               r.temp    <- r.corr[i]
               Q.skat    <- (1-r.temp) * rowSums(temp^2)
               Q.burden  <- r.temp * p.m^2 * rowMeans(temp)^2 
               Q.sim[,i] <- Q.skat + Q.burden
          }
     }
     
     return(list(Q.r=Q.r,Q.sim=Q.sim))
     
}


IPWSKAT_Get_Pvalue <- function(Q.r, Z1, obj, r.corr, Q.sim) {
     
     n.r <- length(r.corr)
     p.m <- dim(Z1)[2]
     
     Z2.r      <- list()
     for (i in 1:n.r) {
          r.temp    <- r.corr[i]
          R.M       <- diag(rep(1-r.temp,p.m)) + matrix(rep(r.temp,p.m*p.m),ncol=p.m)
          L         <- chol(R.M, pivot=TRUE)
          Z2.r[[i]] <- Z1 %*% t(L)
     }
     
     # get mixture parameters
     param.m <- IPWSKAT_Param(Z1, obj$weights.ipw, r.corr, obj$res.moments, obj$out_type, obj$mu)
     muQ  <- param.m$param$muQ
     varQ <- param.m$param$varQ + param.m$VarRemain
     df   <- param.m$param$df
     tau  <- param.m$tau
     
     # Q.r specific statistics
     Each_Info <- IPWSKAT_Each_Q(param.m, Q.r, r.corr, Z2.r, Q.sim, obj$out_type, obj$mu)
     pmin.q    <- Each_Info$pmin.q
     p.val.each<- Each_Info$pval
     
     # calculate SKAT-O p-value
     pval <- SKAT:::SKAT_Optimal_PValue_VarMatching(pmin.q,muQ,varQ,df,tau,r.corr)
     
     return(list(p.value=pval, p.val.each=p.val.each))
     
}


IPWSKAT_Param <- function(Z1, weights.ipw, r.corr, res.moments, out_type, mu) {
     
     n    <- dim(Z1)[1]
     p.m  <- dim(Z1)[2]
     r.n  <- length(r.corr)
     
     z_mean    <- rowMeans(Z1)
     Z_mean    <- matrix(rep(z_mean,p.m),ncol=p.m,byrow=FALSE)
     cof1 <- (t(z_mean) %*% Z1)[1,] / sum(z_mean^2)
     
     Z.item1   <- Z_mean %*% diag(cof1)     # MZ
     Z.item2   <- Z1 - Z.item1 # (I-M)Z
     
     Q.sim = NULL
     if(!is.null(res.moments)) {
          Q.temp = t(sign(res.moments))%*%Z.item2
          Q.sim  = rowSums(Q.temp^2)
     }
     
     re.param  <- IPWSKAT_Param2(Z.item2, Q.sim)
     
     # W3.3 Term : variance of remaining ...
     W3.3.item <- sum((t(Z.item1) %*% Z.item1) * (t(Z.item2) %*% Z.item2)) * 4
     
     # W3.1 Term : tau1 * chisq_1
     cof2 <- (t(z_mean) %*% Z1)[1,] / sum(z_mean^2)
     tau<-rep(0,r.n)
     for(i in 1:r.n){
          r.temp    <- r.corr[i]
          term1     <- p.m*r.temp + cof2^2 * (1-r.temp)
          tau[i]    <- sum(term1) *  sum(z_mean^2)
     }
     
     out<-list(param=re.param, VarRemain=W3.3.item, tau=tau)
     return(out)
     
}


IPWSKAT_Param2 <- function(A, Q.sim) {
     
     lambda    <- SKAT:::Get_Lambda_U_From_Z(A)$lambda
     muQ       <- sum(lambda)
     varQ.sim  <- var(Q.sim)
     df.sim    <- SKAT:::SKAT_Get_DF_Sim(Q.sim)
     
     # No adjustment
     c1 <- rep(0,4)
     for (i in 1:4) {
          c1[i] <- sum(lambda^i)
     }
     param <- SKAT:::Get_Liu_Params_Mod(c1)
     
     return(list(muQ=muQ, varQ=varQ.sim, df=df.sim, lambda=lambda, param.noadj=param))
     
}


IPWSKAT_Each_Q <- function(param.m, Q.r, r.corr, Z2.r, Q.sim, out_type, mu) {
     
     n.r       <- length(r.corr)
     param.r   <- list()
     pval      <- rep(0,n.r)
     pval.noadj<- rep(0,n.r) 
     pmin.q    <- rep(0,n.r)
     
     # Calculate adjusted and non-adjusted p-values for each Q.r
     for (i in 1:n.r) {
          param          <- IPWSKAT_Param2(Z2.r[[i]], Q.sim[,i])
          param.r[[i]]   <- param
          
          # adjusted p-value
          r.temp    <- r.corr[i]
          muQ       <- param$muQ
          varQ1     <- param$varQ
          varQ      <- (1-r.temp)^2*(param.m$param$varQ+param.m$VarRemain) + param.m$tau[i]^2*2
          df        <- param$df
          Q.Norm    <- (Q.r[i]-muQ) / sqrt(varQ) * sqrt(2*df) + df
          pval[i]   <- 1-pchisq(Q.Norm, df=df, ncp=0)
          
          # non-adjusted p-value
          param.noadj    <- param$param.noadj
          Q.Norm         <- (Q.r[i]-param.noadj$muQ)/param.noadj$sigmaQ
          Q.Norm1        <- Q.Norm * param.noadj$sigmaX + param.noadj$muX
          pval.noadj[i]  <- 1-pchisq(Q.Norm1, df=param.noadj$l, ncp=0)
     }
     
     pmin <- min(pval)
     
     for (i in 1:n.r) {
          r.temp <- r.corr[i]
          muQ  <- param.r[[i]]$muQ
          varQ1<- param.r[[i]]$varQ
          varQ <- (1-r.temp)^2*(param.m$param$varQ+param.m$VarRemain) + param.m$tau[i]^2*2
          df   <- param.r[[i]]$df
          q.org<- qchisq(1-pmin, df=df)
          q.q  <- (q.org-df)/sqrt(2*df)*sqrt(varQ) + muQ
          pmin.q[i] <- q.q
     }
     
     return(list(pmin=pmin, pval=pval, pmin.q=pmin.q))
     
}
