############################################################
##
## This script implements a fully Bayesian version
## of GAM model selection under high-dimensional ST dependence.
##
## Template for Hipergator.
## AA - seed (100 replications for each)
## BB - sample size (250, 500, 1000)
## CC - SNR (1/4, 1, 4), takes values 1, 2, 3
##
## Hunter R. Merrill, 2017.03.30
##
## CORRELATED SIMULATION STUDY
##
#############################################################

rm(list=ls())
# AA = 1 #for HiperGator (AA = seed)
# BB = 1 #BB = sample size
# CC = 3 #CC = SNR
set.seed(AA)

# setwd("C:/Users/Hunter/Dropbox/Hunter/Research/BGL_Selection_Study/The_Paper_revision_2/")
setwd("/ufrc/bliznyuk/hmerrill/Bayes_GAM_Select")

library(mgcv) #for gam
library(Rcpp) #Rcpp
library(RcppArmadillo) #better matrix algebra
library(spikeSlabGAM)

##### Simulation parameters
n.true = 24 #number of true predictors
n.sim = 250*(2^(BB-1)) #number of observations
n.predictors = 100 #number of total preds
num.knots.sim = 10 #number of knots for simulation
logical.interactions = FALSE #should interactions be considered?
rho.par = BB #AR1 parameter
SNR = CC #overall signal to noise ratio, Gaussian

# if(rho.par == 1){
#   rho = 0.25
# } else if(rho.par == 2){
#   rho = .5
# } else if(rho.par == 3){
#   rho = .75
# } else if(rho.par == 4){
#   rho = 0
# }

n.burnin = 25000
n.samples = 5000 + n.burnin
thin = 10
# n.burnin = 1000; n.samples = 2000; thin = 1
keepers = seq(from=n.burnin+1, to= n.samples, by=thin) #burnin and thin the chain

# 
# if(n.sim == 250){
#   n.burnin = 50000
#   n.samples = 50000 + n.burnin
#   thin = 10
#   keepers = seq(from=n.burnin+1, to= n.samples, by=thin) #burnin and thin the chain
# } else {
#   n.burnin = 25000
#   n.samples = 25000 + n.burnin
#   thin = 10
#   keepers = seq(from=n.burnin+1, to= n.samples, by=thin) #burnin and thin the chain
# }

##################
### Data Setup ###
##################

# sourceCpp("20160901_BGL_SS_MCMC_1lam.cpp") #get AR(1) version
# sourceCpp("M_BGL_SS_AR1_MCMC_wip_norho.cpp") #get AR(1) version
# sourceCpp("M_BGL_SS_MCMC_wip.cpp")
sourceCpp("SRC/BGL_GAM_MCMC_ST.cpp")
#source("C:/Users/Hunter/Dropbox/Hunter/Research/My_R_Packages/hunterr/R/Basis_expansion_function.R") #basis expansion
source("SRC/Basis_expansion_function_new.R") #basis expansion
# source("Basis_expansion_function_2.R") #basis expansion

#Simulate data
BCO.dat = data.frame(rep(1,n.sim))
names(BCO.dat) = "Int"

covars = matrix(runif(n.predictors*n.sim), ncol=n.predictors, nrow=n.sim)

for(i in 1:n.predictors){
  BCO.dat[[i+1]] = covars[,i]
  names(BCO.dat)[i+1] = paste("X",i,sep="")
}

BCO.dat[,2] = seq(0,1,length=nrow(BCO.dat))

############################################################################
### get basis expansion matrix
num.knots = rep(num.knots.sim,n.predictors)
bases = rep("cr",n.predictors)
basis.expansion = get_basis_expansion_new(df=BCO.dat,
                                          num.knots = num.knots,
                                          bases = bases,
                                          interactions=logical.interactions)

C.full = basis.expansion$C.full
partition = basis.expansion$partition
nvar = basis.expansion$nvar
pc = basis.expansion$pc
pz = basis.expansion$pz
px = basis.expansion$px

#cbind(matrix(colnames(C.full)),partition)

############################################################################
### create trends

# R.list = list() #QR decomposition
# w.list = list() #coefficients
# for(i in 1:max(unique(partition))){
#   R.list[[i]] = qr.R(qr(C.full[,partition == i]))
#   w.list[[i]] = solve(R.list[[i]]) %*% rep(7, dim(R.list[[i]])[1])
# }

trends = matrix(numeric(n.sim*n.true),nrow=n.sim,ncol=n.true)
# for(i in 1:n.true){ #only the first n.true predictors count.
#   trends[,i] = C.full[,partition == i] %*% w.list[[i]]
# }
trends[,1] = 2*sin(pi*BCO.dat[,2])
trends[,2] = exp(BCO.dat[,3]^2)
trends[,3] = -BCO.dat[,4]
trends[,4] = BCO.dat[,5]^11*(10*(1-BCO.dat[,5]))^6 + 10*(10*BCO.dat[,5])^3*(1-BCO.dat[,5])^10
trends[,5] = (BCO.dat[,6]^3 + sin(pi*BCO.dat[,6]^3))
trends[,6] = cos(2*pi*BCO.dat[,7]) + sin(pi*BCO.dat[,7])

for(i in 1:6){
  trends[,i] = (trends[,i] - min(trends[,i]))/diff(range(trends[,i]))
}

trends[,7:12] <- trends[,1:6]/2
trends[,13:18] <- trends[,1:6]/4
trends[,19:24] <- trends[,1:6]/10

# par(mfrow=c(2,3), mar=c(3,3,1,1), mgp=c(2,1,0))
# for(i in 1:n.true){
#   plot(BCO.dat[,i+1], trends[,i], ylim=range(trends))
# }

mean.trend = apply(trends, 1, sum)
signal.var = var(mean.trend)

if(SNR == 3){
  sd.noise = sqrt(signal.var)/(.7/(1-.7)) #SNR = 3
} else if(SNR == 2){
  sd.noise = sqrt(signal.var)/(.55/(1-.55))   #SNR = 1
} else if(SNR == 1){
  sd.noise = sqrt(signal.var)/(.4/(1-.4)) #SNR = 1/3
}

# if(SNR == 3){
#   sd.noise = .1 #SNR = 4
# } else if(SNR == 2){
#   sd.noise = .2   #SNR = 1
# } else if(SNR == 1){
#   sd.noise = .4 #SNR = 1/4
# }

# x.s = BCO.dat$X2
# z.s = BCO.dat$X3
x.s = rep(sort(runif(n.sim/10)), 10)
z.s = rep(runif(n.sim/10), 10)

Dist = as.matrix(dist(cbind(x.s, z.s)))
# cstr <- corExp(.03,form = ~x.s+z.s)  
# cstr <- Initialize(cstr,data.frame(x=z.s,z=z.s))
# V <- corMatrix(cstr) ## correlation matrix for data
V <- exp(-Dist/0.1)*exp(as.matrix(dist(rep(1:(10), each = n.sim/10)))*log(.75))
Cv <- chol(V)*sd.noise

# plot(BCO.dat$X2, BCO.dat$X3, pch = 16,
#      col = fields::tim.colors(n.sim)[rank(trends[,2] + trends[,3])])
# plot(x.s, z.s, pch = 16,
#      col = fields::tim.colors(n.sim)[rank(t(Cv)%*%rnorm(n.sim))])
# plot(BCO.dat$X2, BCO.dat$X3, pch = 16,
#      col = fields::tim.colors(n.sim)[rank(t(Cv)%*%rnorm(n.sim) + trends[,2] + trends[,3])])

# pres.mat = matrix(0, nrow=n.sim, ncol=n.sim)
# diag(pres.mat) = c(1, rep(1+rho^2, n.sim-2) ,1)
# for(i in 1:(n.sim-1)){
#   pres.mat[i,i+1] = -rho
#   pres.mat[i+1,i] = -rho
# }
# cor.mat = solve(pres.mat)
# cov.mat = (sd.noise^2)*cor.mat

# BCO.dat$Y.norm = 0
# BCO.dat$Y.norm[1] = mean.trend[1] + rnorm(1, 0, sd=sd.noise)
# for(i in 2:n.sim){
#   BCO.dat$Y.norm[i] = mean.trend[i] + rho*(BCO.dat$Y.norm[i-1] - mean.trend[i-1]) + rnorm(1, 0, sd=sd.noise)
# }

BCO.dat$Y.norm = c(mean.trend + t(Cv)%*%rnorm(n.sim))
# BCO.dat$Y.norm = c(mean.trend + sd.noise*rnorm(n.sim))
############################################################################



################
#### MCMCMC ####
################

##test: Gaussian
Y = BCO.dat$Y.norm
sig2.alpha = .1
sig2.gamma = .1
pi0.a = 1
pi0.b = 1
lamsep = FALSE
rhomethod = "empirical" #depricated

MCMC_Results_Rcpp = M_BGL_SS_AR1_MCMC(C_full = C.full, #design matrix
                                      Y = Y, #data
                                      partition = partition, #groups definition 
                                      sig2_alpha = sig2.alpha, #hyperparameters for sig^2
                                      sig2_gamma = sig2.gamma,
                                      pi0_a = pi0.a, #hyperparameters for pi_0
                                      pi0_b = pi0.b,
                                      n_samples = n.samples, #how many iterations
                                      sig2_beta = 1e6,
                                      print_step = 100,
                                      tune_rho = 0.1,
                                      DistS = Dist,
                                      DistT = as.matrix(dist(rep(1:(n.sim/10), each = 10))),
                                      rhomethod = rhomethod,
                                      lamsep = lamsep) 

## results
Z.samples.C = MCMC_Results_Rcpp$Z.samples
Zeta.samples.C = MCMC_Results_Rcpp$Zeta.samples
w.samples.C = MCMC_Results_Rcpp$w.samples
sig2.samples.C = MCMC_Results_Rcpp$sig2.samples
tau2.samples.C = MCMC_Results_Rcpp$tau2.samples
pi0.samples.C = MCMC_Results_Rcpp$pi0.samples
# lam.samples.C = MCMC_Results_Rcpp$lam.samples
# rho.samples.C = MCMC_Results_Rcpp$rho.samples

true.model.ind = rep(0, length(unique(partition)))
true.model.ind[(1:n.true+1)] = 1
true.model.ind = true.model.ind[-1] #don't care about intercept.

Selection.Matrix = ifelse(Z.samples.C + Zeta.samples.C > 0, 1L, 0)
est.model.Gaus = apply(Selection.Matrix[,keepers], 1, mean)

HPPM <- function(x) {
  ux <- unique(x, MARGIN=2)#find unique columns in matrix of indicator variables
  nmatches = numeric(ncol(ux))
  for(i in 1:ncol(ux)){
    nmatches[i] = sum(apply(x, 2, identical, ux[,i]))
  }
  ux[,which.max(nmatches)]
}

HPPM.Est = HPPM(Selection.Matrix[,keepers])
which.HPPM = which(apply(Selection.Matrix[,keepers], 2, identical, HPPM.Est))
# plot(rho.samples.C[keepers][which.HPPM], type="l")
MTM.Est = as.integer(est.model.Gaus > 0.5)
which.MTM = which(apply(Selection.Matrix[,keepers], 2, identical, as.numeric(MTM.Est)))

which.true = which(apply(Selection.Matrix[,keepers], 2, identical, c(0,true.model.ind)))
which.true2 = which(apply(Selection.Matrix[,keepers], 2, identical, c(1,true.model.ind)))
which.true = unique(c(which.true, which.true2))

# par(mfrow=c(3,1), mar=c(3,3,1,1), mgp=c(2,1,0))
# plot(rho.samples.C[keepers][which.HPPM], type="l")
# plot(rho.samples.C[keepers][which.MTM], type="l")
# plot(rho.samples.C[keepers][which.true], type="l")

w.est.Gaus = apply(w.samples.C[,keepers], 1, median)

# par(mfrow=c(1,1), mar=c(3,3,1,1), mgp = c(2,1,0))
# plot(BCO.dat$X1, C.full[,partition==1] %*% w.est.Gaus[partition==1],
#      type="l",
#      xlab="Time",
#      ylab="f(Time)",
#      ylim=range(c(trends[,1], C.full[,partition==1] %*% w.est.Gaus[partition==1])))
# lines(BCO.dat$X1, trends[,1], col="blue", lty=2)

count = 0
Dev = numeric(length(keepers))
for(i in keepers){
  count = count + 1
  Dev[count] = -2*sum(dnorm(BCO.dat$Y.norm, mean=C.full %*% w.samples.C[,i], sd=sqrt(sig2.samples.C[i]), log = TRUE))
}

Dbar.Gaus = mean(Dev)
Dhat.Gaus = -2*sum(dnorm(BCO.dat$Y.norm, mean = C.full%*%w.est.Gaus, sd=sqrt(mean(sig2.samples.C[keepers])), log = TRUE))

pD.Gaus = Dbar.Gaus - Dhat.Gaus
DIC.Gaus = Dhat.Gaus + 2*pD.Gaus
###################################################################################

### estimated linear predictors
est.gaus.mean = C.full %*% w.est.Gaus
MSE.gaus = sqrt(mean((est.gaus.mean - mean.trend)^2))

sel.model.Gaus = ifelse(est.model.Gaus < 0.5, 0, 1)[-1]
sel.lin.Gaus = ifelse(apply(Z.samples.C[,keepers], 1, mean) < 0.5, 0, 1)[-1]
sel.nlin.Gaus = ifelse(apply(Zeta.samples.C[,keepers], 1, mean) < 0.5, 0, 1)[-1]
# sel.model.Gaus = HPPM.Est[-1]

true.lin.ind = c( rep(c(0,1,1,1,1,0), 4), rep(0, 100 - 24))
true.nlin.ind = c(rep(c(1,1,0,1,1,1), 4), rep(0, 100 - 24))
FPl.Gaus = sum(true.lin.ind - sel.lin.Gaus == -1)
FNl.Gaus = sum(true.lin.ind - sel.lin.Gaus == 1)
FPn.Gaus = sum(true.nlin.ind - sel.nlin.Gaus == -1)
FNn.Gaus = sum(true.nlin.ind - sel.nlin.Gaus == 1)
FP.Gaus = sum(true.model.ind - sel.model.Gaus == -1)
FN.Gaus = sum(true.model.ind - sel.model.Gaus == 1)

# ############################################################
# ### for ROC curve
# 
# prob.vals = seq(-0.001,1.001,len=100)
# TPR = numeric(length(prob.vals))
# TNR = numeric(length(prob.vals))
# for(i in 1:length(TPR)){
#   this.model = ifelse(est.model.Gaus < prob.vals[i], 0, 1)[-1]
#   TPR[i] = 1 - sum(true.model.ind - this.model == -1)/6
#   TNR[i] = 1 - sum(true.model.ind - this.model == 1)/6
# }
# 
# res.ROC = data.frame("TPR" = TPR,
#                      "TNR" = TNR,
#                      "seed" = AA)
# 
# M = "20161109_ROC_BB_CC.csv"
# does.M.exist = file.exists(M)
# if (!does.M.exist) file.create(M)
# write.table(res.ROC, 
#             file = M, 
#             sep=",",
#             row.names=FALSE,
#             append=TRUE,
#             col.names=!does.M.exist)

############################################################################
### spikeSlabGAM version
BCO.dat$x.s = x.s
BCO.dat$z.s = z.s
options(mc.cores = 1)
mcmc = list(nChains = 1, chainLength = n.samples-n.burnin, burnin = n.burnin, thin = thin, sampleY = FALSE)


m = spikeSlabGAM(formula = as.formula(paste("Y.norm ~", 
                                            paste(names(BCO.dat)[grep("X", names(BCO.dat))], 
                                                  collapse = " + "),
                                            "+ x.s + z.s + X1:x.s:z.s")),
                 data = BCO.dat,
                 mcmc = mcmc)

which.selected = function(x){
  marg.p = summary(x)$trmSummary[-1,1][-(201:212)]
  marg.p.mat = matrix(marg.p, nrow = 2)
  did.select = ifelse(marg.p.mat < 0.5, 0, 1)
  did.sel.term = apply(did.select, 2, sum)
  ifelse(did.sel.term > 0, 1, 0)
}

mny = predict(m, type="link")

sel.model.SSG = which.selected(m)
FP.SSG = sum(true.model.ind - sel.model.SSG == -1)
FN.SSG = sum(true.model.ind - sel.model.SSG == 1)

sel.lin.SSG = as.numeric(summary(m)$trmSummary[-1,1][-(201:212)][grep('lin',names(summary(m)$trmSummary[-1,1][-(201:212)]))] > .5)
sel.nlin.SSG = as.numeric(summary(m)$trmSummary[-1,1][-(201:212)][grep('sm',names(summary(m)$trmSummary[-1,1][-(201:212)]))] > .5)
FPl.SSG = sum(true.lin.ind - sel.lin.SSG == -1)
FNl.SSG = sum(true.lin.ind - sel.lin.SSG == 1)
FPn.SSG = sum(true.nlin.ind - sel.nlin.SSG == -1)
FNn.SSG = sum(true.nlin.ind - sel.nlin.SSG == 1)

Dbar.SSG = mean(-2*m$samples$logLik[[1]])
Dhat.SSG = (-2*sum(dnorm(BCO.dat$Y.norm, mean=mny, sd=sqrt(m$postMeans$sigma2), log=TRUE)))
pD.SSG = Dbar.SSG - Dhat.SSG

DIC.SSG = Dhat.SSG + 2*pD.SSG

sel.terms = ifelse(summary(m)$trmSummary[-1,1] < 0.5, 0, 1)

mny.terms = predict(m, type="terms")[-(201:212)]
for(i in 1:(length(mny.terms)-1)){
  mny.terms[[i]] = mny.terms[[i]]*sel.terms[i]
}
mny.sel = Reduce("+", mny.terms)

# MSE.SSG = sqrt(mean((mny - mean.trend)^2))
MSE.SSG = sqrt(mean((mny.sel - mean.trend)^2))

# # ############################################################################
# ### Marra version: Gaussian
# ctrl = lmeControl(opt='optim')
# fit.gam.dp = gamm(Y.norm ~ ti(X1, k=num.knots.sim) + ti(X2, k=num.knots.sim) + ti(X3, k=num.knots.sim) + ti(X4, k=num.knots.sim) +
#                    ti(X5, k=num.knots.sim) + ti(X6, k=num.knots.sim) + ti(X7, k=num.knots.sim) + ti(X8, k=num.knots.sim) +
#                    ti(X9, k=num.knots.sim) + ti(X10, k=num.knots.sim) + ti(X11, k=num.knots.sim) + ti(X12, k=num.knots.sim) +
#                     te(x.s,z.s,X1,
#                        k=c(25,10),
#                        d=c(2,1),
#                        bs=c("tp","cr")),
#                  dat = BCO.dat, select = TRUE, method="REML", 
#                  correlation = corAR1(form = ~ 1),
#                  control = ctrl)
# 
# sel.model.Gaus.dp = ifelse(summary(fit.gam.dp$gam)$s.pv < 0.05, 1, 0)[-13]
# FP.Gaus.dp = sum(true.model.ind - sel.model.Gaus.dp == -1)
# FN.Gaus.dp = sum(true.model.ind - sel.model.Gaus.dp == 1)
# 
# est.terms = predict(fit.gam.dp$gam, type = "terms")
# est.gaus.mean.dp = est.terms[,1:12]%*%sel.model.Gaus.dp + attr(est.terms, "constant")
# MSE.gaus.dp = sqrt(mean((est.gaus.mean.dp - mean.trend)^2))
# 
# fit.gam.sh = gamm(Y.norm ~ ti(X1, k=num.knots.sim, bs="cs") + ti(X2, k=num.knots.sim, bs="cs") + ti(X3, k=num.knots.sim, bs="cs") + ti(X4, k=num.knots.sim, bs="cs") +
#                    ti(X5, k=num.knots.sim, bs="cs") + ti(X6, k=num.knots.sim, bs="cs") + ti(X7, k=num.knots.sim, bs="cs") + ti(X8, k=num.knots.sim, bs="cs") +
#                    ti(X9, k=num.knots.sim, bs="cs") + ti(X10, k=num.knots.sim, bs="cs") + ti(X11, k=num.knots.sim, bs="cs") + ti(X12, k=num.knots.sim, bs="cs") +
#                     te(x.s,z.s,X1,
#                        k=c(25,10),
#                        d=c(2,1),
#                        bs=c("tp","cr")),
#                  dat = BCO.dat, select = FALSE, method="REML",
#                  correlation = corAR1(form = ~ 1),
#                  control = ctrl)
# 
# sel.model.Gaus.sh = ifelse(summary(fit.gam.sh$gam)$s.pv < 0.05, 1, 0)[-13]
# FP.Gaus.sh = sum(true.model.ind - sel.model.Gaus.sh == -1)
# FN.Gaus.sh = sum(true.model.ind - sel.model.Gaus.sh == 1)
# 
# est.terms.sh = predict(fit.gam.sh$gam, type = "terms")
# est.gaus.mean.sh = est.terms.sh[,1:12]%*%sel.model.Gaus.sh + attr(est.terms.sh, "constant")
# MSE.gaus.sh = sqrt(mean((est.gaus.mean.sh - mean.trend)^2))
# 
# #############################################################


#############################################################

results_AA_BB_CC = data.frame("MSE" = c(MSE.gaus,
                                        # MSE.gaus.sh,
                                        # MSE.gaus.dp,
                                        MSE.SSG),
                              "FP" = c(FP.Gaus,
                                       # FP.Gaus.sh,
                                       # FP.Gaus.dp,
                                       FP.SSG),
                              "FN" = c(FN.Gaus,
                                       # FN.Gaus.sh,
                                       # FN.Gaus.dp,
                                       FN.SSG),
                              "FPl" = c(FPl.Gaus,
                                       # FP.Gaus.sh,
                                       # FP.Gaus.dp,
                                       FPl.SSG),
                              "FNl" = c(FNl.Gaus,
                                       # FN.Gaus.sh,
                                       # FN.Gaus.dp,
                                       FNl.SSG),
                              "FPn" = c(FPn.Gaus,
                                       # FP.Gaus.sh,
                                       # FP.Gaus.dp,
                                       FPn.SSG),
                              "FNn" = c(FNn.Gaus,
                                       # FN.Gaus.sh,
                                       # FN.Gaus.dp,
                                       FNn.SSG),
                              "DIC" = c(DIC.Gaus,
                                        # NA,
                                        # NA,
                                        DIC.SSG),
                              "Method" = rep(c("Bayesian","SSGAM"),1),
                              "SNR" = rep(as.character(SNR), 2),
                              "SampleSize250*2^i-1" = rep(as.character(BB), 2))
results_AA_BB_CC

M = "20171008_results_BGLGAM_STHD.csv"
does.M.exist = file.exists(M)
if (!does.M.exist) file.create(M)
write.table(results_AA_BB_CC, 
            file = M, 
            sep=",",
            row.names=FALSE,
            append=TRUE,
            col.names=!does.M.exist)






#print(object.size(x=lapply(ls(), get)), units="Mb") #for reference: Hipergator
