
rm(list=ls())
library(hisse)
library(tidyverse)
library(Metrics)
library(utils)
library(parallel)

test_fn <- function(rep_index, true.lambda, true.mu, true.psi) {
    
    #Single rate
    all.dead <- TRUE
    while(all.dead){
        try(sim.tab <- hisse::SimulateHisse(turnover=rep(true.lambda+true.mu,2), eps=rep(true.mu/true.lambda,2), max.taxa=200, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
        if(sim.tab$n.surviving > 0){
            try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
            if(Ntip(phy)>200){
                f <- hisse:::GetFossils(phy, psi=true.psi)
                if(sum(f$has_sampled_descendant) > 0){
                    try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
                    if(!is.null(test.tol)){
                        if(Ntip(test.tol) == 200){
                            pp <- hisse:::ProcessSimSample(phy, f)
                            pp.strat <- hisse:::ProcessSimStrat(phy, f)
                            all.dead <- FALSE
                        }
                    }
                }
            }
        }
    }
    
    misse.mk.full <- MiSSE(pp$phy, turnover=1, eps=1, includes.fossils=TRUE, k.sample=pp$k.samples, sann=FALSE)
    k.iffy <- pp$k.samples[sample.int(n=nrow(pp$k.samples), size=max(1, round(nrow(pp$k.samples)/2))),]
    misse.mk.missinghalf <- MiSSE(pp$phy, turnover=1, eps=1, includes.fossils=TRUE, k.sample=k.iffy, sann=FALSE)
    misse.zero_k <- MiSSE(pp$phy, turnover=1, eps=1, includes.fossils=TRUE, sann=FALSE)
    phy.extant <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50)
    misse.strat <- MiSSE(pp.strat$phy, turnover=1, eps=1, includes.fossils=TRUE, strat.intervals=pp.strat$strat.intervals, sann=FALSE)
    misse.extant <- MiSSE(phy.extant, turnover=1, eps=1, includes.fossils=FALSE, sann=FALSE)
    m = ape::Ntip(phy) - ape::Ntip(phy.extant)
    m.used = ape::Ntip(pp$phy) - ape::Ntip(phy.extant)
    k = nrow(pp$k.samples)
    n.extant <- ape::Ntip(phy.extant)
    true.turnover <- true.lambda + true.mu
    true.eps <- true.mu/true.lambda
    
    local_results <- rbind(
    data.frame(t(c(misse="extant",  rep=rep_index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi, m=m, m.used=NA, k=k, k.used=NA, n.extant=n.extant, true.turnover=true.turnover, true.eps=true.eps, misse.extant$solution[c("turnover0A", "eps0A", "psi")])), stringsAsFactors=FALSE),
    data.frame(t(c(misse="full_k",  rep=rep_index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=k, n.extant=n.extant, true.turnover=true.turnover, true.eps=true.eps, misse.mk.full$solution[c("turnover0A", "eps0A", "psi")])), stringsAsFactors=FALSE),
    data.frame(t(c(misse="half_k",  rep=rep_index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=nrow(k.iffy), n.extant=n.extant, true.turnover=true.turnover, true.eps=true.eps, misse.mk.missinghalf$solution[c("turnover0A", "eps0A", "psi")])), stringsAsFactors=FALSE),
    data.frame(t(c(misse="zero_k",  rep=rep_index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=0, n.extant=n.extant, true.turnover=true.turnover, true.eps=true.eps,misse.zero_k$solution[c("turnover0A", "eps0A", "psi")])), stringsAsFactors=FALSE),
    data.frame(t(c(misse="strat_int",  rep=rep_index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=dim(pp.strat$strat.intervals)[1], n.extant=n.extant, true.turnover=true.turnover, true.eps=true.eps,misse.strat$solution[c("turnover0A", "eps0A", "psi")])), stringsAsFactors=FALSE),
    stringsAsFactors=FALSE)

    session.stuff <- sessionInfo()
    save(session.stuff, phy, f, pp, local_results, file=paste("scenario", true.lambda, true.mu, true.psi, rep_index, "Rsave", sep="."))
    
    return(local_results)
}



DoSim <- function(rep.index, true.lambda, true.mu, true.psi){
    failed <- TRUE
    while(failed == TRUE){
        tmp <- NA
        try(tmp <- test_fn(rep_index=rep.index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi))
        print("you have done it!")
        if(length(tmp)>1){
            failed <- FALSE
        }else{
            print("oops!")
        }
    }
}


set1 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.3, true.psi=0.005, mc.cores=40)
set2 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.2, true.psi=0.005, mc.cores=40)
set3 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.1, true.psi=0.005, mc.cores=40)

set4 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.3, true.psi=0.01, mc.cores=40)
set5 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.2, true.psi=0.01, mc.cores=40)
set6 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.1, true.psi=0.01, mc.cores=40)

set7 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.3, true.psi=0.05, mc.cores=40)
set8 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.2, true.psi=0.05, mc.cores=40)
set9 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.1, true.psi=0.05, mc.cores=40)

set10 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.3, true.psi=0.1, mc.cores=40)
set11 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.2, true.psi=0.1, mc.cores=40)
set12 <- mclapply(1:100, DoSim, true.lambda=0.4, true.mu=0.1, true.psi=0.1, mc.cores=40)



#results <- data.frame()
#true.lambda <- 0.4
#true.mu.vector <- c(0.3, 0.2, 0.1)
#true.psi.vector <- c(0.005, 0.01, 0.05, 0.1)
#for (mu_index in seq_along(true.mu.vector)) {
#    for (psi_index in seq_along(true.psi.vector)) {
#        true.mu <- true.mu.vector[mu_index]
#        true.psi <- true.psi.vector[psi_index]
#        for (rep_index in sequence(50)) {
#            true.turnover = true.lambda+true.mu
#            true.eps=true.mu/true.lambda
#            local_results <- data.frame()
#            try(local_results <- test_fn(true.lambda, true.mu, true.psi))
#            if(nrow(local_results)>0) {
#                try(results <- rbind(results, local_results, stringsAsFactors=FALSE))
#                write.csv(results, file="results.csv", quote=FALSE, row.names=FALSE)
#            }
#        }
#    }
#}

#save(results, file="OutputBeforeSummarizing.Rsave")

#load("OutputBeforeSummarizing.Rsave")
#results <- read.csv("/Users/jmbeauli/Desktop/MiSSESims/results.csv")
#for (i in 2:ncol(results)) {
#    results[,i] <- as.numeric(results[,i])
#}
#results$error.eps <- as.numeric(results$eps0A)-as.numeric(results$true.eps)
#results$error.turnover <- as.numeric(results$turnover0A)-as.numeric(results$true.turnover)
#results$m.prop <- as.numeric(results$m.used) / as.numeric(results$m)
#results$k.prop <- as.numeric(results$k.used) / as.numeric(results$k)
#results$misse <- as.factor(results$misse)
#write.csv(results, file="results.csv")

#results_summary <- results %>% group_by(true.mu, true.psi, misse) %>% summarise(median_eps_error=median(error.eps), rmse_eps=rmse(true.eps, eps0A), median_turnover_error=median(error.turnover), rmse_turnover=rmse(true.turnover, turnover0A), m_prop_median=median(m.prop), m_median=median(m.used), k_prop_median=median(k.prop), k_median=median(k.used))
#write.csv(results_summary, file="results_summary.csv")
#results_highlevel <- results %>% group_by(misse) %>% summarise(median_eps_error=median(error.eps), rmse_eps=rmse(true.eps, eps0A), median_turnover_error=median(error.turnover), rmse_turnover=rmse(true.turnover, turnover0A), m_prop_median=median(m.prop), m_median=median(m.used),  k_prop_median=median(k.prop), k_median=median(k.used))
#write.csv(results_highlevel, file="results_highlevel.csv")









