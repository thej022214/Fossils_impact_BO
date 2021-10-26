
rm(list=ls())

library(diversitree)
library(hisse)
library(parallel)
library(TreeSim)
library(utils)

#Need this function because A and B are arbitrary
fixParOrder <- function(misse.fit, to.order="turn"){
    if(to.order == "turn"){
        if(misse.fit$solution["turnover0A"] > misse.fit$solution["turnover0B"]){
            tmp <- misse.fit$solution[1:4]
            misse.fit$solution["turnover0A"] <- tmp["turnover0B"]
            misse.fit$solution["turnover0B"] <- tmp["turnover0A"]
            misse.fit$solution["eps0A"] <- tmp["eps0B"]
            misse.fit$solution["eps0B"] <- tmp["eps0A"]
        }
        return(misse.fit)
    }
    if(to.order == "eps"){
        if(misse.fit$solution["eps0A"] > misse.fit$solution["eps0B"]){
            tmp <- misse.fit$solution[1:4]
            misse.fit$solution["turnover0A"] <- tmp["turnover0B"]
            misse.fit$solution["turnover0B"] <- tmp["turnover0A"]
            misse.fit$solution["eps0A"] <- tmp["eps0B"]
            misse.fit$solution["eps0B"] <- tmp["eps0A"]
        }
        return(misse.fit)
    }
}


GetSimModelMat <- function(model.list){
    model.pars.mat <- c()
    for(model.index in 1:length(model.list)){
        if(model.index == 1){
            model.pars.mat <- rbind(model.pars.mat, c(hisse:::ParameterTransformMiSSESpecial(model.list[[1]]$solution)[1:5], hisse:::ParameterTransformMiSSESpecial(model.list[[1]]$solution)[1:5], model.list[[1]]$solution["psi"]))
        }else{
            model.pars.mat <- rbind(model.pars.mat, c(hisse:::ParameterTransformMiSSESpecial(model.list[[model.index]]$solution)[1:10], model.list[[model.index]]$solution["psi"]))
        }
    }
    return(model.pars.mat)
}


DoReconOnSet <- function(phy, model.list, strat.intervals=FALSE){
    recon.list <- as.list(1:length(model.list))
    for(model.index in 1:length(model.list)){
        if(strat.intervals == TRUE){
            recon.list[[model.index]] <- MarginReconMiSSE(phy=phy, f=model.list[[model.index]]$f, pars=model.list[[model.index]]$solution, hidden.states=model.list[[model.index]]$hidden.states, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils=model.list[[model.index]]$includes.fossils, strat.intervals=model.list[[model.index]]$strat.intervals, AIC=model.list[[model.index]]$AIC, get.tips.only=TRUE)
        }else{
            recon.list[[model.index]] <- MarginReconMiSSE(phy=phy, f=model.list[[model.index]]$f, pars=model.list[[model.index]]$solution, hidden.states=model.list[[model.index]]$hidden.states, condition.on.survival=TRUE, root.type="madfitz", root.p=NULL, includes.fossils=model.list[[model.index]]$includes.fossils, k.samples=model.list[[model.index]]$k.samples, AIC=model.list[[model.index]]$AIC, get.tips.only=TRUE)
        }
    }
    return(recon.list)
}


GetSimTree <- function(nreps, ntax, scenario){
    trees <- as.list(1:nreps)
    for(nrep.index in 1:nreps){
        print(nrep.index)
        if(scenario == "scenario1"){
            #Net div and turnover is higher by ~2x in rate class 2
            pars <- c(0.25,0.5, 0.75, 0.75)
            all.dead <- TRUE
            phy <- NULL
            while(all.dead){
                try(sim.tab <- hisse::SimulateHisse(turnover=c(0.25,0.5), eps=rep(0.75,2), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
                if(sim.tab$n.surviving > 0){
                    try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
                    if(Ntip(phy) > ntax){
                        try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
                        if(!is.null(test.tol)){
                            if(Ntip(test.tol) == ntax){
                                all.dead <- FALSE
                            }
                        }
                    }
                }
            }
        }
        
        if(scenario == "scenario2"){
            #Net div the same, but turnover varies ~2x, and eps varies
            pars <- c(0.25, 0.50, 0.55, 0.75)
            all.dead <- TRUE
            phy <- NULL
            while(all.dead){
                try(sim.tab <- hisse::SimulateHisse(turnover=c(0.25,0.50), eps=c(0.55,0.75), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
                if(sim.tab$n.surviving > 0){
                    try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
                    if(Ntip(phy) > ntax){
                        try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
                        if(!is.null(test.tol)){
                            if(Ntip(test.tol) == ntax){
                                all.dead <- FALSE
                            }
                        }
                    }
                }
            }
        }
        
        if(scenario == "scenario3"){
            #Single rate
            pars <- c(0.375, 0.375, 0.75, 0.75)
            all.dead <- TRUE
            while(all.dead){
                try(sim.tab <- hisse::SimulateHisse(turnover=c(0.30,0.30), eps=c(0.75,0.75), max.t=50, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
                if(sim.tab$n.surviving > 0){
                    try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
                    print(phy)
                    #if(Ntip(phy) > ntax){
                        test.tol <- NULL
                        try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
                        if(!is.null(test.tol)){
                            #if(Ntip(test.tol) == ntax){
                                all.dead <- FALSE
                            #}
                        }
                    #}
                }
            }
        }
        
        if(scenario == "scenario4") {
            #Realistic with mass extinction and heritable sampling
            pars <- c(0.5, 0.2, 0.7, 0.9)
            turnover <- c(0.3, 0.6, 1.5, 1.0)
            eps <- c(0.7, 0.9, 0.95, 0.8)
            transition <- 0.1*matrix(c(NA, 0.5, 0.2, 0.6, 0.1, NA, 0.3, 0.2, 0.1, 0.4, NA, 0.8, 0.1, 0.2, 0.3, NA),nrow=4)
            mass.extinction.heights <- c(40, 70)
            mass.extinction.magnitudes <- c(0.7, 0.9)
            
            sim.tab <- list()
            sim.tab$n.surviving = 0
            while(sim.tab$n.surviving==0) {
                sim.tab <- hisse::SimulateHisse(turnover.rates=turnover, eps.values=eps, transition.rates=transition, max.taxa=ntax/.75, x0=0, mass.extinction.heights=mass.extinction.heights, mass.extinction.magnitudes=mass.extinction.magnitudes, nstart=2)
            }
            phy <- hisse::SimToPhylo(sim.tab, include.extinct=TRUE)
            phy.full <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50)
            # now biased samples down to ntax extant taxa
            sample.prob.raw <- abs(geiger::sim.char(phy.full, par=1, model="BM"))[,1,1]
            sample.prob.raw[which(is.na(sample.prob.raw))] <- median(sample.prob.raw, na.rm=TRUE )
            sample.prob <- sample.prob.raw/max(sample.prob.raw)
            while(ape::Ntip(geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))>ntax) {
                to.kill <- sample(names(sample.prob), size=1, prob=sample.prob)
                phy <- ape::drop.tip(phy, to.kill)
                phy$tip.state <- phy$tip.state[-which(names(phy$tip.state)==to.kill)]
                sample.prob <- sample.prob[-which(names(sample.prob)==to.kill)]
            }
        }
        trees[[nrep.index]] <- phy
    }
    return(trees)
}



test_fn <- function(rep_index, trees, scenario, ntax, true.psi) {
    
    #Process tree
    phy <- trees[[rep_index]]
    failed <- TRUE
    while(failed == TRUE){
        f <- hisse:::GetFossils(phy, psi=true.psi)
        if(sum(f$has_sampled_descendant) > 0){
            pp <- ProcessSimSample(phy, f)
            pp.strat <- ProcessSimStrat(phy, f)
            strat.cache <- hisse:::GetStratInfo(pp.strat$strat.intervals)
            if(strat.cache$l_s > 0){
                if(scenario == "scenario4"){
                    sf <- 0.75
                    failed <- FALSE
                }else{
                    sf <- 1
                    failed <- FALSE
                }
            }
        }
    }

    if(scenario == "scenario1"){
        pars <- c(0.25,0.5, 0.75, 0.75)
    }
    if(scenario == "scenario2"){
        pars <- c(0.25, 0.50, 0.55, 0.75)
    }
    if(scenario == "scenario3"){
        pars <- c(0.375, 0.375, 0.75, 0.75)
    }
    if(scenario == "scenario4"){
        pars <- c(0.5, 0.2, 0.7, 0.9)
    }
    
    #Fit full:
    misse.mk.full1 <- MiSSE(pp$phy, f=sf, turnover=1, eps=1, includes.fossils=TRUE, k.sample=pp$k.samples, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.full2 <- MiSSE(pp$phy, f=sf, turnover=c(1,2), eps=c(1,1), includes.fossils=TRUE, k.sample=pp$k.samples, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.full2 <- fixParOrder(misse.mk.full2, to.order="turn")
    misse.mk.full3 <- MiSSE(pp$phy, f=sf, turnover=c(1,1), eps=c(1,2), includes.fossils=TRUE, k.sample=pp$k.samples, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.full3 <- fixParOrder(misse.mk.full3, to.order="eps")
    misse.mk.full4 <- MiSSE(pp$phy, f=sf, turnover=c(1,2), eps=c(1,2), includes.fossils=TRUE, k.sample=pp$k.samples, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.full4 <- fixParOrder(misse.mk.full4, to.order="turn")
    model.list <- as.list(1:4)
    model.list[[1]] <- misse.mk.full1
    model.list[[2]] <- misse.mk.full2
    model.list[[3]] <- misse.mk.full3
    model.list[[4]] <- misse.mk.full4
    model.list.red <- PruneRedundantModels(model.list)
    
    aic.weights <- GetAICWeights(model.list.red, criterion = "AIC")
    model.pars.mat <- GetSimModelMat(model.list.red)
    
    model.list.full <- model.list
    model.ave.pars.full <- colSums(model.pars.mat * aic.weights)
    names(model.ave.pars.full) <- c("turnoverA", "epsA", "net_divA", "lambdaA", "muA", "turnoverB", "epsB", "net_divB", "lambdaB", "muB", "psi")
    
    ##Now do reconstructions:
    recon.list.full <- DoReconOnSet(pp$phy, model.list.red)
    misse.mk.full.tip <- GetModelAveRates(recon.list.full, type="tips")
    
    #Fit k.iffy:
    k.iffy <- pp$k.samples[sample.int(n=nrow(pp$k.samples), size=max(1, round(nrow(pp$k.samples)/2))),]
    misse.mk.missinghalf1 <- MiSSE(pp$phy, f=sf, turnover=1, eps=1, includes.fossils=TRUE, k.sample=k.iffy, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.missinghalf2 <- MiSSE(pp$phy, f=sf, turnover=c(1,2), eps=c(1,1), includes.fossils=TRUE, k.sample=k.iffy, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.missinghalf2 <- fixParOrder(misse.mk.missinghalf2, to.order="turn")
    misse.mk.missinghalf3 <- MiSSE(pp$phy, f=sf, turnover=c(1,1), eps=c(1,2), includes.fossils=TRUE, k.sample=k.iffy, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.missinghalf3 <- fixParOrder(misse.mk.missinghalf3, to.order="eps")
    misse.mk.missinghalf4 <- MiSSE(pp$phy, f=sf, turnover=c(1,2), eps=c(1,2), includes.fossils=TRUE, k.sample=k.iffy, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.missinghalf4 <- fixParOrder(misse.mk.missinghalf4, to.order="turn")
    model.list <- as.list(1:4)
    model.list[[1]] <- misse.mk.missinghalf1
    model.list[[2]] <- misse.mk.missinghalf2
    model.list[[3]] <- misse.mk.missinghalf3
    model.list[[4]] <- misse.mk.missinghalf4
    model.list.red <- PruneRedundantModels(model.list)
    
    aic.weights <- GetAICWeights(model.list.red, criterion = "AIC")
    model.pars.mat <- GetSimModelMat(model.list.red)
    
    model.list.kiffy <- model.list
    model.ave.pars.kiffy <- colSums(model.pars.mat * aic.weights)
    names(model.ave.pars.kiffy) <- c("turnoverA", "epsA", "net_divA", "lambdaA", "muA", "turnoverB", "epsB", "net_divB", "lambdaB", "muB", "psi")
    
    ##Now do reconstructions:
    recon.list.kiffy <- DoReconOnSet(pp$phy, model.list.red)
    misse.mk.missinghalf.tip <- GetModelAveRates(recon.list.kiffy, type="tips")
    
    #k.zero
    misse.zero_k1 <- MiSSE(pp$phy, f=sf, turnover=1, eps=1, includes.fossils=TRUE, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.zero_k2 <- MiSSE(pp$phy, f=sf, turnover=c(1,2), eps=c(1,1), includes.fossils=TRUE, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.zero_k2 <- fixParOrder(misse.zero_k2, to.order="turn")
    misse.zero_k3 <- MiSSE(pp$phy, f=sf, turnover=c(1,1), eps=c(1,2), includes.fossils=TRUE, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.zero_k3 <- fixParOrder(misse.zero_k3, to.order="eps")
    misse.zero_k4 <- MiSSE(pp$phy, f=sf, turnover=c(1,2), eps=c(1,2), includes.fossils=TRUE, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.zero_k4 <- fixParOrder(misse.zero_k4, to.order="turn")
    model.list[[1]] <- misse.zero_k1
    model.list[[2]] <- misse.zero_k2
    model.list[[3]] <- misse.zero_k3
    model.list[[4]] <- misse.zero_k4
    model.list.red <- PruneRedundantModels(model.list)
    
    aic.weights <- GetAICWeights(model.list.red, criterion = "AIC")
    model.pars.mat <- GetSimModelMat(model.list.red)
    
    model.list.zero <- model.list
    model.ave.pars.zero <- colSums(model.pars.mat * aic.weights)
    names(model.ave.pars.zero) <- c("turnoverA", "epsA", "net_divA", "lambdaA", "muA", "turnoverB", "epsB", "net_divB", "lambdaB", "muB", "psi")
    
    ##Now do reconstructions:
    recon.list.zero_k <- DoReconOnSet(pp$phy, model.list.red)
    misse.zero_k.tip <- GetModelAveRates(recon.list.zero_k, type="tips")
    

    #Fit to strat intervals
    misse.mk.strat1 <- MiSSE(pp.strat$phy, f=sf, turnover=1, eps=1, includes.fossils=TRUE, strat.intervals=pp.strat$strat.intervals, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.strat2 <- MiSSE(pp.strat$phy, f=sf, turnover=c(1,2), eps=c(1,1), includes.fossils=TRUE, strat.intervals=pp.strat$strat.intervals, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.strat2 <- fixParOrder(misse.mk.strat2, to.order="turn")
    misse.mk.strat3 <- MiSSE(pp.strat$phy, f=sf, turnover=c(1,1), eps=c(1,2), includes.fossils=TRUE, strat.intervals=pp.strat$strat.intervals, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.strat3 <- fixParOrder(misse.mk.strat3, to.order="eps")
    misse.mk.strat4 <- MiSSE(pp.strat$phy, f=sf, turnover=c(1,2), eps=c(1,2), includes.fossils=TRUE, strat.intervals=pp.strat$strat.intervals, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.mk.strat4 <- fixParOrder(misse.mk.strat4, to.order="turn")
    model.list <- as.list(1:4)
    model.list[[1]] <- misse.mk.strat1
    model.list[[2]] <- misse.mk.strat2
    model.list[[3]] <- misse.mk.strat3
    model.list[[4]] <- misse.mk.strat4
    model.list.red <- PruneRedundantModels(model.list)
    
    aic.weights <- GetAICWeights(model.list.red, criterion = "AIC")
    model.pars.mat <- GetSimModelMat(model.list.red)
    
    model.list.strat <- model.list
    model.ave.pars.strat <- colSums(model.pars.mat * aic.weights)
    names(model.ave.pars.strat) <- c("turnoverA", "epsA", "net_divA", "lambdaA", "muA", "turnoverB", "epsB", "net_divB", "lambdaB", "muB", "psi")
    
    ##Now do reconstructions:
    recon.list.strat <- DoReconOnSet(pp$phy, model.list.red, strat.intervals=TRUE)
    misse.mk.strat.tip <- GetModelAveRates(recon.list.strat, type="tips")
    
    #Extant only
    phy.extant <- geiger::drop.extinct(pp$phy, tol=.Machine$double.eps^.50)
    misse.extant1 <- MiSSE(phy.extant, f=sf, turnover=1, eps=1, includes.fossils=FALSE, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.extant2 <- MiSSE(phy.extant, f=sf, turnover=c(1,2), eps=c(1,1), includes.fossils=FALSE, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.extant2 <- fixParOrder(misse.extant2, to.order="turn")
    misse.extant3 <- MiSSE(phy.extant, f=sf, turnover=c(1,1), eps=c(1,2), includes.fossils=FALSE, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.extant3 <- fixParOrder(misse.extant3, to.order="eps")
    misse.extant4 <- MiSSE(phy.extant, f=sf, turnover=c(1,2), eps=c(1,2), includes.fossils=FALSE, trans.upper=10, turnover.upper=100, eps.upper=1.5)
    misse.extant4 <- fixParOrder(misse.extant4, to.order="turn")
    model.list[[1]] <- misse.extant1
    model.list[[2]] <- misse.extant2
    model.list[[3]] <- misse.extant3
    model.list[[4]] <- misse.extant4
    model.list.red <- PruneRedundantModels(model.list)
    
    aic.weights <- GetAICWeights(model.list.red, criterion = "AIC")
    model.pars.mat <- GetSimModelMat(model.list.red)
    
    model.list.extant <- model.list
    model.ave.pars.extant <- colSums(model.pars.mat * aic.weights)
    names(model.ave.pars.extant) <- c("turnoverA", "epsA", "net_divA", "lambdaA", "muA", "turnoverB", "epsB", "net_divB", "lambdaB", "muB", "psi")
    
    ##Now do reconstructions:
    recon.list.extant <- DoReconOnSet(phy.extant, model.list.red)
    misse.extant.tip <- GetModelAveRates(recon.list.extant, type="tips")
    
    m = ape::Ntip(phy) - ape::Ntip(phy.extant)
    m.used = ape::Ntip(pp$phy) - ape::Ntip(phy.extant)
    k = nrow(pp$k.samples)
    n.extant <- ape::Ntip(phy.extant)
    true.tip.rates <- phy$tip.state[names(phy$tip.state) %in% phy.extant$tip.label]
    
    local_results <- rbind(
    data.frame(t(c(misse="extant", rep=rep_index, true.lambda=pars[1:2]/(1+pars[3:4]), true.mu=(pars[1:2]*pars[3:4])/(1+pars[3:4]), true.psi=true.psi, m=m, m.used=NA, k=k, k.used=NA, n.extant=n.extant, true.turnover=pars[1:2], true.eps=pars[3:4], model.ave.pars.extant)), stringsAsFactors=FALSE),
    data.frame(t(c(misse="full_k", rep=rep_index, true.lambda=pars[1:2]/(1+pars[3:4]), true.mu=(pars[1:2]*pars[3:4])/(1+pars[3:4]), true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=k, n.extant=n.extant, true.turnover=pars[1:2], true.eps=pars[3:4], model.ave.pars.full)), stringsAsFactors=FALSE),
    data.frame(t(c(misse="half_k", rep=rep_index, true.lambda=pars[1:2]/(1+pars[3:4]), true.mu=(pars[1:2]*pars[3:4])/(1+pars[3:4]), true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=nrow(k.iffy), n.extant=n.extant, true.turnover=pars[1:2], true.eps=pars[3:4], model.ave.pars.kiffy)), stringsAsFactors=FALSE),
    data.frame(t(c(misse="zero_k", rep=rep_index, true.lambda=pars[1:2]/(1+pars[3:4]), true.mu=(pars[1:2]*pars[3:4])/(1+pars[3:4]), true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=0, n.extant=n.extant, true.turnover=pars[1:2], true.eps=pars[3:4], model.ave.pars.zero)), stringsAsFactors=FALSE),
    data.frame(t(c(misse="strat_int", rep=rep_index, true.lambda=pars[1:2]/(1+pars[3:4]), true.mu=(pars[1:2]*pars[3:4])/(1+pars[3:4]), true.psi=true.psi, m=m, m.used=m.used, k=0, k.used=0, n.extant=n.extant, true.turnover=pars[1:2], true.eps=pars[3:4], model.ave.pars.strat)), stringsAsFactors=FALSE),
    stringsAsFactors=FALSE)
    
    session.stuff <- sessionInfo()
    save(session.stuff, model.list.extant, model.list.zero, model.list.kiffy, model.list.full, model.list.strat, phy, f, pp, local_results, misse.mk.full.tip, misse.mk.missinghalf.tip, misse.zero_k.tip, misse.mk.strat.tip, misse.extant.tip, true.tip.rates, recon.list.full, recon.list.kiffy, recon.list.zero_k, recon.list.extant, recon.list.strat, file=paste(scenario, ntax, true.psi, rep_index, "Rsave", sep="."))
    return(local_results)
}


DoSim <- function(rep.index, trees, scenario, ntax, true.psi){
    failed <- TRUE
    while(failed == TRUE){
        tmp <- NA
        try(tmp <- test_fn(rep_index=rep.index, trees=trees, scenario=scenario, ntax=ntax, true.psi=true.psi))
        if(length(tmp)>1){
            print("you have done it!")
            failed <- FALSE
        }else{
            print("oops!")
        }
    }
}


#load("Scenario1trees.Rsave")
#set1 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario1", ntax=200, true.psi=0.01, mc.cores=20)
#set2 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario1", ntax=200, true.psi=0.05, mc.cores=20)
#set3 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario1", ntax=200, true.psi=0.1, mc.cores=20)

#load("Scenario2trees.Rsave")
#set1 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario2", ntax=200, true.psi=0.01, mc.cores=20)
#set2 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario2", ntax=200, true.psi=0.05, mc.cores=20)
#set3 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario2", ntax=200, true.psi=0.1, mc.cores=20)

#load("Scenario3trees.Rsave")
#set1 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario3", ntax=200, true.psi=0.01, mc.cores=20)
#set2 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario3", ntax=200, true.psi=0.05, mc.cores=20)
#set3 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario3", ntax=200, true.psi=0.1, mc.cores=20)

load("Scenario4trees.Rsave")
#set1 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario4", ntax=200, true.psi=0.01, mc.cores=20)
#set2 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario4", ntax=200, true.psi=0.05, mc.cores=20)
#set3 <- mclapply(1:100, DoSim, trees=trees, scenario="scenario4", ntax=200, true.psi=0.1, mc.cores=20)




