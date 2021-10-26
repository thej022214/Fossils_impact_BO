
rm(list=ls())

library(diversitree)
library(hisse)
library(parallel)

set.seed(42)

test_fn <- function(rep_index, scenario, true.psi) {
    #phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = true.lambda, mu = true.mu)[[1]]
    if(scenario == "scenario1"){
        #Net div and turnover vary by ~3x in rate class 2
        pars <- c(0.1677143, 0.5031429, 0.1677143*.75, 0.5031429*.75, 0.005, 0.005)
        all.dead <- TRUE
        phy <- NULL
        while(all.dead){
            try(phy <- tree.bisse(pars, max.taxa=300, x0=0, include.extinct=TRUE))
            if(!is.null(phy)){
                if(Ntip(phy)>300){
                    all.dead <- FALSE
                }
            }
        }
    }
    if(scenario == "scenario2"){
        #Net div and turnover is higher by ~1.5x in rate class 2
        pars <- c(0.2683429, 0.4025143, 0.2683429*.75, 0.4025143*.75, 0.005, 0.005)
        all.dead <- TRUE
        phy <- NULL
        while(all.dead){
            try(phy <- tree.bisse(pars, max.taxa=300, x0=0, include.extinct=TRUE))
            if(!is.null(phy)){
                if(Ntip(phy)>300){
                    all.dead <- FALSE
                }
            }
        }
    }
    if(scenario == "scenario3"){
        #Net div the same, but turnover varies ~2x, and eps varies
        pars <- c(0.1417773, 0.2835546, 0.1417773*0.5, 0.2835546*.75, 0.005, 0.005)
        all.dead <- TRUE
        phy <- NULL
        while(all.dead){
            try(phy <- tree.bisse(pars, max.taxa=300, x0=0, include.extinct=TRUE))
            if(!is.null(phy)){
                if(Ntip(phy)>300){
                    all.dead <- FALSE
                }
            }
        }
    }
    if(scenario == "scenario4"){
        #Single rate
        pars <- c(0.3354962, 0.3354962, 0.3354962*0.75, 0.3354962*0.75, 0.005, 0.005)
        phy <- TreeSim::sim.bd.taxa(n = 300, numbsim = 1, lambda = 0.3354962, mu = 0.3354962*0.75)[[1]]
    }
    
    f <- hisse:::GetFossils(phy, psi=true.psi)
    pp <- hisse:::ProcessSimSample(phy, f)

    #full:
    misse.mk.full1 <- MiSSE(pp$phy, turnover=1, eps=1, includes.fossils=TRUE, k.sample=pp$k.samples, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.mk.full2 <- MiSSE(pp$phy, turnover=c(1,2), eps=c(1,1), includes.fossils=TRUE, k.sample=pp$k.samples, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.mk.full3 <- MiSSE(pp$phy, turnover=c(1,1), eps=c(1,2), includes.fossils=TRUE, k.sample=pp$k.samples, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.mk.full4 <- MiSSE(pp$phy, turnover=c(1,2), eps=c(1,2), includes.fossils=TRUE, k.sample=pp$k.samples, sann=TRUE, trans.upper=10, turnover.upper=100)
    model.list <- as.list(1:4)
    model.list[[1]] <- misse.mk.full1
    model.list[[2]] <- misse.mk.full2
    model.list[[3]] <- misse.mk.full3
    model.list[[4]] <- misse.mk.full4
    aic.weights <- GetAICWeights(model.list)
    model.pars.mat <- rbind(
    c(hisse:::ParameterTransformMiSSESpecial(misse.mk.full1$solution)[1:5], hisse:::ParameterTransformMiSSESpecial(misse.mk.full1$solution)[1:5], misse.mk.full1$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.mk.full2$solution)[1:10], misse.mk.full2$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.mk.full3$solution)[1:10], misse.mk.full3$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.mk.full4$solution)[1:10], misse.mk.full4$solution["psi"])
    )
    model.ave.pars.full <- colSums(model.pars.mat * aic.weights)
    names(model.ave.pars.full) <- c("turnoverA", "epsA", "net_divA", "lambdaA", "muA", "turnoverB", "epsB", "net_divB", "lambdaB", "muB", "psi")

    #k.iffy:
    k.iffy <- pp$k.samples[sample.int(n=nrow(pp$k.samples), size=max(1, round(nrow(pp$k.samples)/2))),]
    misse.mk.missinghalf1 <- MiSSE(pp$phy, turnover=1, eps=1, includes.fossils=TRUE, k.sample=k.iffy,  sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.mk.missinghalf2 <- MiSSE(pp$phy, turnover=c(1,2), eps=c(1,1), includes.fossils=TRUE, k.sample=k.iffy, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.mk.missinghalf3 <- MiSSE(pp$phy, turnover=c(1,1), eps=c(1,2), includes.fossils=TRUE, k.sample=k.iffy, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.mk.missinghalf4 <- MiSSE(pp$phy, turnover=c(1,2), eps=c(1,2), includes.fossils=TRUE, k.sample=k.iffy, sann=TRUE, trans.upper=10, turnover.upper=100)
    model.list <- as.list(1:4)
    model.list[[1]] <- misse.mk.missinghalf1
    model.list[[2]] <- misse.mk.missinghalf2
    model.list[[3]] <- misse.mk.missinghalf3
    model.list[[4]] <- misse.mk.missinghalf4
    aic.weights <- GetAICWeights(model.list)
    model.pars.mat <- rbind(
    c(hisse:::ParameterTransformMiSSESpecial(misse.mk.missinghalf1$solution)[1:5], hisse:::ParameterTransformMiSSESpecial(misse.mk.missinghalf1$solution)[1:5], misse.mk.missinghalf1$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.mk.missinghalf2$solution)[1:10], misse.mk.missinghalf2$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.mk.missinghalf3$solution)[1:10], misse.mk.missinghalf3$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.mk.missinghalf4$solution)[1:10], misse.mk.missinghalf4$solution["psi"])
    )
    model.ave.pars.kiffy <- colSums(model.pars.mat * aic.weights)
    names(model.ave.pars.kiffy) <- c("turnoverA", "epsA", "net_divA", "lambdaA", "muA", "turnoverB", "epsB", "net_divB", "lambdaB", "muB", "psi")


    #k.zero
    misse.zero_k1 <- MiSSE(pp$phy, turnover=1, eps=1, includes.fossils=TRUE, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.zero_k2 <- MiSSE(pp$phy, turnover=c(1,2), eps=c(1,1), includes.fossils=TRUE, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.zero_k3 <- MiSSE(pp$phy, turnover=c(1,1), eps=c(1,2), includes.fossils=TRUE, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.zero_k4 <- MiSSE(pp$phy, turnover=c(1,2), eps=c(1,2), includes.fossils=TRUE, sann=TRUE, trans.upper=10, turnover.upper=100)
    model.list[[1]] <- misse.zero_k1
    model.list[[2]] <- misse.zero_k2
    model.list[[3]] <- misse.zero_k3
    model.list[[4]] <- misse.zero_k4
    aic.weights <- GetAICWeights(model.list)
    model.pars.mat <- rbind(
    c(hisse:::ParameterTransformMiSSESpecial(misse.zero_k1$solution)[1:5], hisse:::ParameterTransformMiSSESpecial(misse.zero_k1$solution)[1:5], misse.zero_k1$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.zero_k2$solution)[1:10], misse.zero_k2$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.zero_k3$solution)[1:10], misse.zero_k3$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.zero_k4$solution)[1:10], misse.zero_k4$solution["psi"])
    )
    model.ave.pars.zero <- colSums(model.pars.mat * aic.weights)
    names(model.ave.pars.zero) <- c("turnoverA", "epsA", "net_divA", "lambdaA", "muA", "turnoverB", "epsB", "net_divB", "lambdaB", "muB", "psi")


    phy.extant <- geiger::drop.extinct(phy)
    misse.extant1 <- MiSSE(phy.extant, turnover=1, eps=1, includes.fossils=FALSE, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.extant2 <- MiSSE(phy.extant, turnover=c(1,2), eps=c(1,1), includes.fossils=FALSE, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.extant3 <- MiSSE(phy.extant, turnover=c(1,1), eps=c(1,2), includes.fossils=FALSE, sann=TRUE, trans.upper=10, turnover.upper=100)
    misse.extant4 <- MiSSE(phy.extant, turnover=c(1,2), eps=c(1,2), includes.fossils=FALSE, sann=TRUE, trans.upper=10, turnover.upper=100)
    model.list[[1]] <- misse.extant1
    model.list[[2]] <- misse.extant2
    model.list[[3]] <- misse.extant3
    model.list[[4]] <- misse.extant4
    aic.weights <- GetAICWeights(model.list)
    model.pars.mat <- rbind(
    c(hisse:::ParameterTransformMiSSESpecial(misse.extant1$solution)[1:5], hisse:::ParameterTransformMiSSESpecial(misse.extant1$solution)[1:5], misse.extant1$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.extant2$solution)[1:10], misse.extant2$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.extant3$solution)[1:10], misse.extant3$solution["psi"]),
    c(hisse:::ParameterTransformMiSSESpecial(misse.extant4$solution)[1:10], misse.extant4$solution["psi"])
    )
    model.ave.pars.extant <- colSums(model.pars.mat * aic.weights)
    names(model.ave.pars.extant) <- c("turnoverA", "epsA", "net_divA", "lambdaA", "muA", "turnoverB", "epsB", "net_divB", "lambdaB", "muB", "psi")

    m = ape::Ntip(phy) - ape::Ntip(phy.extant)
    m.used = ape::Ntip(pp$phy) - ape::Ntip(phy.extant)
    k = nrow(pp$k.samples)
    n.extant <- ape::Ntip(phy.extant)
    
    local_results <- rbind(
    data.frame(t(c(misse="extant", rep=rep_index, true.lambda=pars[1:2], true.mu=pars[3:4], true.psi=true.psi, m=m, m.used=NA, k=k, k.used=NA, n.extant=n.extant, true.turnover=pars[1:2]+pars[3:4], true.eps=pars[3:4]/pars[1:2], model.ave.pars.extant)), stringsAsFactors=FALSE),
    data.frame(t(c(misse="full_k", rep=rep_index, true.lambda=pars[1:2], true.mu=pars[3:4], true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=k, n.extant=n.extant, true.turnover=pars[1:2]+pars[3:4], true.eps=pars[3:4]/pars[1:2], model.ave.pars.full)), stringsAsFactors=FALSE),
    data.frame(t(c(misse="half_k", rep=rep_index, true.lambda=pars[1:2], true.mu=pars[3:4], true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=nrow(k.iffy), n.extant=n.extant, true.turnover=pars[1:2]+pars[3:4], true.eps=pars[3:4]/pars[1:2], model.ave.pars.kiffy)), stringsAsFactors=FALSE),
    data.frame(t(c(misse="zero_k", rep=rep_index, true.lambda=pars[1:2], true.mu=pars[3:4], true.psi=true.psi, m=m, m.used=m.used, k=k, k.used=0, n.extant=n.extant, true.turnover=pars[1:2]+pars[3:4], true.eps=pars[3:4]/pars[1:2], model.ave.pars.zero)), stringsAsFactors=FALSE),
    stringsAsFactors=FALSE)
    save(phy, f, pp, local_results, file=paste(scenario, true.psi, rep_index, "Rsave", sep="."))
    return(local_results)
}


set1 <- mclapply(1:100, test_fn, scenario="scenario1", true.psi=0.01, mc.cores=40)
save(set1, file="scenario1_0.01_all.Rsave)
restart.set <- mclapply(1:100, scenario="scenario1", true.psi=0.05,  mc.cores=40)
save(set2, file="scenario1_0.05_all.Rsave)
restart.set <- mclapply(1:100, scenario="scenario1", true.psi=0.10, mc.cores=40)
save(set3, file="scenario1_0.10_all.Rsave)









