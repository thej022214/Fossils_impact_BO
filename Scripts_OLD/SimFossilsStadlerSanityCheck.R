rm(list=ls())
# set.seed(42)
# phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = 0.3, mu = 0.1)[[1]]
# f <- GetFossils(phy, psi=0.1)
#
library(hisse)

GetStadlerEst <- function(phy, k.samples){
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1)
    fossil.taxa <- which(dat.tab$branch.type == 1)
    #This is for starting values:
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$branch.type == 1)]
    phy.no.k <- phy
    split.times <- paleotree:::dateNodes(phy.no.k, rootAge=max(node.depth.edgelength(phy.no.k)))[-c(1:Ntip(phy.no.k))]
    n <- Ntip(phy.no.k)-length(fossil.taxa)
    m <- length(fossil.taxa)
    x_times <- split.times
    y_times <- fossil.ages
    if(is.null(k.samples)){
        k <- 0
    }else{
        k <- dim(k.samples)[1]
    }
    model.fit <- hisse:::starting.point.generator.fossils(n.tax=n, k=1, samp.freq.tree=1, q.div=5, fossil.taxa=fossil.taxa, fossil.ages=fossil.ages, no.k.samples=k, split.times=split.times, get.likelihood=TRUE)
    res <- c(model.fit[1] +  model.fit[2],  model.fit[2]/model.fit[1], model.fit[3])
    names(res) <- c("turnover0A", "eps0A", "psi")
    return(res)
}


#
# MiSSE(phy, includes.fossils=TRUE)
set.seed(42)
test_fn <- function(true.lambda, true.mu, true.psi) {
    phy <- TreeSim::sim.bd.taxa(n = 100, numbsim = 1, lambda = true.lambda, mu = true.mu)[[1]]
    f <- hisse:::GetFossils(phy, psi=true.psi)
    #save(phy,f, file="checkingSims1.Rsave")
    pp <- hisse:::ProcessSimSample(phy, f)
    #save(pp, file="checkingSims2.Rsave")
    misse.mk.full <- GetStadlerEst(pp$phy, k.sample=pp$k.samples)
    k.iffy <- pp$k.samples[sample.int(n=nrow(pp$k.samples), size=max(1, round(nrow(pp$k.samples)/2))),]
    misse.mk.missinghalf <- GetStadlerEst(pp$phy, k.sample=k.iffy)
    misse.zero_k <- GetStadlerEst(pp$phy, k.sample=NULL)
    phy.extant <- geiger::drop.extinct(phy)
    misse.extant <- GetStadlerEst(phy.extant, k.samples=NULL)
    m = ape::Ntip(pp$phy) - ape::Ntip(phy.extant)
    k = nrow(pp$k.samples)
    n.extant <- ape::Ntip(phy.extant)
    local_results <- rbind(
    data.frame(t(c(misse="extant",  rep=rep_index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi, m=m, k=k, k.used=NA, n.extant=n.extant, true.turnover=true.turnover, true.eps=true.eps, misse.extant[c("turnover0A", "eps0A", "psi")])), stringsAsFactors=FALSE),
    data.frame(t(c(misse="full_k",  rep=rep_index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi, m=m, k=k, k.used=k, n.extant=n.extant, true.turnover=true.turnover, true.eps=true.eps, misse.mk.full[c("turnover0A", "eps0A", "psi")])), stringsAsFactors=FALSE),
    data.frame(t(c(misse="half_k",  rep=rep_index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi, m=m, k=k, k.used=nrow(k.iffy), n.extant=n.extant, true.turnover=true.turnover, true.eps=true.eps, misse.mk.missinghalf[c("turnover0A", "eps0A", "psi")])), stringsAsFactors=FALSE),
    data.frame(t(c(misse="zero_k",  rep=rep_index, true.lambda=true.lambda, true.mu=true.mu, true.psi=true.psi, m=m, k=k, k.used=0, n.extant=n.extant, true.turnover=true.turnover, true.eps=true.eps, misse.zero_k[c("turnover0A", "eps0A", "psi")])), stringsAsFactors=FALSE),
    stringsAsFactors=FALSE
    )
    print(local_results)
    return(local_results)
}


results <- data.frame()
true.lambda = 0.3
true.mu.vector <- c( 0.2, 0.1)
true.psi.vector <- c(0.05, 0.01)
for (mu_index in seq_along(true.mu.vector)) {
    for (psi_index in seq_along(true.psi.vector)) {
        true.mu <- true.mu.vector[mu_index]
        true.psi <- true.psi.vector[psi_index]
        for (rep_index in sequence(50)) {
            true.turnover = true.lambda+true.mu
            true.eps=true.mu/true.lambda
            local_results <- data.frame()
            try(local_results <- test_fn(true.lambda, true.mu, true.psi))
            if(nrow(local_results)>0) {
                try(results <- rbind(results, local_results, stringsAsFactors=FALSE))
                write.csv(results, file="results.csv", quote=FALSE, row.names=FALSE)
            }
        }
    }
}

for (i in 2:ncol(results)) {
    results[,i] <- as.numeric(results[,i])
}
results$error.eps <- as.numeric(results$eps0A)-as.numeric(results$true.eps)
results$error.turnover <- as.numeric(results$turnover0A)-as.numeric(results$true.turnover)
results$misse <- as.factor(results$misse)
write.csv(results, file="results.csv")
library(tidyverse)
library(Metrics)
results_summary <- results %>% group_by(true.mu, true.psi, misse) %>% summarise(median_eps_error=median(error.eps), rmse_eps=rmse(true.eps, eps0A), median_turnover_error=median(error.turnover), rmse_turnover=rmse(true.turnover, turnover0A))
write.csv(results_summary, file="results_summary.csv")#
results_highlevel <- results %>% group_by(misse) %>% summarise(median_eps_error=median(error.eps), rmse_eps=rmse(true.eps, eps0A), median_turnover_error=median(error.turnover), rmse_turnover=rmse(true.turnover, turnover0A))
write.csv(results_highlevel, file="results_highlevel.csv")
