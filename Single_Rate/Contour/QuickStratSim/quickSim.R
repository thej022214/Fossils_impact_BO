

library(hisse)

GetSimRep <- function(phy, f, true.psi){
    
    ntax <- Ntip(geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
    ss <- hisse:::ProcessSimStrat(phy, f)
    #save(phy, f, ss, file="bugfinding.Rsave")
    split.times.plus.tips <- paleotree::dateNodes(ss$phy, rootAge=max(node.depth.edgelength(ss$phy)))
    split.times <- split.times.plus.tips[-c(1:Ntip(ss$phy))]
    strat.cache <- hisse:::GetStratInfo(ss$strat.intervals)
    print(strat.cache)
    if(strat.cache$l_s == 0){
        f <- hisse:::GetFossils(phy, psi=true.psi)
        ss <- hisse:::ProcessSimStrat(phy, f)
        #save(phy, f, ss, file="bugfindingAfterRedo.Rsave")
        split.times.plus.tips <- paleotree::dateNodes(ss$phy, rootAge=max(node.depth.edgelength(ss$phy)))
        split.times <- split.times.plus.tips[-c(1:Ntip(ss$phy))]
        strat.cache <- hisse:::GetStratInfo(ss$strat.intervals)
        #print(strat.cache)
    }
    k.samples <- hisse:::GetIntervalToK(ss$strat.intervals, intervening.intervals=strat.cache$intervening.intervals)
    extinct.tips <- which(round(as.numeric(k.samples$timefrompresent),8) %in% round(as.numeric(split.times.plus.tips[c(1:Ntip(ss$phy))]),8))
    if(length(extinct.tips > 0)){
        k.samples <- k.samples[-extinct.tips,]
    }
    phy <- hisse:::AddKNodes(ss$phy, k.samples)
    fix.type <- hisse:::GetKSampleMRCA(phy, k.samples, strat.intervals=TRUE)
    edge_details <- hisse:::GetEdgeDetails(phy, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    
    gen <- hisse:::FindGenerations(phy)
    dat.tab <- hisse:::OrganizeDataMiSSE(phy=phy, f=1, hidden.states=1, includes.intervals=TRUE, intervening.intervals=strat.cache$intervening.intervals)
    fossil.taxa <- edge_details$tipward_node[which(edge_details$type == "extinct_tip" | edge_details$type == "k_extinct_interval")]
    fossil.ages <- dat.tab$TipwardAge[which(dat.tab$DesNode %in% fossil.taxa)]
    
    cols <- c("FocalNode","DesNode", "RootwardAge", "TipwardAge", "branch.type")
    seg.map <- dat.tab[, cols, with=FALSE]
    #remove k tips -- we do not do anything with them.
    data.table::setkey(seg.map, branch.type)
    #drop the k.tips because we do not do calculation on these zero length edges:
    seg.map <- seg.map[branch.type != 2]
    #######################
    
    prop.strat.length <- strat.cache$l_s / sum(phy$edge.length)
    
    start.points <- hisse:::starting.point.generator.intervals(k=1, n=ntax, samp.freq=1, seg_map=seg.map, split.times=split.times, fossil.ages=fossil.ages, strat.cache=strat.cache, get.likelihood=TRUE)
    
    ests <- c(start.points[1]+start.points[2], start.points[2]/start.points[1], start.points[3])
    
    return(c(ests,prop.strat.length))
    
}


DoReps <- function(nreps, trees=NULL){
    true.psi <- seq(0.05, 0.15, 0.01)
    ntax=200
    if(is.null(trees[[1]])){
        set.seed(4)
        try(sim.tab <- hisse::SimulateHisse(turnover=c(0.7,0.7), eps=rep(0.75,2), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
        phy <- hisse:::SimToPhylo(sim.tab, include.extinct=TRUE)
        res <- c()
        for(psi.index in 1:length(true.psi)){
            print(psi.index)
            for(rep.index in 1:nreps){
                print(rep.index)
                f <- hisse:::GetFossils(phy, psi=true.psi[psi.index])
                est.stad <- GetSimRep(phy, f, true.psi[psi.index])
                res <- rbind(res, c(est.stad,true.psi[psi.index]))
            }
        }
    }else{
        res <- c()
        for(psi.index in 1:length(true.psi)){
            print(psi.index)
            print(true.psi[psi.index])
            for(rep.index in 1:nreps){
                print(rep.index)
                phy <- trees[[rep.index]]
                try(f <- hisse:::GetFossils(phy, psi=true.psi[psi.index]))
                print("got the fossil sample")
                est.stad <- NULL
                try(est.stad <- GetSimRep(phy, f, true.psi[psi.index]))
                print("est.stad done")
                print(est.stad)
                if(!is.null(est.stad[1])){
                    res <- rbind(res, c(est.stad,true.psi[psi.index]))
                }
            }
        }
    }
    colnames(res) <- c("turn", "eps", "psi", "prop.strat", "true.psi")
    return(res)
}

#load("/Users/jeremybeaulieu/Desktop/randomSet.quicksim.Rsave")
#pp <- DoReps(100, trees=trees)
#save(pp, file="quickSimRandom.Rsave")


#load("quickSimCHECK.Rsave")
#turn.bias <- pp[,1] - 0.30
#eps.bias <- pp[,2] - 0.75
#psi.bias <- pp[,3] - pp[,5]

#turn.bias <- turn.bias[-which(eps.bias > .2)]
#psi.bias <- psi.bias[-which(eps.bias > .2)]
#prop.int <- pp[-which(eps.bias > .2),4]
#eps.bias <- eps.bias[-which(eps.bias > .2)]


#layout(matrix(c(1:3), nrow=1, byrow=TRUE))
#plot(prop.int,turn.bias, xlab="Prop. Interval of Tree length", ylab="bias in turn")
#abline(h=0, lty=2, col="red")
#plot(prop.int,eps.bias, xlab="Prop. Interval of Tree length", ylab="bias in eps")
#abline(h=0, lty=2, col="red")
#plot(prop.int,psi.bias, xlab="Prop. Interval of Tree length", ylab="bias in psi")
#abline(h=0, lty=2, col="red")


#layout(matrix(c(1:2), nrow=1, byrow=TRUE))
#boxplot(turn.bias~pp[,5], xlab="true psi", ylab="bias in turn")
#abline(h=0, lty=2, col="red")
#boxplot(eps.bias~pp[,5], xlab="true psi", ylab="bias in eps")
#abline(h=0, lty=2, col="red")

