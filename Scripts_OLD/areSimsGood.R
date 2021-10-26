
library(hisse)
library(paleotree)

#TESTING SCENARIO 1:

prefix <- paste("/Users/jmbeauli/hisse_fossil_sims/Two_Rate/Scenario1_redone_recons/1000_taxa/scenario1.1000.*")
files <- system(paste("ls -1 ", prefix, ".*", sep=""), intern=TRUE)

max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()

for(index in 1:length(files)){
    load(files[index])
    max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
    max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
    max.age.extant <- c(max.age.extant, max(dateNodes(model.list.extant[[1]]$phy, rootAge=max(node.depth.edgelength(model.list.extant[[1]]$phy)))[-c(1:Ntip(model.list.extant[[1]]$phy))]))
}

#> mean(max.age.true)
#[1] 130.4624
#> mean(max.age.sampled)
#[1] 117.812
#> mean(max.age.extant)
#[1] 107.9826


max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()
transitions.raw <- c()
ten.percent <- c()
state.freqs <- c()
ntax=200
true.psi=0.01

for(index in 1:20){
    print(index)
    #Net div and turnover is higher by ~2x in rate class 2
    pars <- c(0.25, 0.5, 0.75, 0.75)
    all.dead <- TRUE
    phy <- NULL
    while(all.dead){
        try(sim.tab <- hisse::SimulateHisse(turnover=c(0.25,0.25), eps=rep(0.25,2), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
        if(sim.tab$n.surviving > 0){
            try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
            if(Ntip(phy) > ntax){
                f <- hisse:::GetFossils(phy, psi=true.psi)
                if(sum(f$has_sampled_descendant) > 0){
                    try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
                    if(!is.null(test.tol)){
                        if(Ntip(test.tol) == ntax){
                            pp <- hisse:::ProcessSimSample(phy, f)
                            sf <- 1
                            all.dead <- FALSE
                            max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
                            max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
                            max.age.extant <- c(max.age.extant, max(dateNodes(test.tol, rootAge=max(node.depth.edgelength(test.tol)))[-c(1:Ntip(test.tol))]))
                            transitions.raw <- rbind(transitions.raw, sim.tab$transition.counts[!is.na(sim.tab$transition.counts)])
                            ten.percent <- c(ten.percent, sum(sim.tab$transition.counts[!is.na(sim.tab$transition.counts)]) / length(phy$edge.length))
                            true.tip.states <- phy$tip.state[names(phy$tip.state) %in% test.tol$tip.label]
                            state.freqs <- rbind(state.freqs, table(true.tip.states))
                        }
                    }
                }
            }
        }
    }
}

#> mean(max.age.true)
#[1] 132.6608
#> mean(max.age.sampled)
#[1] 107.552
#> mean(max.age.extant)
#[1] 105.0721
#> mean(ten.percent)
#[1] 0.01144474
#> colMeans(state.freqs)
#0     1
#304.5 695.5

max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()
transitions.raw <- c()
ten.percent <- c()
state.freqs <- c()
k.foss <- c()
m.foss <- c()
ntax=200
true.psi=0.1

for(index in 1:20){
    print(index)
    #Net div and turnover is higher by ~2x in rate class 2
    pars <- c(0.25, 0.50, 0.55, 0.75)
    all.dead <- TRUE
    phy <- NULL
    while(all.dead){
        try(sim.tab <- hisse::SimulateHisse(turnover=c(0.25,0.5), eps=rep(0.75,2), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
        if(sim.tab$n.surviving > 0){
            try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
            if(Ntip(phy) > ntax){
                f <- hisse:::GetFossils(phy, psi=true.psi)
                if(sum(f$has_sampled_descendant) > 0){
                    try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
                    if(!is.null(test.tol)){
                        if(Ntip(test.tol) == ntax){
                            pp <- hisse:::ProcessSimSample(phy, f)
                            sf <- 1
                            all.dead <- FALSE
                            max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
                            max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
                            max.age.extant <- c(max.age.extant, max(dateNodes(test.tol, rootAge=max(node.depth.edgelength(test.tol)))[-c(1:Ntip(test.tol))]))
                            transitions.raw <- rbind(transitions.raw, sim.tab$transition.counts[!is.na(sim.tab$transition.counts)])
                            ten.percent <- c(ten.percent, sum(sim.tab$transition.counts[!is.na(sim.tab$transition.counts)]) / length(phy$edge.length))
                            true.tip.states <- phy$tip.state[names(phy$tip.state) %in% test.tol$tip.label]
                            k.foss <- c(k.foss, dim(pp$k.samples)[1])
                            m.foss <- c(m.foss, Ntip(pp$phy)-ntax)
                            state.freqs <- rbind(state.freqs, table(true.tip.states))
                        }
                    }
                }
            }
        }
    }
}

#> mean(k.foss)
#[1] 17.65

#0.01
#> mean(m.foss)
#[1] 21.9
#0.05
#> mean(k.foss)
#[1] 106.65
#> mean(m.foss)
#[1] 102.6
#0.10
#> mean(k.foss)
#[1] 248.2
#> mean(m.foss)
#[1] 175.9

#> mean(max.age.true)
#[1] 111.563
#> mean(max.age.sampled)
#[1] 100.4082
#> mean(max.age.extant)
#[1] 83.30586
#> mean(ten.percent)
#[1] 0.01367862
#> colMeans(state.freqs)
#    0     1
# 94.7 105.3


#TESTING SCENARIO 2:

prefix <- paste("/Users/jmbeauli/hisse_fossil_sims/Two_Rate/Scenario2_redone_recons/1000_taxa/scenario2.1000.*")
files <- system(paste("ls -1 ", prefix, ".*", sep=""), intern=TRUE)

max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()

for(index in 1:length(files)){
    load(files[index])
    max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
    max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
    max.age.extant <- c(max.age.extant, max(dateNodes(model.list.extant[[1]]$phy, rootAge=max(node.depth.edgelength(model.list.extant[[1]]$phy)))[-c(1:Ntip(model.list.extant[[1]]$phy))]))
}

#> mean(max.age.true)
#[1] 88.54889
#> mean(max.age.sampled)
#[1] 80.20996
#> mean(max.age.extant)
#[1] 78.05654

max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()
transitions.raw <- c()
ten.percent <- c()
state.freqs <- c()
k.foss <- c()
m.foss <- c()
ntax=200
true.psi=0.1

for(index in 1:20){
    print(index)
    #Net div and turnover is higher by ~2x in rate class 2
    pars <- c(0.25, 0.50, 0.55, 0.75)
    all.dead <- TRUE
    phy <- NULL
    while(all.dead){
        try(sim.tab <- hisse::SimulateHisse(turnover=c(0.25,0.50), eps=c(0.55,0.75), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0, 0.01, NA), nrow=2), nstart=2))
        if(sim.tab$n.surviving > 0){
            try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
            if(Ntip(phy) > ntax){
                f <- hisse:::GetFossils(phy, psi=true.psi)
                if(sum(f$has_sampled_descendant) > 0){
                    try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
                    if(!is.null(test.tol)){
                        if(Ntip(test.tol) == ntax){
                            pp <- hisse:::ProcessSimSample(phy, f)
                            sf <- 1
                            all.dead <- FALSE
                            max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
                            max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
                            max.age.extant <- c(max.age.extant, max(dateNodes(test.tol, rootAge=max(node.depth.edgelength(test.tol)))[-c(1:Ntip(test.tol))]))
                            transitions.raw <- rbind(transitions.raw, sim.tab$transition.counts[!is.na(sim.tab$transition.counts)])
                            ten.percent <- c(ten.percent, sum(sim.tab$transition.counts[!is.na(sim.tab$transition.counts)]) / length(phy$edge.length))
                            true.tip.states <- phy$tip.state[names(phy$tip.state) %in% test.tol$tip.label]
                            k.foss <- c(k.foss, dim(pp$k.samples)[1])
                            m.foss <- c(m.foss, Ntip(pp$phy)-ntax)
                            state.freqs <- rbind(state.freqs, table(true.tip.states))
                        }
                    }
                }
            }
        }
    }
}


library(diversitree)

for(index in 1:20){
    print(index)
    pars <- c(0.25/(1+.55), 0.5/(1+.75), (0.25*.55)/(1+.55), (0.5*.75)/(1+.75), 0.01, 0)
    phy <- NULL
    while( is.null( phy ) ){
        phy <- tree.bisse(pars, max.taxa=1000, x0=0)
    }
    quartz()
    h <- history.from.sim.discrete(phy, 0:1)
    plot(h, phy)
}

#0.1
> mean(k.foss)
[1] 172.3
> mean(m.foss)
[1] 89.75

#> mean(max.age.true)
#[1] 88.32575
#> mean(max.age.sampled)
#[1] 83.5875
#> mean(max.age.extant)
#[1] 82.98902
#> mean(ten.percent)
#[1] 0.01356247
#> colMeans(state.freqs)
#0       1
#750.375 249.625


#TESTING SCENARIO 3:

prefix <- paste("/Users/jmbeauli/hisse_fossil_sims/Two_Rate/Scenario3_redone_recons/1000_taxa/scenario3.1000.*")
files <- system(paste("ls -1 ", prefix, ".*", sep=""), intern=TRUE)

max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()

for(index in 1:length(files)){
    load(files[index])
    max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
    max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
    max.age.extant <- c(max.age.extant, max(dateNodes(model.list.extant[[1]]$phy, rootAge=max(node.depth.edgelength(model.list.extant[[1]]$phy)))[-c(1:Ntip(model.list.extant[[1]]$phy))]))
}

#> mean(max.age.true)
#[1] 110.8972
#>  mean(max.age.sampled)
#[1] 99.91218
#> mean(max.age.extant)
#[1] 95.04251

max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()
transitions.raw <- c()
ten.percent <- c()
state.freqs <- c()
ntax=1000
true.psi=0.005

for(index in 1:20){
    print(index)
    #Net div and turnover is higher by ~2x in rate class 2
    pars <- c(0.375, 0.375, 0.75, 0.75)
    all.dead <- TRUE
    while(all.dead){
        try(sim.tab <- hisse::SimulateHisse(turnover=c(0.375,0.375), eps=c(0.75,0.75), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
        if(sim.tab$n.surviving > 0){
            try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
            if(Ntip(phy) > ntax){
                f <- hisse:::GetFossils(phy, psi=true.psi)
                if(sum(f$has_sampled_descendant) > 0){
                    try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
                    if(!is.null(test.tol)){
                        if(Ntip(test.tol) == ntax){
                            pp <- hisse:::ProcessSimSample(phy, f)
                            sf <- 1
                            all.dead <- FALSE
                            max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
                            max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
                            max.age.extant <- c(max.age.extant, max(dateNodes(test.tol, rootAge=max(node.depth.edgelength(test.tol)))[-c(1:Ntip(test.tol))]))
                            transitions.raw <- rbind(transitions.raw, sim.tab$transition.counts[!is.na(sim.tab$transition.counts)])
                            ten.percent <- c(ten.percent, sum(sim.tab$transition.counts[!is.na(sim.tab$transition.counts)]) / length(phy$edge.length))
                            true.tip.states <- phy$tip.state[names(phy$tip.state) %in% test.tol$tip.label]
                            state.freqs <- rbind(state.freqs, table(true.tip.states))
                        }
                    }
                }
            }
        }
    }
}

#> mean(max.age.true)
#[1] 110.5973
#> mean(max.age.sampled)
#[1] 100.2246
#> mean(max.age.extant)
#[1] 97.05678
#> mean(ten.percent)
#[1] 0.01153843
#> colMeans(state.freqs)
#0     1
#648.5 351.5


max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()
transitions.raw <- c()
ten.percent <- c()
state.freqs <- c()
k.foss <- c()
m.foss <- c()
ntax=200
true.psi=0.1

for(index in 1:20){
    print(index)
    #Net div and turnover is higher by ~2x in rate class 2
    pars <- c(0.25, 0.50, 0.55, 0.75)
    all.dead <- TRUE
    phy <- NULL
    while(all.dead){
        try(sim.tab <- hisse::SimulateHisse(turnover=c(0.375,0.375), eps=c(0.75,0.75), max.taxa=ntax, x0=0, transition.rates=matrix(c(NA, 0.005, 0.005, NA), nrow=2), nstart=2))
        if(sim.tab$n.surviving > 0){
            try(phy <- SimToPhylo(sim.tab, include.extinct=TRUE))
            if(Ntip(phy) > ntax){
                f <- hisse:::GetFossils(phy, psi=true.psi)
                if(sum(f$has_sampled_descendant) > 0){
                    try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
                    if(!is.null(test.tol)){
                        if(Ntip(test.tol) == ntax){
                            pp <- hisse:::ProcessSimSample(phy, f)
                            sf <- 1
                            all.dead <- FALSE
                            max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
                            max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
                            max.age.extant <- c(max.age.extant, max(dateNodes(test.tol, rootAge=max(node.depth.edgelength(test.tol)))[-c(1:Ntip(test.tol))]))
                            transitions.raw <- rbind(transitions.raw, sim.tab$transition.counts[!is.na(sim.tab$transition.counts)])
                            ten.percent <- c(ten.percent, sum(sim.tab$transition.counts[!is.na(sim.tab$transition.counts)]) / length(phy$edge.length))
                            true.tip.states <- phy$tip.state[names(phy$tip.state) %in% test.tol$tip.label]
                            k.foss <- c(k.foss, dim(pp$k.samples)[1])
                            m.foss <- c(m.foss, Ntip(pp$phy)-ntax)
                            state.freqs <- rbind(state.freqs, table(true.tip.states))
                        }
                    }
                }
            }
        }
    }
}

#0.1
#> mean(k.foss)
#[1] 201.95
#> mean(m.foss)
#[1] 140.75
#> mean(max.age.true)
#[1] 74.58075
#> mean(max.age.extant)
#[1] 60.63529
#> colMeans(state.freqs)
#    0     1
#149.3  50.7
#> mean(ten.percent)
#[1] 0.01147392



#TESTING SCENARIO 4:

prefix <- paste("/Users/jmbeauli/hisse_fossil_sims/Two_Rate/Scenario4_redone_recons/1000_taxa/scenario4.1000.*")
files <- system(paste("ls -1 ", prefix, ".*", sep=""), intern=TRUE)

max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()

for(index in 1:length(files)){
    load(files[index])
    max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
    max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
    max.age.extant <- c(max.age.extant, max(dateNodes(model.list.extant[[1]]$phy, rootAge=max(node.depth.edgelength(model.list.extant[[1]]$phy)))[-c(1:Ntip(model.list.extant[[1]]$phy))]))
}

#> mean(max.age.true)
#[1] 288.7356
#> mean(max.age.sampled)
#[1] 241.7721
#> mean(max.age.extant)
#[1] 206.3047



#TESTING SCENARIO 4 no mass extinction:

max.age.true <- c()
max.age.sampled <- c()
max.age.extant <- c()
transitions.raw <- c()
ten.percent <- c()
state.freqs <- c()
ntax=1000
true.psi=0.005

for(index in 1:20){
    print(index)
    #Realistic with mass extinction and heritable sampling
    pars <- c(0.5, 0.2, 0.7, 0.9)
    turnover <- c(0.1, 0.3, 1.0, 0.6)
    eps <- c(0.7, 0.9, 0.95, 0.8)
    transition <- 0.1*matrix(c(NA, 0.5, 0.2, 0.6, 0.1, NA, 0.3, 0.2, 0.1, 0.4, NA, 0.8, 0.1, 0.2, 0.3, NA),nrow=4)
    
    sim.tab <- list()
    sim.tab$n.surviving = 0
    while(sim.tab$n.surviving==0) {
        sim.tab <- hisse::SimulateHisse(turnover.rates=turnover, eps.values=eps, transition.rates=transition, max.taxa=ntax/.75, x0=0, nstart=2)
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
    if(Ntip(phy)>ntax){
        f <- hisse:::GetFossils(phy, psi=true.psi)
        if(sum(f$has_sampled_descendant) > 0){
            try(test.tol <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))
            if(!is.null(test.tol)){
                if(Ntip(test.tol) == ntax){
                    pp <- hisse:::ProcessSimSample(phy, f)
                    sf <- 0.75
                    all.dead <- FALSE
                    max.age.true <- c(max.age.true, max(dateNodes(phy, rootAge=max(node.depth.edgelength(phy)))[-c(1:Ntip(phy))]))
                    max.age.sampled <- c(max.age.sampled, max(dateNodes(pp$phy, rootAge=max(node.depth.edgelength(pp$phy)))[-c(1:Ntip(pp$phy))]))
                    max.age.extant <- c(max.age.extant, max(dateNodes(test.tol, rootAge=max(node.depth.edgelength(test.tol)))[-c(1:Ntip(test.tol))]))
                    transitions.raw <- rbind(transitions.raw, sim.tab$transition.counts[!is.na(sim.tab$transition.counts)])
                    ten.percent <- c(ten.percent, sum(sim.tab$transition.counts[!is.na(sim.tab$transition.counts)]) / length(phy$edge.length))
                    true.tip.states <- phy$tip.state[names(phy$tip.state) %in% test.tol$tip.label]
                    state.freqs <- rbind(state.freqs, table(true.tip.states))
                }
            }
        }
    }
}

# 17 reps only:
#> mean(max.age.true)
#[1] 235.1108
#> mean(max.age.sampled)
#[1] 219.7112
#> mean(max.age.extant)
#[1] 201.7819



#TESTING SCENARIO 4 with mass extinction:


# 17 reps only:
#> mean(max.age.true)
#[1] 235.1108
#> mean(max.age.sampled)
#[1] 219.7112
#> mean(max.age.extant)
#[1] 201.7819



### Sim tree figure:
phy <- tree.set[[2]]
f <- hisse:::GetFossils(phy, psi=0.05)

phy.full <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50)
sample.prob.raw <- abs(geiger::sim.char(phy.full, par=1, model="BM"))[,1,1]
sample.prob.raw[which(is.na(sample.prob.raw))] <- median(sample.prob.raw, na.rm=TRUE )
sample.prob <- sample.prob.raw/max(sample.prob.raw)
while(ape::Ntip(geiger::drop.extinct(phy, tol=.Machine$double.eps^.50))>ntax) {
    to.kill <- sample(names(sample.prob), size=1, prob=sample.prob)
    phy <- ape::drop.tip(phy, to.kill)
    phy$tip.state <- phy$tip.state[-which(names(phy$tip.state)==to.kill)]
    sample.prob <- sample.prob[-which(names(sample.prob)==to.kill)]
}


phy.extant.only <- geiger::drop.extinct(phy, tol=.Machine$double.eps^.50)
phy <- tree.set[[2]]
trait <- cbind(phy$tip.label, 0)
trait[which(phy$tip.label %in% phy.extant.only$tip.label),2] <- 1
trait.final <- data.frame(taxon=trait[,1], trait <- as.numeric(trait[,2]))

library(picante)
rownames(trait.final) <- trait.final[,1]
phy0 <- picante::match.phylo.data(phy, trait.final)$phy
trait0 <- picante::match.phylo.data(phy, trait.final)$data
mode <- trait0[,2]
names(mode) <- trait0[,1]
library(viridis)


getphylo_x <- function(tree, node) {
    if(is.character(node)) {
        node <- which(c(tree$tip.label, tree$node.label)==node)
    }
    pi <- tree$edge[tree$edge[,2]==node, 1]
    if (length(pi)) {
        ei <- which(tree$edge[,1]==pi & tree$edge[,2]==node)
        tree$edge.length[ei] + Recall(tree, pi)
    } else {
        if(!is.null(tree$root.edge)) {
            tree$root.edge
        } else {
            0
        }
    }
}


getphylo_y <- function(tree, node) {
    if(is.character(node)) {
        node <- which(c(tree$tip.label, tree$node.label)==node)
    }
    ci <- tree$edge[tree$edge[,1]==node, 2]
    if (length(ci)==2) {
        mean(c(Recall(tree, ci[1]), Recall(tree, ci[2])))
    } else if (length(ci)==0) {
        Ntip <- length(tree$tip.label)
        which(tree$edge[tree$edge[, 2] <= Ntip, 2] == node)
    } else {
        stop(paste("error", length(ci)))
    }
}


library(phytools)

pdf("TEST2.pdf", width=6, height=12)
colors_states <- c("transparent", viridis(4)[1])
# Plot states at tips
plot(ladderize(phy, right=FALSE), show.tip.label=F, edge.width=0.3, adj=1, cex=0.08)
abline(v=40, col="blue", lty=3, lwd=1)
abline(v=70, col="blue", lty=3, lwd=1)
par(fg="transparent")
AddFossilPoints(f, ladderize(phy, right=FALSE), pch=19, cex=0.8, col=viridis(4)[2])
tiplabels(pie=to.matrix(mode, sort(unique(mode))),piecol=colors_states, cex=0.25, lwd=0.3, frame = "n")
par(fg="black")
#tiplabels(pch=1, bg=tip.cols, adj=1, cex=0.1, width = 0.1)
axisPhylo(cex.axis=.75, lwd=.75)
dev.off()


files.fullset <- list.files(pattern="Rsave")
for(i in 1:length(files.fullset)){
    load(files.fullset[i])
    loglik <- misse.list[[1]]$loglik
    if(round(loglik,2) == -6220.64){
        print(files.fullset[i])
    }
}


pdf("TEST3.pdf")
ss <-ape::ltt.plot.coords(phy)
plot(ss, axes=FALSE, ylab="", xlab="", ylim=c(0,300), xlim=c(-120,0), col="white")
lines(ss)
par(tck=.01)
axis(2, at = seq(0, 300, by = 50), las=1, lwd=1, cex.axis=1, labels=TRUE, mgp=c(.75,.5,0))
axis(1, at = seq(-120, 0, by = 20), las=1, lwd=1, cex.axis=1, labels=c(120,100,80,60,40,20,0), mgp=c(.75,.5,0))
title(ylab="N", line=2.5, cex=2)
title(xlab="Time", line=2.5)
mtext("A",side=3, line=0, adj=0, cex=1.5)
dev.off()


